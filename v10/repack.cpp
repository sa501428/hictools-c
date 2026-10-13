#include "repack.h"
#include "common/thread_pool.h"
#include <cerrno>
#include <deque>
#include <unistd.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <sys/stat.h>

namespace hic10 {
namespace {
void write_bytes(FILE *output, const Bytes &bytes) {
    check(std::fwrite(bytes.data(), 1, bytes.size(), output) == bytes.size(),
          "cannot write repacked V10 file");
}
uint64_t position(FILE *output) {
    auto result = ftello(output);
    check(result >= 0, "cannot determine repacked V10 position");
    return result;
}
uint32_t get_word(const Bytes &bytes, size_t at) {
    check(at <= bytes.size() && 4 <= bytes.size() - at, "truncated V10 u32 field");
    uint32_t value = 0;
    for (unsigned i = 0; i < 4; ++i)
        value |= uint32_t(bytes[at + i]) << (8 * i);
    return value;
}
uint64_t get_wide(const Bytes &bytes, size_t at) {
    check(at <= bytes.size() && 8 <= bytes.size() - at, "truncated V10 u64 field");
    uint64_t value = 0;
    for (unsigned i = 0; i < 8; ++i)
        value |= uint64_t(bytes[at + i]) << (8 * i);
    return value;
}
void set_wide(Bytes &bytes, size_t at, uint64_t value) {
    check(at <= bytes.size() && 8 <= bytes.size() - at, "truncated V10 u64 field");
    for (unsigned i = 0; i < 8; ++i)
        bytes[at + i] = static_cast<uint8_t>(value >> (8 * i));
}
Bytes serialize_header(const Header &header) {
    Bytes bytes;
    magic(bytes, "HIC\0");
    put(bytes, 10, 4);
    bytes.resize(88, 0); // Includes reserved footer/NVI/EVI/NEVI locator space.
    str(bytes, header.genome);
    put(bytes, header.attributes.size(), 4);
    for (const auto &attribute : header.attributes) {
        str(bytes, attribute.first);
        str(bytes, attribute.second);
    }
    put(bytes, header.chromosomes.size(), 4);
    for (const auto &chromosome : header.chromosomes) {
        str(bytes, chromosome.name);
        put(bytes, chromosome.length, 8);
    }
    for (const auto &resolutions : header.resolutions) {
        put(bytes, resolutions.size(), 4);
        for (const auto &resolution : resolutions) {
            put(bytes, resolution.bin, 4);
            put(bytes, resolution.mode, 1);
            put(bytes, resolution.aggregation, 1);
            put(bytes, 0, 2);
            put(bytes, resolution.source, 4);
        }
    }
    if (!header.resolutions[1].empty())
        for (const auto &chromosome : header.chromosomes) {
            put(bytes, chromosome.sites.size(), 4);
            for (uint64_t site : chromosome.sites)
                put(bytes, site, 8);
        }
    put(bytes, header.norms.size(), 4);
    for (const auto &name : header.norms)
        str(bytes, name);
    set_wide(bytes, 8, bytes.size());
    return bytes;
}
} // namespace

Bytes copy_matrix_storage(FILE *output, Reader &reader, ThreadPool *pool) {

    // Copy only referenced matrix records, indexes, and blocks. V10 permits
    // arbitrary physical ordering, including vectors interleaved with matrices.
    std::vector<FileLocator> intervals;
    auto fields = reader.matrix_relocation_fields(&intervals);
    std::sort(intervals.begin(), intervals.end(),
              [](auto a, auto b) { return a.position < b.position; });
    std::map<uint64_t, uint64_t> relocated;
    uint64_t previous_end = reader.header_length(), next_output = position(output);
    for (auto interval : intervals) {
        check(interval.position >= previous_end && interval.position <= reader.file_size() &&
                  interval.length <= reader.file_size() - interval.position,
              "overlapping or invalid V10 matrix interval");
        relocated[interval.position] = next_output;
        next_output = plus(next_output, interval.length);
        previous_end = plus(interval.position, interval.length);
    }
    // Patch pointers while copying records in batches rather than issuing an
    // independent disk seek/read/write for every block locator in a large file.
    constexpr uint64_t chunk = 8 * 1024 * 1024;
    size_t field_index = 0;
    std::deque<std::future<void>> pending;
    if (pool)
        check(std::fflush(output) == 0, "cannot flush repack header");
    const int fd = fileno(output);
    try {
        for (auto interval : intervals) {
            for (uint64_t offset = 0; offset < interval.length;) {
                uint64_t at = interval.position + offset;
                uint64_t count = std::min(chunk, interval.length - offset);
                auto end = std::lower_bound(fields.begin() + field_index, fields.end(), at + count);
                if (end != fields.begin() + field_index && *(end - 1) + 8 > at + count)
                    count = *(end - 1) - at;
                check(count > 0, "split relocation field outside matrix interval");
                size_t first = field_index;
                while (field_index < fields.size() && fields[field_index] < at + count)
                    ++field_index;
                size_t last = field_index;
                uint64_t output_at = relocated.at(interval.position) + offset;
                auto copy = [&, at, count, first, last, output_at, fd]() {
                    auto bytes =
                        pool ? reader.read_bytes_at(at, count) : reader.read_bytes(at, count);
                    for (size_t j = first; j < last; ++j) {
                        uint64_t field = fields[j];
                        check(field >= at && field - at + 8 <= count,
                              "V10 relocation field outside copied matrix record");
                        const auto &positions = relocated;
                        auto found = positions.find(get_wide(bytes, field - at));
                        check(found != positions.end(),
                              "V10 matrix locator outside copied storage");
                        set_wide(bytes, field - at, found->second);
                    }
                    if (!pool)
                        write_bytes(output, bytes);
                    else
                        for (size_t done = 0; done < bytes.size();) {
                            auto n = pwrite(fd, bytes.data() + done, bytes.size() - done,
                                            output_at + done);
                            if (n < 0 && errno == EINTR)
                                continue;
                            check(n > 0, "cannot copy repacked V10 matrix range");
                            done += n;
                        }
                };
                if (!pool)
                    copy();
                else {
                    pending.push_back(pool->submit(copy));
                    if (pending.size() >= pool->num_threads()) {
                        pending.front().get();
                        pending.pop_front();
                    }
                }
                offset += count;
            }
        }
        while (!pending.empty()) {
            pending.front().get();
            pending.pop_front();
        }
    } catch (...) {
        // Reader, relocation map and output must outlive every scheduled copy.
        for (auto &task : pending)
            if (task.valid())
                task.wait();
        throw;
    }
    if (pool)
        check(fseeko(output, next_output, SEEK_SET) == 0, "cannot seek repacked matrix end");
    check(field_index == fields.size(), "uncopied V10 matrix relocation fields");

    Bytes footer = reader.read_bytes(reader.footer().position, reader.footer().length);
    check(footer.size() >= 24 && std::equal(footer.begin(), footer.begin() + 4, "H10F") &&
              get_wide(footer, 8) == footer.size(),
          "invalid V10 footer during repack");
    uint32_t matrices = get_word(footer, 16);
    check(matrices == reader.matrices().size() && footer.size() == 24 + uint64_t(matrices) * 24,
          "invalid V10 footer matrix count");
    for (uint32_t i = 0; i < matrices; ++i) {
        size_t field = 24 + uint64_t(i) * 24 + 8;
        auto found = relocated.find(get_wide(footer, field));
        check(found != relocated.end(), "V10 footer locator outside copied storage");
        set_wide(footer, field, found->second);
    }
    return footer;
}
Bytes repack_matrix_prefix(FILE *output, Reader &reader, const Header &header, ThreadPool *pool,
                           Reader *extra) {
    write_bytes(output, serialize_header(header));
    Bytes footer = copy_matrix_storage(output, reader, pool);
    if (extra) {
        Bytes added = copy_matrix_storage(output, *extra, pool);
        uint32_t old_count = get_word(footer, 16), new_count = get_word(added, 16);
        // Extra pairs are appended after all original chromosome IDs.
        if (old_count && new_count)
            check(get_word(footer, footer.size() - 24) < get_word(added, 24),
                  "appended matrix directory is not sorted");
        footer.insert(footer.end(), added.begin() + 24, added.end());
        set_wide(footer, 8, footer.size());
        uint32_t count = narrow(uint64_t(old_count) + new_count);
        for (unsigned i = 0; i < 4; ++i) footer[16 + i] = count >> (8 * i);
    }
    return footer;
}
void append_matrix_atomic(Reader &reader, const Header &header, const std::string &path,
                          uint32_t threads, Reader *extra,
                          const std::function<void(const std::string &)> &validate) {
    check(threads > 0 && threads <= 256, "All repair copy threads must be 1..256");
    std::string pattern = path + ".add-all.XXXXXX";
    std::vector<char> name(pattern.begin(), pattern.end()); name.push_back(0);
    int fd = mkstemp(name.data());
    check(fd >= 0, "cannot create All repair temporary file");
    struct Cleanup {
        std::string path;
        FILE *file = nullptr;
        ~Cleanup() { if (file) std::fclose(file); if (!path.empty()) std::remove(path.c_str()); }
    } cleanup{name.data(), nullptr};
    struct stat target{};
    if (stat(path.c_str(), &target) == 0 && fchmod(fd, target.st_mode & 07777) != 0) {
        close(fd); check(false, "cannot preserve All repair output permissions");
    }
    cleanup.file = fdopen(fd, "w+b");
    if (!cleanup.file) { close(fd); check(false, "cannot open All repair temporary file"); }
    FILE *output = cleanup.file;
    auto bytes = serialize_header(header);
    check(bytes.size() >= reader.header_length(), "All repair cannot shrink the header");
    uint64_t delta = bytes.size() - reader.header_length();
    uint64_t end = plus(reader.file_size(), delta);
    check(end <= uint64_t(INT64_MAX), "All repair output exceeds signed 64-bit file offsets");
    auto write_at = [&](uint64_t at, const Bytes &data) {
        check(at <= uint64_t(INT64_MAX) && data.size() <= uint64_t(INT64_MAX) - at,
              "All repair output offset overflow");
        for (size_t done = 0; done < data.size();) {
            auto n = pwrite(fd, data.data() + done, data.size() - done, at + done);
            if (n < 0 && errno == EINTR) continue;
            check(n > 0, "cannot write All repair output (check space/quota)");
            done += n;
        }
    };
    write_at(0, bytes);
    constexpr uint64_t chunk = 8 * 1024 * 1024;
    uint64_t start = reader.header_length(), copied = 0;
    auto begun = std::chrono::steady_clock::now(), last = begun;
    std::cerr << "All repair: copying " << reader.file_size() - start
              << " bytes with " << threads << " workers; no per-block relocation map\n";
    {
        ThreadPool pool(threads);
        std::deque<std::pair<std::future<void>, uint64_t>> pending;
        auto finish_one = [&]() {
            pending.front().first.get(); copied += pending.front().second; pending.pop_front();
            auto now = std::chrono::steady_clock::now();
            if (now - last >= std::chrono::seconds(30)) {
                double seconds = std::chrono::duration<double>(now - begun).count();
                std::cerr << "All repair: copied " << copied << '/' << reader.file_size() - start
                          << " bytes, " << copied / seconds / 1e6 << " MB/s\n";
                last = now;
            }
        };
        // pool is destroyed (joining all workers) before captures/output on any exception.
        for (uint64_t at = start; at < reader.file_size();) {
            uint64_t count = std::min(chunk, reader.file_size() - at);
            pending.emplace_back(pool.submit([&, at, count]() {
                auto data = reader.read_bytes_at(at, count);
                write_at(plus(at, delta), data);
            }), count);
            if (pending.size() >= threads) finish_one();
            at += count;
        }
        while (!pending.empty()) finish_one();
    }
    std::cerr << "All repair: streaming matrix and vector pointer updates\n";
    reader.relocate_matrix_pointers(delta, write_at);
    reader.for_each_vector_entry([&](uint64_t at, VectorEntry v) {
        for (uint32_t j = 0; j < v.chunks; ++j) {
            size_t field = v.descriptors + uint64_t(j) * 32 + 16;
            set_wide(v.bytes, field, plus(get_wide(v.bytes, field), delta));
        }
        write_at(plus(at, delta), v.bytes);
    });
    auto footer = reader.read_bytes(reader.footer().position, reader.footer().length);
    for (size_t field = 32; field < footer.size(); field += 24)
        set_wide(footer, field, plus(get_wide(footer, field), delta));
    check(fseeko(output, end, SEEK_SET) == 0, "cannot seek appended All matrix");
    if (extra) {
        auto added = copy_matrix_storage(output, *extra, nullptr);
        uint32_t old_count = get_word(footer, 16), new_count = get_word(added, 16);
        if (old_count && new_count)
            check(get_word(footer, footer.size() - 24) < get_word(added, 24),
                  "appended All matrix directory is not sorted");
        footer.insert(footer.end(), added.begin() + 24, added.end());
        set_wide(footer, 8, footer.size());
        uint32_t count = narrow(uint64_t(old_count) + new_count);
        for (unsigned j = 0; j < 4; ++j) footer[16 + j] = count >> (8 * j);
    }
    uint64_t footer_at = position(output);
    write_bytes(output, footer);
    set_wide(bytes, 16, footer_at); set_wide(bytes, 24, footer.size());
    for (unsigned kind = 0; kind < 3; ++kind) {
        auto loc = reader.vector_indexes()[kind];
        set_wide(bytes, 32 + kind * 16, loc.length ? plus(loc.position, delta) : 0);
        set_wide(bytes, 40 + kind * 16, loc.length);
    }
    check(std::fflush(output) == 0, "cannot flush All repair output");
    write_at(0, bytes);
    std::cerr << "All repair: syncing and validating before atomic publication\n";
    check(fsync(fd) == 0, "cannot sync All repair output");
    int closed = std::fclose(output); cleanup.file = nullptr;
    check(closed == 0, "cannot close All repair output");
    validate(cleanup.path);
    reader.assert_unchanged();
    check(std::rename(cleanup.path.c_str(), path.c_str()) == 0,
          "cannot atomically publish All repair output");
    cleanup.path.clear();
    double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begun).count();
    std::cerr << "All repair: published " << path << " in " << seconds << " seconds\n";
}
} // namespace hic10
