#include "all.h"
#include "repack.h"
#include <iostream>
#include <cctype>
#include <cmath>
#include <sys/stat.h>
#include <unistd.h>

namespace hic10 {
bool is_all(const std::string &name) {
    return name.size() == 3 && std::tolower(static_cast<unsigned char>(name[0])) == 'a' &&
        std::tolower(static_cast<unsigned char>(name[1])) == 'l' &&
        std::tolower(static_cast<unsigned char>(name[2])) == 'l';
}
AllBuilder::AllBuilder(Header &h, const AllOptions &o) {
    check(o.bins > 0 && o.bins <= 4096 && o.source > 0, "invalid All overview options");
    check(!h.resolutions[0].empty(), "All requires a BP resolution");
    for (const auto &c : h.chromosomes) check(!is_all(c.name), "All already exists");
    finest_ = h.resolutions[0].front().bin;
    // Select a stored anchor no coarser than the requested source. Never decode
    // a derived matrix into memory just to build the overview.
    source_ = finest_;
    for (const auto &r : h.resolutions[0])
        if (!r.mode && r.bin <= o.source) source_ = r.bin;
    index_ = narrow(h.chromosomes.size());
    uint64_t total = 0, maximum_length = 0;
    std::string lengths;
    for (const auto &c : h.chromosomes) {
        offsets_.push_back(total);
        uint64_t n = ceil_div(c.length, source_);
        total = plus(total, n);
        maximum_length = std::max(maximum_length, c.length);
        if (!lengths.empty()) lengths += ',';
        lengths += std::to_string(n);
    }
    check(total > 0, "All requires real chromosomes");
    // Keep the largest expected-vector length unchanged, even for tiny genomes.
    uint64_t limit = std::min<uint64_t>(o.bins, std::max<uint64_t>(1, maximum_length / finest_));
    // Desktop 1.x stores All boundary coordinates in signed 32-bit arrays.
    // This is a producer sizing policy, not a new wire-format restriction.
    limit = std::min<uint64_t>(limit, std::max<uint64_t>(1, uint64_t(INT32_MAX) / finest_));
    group_ = ceil_div(total, limit);
    check(group_ <= UINT64_MAX / source_, "All physical scale overflow");
    uint64_t bins = ceil_div(total, group_);
    check(bins <= UINT64_MAX / finest_, "All chromosome length overflow");
    h.chromosomes.push_back({"All", bins * finest_, {}});
    h.attributes.insert(h.attributes.end(), {
        {"hicAll.version", "1"}, {"hicAll.sourceResolution", std::to_string(source_)},
        {"hicAll.groupFactor", std::to_string(group_)},
        {"hicAll.finestResolution", std::to_string(finest_)},
        {"hicAll.binCount", std::to_string(bins)}, {"hicAll.chromosomeBins", lengths}});
    header_ = h;
}
void AllBuilder::add(uint32_t a, uint32_t b, uint64_t x, uint64_t y,
                     uint32_t bin, uint64_t value, bool score) {
    check(a < index_ && b < index_ && bin > 0, "invalid All source chromosome/bin");
    auto coordinate = [&](uint32_t chr, uint64_t p) {
        check(p <= UINT64_MAX / bin, "All source coordinate overflow");
        uint64_t bp = p * bin;
        check(bp < header_.chromosomes[chr].length, "All source coordinate out of bounds");
        return narrow(plus(offsets_[chr], bp / source_) / group_);
    };
    uint32_t i = coordinate(a, x), j = coordinate(b, y);
    if (i > j) std::swap(i, j);
    auto key = std::make_pair(i, j);
    if (score) {
        if (!scored_) {
            for (auto c : counts_) scores_[c.first] = static_cast<double>(c.second);
            counts_.clear();
            scored_ = true;
        }
        scores_[key] += floating(narrow(value));
    } else if (scored_) scores_[key] += static_cast<double>(value);
    else counts_[key] = plus(counts_[key], value);
}
void AllBuilder::collect(Reader &reader) {
    uint32_t ri = reader.header().resolution(0, source_);
    for (auto key : reader.matrices()) {
        if (key.chr1 >= index_ || key.chr2 >= index_) continue;
        // stream_materialized emits canonical raw cells and preserves u64 counts.
        // The type is a matrix property; Reader exposes it through this small
        // descriptor read rather than materializing the source matrix.
        auto loc = reader.matrix_location(key);
        auto desc = reader.read_bytes(loc.position + 24 + uint64_t(ri) * 76, 4);
        bool scored = desc.at(3) == 1;
        reader.stream_materialized(key.chr1, key.chr2, 0, ri, [&](const Cell &c) {
            add(key.chr1, key.chr2, c.x, c.y, source_, c.value, scored);
        }, true);
    }
}
Matrix AllBuilder::matrix(uint8_t unit, uint32_t ri) const {
    Matrix result;
    if (unit) return result; // No restriction-fragment interpretation for All.
    uint32_t bin = header_.resolutions[0].at(ri).bin;
    auto key = [&](std::pair<uint32_t, uint32_t> c) {
        return std::make_pair(narrow(uint64_t(c.first) * finest_ / bin),
                              narrow(uint64_t(c.second) * finest_ / bin));
    };
    result.scores = scored_;
    if (scored_) {
        std::map<std::pair<uint32_t, uint32_t>, double> sums;
        for (auto c : scores_) sums[key(c.first)] += c.second;
        for (auto c : sums) result.cells.push_back({c.first.first, c.first.second,
                                                  bits(static_cast<float>(c.second))});
    } else {
        std::map<std::pair<uint32_t, uint32_t>, uint64_t> sums;
        for (auto c : counts_) sums[key(c.first)] = plus(sums[key(c.first)], c.second);
        for (auto c : sums) if (c.second)
            result.cells.push_back({c.first.first, c.first.second, c.second});
    }
    return result;
}
void AllBuilder::write(Writer &writer) const {
    writer.matrix(index_, index_, [&](uint8_t u, uint32_t ri) { return matrix(u, ri); });
}
void add_all(const std::string &input, const std::string &output,
             const AllOptions &options, int level, uint32_t threads) {
    check(options.bins > 0 && options.bins <= 4096 && options.source > 0, "invalid All overview options");
    check(threads > 0 && threads <= 256, "All repair copy threads must be 1..256");
    Reader reader(input);
    for (const auto &c : reader.header().chromosomes) if (is_all(c.name)) {
        // Refuse a misleading successful repair of an existing broken All.
        uint32_t chr = narrow(&c - reader.header().chromosomes.data());
        auto loc = reader.matrix_location({chr, chr});
        check(loc.length > 0, "All header exists without an All/All matrix");
        if (input != output) {
            append_matrix_atomic(reader, reader.header(), output, threads, nullptr,
                [](const std::string &path) {
                    Reader validated(path); validated.validate_matrix_storage();
                    validated.for_each_vector_entry([](uint64_t, VectorEntry) {});
                });
        }
        return;
    }
    Header header = reader.header();
    AllBuilder all(header, options);
    std::cerr << "All repair: aggregating stored " << all.source() << " bp matrices\n";
    all.collect(reader);
    std::string pattern = output + ".all-matrix.XXXXXX";
    std::vector<char> name(pattern.begin(), pattern.end()); name.push_back(0);
    int fd = mkstemp(name.data()); check(fd >= 0, "cannot create All matrix temporary");
    close(fd);
    struct Cleanup { std::string path; ~Cleanup() { std::remove(path.c_str()); } } cleanup{name.data()};
    Options writer_options; writer_options.level = level; writer_options.threads = threads;
    writer_options.verifyDerived = false;
    // Preserve the original global storage policy, including explicit overrides.
    for (const auto &r : header.resolutions[0]) {
        if (r.mode) writer_options.derived.emplace_back(r.bin, header.resolutions[0][r.source].bin);
        else writer_options.materialized.push_back(r.bin);
    }
    { Writer writer(cleanup.path, header, writer_options); all.write(writer); writer.finish({}); }
    Reader extra(cleanup.path);
    append_matrix_atomic(reader, header, output, threads, &extra, [&](const std::string &path) {
        Reader validated(path);
        check(validated.header().chromosomes.size() == header.chromosomes.size(), "All repair header mismatch");
        validated.validate_matrix_storage();
        validated.for_each_vector_entry([](uint64_t, VectorEntry) {});
        for (uint32_t ri = 0; ri < header.resolutions[0].size(); ++ri)
            validated.matrix(all.index(), all.index(), 0, ri);
    });
}
} // namespace hic10
