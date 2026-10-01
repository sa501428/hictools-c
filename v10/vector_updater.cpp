#include "vector_updater.h"
#include "repack.h"
#include <sys/stat.h>
#include <unistd.h>
#include <zstd.h>

namespace hic10 {
namespace {
uint64_t get(const Bytes &b, size_t at, unsigned n) {
    check(at <= b.size() && n <= b.size() - at, "truncated V10 vector field");
    uint64_t value = 0;
    for (unsigned i = 0; i < n; ++i) value |= uint64_t(b[at + i]) << (8 * i);
    return value;
}
void set(Bytes &b, size_t at, uint64_t value, unsigned n) {
    check(at <= b.size() && n <= b.size() - at, "truncated V10 vector field");
    for (unsigned i = 0; i < n; ++i) b[at + i] = value >> (8 * i);
}
}
class VectorOutputImpl {
    FILE *file_ = nullptr;
    std::string target_, temporary_;
    int level_;
    Header header_;
    Bytes footer_bytes_;
    struct Entry {
        std::tuple<uint32_t, uint32_t, uint8_t, uint32_t> key;
        Bytes bytes;
    };
    std::array<std::vector<Entry>, 3> entries_;
    uint64_t position() {
        auto p = ftello(file_);
        check(p >= 0, "cannot determine V10 output position");
        return p;
    }
    void write(const Bytes &b) {
        check(std::fwrite(b.data(), 1, b.size(), file_) == b.size(), "cannot write V10 output");
    }
    Bytes compress(const Bytes &raw) {
        Bytes out(ZSTD_compressBound(raw.size()));
        size_t n = ZSTD_compress(out.data(), out.size(), raw.data(), raw.size(), level_);
        check(!ZSTD_isError(n), "Zstandard compression failure");
        out.resize(n);
        return out;
    }

  public:
    VectorOutputImpl(Reader &reader, const Header &header, const std::string &path, int level)
        : target_(path), level_(level), header_(header) {
        temporary_ = path + ".addnorm.XXXXXX";
        std::vector<char> name(temporary_.begin(), temporary_.end());
        name.push_back(0);
        int fd = mkstemp(name.data());
        check(fd >= 0, "cannot create addnorm temporary file");
        temporary_ = name.data();
        struct stat metadata {};
        if (stat(path.c_str(), &metadata) == 0 && fchmod(fd, metadata.st_mode & 07777) != 0) {
            close(fd);
            std::remove(temporary_.c_str());
            temporary_.clear();
            check(false, "cannot preserve V10 file permissions");
        }
        file_ = fdopen(fd, "w+b");
        if (!file_) {
            close(fd);
            std::remove(temporary_.c_str());
            temporary_.clear();
            check(false, "cannot open addnorm temporary file");
        }
        try {
            footer_bytes_ = repack_matrix_prefix(file_, reader, header);
            for (auto v : reader.vector_entries()) {
                for (uint32_t j = 0; j < v.chunks; ++j) {
                    size_t at = v.descriptors + uint64_t(j) * 32;
                    uint64_t old = get(v.bytes, at + 16, 8);
                    uint32_t length = narrow(get(v.bytes, at + 24, 4));
                    uint64_t pos = position();
                    // Stored chunk bytes and all scale-factor bits stay exact.
                    constexpr uint64_t batch = 8 * 1024 * 1024;
                    for (uint64_t offset = 0; offset < length; offset += batch)
                        write(reader.read_bytes(old + offset, std::min(batch, uint64_t(length) - offset)));
                    set(v.bytes, at + 16, pos, 8);
                }
                entries_[v.kind].push_back({{v.norm, v.chr, v.unit, v.ri}, std::move(v.bytes)});
            }
        } catch (...) {
            std::fclose(file_); file_ = nullptr;
            std::remove(temporary_.c_str()); temporary_.clear();
            throw;
        }
    }
    ~VectorOutputImpl() {
        if (file_)
            std::fclose(file_);
        if (!temporary_.empty())
            std::remove(temporary_.c_str());
    }
    bool contains(uint8_t kind, uint32_t norm, uint32_t chr, uint8_t unit, uint32_t ri) const {
        auto key = std::make_tuple(kind == 1 ? 0 : norm, kind == 0 ? chr : 0, unit, ri);
        for (const auto &e : entries_[kind]) if (e.key == key) return true;
        return false;
    }
    void add(const Vector &v) {
        check(v.kind <= 2 && v.unit <= 1, "invalid new V10 vector kind/unit");
        check(v.ri < header_.resolutions[v.unit].size() &&
                  (v.kind == 1 || v.norm < header_.norms.size()) &&
                  (v.kind != 0 || v.chr < header_.chromosomes.size()),
              "new V10 vector key outside header");
        uint64_t expected = 0;
        if (v.kind == 0) expected = header_.bins(v.chr, v.unit, v.ri);
        else for (uint32_t chr = 0; chr < header_.chromosomes.size(); ++chr)
            expected = std::max(expected, header_.bins(chr, v.unit, v.ri));
        check(v.value_count() == expected, "new V10 vector length mismatch");
        for (auto scale : v.scales)
            check(scale.first < header_.chromosomes.size(), "new V10 expected scale outside header");
        if (contains(v.kind, v.norm, v.chr, v.unit, v.ri)) return;
        Bytes descriptors;
        constexpr uint32_t nominal = 65536;
        for (uint64_t begin = 0; begin < v.value_count(); begin += nominal) {
            uint32_t n = narrow(std::min<uint64_t>(nominal, v.value_count() - begin));
            auto words = v.loader ? v.loader(begin, n)
                                  : std::vector<uint32_t>(v.values.begin() + begin, v.values.begin() + begin + n);
            check(words.size() == n, "short V10 vector loader read");
            Bytes best;
            uint8_t transform = 0;
            for (uint8_t t = 0; t < 3; ++t) {
                Bytes raw;
                raw.reserve(uint64_t(n) * 4);
                if (t == 1)
                    for (unsigned lane = 0; lane < 4; ++lane)
                        for (uint32_t j = 0; j < n; ++j)
                            put(raw, words[j] >> (8 * lane), 1);
                else
                    for (uint32_t j = 0; j < n; ++j)
                        put(raw, words[j] ^ (t == 2 && j ? words[j - 1] : 0),
                            4);
                auto frame = compress(raw);
                if (best.empty() || frame.size() < best.size()) {
                    best = std::move(frame);
                    transform = t;
                }
            }
            Bytes stored;
            magic(stored, "H10V");
            put(stored, 1, 1);
            put(stored, transform, 1);
            put(stored, 0, 2);
            put(stored, uint64_t(n) * 4, 4);
            put(stored, n, 4);
            append(stored, best);
            uint64_t pos = position();
            write(stored);
            put(descriptors, begin, 8);
            put(descriptors, n, 4);
            put(descriptors, transform, 1);
            put(descriptors, 1, 1);
            put(descriptors, 0, 2);
            put(descriptors, pos, 8);
            put(descriptors, stored.size(), 4);
            put(descriptors, uint64_t(n) * 4, 4);
        }
        Bytes entry;
        put(entry, 0, 4);
        if (v.kind != 1)
            put(entry, v.norm, 4);
        if (v.kind == 0)
            put(entry, v.chr, 4);
        put(entry, v.unit, 1);
        put(entry, 0, 3);
        put(entry, v.ri, 4);
        put(entry, 0, 4);
        put(entry, v.value_count(), 8);
        put(entry, nominal, 4);
        put(entry, descriptors.size() / 32, 4);
        if (v.kind) {
            put(entry, v.scales.size(), 4);
            put(entry, 0, 4);
            for (auto s : v.scales) {
                put(entry, s.first, 4);
                put(entry, s.second, 4);
            }
        }
        append(entry, descriptors);
        entries_[v.kind].push_back(
            {{v.kind == 1 ? 0 : v.norm, v.kind == 0 ? v.chr : 0, v.unit, v.ri}, std::move(entry)});
    }
    void finish(const Header &header) {
        for (uint8_t kind = 0; kind < 3; ++kind) {
            auto &list = entries_[kind];
            if (list.empty())
                continue;
            std::sort(list.begin(), list.end(), [](auto &a, auto &b) { return a.key < b.key; });
            Bytes index;
            magic(index, kind == 0 ? "NVI0" : kind == 1 ? "EVI0" : "NEVI");
            put(index, 1, 4);
            put(index, list.size(), 4);
            put(index, 0, 4);
            for (auto &item : list) {
                auto &entry = item.bytes;
                uint8_t unit = std::get<2>(item.key);
                uint32_t ri = std::get<3>(item.key);
                uint32_t bin = header.resolutions[unit][ri].bin;
                size_t binOffset = 4 + (kind != 1 ? 4 : 0) + (kind == 0 ? 4 : 0) + 4 + 4;
                for (unsigned j = 0; j < 4; ++j)
                    entry[binOffset + j] = static_cast<uint8_t>(bin >> (8 * j));
                uint32_t length = narrow(entry.size());
                for (unsigned j = 0; j < 4; ++j)
                    entry[j] = static_cast<uint8_t>(length >> (8 * j));
                append(index, entry);
            }
            uint64_t pos = position();
            write(index);
            Bytes loc;
            put(loc, pos, 8);
            put(loc, index.size(), 8);
            auto saved = position();
            check(fseeko(file_, 32 + 16 * kind, SEEK_SET) == 0, "cannot patch V10 vector locator");
            write(loc);
            check(fseeko(file_, saved, SEEK_SET) == 0, "cannot restore V10 output position");
        }
        uint64_t footer_position = position();
        write(footer_bytes_);
        Bytes footer_locator;
        put(footer_locator, footer_position, 8);
        put(footer_locator, footer_bytes_.size(), 8);
        auto saved = position();
        check(fseeko(file_, 16, SEEK_SET) == 0, "cannot patch V10 footer locator");
        write(footer_locator);
        check(fseeko(file_, saved, SEEK_SET) == 0, "cannot restore V10 output position");
        check(std::fflush(file_) == 0 && fsync(fileno(file_)) == 0,
              "cannot sync normalized V10 file");
        auto f = file_;
        file_ = nullptr;
        check(std::fclose(f) == 0, "cannot close normalized V10 file");
        check(std::rename(temporary_.c_str(), target_.c_str()) == 0,
              "cannot replace normalized V10 file");
        temporary_.clear();
    }
};

VectorOutput::VectorOutput(Reader &r, const Header &h, const std::string &path, int level)
    : impl_(new VectorOutputImpl(r, h, path, level)) {}
VectorOutput::~VectorOutput() = default;
bool VectorOutput::contains(uint8_t kind, uint32_t norm, uint32_t chr, uint8_t unit, uint32_t ri) const {
    return impl_->contains(kind, norm, chr, unit, ri);
}
void VectorOutput::add(const Vector &v) { impl_->add(v); }
void VectorOutput::finish(const Header &h) { impl_->finish(h); }
} // namespace hic10
