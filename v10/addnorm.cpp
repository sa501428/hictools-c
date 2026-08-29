#include "addnorm.h"
#include "hic_addnorm/scale_norm.h"
#include "reader.h"
#include "repack.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <map>
#include <set>
#include <sys/stat.h>
#include <tuple>
#include <unistd.h>
#include <zstd.h>

namespace hic10 {
namespace {
struct Sparse {
    std::vector<uint32_t> row, col;
    std::vector<double> value;
    uint32_t bins = 0;
};
Sparse sparse_matrix(const Matrix &matrix, uint32_t bins) {
    Sparse s;
    s.bins = bins;
    s.row.reserve(matrix.cells.size());
    s.col.reserve(matrix.cells.size());
    s.value.reserve(matrix.cells.size());
    for (auto cell : matrix.cells) {
        double value = matrix.scores ? floating(static_cast<uint32_t>(cell.value))
                                     : static_cast<double>(cell.value);
        if (!std::isfinite(value) || value <= 0)
            continue;
        check(cell.x < bins && cell.y < bins, "normalization cell outside chromosome");
        s.row.push_back(cell.x);
        s.col.push_back(cell.y);
        s.value.push_back(value);
    }
    return s;
}
std::vector<float> raw_vc(const Sparse &s) {
    std::vector<double> sums(s.bins, 0);
    for (size_t i = 0; i < s.row.size(); ++i) {
        sums[s.row[i]] += s.value[i];
        if (s.row[i] != s.col[i])
            sums[s.col[i]] += s.value[i];
    }
    std::vector<float> result(s.bins);
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = static_cast<float>(sums[i]);
    return result;
}
void fix_sum(const Sparse &s, std::vector<float> &norm) {
    double raw = 0, normalized = 0;
    for (size_t i = 0; i < s.row.size(); ++i) {
        uint32_t x = s.row[i], y = s.col[i];
        float nx = norm[x], ny = norm[y];
        if (!(nx > 0) || !(ny > 0) || !std::isfinite(nx) || !std::isfinite(ny))
            continue;
        double multiple = x == y ? 1 : 2;
        raw += multiple * s.value[i];
        normalized += multiple * s.value[i] / (double(nx) * ny);
    }
    if (raw <= 0 || normalized <= 0)
        return;
    float factor = static_cast<float>(std::sqrt(normalized / raw));
    for (float &v : norm)
        if (v > 0 && std::isfinite(v))
            v *= factor;
}
struct ExpectedAccumulator {
    std::vector<double> actual;
    std::map<uint32_t, double> observed;
    std::map<uint32_t, uint32_t> chromosome_bins;
    explicit ExpectedAccumulator(uint32_t size) : actual(size, 0) {}
    void add(uint32_t chr, const Sparse &s, const std::vector<float> *norm = nullptr) {
        double total = 0;
        chromosome_bins[chr] = s.bins;
        for (size_t i = 0; i < s.row.size(); ++i) {
            double value = s.value[i];
            if (norm) {
                float a = (*norm)[s.row[i]], b = (*norm)[s.col[i]];
                if (!(a > 0) || !(b > 0) || !std::isfinite(a) || !std::isfinite(b))
                    continue;
                value /= double(a) * b;
            }
            uint32_t distance = s.col[i] - s.row[i];
            actual[distance] += value;
            total += value;
        }
        if (total > 0)
            observed[chr] += total;
    }
    Vector finish(uint8_t kind, uint32_t norm, uint8_t unit, uint32_t ri, bool smooth) {
        std::vector<uint32_t> lengths;
        for (auto e : chromosome_bins)
            if (observed.count(e.first))
                lengths.push_back(e.second);
        std::sort(lengths.begin(), lengths.end());
        std::vector<uint64_t> suffix(lengths.size() + 1, 0);
        for (size_t i = lengths.size(); i-- > 0;)
            suffix[i] = suffix[i + 1] + lengths[i];
        auto possible = [&](uint64_t distance) {
            size_t first = std::upper_bound(lengths.begin(), lengths.end(), distance) -
                           lengths.begin();
            return double(suffix[first]) - double(lengths.size() - first) * distance;
        };
        Vector v;
        v.kind = kind;
        v.norm = norm;
        v.unit = unit;
        v.ri = ri;
        v.values.reserve(actual.size());
        uint64_t support = lengths.empty() ? 0 : lengths.back();
        if (smooth && support) {
            const double minimum = 400;
            double numerator = actual[0], denominator = possible(0);
            uint64_t lo = 0, hi = 0;
            for (uint64_t i = 0; i < support; ++i) {
                if (numerator < minimum) {
                    while (numerator < minimum && hi + 1 < support) {
                        ++hi;
                        numerator += actual[hi];
                        denominator += possible(hi);
                    }
                } else
                    while (hi > lo && numerator - actual[lo] - actual[hi] >= minimum) {
                        numerator -= actual[lo] + actual[hi];
                        denominator -= possible(lo) + possible(hi);
                        ++lo;
                        --hi;
                    }
                v.values.push_back(
                    bits(denominator > 0 ? static_cast<float>(numerator / denominator) : 0));
                if (hi + 2 < support) {
                    numerator += actual[hi + 1] + actual[hi + 2];
                    denominator += possible(hi + 1) + possible(hi + 2);
                    hi += 2;
                } else if (hi + 1 < support) {
                    ++hi;
                    numerator += actual[hi];
                    denominator += possible(hi);
                }
            }
        } else
            for (uint64_t i = 0; i < support; ++i) {
                double denominator = possible(i);
                v.values.push_back(
                    bits(denominator > 0 ? static_cast<float>(actual[i] / denominator) : 0));
            }
        v.values.resize(actual.size(), 0x7fc00000);
        for (auto e : observed) {
            double total = 0;
            uint32_t n = chromosome_bins[e.first];
            for (uint32_t d = 0; d < n && d < v.values.size(); ++d)
                total += double(n - d) * floating(v.values[d]);
            if (total > 0 && e.second > 0)
                v.scales[e.first] = bits(static_cast<float>(total / e.second));
        }
        return v;
    }
};
class VectorOutput {
    FILE *file_ = nullptr;
    std::string target_, temporary_;
    int level_;
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
    VectorOutput(Reader &reader, const Header &header, const std::string &path, int level)
        : target_(path), level_(level) {
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
        footer_bytes_ = repack_matrix_prefix(file_, reader, header);
    }
    ~VectorOutput() {
        if (file_)
            std::fclose(file_);
        if (!temporary_.empty())
            std::remove(temporary_.c_str());
    }
    void add(const Vector &v) {
        Bytes descriptors;
        constexpr uint32_t nominal = 65536;
        for (uint64_t begin = 0; begin < v.values.size(); begin += nominal) {
            uint32_t n = narrow(std::min<uint64_t>(nominal, v.values.size() - begin));
            Bytes best;
            uint8_t transform = 0;
            for (uint8_t t = 0; t < 3; ++t) {
                Bytes raw;
                raw.reserve(uint64_t(n) * 4);
                if (t == 1)
                    for (unsigned lane = 0; lane < 4; ++lane)
                        for (uint32_t j = 0; j < n; ++j)
                            put(raw, v.values[begin + j] >> (8 * lane), 1);
                else
                    for (uint32_t j = 0; j < n; ++j)
                        put(raw, v.values[begin + j] ^ (t == 2 && j ? v.values[begin + j - 1] : 0),
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
        put(entry, v.values.size(), 8);
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
uint32_t norm_id(const Header &h, const std::string &name) {
    auto it = std::find(h.norms.begin(), h.norms.end(), name);
    check(it != h.norms.end(), "normalization dictionary lacks requested type " + name);
    return static_cast<uint32_t>(it - h.norms.begin());
}
bool real_chromosome(const Chromosome &c) {
    return c.name != "ALL" && c.name != "All" && c.name != "all";
}
} // namespace
void add_norm_v10(const std::string &path, const AddNormOptions &options) {
    Reader reader(path);
    Header h = reader.header();
    std::set<std::string> standard{"VC", "VC_SQRT", "SCALE"};
    for (auto name : h.norms)
        if (!standard.count(name) &&
            (reader.vector_indexes()[0].length || reader.vector_indexes()[2].length))
            throw std::runtime_error(
                "V10: refusing to discard existing unsupported normalization " + name);
    bool build_scale = false;
    if (options.scale)
        for (uint8_t unit = 0; unit < 2; ++unit)
            for (const auto &resolution : h.resolutions[unit])
                if (unit != 0 || options.minimum_scale_resolution == 0 ||
                    int(resolution.bin) >= options.minimum_scale_resolution)
                    build_scale = true;
    h.norms.clear();
    if (options.vc)
        h.norms.push_back("VC");
    if (options.vc_sqrt)
        h.norms.push_back("VC_SQRT");
    if (build_scale)
        h.norms.push_back("SCALE");
    uint32_t vc = options.vc ? norm_id(h, "VC") : UINT32_MAX;
    uint32_t vcs = options.vc_sqrt ? norm_id(h, "VC_SQRT") : UINT32_MAX;
    uint32_t scale = build_scale ? norm_id(h, "SCALE") : UINT32_MAX;
    VectorOutput output(reader, h, path, options.compression_level);
    for (uint8_t unit = 0; unit < 2; ++unit)
        for (uint32_t ri = 0; ri < h.resolutions[unit].size(); ++ri) {
            uint32_t bin = h.resolutions[unit][ri].bin;
            bool derived = h.resolutions[unit][ri].mode != 0;
            std::fprintf(stderr, "\n%s resolution %u %s%s\n", unit ? "FRAG" : "BP", bin,
                         derived ? "(derived on the fly) " : "", "");
            uint32_t maximum = 0;
            for (uint32_t chr = 0; chr < h.chromosomes.size(); ++chr)
                maximum = std::max(maximum, narrow(h.bins(chr, unit, ri)));
            check(maximum != 0, "cannot create an empty expected vector");
            ExpectedAccumulator raw(maximum);
            for (uint32_t chr = 0; chr < h.chromosomes.size(); ++chr)
                if (real_chromosome(h.chromosomes[chr])) {
                    auto s = sparse_matrix(reader.matrix(chr, chr, unit, ri),
                                           narrow(h.bins(chr, unit, ri)));
                    raw.add(chr, s);
                }
            output.add(raw.finish(1, 0, unit, ri, true));
            struct Work {
                const char *name;
                uint32_t id;
                bool enabled;
            };
            std::vector<Work> work{
                {"VC", vc, options.vc},
                {"VC_SQRT", vcs, options.vc_sqrt},
                {"SCALE", scale,
                 build_scale && (unit != 0 || options.minimum_scale_resolution == 0 ||
                                 int(bin) >= options.minimum_scale_resolution)}};
            for (auto item : work) {
                if (!item.enabled)
                    continue;
                std::fprintf(stderr, "  %s", item.name);
                std::fflush(stderr);
                ExpectedAccumulator expected(maximum);
                uint32_t written = 0;
                for (uint32_t chr = 0; chr < h.chromosomes.size(); ++chr) {
                    if (!real_chromosome(h.chromosomes[chr]))
                        continue;
                    auto s = sparse_matrix(reader.matrix(chr, chr, unit, ri),
                                           narrow(h.bins(chr, unit, ri)));
                    if (s.row.empty())
                        continue;
                    std::vector<float> norm;
                    if (item.id == scale) {
                        ScaleParams p;
                        p.tolerance = options.tolerance;
                        p.total_max_iter = options.max_iterations;
                        p.num_threads = options.threads;
                        std::vector<float> values;
                        values.reserve(s.value.size());
                        for (double value : s.value)
                            values.push_back(static_cast<float>(value));
                        std::vector<double> b(s.bins, std::numeric_limits<double>::quiet_NaN());
                        scale_balance(s.row.size(), s.row, s.col, values, s.bins, b, p);
                        pp_norm_vector(s.row.size(), s.row, s.col, values, s.bins, b,
                                       options.threads);
                        norm.resize(s.bins);
                        for (size_t i = 0; i < b.size(); ++i)
                            norm[i] = static_cast<float>(b[i]);
                    } else {
                        norm = raw_vc(s);
                        if (item.id == vcs)
                            for (float &value : norm)
                                value = std::sqrt(value);
                    }
                    fix_sum(s, norm);
                    bool valid = false;
                    for (float value : norm)
                        if (value > 0 && std::isfinite(value)) {
                            valid = true;
                            break;
                        }
                    if (!valid)
                        continue;
                    Vector v;
                    v.kind = 0;
                    v.norm = item.id;
                    v.chr = chr;
                    v.unit = unit;
                    v.ri = ri;
                    for (float value : norm)
                        v.values.push_back(bits(value));
                    output.add(v);
                    expected.add(chr, s, &norm);
                    ++written;
                }
                if (written)
                    output.add(expected.finish(2, item.id, unit, ri, true));
                std::fprintf(stderr, " (%u chromosomes)\n", written);
            }
        }
    output.finish(h);
    std::fprintf(stderr, "\nV10 normalization complete: %s\n", path.c_str());
}
} // namespace hic10
