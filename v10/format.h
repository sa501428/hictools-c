#pragma once
// Internal V10 writer model. No dependency on straw or the V9 writer.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

class ThreadPool;

namespace hic10 {
using Bytes = std::vector<uint8_t>;
inline void check(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error("V10: " + message);
}
inline uint64_t plus(uint64_t a, uint64_t b) {
    check(b <= UINT64_MAX - a, "uint64 overflow");
    return a + b;
}
inline uint32_t narrow(uint64_t x) {
    check(x <= UINT32_MAX, "uint32 overflow");
    return static_cast<uint32_t>(x);
}
inline uint64_t ceil_div(uint64_t value, uint64_t divisor) {
    check(divisor != 0, "division by zero");
    return value / divisor + (value % divisor != 0);
}
// V10 block numbers are u32. A fixed block scale can overflow at very fine
// resolutions: hg38 chr1 has about 24.9 million bins at 10 bp, which produces
// more than 2^32 positions with 256-bin blocks. Keep the requested scale as a
// minimum, but enlarge it until each grid axis has at most 65,536 blocks.
inline uint32_t safe_block_bin_count(uint64_t column_bins, uint64_t row_bins, uint32_t requested) {
    check(requested != 0 && column_bins != 0 && row_bins != 0, "invalid block geometry");
    constexpr uint64_t MAX_BLOCKS_PER_AXIS = uint64_t{1} << 16;
    uint64_t block_bins =
        std::max<uint64_t>(requested, std::max(ceil_div(column_bins, MAX_BLOCKS_PER_AXIS),
                                               ceil_div(row_bins, MAX_BLOCKS_PER_AXIS)));
    uint64_t columns = ceil_div(column_bins, block_bins);
    uint64_t rows = ceil_div(row_bins, block_bins);
    check(columns * rows <= uint64_t{UINT32_MAX} + 1,
          "matrix grid cannot be represented by u32 block numbers");
    return narrow(block_bins);
}
// V10 keeps the V9 high-resolution block-sizing policy: as bin sizes become
// finer, grow logical blocks so a chromosome axis has hundreds (cis) or tens
// (trans) of blocks rather than millions. --block-bins remains an additional
// minimum; it cannot opt back into the pathological tiny-block layout.
inline uint32_t adaptive_block_bin_count(uint64_t column_bins, uint64_t row_bins,
                                         uint32_t bin_size, uint32_t requested,
                                         bool cis) {
    check(requested && column_bins && row_bins && bin_size, "invalid adaptive block geometry");
    uint64_t block_bins = requested;
    constexpr uint64_t BLOCK_CAPACITY = 1000;
    const uint64_t cutoff = cis ? 500 : 5000;
    const uint64_t maximum_bins = std::max(column_bins, row_bins);
    uint64_t target_columns;
    if (bin_size < cutoff) {
        using Wide = unsigned __int128;
        target_columns = static_cast<uint64_t>(
            Wide(maximum_bins) * bin_size / (BLOCK_CAPACITY * cutoff)) + 1;
    } else {
        target_columns = maximum_bins / BLOCK_CAPACITY + 1;
    }
    block_bins = std::max<uint64_t>(block_bins, maximum_bins / target_columns + 1);
    // This also makes rectangular trans block numbers fit u32. Rotated cis
    // numbering has far fewer rows; validate its exact maximum at the call site.
    return safe_block_bin_count(column_bins, row_bins, narrow(block_bins));
}
inline uint32_t rotated_depth(uint64_t distance, uint32_t block_bins) {
    check(block_bins, "zero rotated block scale");
    using Wide = unsigned __int128;
    const Wide lhs = Wide(distance) * distance;
    const Wide scale = Wide(2) * block_bins * block_bins;
    uint32_t depth = 0;
    while (depth < 32) {
        const Wide threshold = (Wide(1) << (depth + 1)) - 1;
        if (threshold * threshold > lhs / scale)
            break;
        ++depth;
    }
    return depth;
}
inline uint32_t rotated_block_number(uint32_t x, uint32_t y, uint32_t block_bins,
                                     uint32_t columns) {
    check(x <= y && columns, "invalid rotated cis coordinate");
    uint64_t depth = rotated_depth(uint64_t(y) - x, block_bins);
    uint64_t along = (uint64_t(x) + y) / (uint64_t(2) * block_bins);
    return narrow(depth * columns + along);
}
inline uint32_t bits(float x) {
    uint32_t b;
    std::memcpy(&b, &x, 4);
    return b;
}
inline float floating(uint32_t b) {
    float x;
    std::memcpy(&x, &b, 4);
    return x;
}
inline uint64_t bits(double x) {
    uint64_t b;
    std::memcpy(&b, &x, 8);
    return b;
}
inline void put(Bytes &b, uint64_t v, unsigned n) {
    for (unsigned i = 0; i < n; ++i)
        b.push_back(static_cast<uint8_t>(v >> (i * 8)));
}
inline void var(Bytes &b, uint64_t v) {
    do {
        uint8_t x = v & 127;
        v >>= 7;
        b.push_back(x | (v ? 128 : 0));
    } while (v);
}
inline void magic(Bytes &b, const char *s) {
    b.insert(b.end(), s, s + 4);
}
inline void str(Bytes &b, const std::string &s) {
    check(s.find('\0') == std::string::npos, "embedded NUL in string");
    b.insert(b.end(), s.begin(), s.end());
    b.push_back(0);
}
inline void append(Bytes &a, const Bytes &b) {
    a.insert(a.end(), b.begin(), b.end());
}
struct Chromosome {
    std::string name;
    uint64_t length;
    std::vector<uint64_t> sites;
};
struct Resolution {
    uint32_t bin;
    uint8_t mode = 0, aggregation = 1;
    uint32_t source = UINT32_MAX;
};
inline uint32_t required_derived_source(uint32_t bin) {
    switch (bin) {
    case 20:
    case 50:
        return 10;
    case 200:
    case 500:
        return 100;
    case 2000:
        return 1000;
    default:
        return 0;
    }
}
inline bool required_materialized_resolution(uint32_t bin) {
    return bin == 500000;
}
inline bool required_bp_resolution_policy(const std::vector<Resolution> &list) {
    for (uint32_t i = 0; i < list.size(); ++i) {
        const auto &r = list[i];
        if (uint32_t source_bin = required_derived_source(r.bin)) {
            auto source = std::find_if(list.begin(), list.end(),
                                       [&](const Resolution &s) { return s.bin == source_bin; });
            if (source == list.end() || !r.mode || r.source != uint32_t(source - list.begin()) ||
                source->mode)
                return false;
        }
        if (required_materialized_resolution(r.bin) && r.mode)
            return false;
    }
    return true;
}
struct Header {
    std::string genome;
    std::vector<std::pair<std::string, std::string>> attributes;
    std::vector<Chromosome> chromosomes;
    std::array<std::vector<Resolution>, 2> resolutions;
    std::vector<std::string> norms;
    uint64_t bins(uint32_t chr, uint8_t unit, uint32_t ri) const {
        uint64_t length = unit ? chromosomes.at(chr).sites.size() + 1 : chromosomes.at(chr).length;
        uint32_t bin = resolutions.at(unit).at(ri).bin;
        return length / bin + (length % bin != 0);
    }
    uint32_t resolution(uint8_t unit, uint32_t bin) const {
        const auto &list = resolutions.at(unit);
        for (uint32_t i = 0; i < list.size(); ++i)
            if (list[i].bin == bin)
                return i;
        throw std::runtime_error("V10: unknown resolution " + std::to_string(bin));
    }
};
struct Cell {
    uint32_t x, y;
    uint64_t value;
}; // raw u64 count or f32 bits
struct Matrix {
    bool scores = false;
    std::vector<Cell> cells;
};
struct Vector {
    uint8_t kind = 0, unit = 0; // 0 normalization, 1 expected, 2 normalized expected
    uint32_t norm = 0, chr = 0, ri = 0;
    std::vector<uint32_t> values; // preserve every f32 bit, including NaN payloads
    std::map<uint32_t, uint32_t> scales;
};
struct Options {
    int level = 6;
    uint32_t threads = 4;
    uint32_t blockBins = 256, pageBytes = 512 * 1024;
    std::vector<std::pair<uint32_t, uint32_t>> derived;
    bool scores = false, verifyDerived = true;
};
class Writer {
  public:
    Writer(const std::string &output, Header header, const Options &options,
           std::shared_ptr<ThreadPool> pool = {});
    ~Writer();
    Writer(const Writer &) = delete;
    Writer &operator=(const Writer &) = delete;
    // Callback loads one resolution at a time; matrix storage is bounded by one
    // chromosome-pair resolution plus the currently assembled page.
    void matrix(uint32_t chr1, uint32_t chr2, const std::function<Matrix(uint8_t, uint32_t)> &load);
    void finish(const std::vector<Vector> &vectors);
    const Header &header() const {
        return header_;
    }

  private:
    FILE *file_ = nullptr;
    std::string output_, temporary_;
    Header header_;
    Options options_;
    std::shared_ptr<ThreadPool> pool_;
    struct Entry {
        uint32_t a, b;
        uint64_t pos, len;
    };
    std::vector<Entry> entries_;
    uint64_t position() const;
    uint64_t write(const Bytes &bytes);
    void patch(uint64_t pos, const Bytes &bytes);
};
void convert(const std::string &input, const std::string &output, const Options &options);
} // namespace hic10
