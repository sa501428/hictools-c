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
    int level = 3;
    uint32_t threads = 4;
    uint32_t blockBins = 256, pageBytes = 128 * 1024;
    std::vector<std::pair<uint32_t, uint32_t>> derived;
    bool scores = false;
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
