#include "format.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
#include <sys/stat.h>
#include <tuple>
#include <zlib.h>
#include <zstd.h>

namespace hic10 {
namespace {
constexpr uint64_t limit = 512ULL * 1024 * 1024;
class Input {
    FILE *f_;
    uint64_t end_;

  public:
    explicit Input(const std::string &path) : f_(std::fopen(path.c_str(), "rb")) {
        check(f_ != nullptr, "cannot open input");
        if (fseeko(f_, 0, SEEK_END) != 0) {
            std::fclose(f_);
            throw std::runtime_error("V10: input must be seekable");
        }
        auto n = ftello(f_);
        if (n < 0) {
            std::fclose(f_);
            throw std::runtime_error("V10: cannot stat input");
        }
        end_ = n;
        seek(0);
    }
    ~Input() {
        std::fclose(f_);
    }
    uint64_t position() const {
        auto n = ftello(f_);
        check(n >= 0, "input position failed");
        return n;
    }
    uint64_t left() const {
        return end_ - position();
    }
    void interval(uint64_t pos, uint64_t n) const {
        check(pos <= end_ && n <= end_ - pos, "V9 interval outside input");
    }
    void seek(uint64_t pos) {
        interval(pos, 0);
        check(fseeko(f_, pos, SEEK_SET) == 0, "input seek failed");
    }
    Bytes read(uint64_t n) {
        check(n <= left() && n <= limit, "invalid/oversized V9 read");
        Bytes b(n);
        check(std::fread(b.data(), 1, b.size(), f_) == b.size(), "truncated V9 input");
        return b;
    }
    uint64_t integer(unsigned n) {
        auto b = read(n);
        uint64_t v = 0;
        for (unsigned i = 0; i < n; ++i)
            v |= uint64_t(b[i]) << (8 * i);
        return v;
    }
    uint32_t word() {
        return static_cast<uint32_t>(integer(4));
    }
    uint64_t wide() {
        return integer(8);
    }
    uint32_t count(uint32_t minBytes = 1) {
        auto n = word();
        check(n <= INT32_MAX && n <= left() / minBytes, "invalid V9 count");
        return n;
    }
    std::string string() {
        std::string s;
        for (;;) {
            check(s.size() < 1024 * 1024, "V9 string too long");
            auto c = integer(1);
            if (!c)
                return s;
            s.push_back(static_cast<char>(c));
        }
    }
};
struct Cursor {
    const Bytes &b;
    size_t at = 0;
    uint64_t integer(unsigned n) {
        check(n <= b.size() - at, "truncated V9 block");
        uint64_t v = 0;
        for (unsigned i = 0; i < n; ++i)
            v |= uint64_t(b[at++]) << (8 * i);
        return v;
    }
    int32_t i32() {
        return static_cast<int32_t>(integer(4));
    }
    int16_t i16() {
        return static_cast<int16_t>(integer(2));
    }
    uint8_t byte() {
        return static_cast<uint8_t>(integer(1));
    }
};
struct Block {
    uint32_t number;
    uint64_t pos;
    uint32_t bytes;
};
struct Zoom {
    uint8_t unit;
    uint32_t bin;
    std::vector<Block> blocks;
};
struct LegacyMatrix {
    uint32_t a, b;
    std::vector<Zoom> zooms;
};
uint8_t unitId(const std::string &unit) {
    check(unit == "BP" || unit == "FRAG", "unknown V9 unit");
    return unit == "FRAG";
}
uint32_t normId(Header &h, const std::string &name) {
    check(!name.empty() && name != "NONE", "invalid normalization name");
    auto it = std::find(h.norms.begin(), h.norms.end(), name);
    if (it != h.norms.end())
        return static_cast<uint32_t>(it - h.norms.begin());
    h.norms.push_back(name);
    return narrow(h.norms.size() - 1);
}
Bytes inflateBlock(const Bytes &compressed) {
    Bytes result;
    if (compressed.size() >= 4 && compressed[0] == 0x28 && compressed[1] == 0xb5 &&
        compressed[2] == 0x2f && compressed[3] == 0xfd) {
        auto n = ZSTD_getFrameContentSize(compressed.data(), compressed.size());
        check(n <= limit, "invalid/oversized legacy Zstandard frame");
        result.resize(n);
        auto size =
            ZSTD_decompress(result.data(), result.size(), compressed.data(), compressed.size());
        check(!ZSTD_isError(size) && size == n, "legacy Zstandard decompression failure");
        return result;
    }
    z_stream stream{};
    check(inflateInit(&stream) == Z_OK, "zlib initialization failed");
    stream.next_in = const_cast<Bytef *>(compressed.data());
    stream.avail_in = narrow(compressed.size());
    int rc = Z_OK;
    do {
        size_t start = result.size();
        if (start >= limit) {
            inflateEnd(&stream);
            throw std::runtime_error("V10: V9 block exceeds allocation limit");
        }
        result.resize(start + 65536);
        stream.next_out = result.data() + start;
        stream.avail_out = 65536;
        rc = inflate(&stream, Z_NO_FLUSH);
        result.resize(start + 65536 - stream.avail_out);
        if (rc != Z_OK && rc != Z_STREAM_END) {
            inflateEnd(&stream);
            throw std::runtime_error("V10: corrupt V9 zlib block");
        }
    } while (rc != Z_STREAM_END);
    bool consumed = stream.avail_in == 0;
    inflateEnd(&stream);
    check(consumed, "trailing bytes in V9 zlib block");
    return result;
}
std::vector<Cell> decode(Input &f, const Block &block) {
    f.seek(block.pos);
    auto compressed = f.read(block.bytes);
    auto bytes = inflateBlock(compressed);
    Cursor c{bytes};
    int32_t n = c.i32(), x = c.i32(), y = c.i32();
    uint8_t floatFlag = c.byte(), xFlag = c.byte(), yFlag = c.byte(), type = c.byte();
    check(n >= 0 && uint64_t(n) <= bytes.size() / 2 && uint64_t(n) <= limit / sizeof(Cell) &&
              x >= 0 && y >= 0 && floatFlag <= 1 && xFlag <= 1 && yFlag <= 1 &&
              (type == 1 || type == 2),
          "invalid V9 block header");
    std::vector<Cell> cells;
    cells.reserve(n);
    auto count = [&]() {
        return floatFlag ? static_cast<uint32_t>(c.integer(4)) : bits(static_cast<float>(c.i16()));
    };
    auto emit = [&](int64_t bx, int64_t by, uint32_t value) {
        check(bx >= 0 && by >= 0 && bx <= UINT32_MAX && by <= UINT32_MAX, "invalid V9 coordinates");
        check(cells.size() < size_t(n), "V9 record count mismatch");
        cells.push_back({static_cast<uint32_t>(bx), static_cast<uint32_t>(by), value});
    };
    if (type == 1) {
        int32_t rows = yFlag ? c.i32() : c.i16();
        check(rows >= 0 && rows <= n, "invalid V9 row count");
        for (int32_t row = 0; row < rows; ++row) {
            int32_t dy = yFlag ? c.i32() : c.i16(), columns = xFlag ? c.i32() : c.i16();
            check(dy >= 0 && columns >= 0 && columns <= n, "invalid V9 row");
            for (int32_t j = 0; j < columns; ++j) {
                int32_t dx = xFlag ? c.i32() : c.i16();
                check(dx >= 0, "negative V9 column");
                uint32_t value = count();
                emit(int64_t(x) + dx, int64_t(y) + dy, value);
            }
        }
    } else {
        int32_t slots = c.i32();
        int16_t width = c.i16();
        check(slots >= 0 && width > 0 && uint64_t(slots) <= limit / 4, "invalid V9 dense size");
        for (int32_t i = 0; i < slots; ++i) {
            uint32_t value = count();
            float v = floating(value);
            if ((floatFlag && std::isnan(v)) || (!floatFlag && v == -32768))
                continue;
            emit(int64_t(x) + i % width, int64_t(y) + i / width, value);
        }
    }
    check(cells.size() == size_t(n) && c.at == bytes.size(),
          "V9 block length/record count mismatch");
    return cells;
}
void fitVector(Header &h, Vector &v) {
    uint64_t n = v.kind == 0 ? h.bins(v.chr, v.unit, v.ri) : 0;
    if (v.kind)
        for (uint32_t c = 0; c < h.chromosomes.size(); ++c)
            n = std::max(n, h.bins(c, v.unit, v.ri));
    check(n <= limit / 4, "vector exceeds allocation limit");
    if (n != v.values.size()) {
        // V9 commonly has a redundant terminal normalization bin or a shorter
        // expected array. Record original length and all surplus words, then
        // fill newly addressable distances with NaNs rather than invent values.
        std::string value = std::to_string(v.values.size()) + ":";
        const char *hex = "0123456789abcdef";
        for (uint64_t i = n; i < v.values.size(); ++i)
            for (unsigned j = 0; j < 8; ++j)
                value.push_back(hex[(v.values[i] >> (28 - 4 * j)) & 15]);
        h.attributes.emplace_back("hictools.v9.vector." + std::to_string(v.kind) + "." +
                                      std::to_string(v.norm) + "." + std::to_string(v.chr) + "." +
                                      std::to_string(v.unit) + "." + std::to_string(v.ri),
                                  value);
        v.values.resize(n, 0x7fc00000);
    }
}
} // namespace
void convert(const std::string &input, const std::string &output, const Options &options) {
    struct stat in{}, out{};
    check(stat(input.c_str(), &in) == 0, "cannot stat input");
    check(stat(output.c_str(), &out) != 0 || in.st_dev != out.st_dev || in.st_ino != out.st_ino,
          "input and output must be different files");
    Input f(input);
    auto magicBytes = f.read(4);
    check(magicBytes == Bytes({'H', 'I', 'C', 0}), "bad HIC magic");
    check(f.word() == 9, "convert accepts only V9 files");
    auto footer = f.wide();
    Header h;
    h.genome = f.string();
    auto nvi = f.wide(), nviLength = f.wide();
    uint32_t n = f.count(2);
    for (uint32_t i = 0; i < n; ++i) {
        auto key = f.string();
        auto value = f.string();
        h.attributes.emplace_back(key, value);
    }
    n = f.count(10);
    check(n > 0, "no V9 chromosomes");
    for (uint32_t i = 0; i < n; ++i) {
        Chromosome c;
        c.name = f.string();
        c.length = f.wide();
        check(c.length > 0 && c.length <= INT64_MAX, "invalid V9 chromosome length");
        h.chromosomes.push_back(c);
    }
    for (auto &list : h.resolutions) {
        n = f.count(4);
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t bin = f.word();
            check(bin > 0 && bin <= INT32_MAX, "invalid V9 bin size");
            list.push_back({bin});
        }
    }
    if (!h.resolutions[1].empty())
        for (auto &chr : h.chromosomes) {
            n = f.count(4);
            for (uint32_t j = 0; j < n; ++j)
                chr.sites.push_back(f.word());
        }
    f.seek(footer);
    auto v5Bytes = f.wide();
    auto v5End = plus(f.position(), v5Bytes);
    f.interval(f.position(), v5Bytes);
    struct Entry {
        std::string key;
        uint64_t pos;
        uint32_t size;
    };
    std::vector<Entry> entries;
    n = f.count(14);
    for (uint32_t i = 0; i < n; ++i) {
        Entry e;
        e.key = f.string();
        e.pos = f.wide();
        e.size = f.word();
        check(e.size >= 12, "invalid V9 matrix size");
        f.interval(e.pos, e.size);
        entries.push_back(e);
    }
    std::vector<Vector> vectors;
    // Initially record bin sizes in ri; remap after collecting legacy All zooms.
    auto expected = [&](bool normalized) {
        uint32_t count = f.count(normalized ? 21 : 19);
        for (uint32_t i = 0; i < count; ++i) {
            Vector v;
            v.kind = normalized ? 2 : 1;
            if (normalized)
                v.norm = normId(h, f.string());
            v.unit = unitId(f.string());
            v.ri = f.word();
            uint64_t nv = f.wide();
            check(nv <= f.left() / 4 && nv <= limit / 4, "invalid V9 expected length");
            for (uint64_t j = 0; j < nv; ++j)
                v.values.push_back(f.word());
            uint32_t ns = f.count(8);
            for (uint32_t j = 0; j < ns; ++j) {
                uint32_t chr = f.word(), value = f.word();
                check(chr < h.chromosomes.size() && v.scales.emplace(chr, value).second,
                      "duplicate/invalid V9 scale");
            }
            vectors.push_back(std::move(v));
        }
    };
    expected(false);
    check(f.position() == v5End, "V9 footer length mismatch");
    expected(true);
    if (!nvi) {
        nvi = f.position();
        nviLength = 0;
    }
    f.seek(nvi);
    n = f.count(26);
    struct Norm {
        Vector vector;
        uint64_t pos, len;
    };
    std::vector<Norm> norms;
    for (uint32_t i = 0; i < n; ++i) {
        Norm nv;
        nv.vector.kind = 0;
        nv.vector.norm = normId(h, f.string());
        nv.vector.chr = f.word();
        nv.vector.unit = unitId(f.string());
        nv.vector.ri = f.word();
        nv.pos = f.wide();
        nv.len = f.wide();
        check(nv.vector.chr < h.chromosomes.size() && nv.len >= 8, "invalid V9 norm vector");
        f.interval(nv.pos, nv.len);
        norms.push_back(std::move(nv));
    }
    if (nviLength)
        check(f.position() == plus(nvi, nviLength), "V9 normalization index length mismatch");
    for (auto &nv : norms) {
        f.seek(nv.pos);
        auto count = f.wide();
        check(count <= limit / 4 && plus(8, count * 4) == nv.len, "V9 norm vector length mismatch");
        for (uint64_t j = 0; j < count; ++j)
            nv.vector.values.push_back(f.word());
        vectors.push_back(std::move(nv.vector));
    }
    std::vector<LegacyMatrix> matrices;
    std::set<std::pair<uint32_t, uint32_t>> keys;
    for (auto entry : entries) {
        f.seek(entry.pos);
        LegacyMatrix m;
        m.a = f.word();
        m.b = f.word();
        n = f.count(39);
        check(m.a <= m.b && m.b < h.chromosomes.size() && keys.insert({m.a, m.b}).second,
              "invalid/duplicate V9 matrix pair");
        check(entry.key == std::to_string(m.a) + "_" + std::to_string(m.b),
              "V9 matrix key mismatch");
        std::set<std::pair<uint8_t, uint32_t>> zoomKeys;
        for (uint32_t i = 0; i < n; ++i) {
            Zoom z;
            z.unit = unitId(f.string());
            f.word();
            for (int j = 0; j < 4; ++j)
                f.word();
            z.bin = f.word();
            f.word();
            f.word();
            uint32_t nb = f.count(16);
            check(z.bin > 0 && z.bin <= INT32_MAX && zoomKeys.insert({z.unit, z.bin}).second,
                  "invalid V9 zoom");
            // The synthetic All matrix can have a resolution absent from the
            // V9 header. Advertise it explicitly rather than dropping the data.
            auto &list = h.resolutions[z.unit];
            if (std::none_of(list.begin(), list.end(), [&](auto r) { return r.bin == z.bin; }))
                list.push_back({z.bin});
            std::set<uint32_t> blockNumbers;
            for (uint32_t j = 0; j < nb; ++j) {
                Block b;
                b.number = f.word();
                b.pos = f.wide();
                b.bytes = f.word();
                check(b.bytes > 0 && blockNumbers.insert(b.number).second,
                      "invalid V9 block index");
                f.interval(b.pos, b.bytes);
                z.blocks.push_back(b);
            }
            m.zooms.push_back(std::move(z));
        }
        check(f.position() == plus(entry.pos, entry.size), "V9 matrix metadata length mismatch");
        matrices.push_back(std::move(m));
    }
    for (auto &list : h.resolutions) {
        std::sort(list.begin(), list.end(), [](auto a, auto b) { return a.bin < b.bin; });
        for (size_t i = 1; i < list.size(); ++i)
            check(list[i].bin != list[i - 1].bin, "duplicate V9 header resolution");
    }
    for (auto &v : vectors) {
        v.ri = h.resolution(v.unit, v.ri);
        fitVector(h, v);
    }
    Writer writer(output, h, options);
    for (const auto &m : matrices) {
        std::fprintf(stderr, "Converting %s x %s\n", h.chromosomes[m.a].name.c_str(),
                     h.chromosomes[m.b].name.c_str());
        writer.matrix(m.a, m.b, [&](uint8_t u, uint32_t ri) {
            Matrix result;
            result.scores = options.scores;
            uint32_t bin = h.resolutions[u][ri].bin;
            for (const auto &z : m.zooms)
                if (z.unit == u && z.bin == bin)
                    for (const auto &b : z.blocks) {
                        auto cells = decode(f, b);
                        for (auto c : cells) {
                            float v = floating(static_cast<uint32_t>(c.value));
                            if (!std::isfinite(v) || v <= 0 || std::floor(v) != v ||
                                double(v) >= std::ldexp(1.0, 64))
                                result.scores = true;
                            result.cells.push_back(c);
                        }
                    }
            if (!result.scores)
                for (auto &c : result.cells)
                    c.value = static_cast<uint64_t>(floating(static_cast<uint32_t>(c.value)));
            return result;
        });
    }
    writer.finish(vectors);
}
} // namespace hic10
