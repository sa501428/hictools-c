#include "reader.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <set>
#include <zstd.h>

namespace hic10 {
namespace {
constexpr uint64_t RECORD_LIMIT = 512ULL * 1024 * 1024;
struct Cursor {
    const uint8_t *data;
    size_t size, at = 0;
    explicit Cursor(const Bytes &b) : data(b.data()), size(b.size()) {}
    Cursor(const uint8_t *p, size_t n) : data(p), size(n) {}
    size_t left() const { return size - at; }
    void need(uint64_t n) const { check(n <= left(), "truncated V10 record"); }
    uint64_t integer(unsigned n) {
        need(n);
        uint64_t v = 0;
        for (unsigned i = 0; i < n; ++i)
            v |= uint64_t(data[at++]) << (8 * i);
        return v;
    }
    uint8_t byte() { return static_cast<uint8_t>(integer(1)); }
    uint32_t word() { return static_cast<uint32_t>(integer(4)); }
    uint64_t wide() { return integer(8); }
    uint64_t var() {
        uint64_t v = 0;
        for (unsigned i = 0; i < 10; ++i) {
            uint8_t b = byte();
            check(i < 9 || b <= 1, "ULEB128 overflow");
            v |= uint64_t(b & 127) << (7 * i);
            if (!(b & 128)) {
                check(i == 0 || b, "noncanonical ULEB128");
                return v;
            }
        }
        throw std::runtime_error("V10: unterminated ULEB128");
    }
    void magic4(const char *s) {
        need(4);
        check(std::memcmp(data + at, s, 4) == 0, "bad V10 record magic");
        at += 4;
    }
    void zero(size_t n) {
        while (n--)
            check(byte() == 0, "nonzero reserved field");
    }
    std::string string() {
        size_t start = at;
        while (at < size && data[at])
            ++at;
        check(at < size && at - start <= 1024 * 1024, "invalid V10 string");
        std::string s(reinterpret_cast<const char *>(data + start), at - start);
        ++at;
        return s;
    }
    Cursor take(uint64_t n) {
        need(n);
        Cursor c(data + at, static_cast<size_t>(n));
        at += n;
        return c;
    }
    void done() const { check(at == size, "trailing V10 record bytes"); }
};
FileLocator locator(Cursor &c) {
    FileLocator l{c.wide(), c.wide()};
    check((l.position == 0) == (l.length == 0), "incomplete V10 locator");
    return l;
}
Bytes decompress(Cursor &c, uint32_t output) {
    check(output <= RECORD_LIMIT, "V10 decompressed record too large");
    size_t frame = ZSTD_findFrameCompressedSize(c.data + c.at, c.left());
    check(!ZSTD_isError(frame) && frame == c.left(), "invalid V10 Zstandard frame length");
    Bytes result(output);
    size_t n = ZSTD_decompress(result.data(), result.size(), c.data + c.at, c.left());
    check(!ZSTD_isError(n) && n == output, "V10 Zstandard decompression failed");
    return result;
}
} // namespace
Reader::Reader(const std::string &path) {
    file_ = std::fopen(path.c_str(), "rb");
    check(file_ != nullptr, "cannot open V10 input: " + path);
    check(fseeko(file_, 0, SEEK_END) == 0, "cannot seek V10 input");
    auto end = ftello(file_);
    check(end >= 0, "cannot size V10 input");
    file_size_ = end;
    auto prefix = read_bytes(0, 88);
    Cursor p(prefix);
    p.magic4("HIC\0");
    check(p.word() == 10, "addnorm requires V10 input");
    uint64_t headerLength = p.wide();
    header_length_ = headerLength;
    footer_ = locator(p);
    vector_indexes_[0] = locator(p);
    vector_indexes_[1] = locator(p);
    vector_indexes_[2] = locator(p);
    p.zero(8);
    auto bytes = read_bytes(0, headerLength);
    Cursor c(bytes);
    c.magic4("HIC\0");
    check(c.word() == 10 && c.wide() == headerLength, "V10 header length mismatch");
    locator(c);
    locator(c);
    locator(c);
    locator(c);
    c.zero(8);
    header_.genome = c.string();
    uint32_t n = c.word();
    for (uint32_t i = 0; i < n; ++i)
        header_.attributes.emplace_back(c.string(), c.string());
    n = c.word();
    check(n > 0, "V10 file has no chromosomes");
    for (uint32_t i = 0; i < n; ++i)
        header_.chromosomes.push_back({c.string(), c.wide(), {}});
    for (auto &list : header_.resolutions) {
        n = c.word();
        for (uint32_t i = 0; i < n; ++i) {
            Resolution r;
            r.bin = c.word();
            r.mode = c.byte();
            r.aggregation = c.byte();
            c.zero(2);
            r.source = c.word();
            list.push_back(r);
        }
        for (uint32_t i = 0; i < list.size(); ++i) {
            const auto &r = list[i];
            check(r.bin && r.mode <= 1 && r.aggregation == 1, "invalid V10 resolution");
            if (!r.mode)
                check(r.source == UINT32_MAX, "materialized V10 resolution has a source");
            else
                check(r.source < list.size() && r.source != i && !list[r.source].mode &&
                          list[r.source].bin < r.bin && r.bin % list[r.source].bin == 0,
                      "invalid derived V10 resolution source");
        }
    }
    if (!header_.resolutions[1].empty())
        for (auto &chr : header_.chromosomes) {
            n = c.word();
            for (uint32_t i = 0; i < n; ++i)
                chr.sites.push_back(c.wide());
        }
    n = c.word();
    for (uint32_t i = 0; i < n; ++i) {
        auto name = c.string();
        check(!name.empty() && std::find(header_.norms.begin(), header_.norms.end(), name) ==
                                   header_.norms.end(),
              "invalid V10 normalization dictionary");
        header_.norms.push_back(std::move(name));
    }
    c.done();
    for (uint32_t chr = 0; chr < header_.chromosomes.size(); ++chr)
        for (uint8_t u = 0; u < 2; ++u)
            for (uint32_t ri = 0; ri < header_.resolutions[u].size(); ++ri)
                check(header_.bins(chr, u, ri) <= UINT32_MAX, "V10 chromosome bin count overflow");
    auto footer = read_bytes(footer_.position, footer_.length);
    Cursor f(footer);
    f.magic4("H10F");
    check(f.word() == 1 && f.wide() == footer_.length, "V10 footer mismatch");
    n = f.word();
    f.zero(4);
    MatrixKey previous{0, 0};
    for (uint32_t i = 0; i < n; ++i) {
        MatrixKey key{f.word(), f.word()};
        auto loc = locator(f);
        check(key.chr1 <= key.chr2 && key.chr2 < header_.chromosomes.size() &&
                  (!i || previous < key),
              "invalid V10 matrix key");
        matrix_locations_[key] = loc;
        matrix_keys_.push_back(key);
        previous = key;
    }
    f.done();
}
Reader::~Reader() {
    if (file_)
        std::fclose(file_);
}
Bytes Reader::read_bytes(uint64_t position, uint64_t length) {
    check(length <= RECORD_LIMIT && position <= file_size_ && length <= file_size_ - position,
          "V10 file interval out of bounds");
    Bytes b(length);
    check(fseeko(file_, position, SEEK_SET) == 0, "V10 seek failed");
    check(std::fread(b.data(), 1, b.size(), file_) == b.size(), "V10 short read");
    return b;
}
uint64_t Reader::vector_data_start() {
    uint64_t first = footer_.position;
    for (uint8_t kind = 0; kind < 3; ++kind) {
        const auto &loc = vector_indexes_[kind];
        if (!loc.length)
            continue;
        first = std::min(first, loc.position);
        auto bytes = read_bytes(loc.position, loc.length);
        Cursor c(bytes);
        c.magic4(kind == 0 ? "NVI0" : kind == 1 ? "EVI0" : "NEVI");
        check(c.word() == 1, "unsupported V10 vector index version");
        uint32_t entries = c.word();
        c.zero(4);
        for (uint32_t i = 0; i < entries; ++i) {
            uint32_t length = c.word();
            check(length >= 4, "invalid V10 vector entry length");
            Cursor e = c.take(length - 4);
            if (kind != 1)
                e.word();
            if (kind == 0)
                e.word();
            check(e.byte() <= 1, "invalid V10 vector unit");
            e.zero(3);
            e.word();
            e.word();
            e.wide();
            e.word();
            uint32_t chunks = e.word();
            if (kind) {
                uint32_t scales = e.word();
                e.zero(4);
                e.take(uint64_t(scales) * 8);
            }
            for (uint32_t j = 0; j < chunks; ++j) {
                e.wide();
                e.word();
                e.byte();
                e.byte();
                e.zero(2);
                uint64_t position = e.wide();
                uint32_t stored = e.word();
                e.word();
                check(position <= file_size_ && stored <= file_size_ - position,
                      "V10 vector chunk outside file");
                first = std::min(first, position);
            }
            e.done();
        }
        c.done();
    }
    check(first >= 88 && first <= footer_.position, "invalid V10 vector section order");
    return first;
}
std::vector<uint64_t> Reader::matrix_relocation_fields() {
    std::vector<uint64_t> fields;
    for (auto key : matrix_keys_) {
        const auto &meta = metadata(key);
        for (const auto &z : meta.zooms) {
            if (!z.page_index.length)
                continue;
            fields.push_back(z.page_index_position_field);
            auto bytes = read_bytes(z.page_index.position, z.page_index.length);
            Cursor c(bytes);
            c.magic4("H10I");
            check(c.word() == 1 && c.word() == z.pages, "V10 page index mismatch");
            uint32_t interval = c.word(), groups = c.word();
            c.zero(4);
            uint64_t blob_length = c.wide();
            check(interval && groups == (uint64_t(z.pages) + interval - 1) / interval,
                  "invalid V10 checkpoint count");
            for (uint32_t i = 0; i < groups; ++i) {
                c.word();
                c.word();
                c.word();
                c.zero(4);
                fields.push_back(z.page_index.position + c.at);
                c.wide();
                c.wide();
            }
            c.take(blob_length);
            c.done();
        }
    }
    std::sort(fields.begin(), fields.end());
    check(std::adjacent_find(fields.begin(), fields.end()) == fields.end(),
          "duplicate V10 relocation field");
    return fields;
}
const Reader::MatrixMeta &Reader::metadata(MatrixKey key) {
    auto cached = matrix_metadata_.find(key);
    if (cached != matrix_metadata_.end())
        return cached->second;
    auto found = matrix_locations_.find(key);
    check(found != matrix_locations_.end(), "missing V10 matrix");
    auto bytes = read_bytes(found->second.position, found->second.length);
    Cursor c(bytes);
    c.magic4("H10M");
    check(c.word() == 1 && c.word() == key.chr1 && c.word() == key.chr2,
          "V10 matrix header mismatch");
    uint32_t n = c.word();
    c.zero(4);
    check(n == header_.resolutions[0].size() + header_.resolutions[1].size(),
          "V10 resolution descriptor count mismatch");
    MatrixMeta meta;
    for (uint32_t i = 0; i < n; ++i) {
        Zoom z;
        z.unit = c.byte();
        z.mode = c.byte();
        z.aggregation = c.byte();
        z.type = c.byte();
        z.resolution = c.word();
        z.bin = c.word();
        z.source = c.word();
        z.grid = c.byte();
        c.zero(3);
        c.wide();
        z.occupied = c.wide();
        c.word();
        c.word();
        z.block_bins = c.word();
        z.columns = c.word();
        z.page_index_position_field = found->second.position + 24 + uint64_t(i) * 76 + 52;
        z.page_index = locator(c);
        z.pages = c.word();
        z.blocks = c.word();
        uint8_t expected_unit = i < header_.resolutions[0].size() ? 0 : 1;
        uint32_t expected_ri = i - (expected_unit ? header_.resolutions[0].size() : 0);
        const auto &r = header_.resolutions[expected_unit][expected_ri];
        check(z.unit == expected_unit && z.resolution == expected_ri && z.bin == r.bin &&
                  z.mode == r.mode && z.aggregation == r.aggregation && z.source == r.source &&
                  z.type <= 1 && z.grid <= 1 && (!z.grid || key.chr1 == key.chr2) &&
                  z.block_bins && z.columns,
              "V10 zoom mismatch");
        if (z.mode)
            check(!z.page_index.length && !z.pages && !z.blocks,
                  "derived V10 resolution has storage");
        else if (z.occupied)
            check(z.page_index.length && z.pages && z.blocks,
                  "materialized V10 resolution lacks storage");
        else
            check(!z.pages && !z.blocks, "empty V10 resolution has pages");
        meta.zooms.push_back(z);
    }
    c.done();
    return matrix_metadata_.emplace(key, std::move(meta)).first->second;
}
Matrix Reader::materialized(MatrixKey key, const Zoom &z) {
    Matrix result;
    result.scores = z.type != 0;
    if (!z.pages)
        return result;
    auto index = read_bytes(z.page_index.position, z.page_index.length);
    Cursor c(index);
    c.magic4("H10I");
    check(c.word() == 1 && c.word() == z.pages, "V10 page index mismatch");
    uint32_t interval = c.word(), groups = c.word();
    c.zero(4);
    uint64_t blobLength = c.wide();
    check(interval && groups == (uint64_t(z.pages) + interval - 1) / interval,
          "invalid V10 checkpoint count");
    struct Checkpoint {
        uint32_t first, n, block;
        uint64_t position, offset;
    };
    std::vector<Checkpoint> checkpoints;
    for (uint32_t i = 0; i < groups; ++i) {
        Checkpoint q{c.word(), c.word(), c.word(), 0, 0};
        c.zero(4);
        q.position = c.wide();
        q.offset = c.wide();
        checkpoints.push_back(q);
    }
    Cursor blob = c.take(blobLength);
    c.done();
    uint32_t pages = 0;
    for (const auto &q : checkpoints) {
        check(q.first == pages && q.offset == blob.at && q.n && q.n <= interval,
              "invalid V10 checkpoint");
        uint64_t position = q.position;
        uint32_t first = q.block, last = 0;
        for (uint32_t i = 0; i < q.n; ++i) {
            if (i)
                first = narrow(uint64_t(last) + 1 + blob.var());
            last = narrow(uint64_t(first) + blob.var());
            uint64_t stored = blob.var();
            uint32_t raw = narrow(blob.var());
            auto page = read_bytes(position, stored);
            position += stored;
            Cursor pc(page);
            pc.magic4("H10P");
            check(pc.byte() == 1 && pc.byte() == 1, "unsupported V10 page codec");
            pc.zero(2);
            check(pc.word() == raw, "V10 page size mismatch");
            uint32_t blockCount = pc.word();
            check(blockCount && blockCount <= raw / 42, "invalid V10 page block count");
            auto payload = decompress(pc, raw);
            Cursor body(payload);
            Cursor directory = body.take(body.word());
            struct BlockRef {
                uint32_t number;
                uint64_t length;
            };
            std::vector<BlockRef> refs;
            uint32_t number = 0;
            for (uint32_t bi = 0; bi < blockCount; ++bi) {
                uint64_t delta = directory.var();
                check(!bi || delta, "duplicate V10 block number");
                number = narrow(bi ? uint64_t(number) + delta : delta);
                uint64_t length = directory.var();
                check(length >= 40, "invalid V10 block length");
                refs.push_back({number, length});
            }
            directory.done();
            check(refs.front().number == first && refs.back().number == last,
                  "V10 page block range mismatch");
            for (const auto &ref : refs) {
                Cursor b = body.take(ref.length);
                check(b.byte() == 1, "unsupported V10 block version");
                uint8_t rep = b.byte(), mode = b.byte(), type = b.byte(), flags = b.byte();
                b.zero(3);
                uint32_t x = b.word(), y = b.word(), w = b.word(), h = b.word();
                uint64_t occupied = b.wide();
                uint32_t np = b.word(), nv = b.word();
                check(type == z.type && rep <= 2 && mode <= 2 && flags <= 1 && w && h && occupied,
                      "invalid V10 block");
                Cursor positions = b.take(np), values = b.take(nv);
                b.done();
                uint64_t cells = uint64_t(w) * h, slots = rep == 2 ? cells : occupied;
                check(occupied <= cells && slots <= RECORD_LIMIT / sizeof(uint64_t),
                      "oversized V10 block");
                std::vector<uint64_t> present;
                if (rep == 0) {
                    check(!flags && occupied <= np, "invalid V10 sparse position stream");
                    uint64_t previous = 0;
                    for (uint64_t k = 0; k < occupied; ++k) {
                        uint64_t d = positions.var();
                        check(!k || d, "duplicate V10 sparse cell");
                        previous = k ? plus(previous, d) : d;
                        check(previous < cells, "V10 sparse cell outside block");
                        present.push_back(previous);
                    }
                } else if (rep == 1 || type) {
                    check(flags == 1 && np == (cells + 7) / 8,
                          "invalid V10 presence bitmap");
                    if (cells % 8)
                        check((positions.data[np - 1] >> (cells % 8)) == 0,
                              "nonzero V10 bitmap padding");
                    for (uint64_t k = 0; k < cells; ++k)
                        if (positions.data[k / 8] & (1u << (k % 8)))
                            present.push_back(k);
                    check(present.size() == occupied, "V10 bitmap population mismatch");
                    positions.at = positions.size;
                } else
                    check(!flags && !np, "dense V10 counts have a presence stream");
                positions.done();
                check(rep != 2 || mode == 2, "dense V10 values must use direct mode");
                auto scalar = [&]() { return type ? uint64_t(values.word()) : values.var(); };
                std::vector<uint64_t> decoded;
                if (mode == 0)
                    decoded.assign(slots, scalar());
                else if (mode == 1) {
                    uint64_t def = scalar(), ne = values.var();
                    check(ne && ne < slots && ne <= values.left(),
                          "invalid V10 exception count");
                    std::vector<uint64_t> ord;
                    uint64_t previous = 0;
                    for (uint64_t k = 0; k < ne; ++k) {
                        uint64_t d = values.var();
                        check(!k || d, "duplicate V10 exception ordinal");
                        previous = k ? plus(previous, d) : d;
                        check(previous < slots, "V10 exception ordinal out of range");
                        ord.push_back(previous);
                    }
                    decoded.assign(slots, def);
                    for (auto o : ord) {
                        uint64_t exception = scalar();
                        check(exception != def, "V10 exception equals default");
                        decoded[o] = exception;
                    }
                } else {
                    decoded.reserve(slots);
                    for (uint64_t k = 0; k < slots; ++k)
                        decoded.push_back(scalar());
                }
                values.done();
                uint64_t emitted = 0;
                size_t pi = 0;
                for (uint64_t k = 0; k < slots; ++k) {
                    bool exists = true;
                    uint64_t local = rep == 2 ? k : present[k], value = decoded[k];
                    if (rep == 2) {
                        if (!type)
                            exists = value != 0;
                        else {
                            exists = pi < present.size() && present[pi] == k;
                            if (exists)
                                ++pi;
                        }
                    }
                    if (!exists)
                        continue;
                    if (!type)
                        check(value > 0, "V10 count cell is zero");
                    uint32_t bx = narrow(uint64_t(x) + local % w),
                             by = narrow(uint64_t(y) + local / w);
                    check(bx < header_.bins(key.chr1, z.unit, z.resolution) &&
                              by < header_.bins(key.chr2, z.unit, z.resolution) &&
                              (key.chr1 != key.chr2 || bx <= by),
                          "V10 block cell outside matrix");
                    result.cells.push_back({bx, by, value});
                    ++emitted;
                }
                check(emitted == occupied, "V10 block occupancy mismatch");
            }
            body.done();
            ++pages;
        }
    }
    blob.done();
    check(pages == z.pages && result.cells.size() == z.occupied, "V10 matrix occupancy mismatch");
    std::sort(result.cells.begin(), result.cells.end(),
              [](auto a, auto b) { return std::tie(a.y, a.x) < std::tie(b.y, b.x); });
    return result;
}
Matrix Reader::matrix(uint32_t chr1, uint32_t chr2, uint8_t unit, uint32_t ri) {
    MatrixKey key{chr1, chr2};
    auto loc = matrix_locations_.find(key);
    if (loc == matrix_locations_.end())
        return {};
    const auto &meta = metadata(key);
    const Zoom *zoom = nullptr;
    for (const auto &z : meta.zooms)
        if (z.unit == unit && z.resolution == ri) {
            zoom = &z;
            break;
        }
    check(zoom != nullptr, "V10 resolution missing from matrix");
    if (!zoom->mode)
        return materialized(key, *zoom);
    Matrix source = matrix(chr1, chr2, unit, zoom->source);
    uint32_t factor = zoom->bin / header_.resolutions[unit][zoom->source].bin;
    struct Acc {
        uint64_t count = 0;
        double score = 0;
    };
    std::map<std::pair<uint32_t, uint32_t>, Acc> sums;
    for (auto cell : source.cells) {
        auto &a = sums[{cell.y / factor, cell.x / factor}];
        if (source.scores) {
            float v = floating(static_cast<uint32_t>(cell.value));
            check(std::isfinite(v), "nonfinite derived score");
            a.score += v;
        } else
            a.count = plus(a.count, cell.value);
    }
    Matrix out;
    out.scores = source.scores;
    for (auto e : sums)
        out.cells.push_back(
            {e.first.second, e.first.first,
             out.scores ? bits(static_cast<float>(e.second.score)) : e.second.count});
    return out;
}
} // namespace hic10
