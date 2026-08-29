#include "format.h"
#include "../common/thread_pool.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <future>
#include <set>
#include <sys/stat.h>
#include <tuple>
#include <unistd.h>
#include <zstd.h>

namespace hic10 {
namespace {
Bytes compressed(const Bytes &bytes, int level) {
    Bytes out(ZSTD_compressBound(bytes.size()));
    size_t n = ZSTD_compress(out.data(), out.size(), bytes.data(), bytes.size(), level);
    check(!ZSTD_isError(n), "Zstandard compression failed");
    out.resize(n);
    return out;
}
void scalar(Bytes &b, uint64_t value, bool scores) {
    if (scores)
        put(b, value, 4);
    else
        var(b, value);
}
Bytes block(const std::vector<Cell> &cells, size_t begin, size_t end, bool scores) {
    check(begin < end && end <= cells.size(), "empty logical block");
    const size_t count = end - begin;
    auto cell = [&](size_t i) -> const Cell & { return cells[begin + i]; };
    uint32_t x = UINT32_MAX, y = UINT32_MAX, xmax = 0, ymax = 0;
    for (size_t i = 0; i < count; ++i) {
        auto c = cell(i);
        x = std::min(x, c.x);
        y = std::min(y, c.y);
        xmax = std::max(xmax, c.x);
        ymax = std::max(ymax, c.y);
    }
    uint32_t w = narrow(uint64_t(xmax) - x + 1), h = narrow(uint64_t(ymax) - y + 1);
    uint64_t area = uint64_t(w) * h;
    Bytes sparse, bitmap((area + 7) / 8, 0), direct;
    uint64_t prev = 0;
    std::map<uint64_t, size_t> frequencies;
    for (size_t i = 0; i < count; ++i) {
        auto c = cell(i);
        uint64_t p = uint64_t(c.y - y) * w + c.x - x;
        check(!i || p > prev, "duplicate/out-of-order cell");
        var(sparse, i ? p - prev : p);
        prev = p;
        bitmap[p / 8] |= 1u << (p % 8);
        scalar(direct, c.value, scores);
        ++frequencies[c.value];
    }
    uint64_t defaultValue = 0;
    size_t best = 0;
    for (auto e : frequencies)
        if (e.second > best) {
            defaultValue = e.first;
            best = e.second;
        }
    uint8_t mode = 2;
    Bytes values = direct;
    if (best == count) {
        values.clear();
        scalar(values, defaultValue, scores);
        mode = 0;
    } else {
        Bytes exceptions;
        scalar(exceptions, defaultValue, scores);
        var(exceptions, count - best);
        uint64_t last = 0;
        bool first = true;
        for (size_t i = 0; i < count; ++i)
            if (cell(i).value != defaultValue) {
                var(exceptions, first ? i : i - last);
                last = i;
                first = false;
            }
        for (size_t i = 0; i < count; ++i)
            if (cell(i).value != defaultValue)
                scalar(exceptions, cell(i).value, scores);
        if (exceptions.size() < values.size()) {
            values = std::move(exceptions);
            mode = 1;
        }
    }
    uint8_t rep = 0, flags = 0;
    Bytes positions = sparse;
    if (bitmap.size() < positions.size()) {
        positions = bitmap;
        rep = 1;
        flags = 1;
    }
    Bytes dense;
    size_t ci = 0;
    for (uint64_t p = 0; p < area; ++p) {
        bool present = ci < count && uint64_t(cell(ci).y - y) * w + cell(ci).x - x == p;
        scalar(dense, present ? cell(ci++).value : 0, scores);
    }
    if (dense.size() + (scores ? bitmap.size() : 0) < positions.size() + values.size()) {
        rep = 2;
        flags = scores ? 1 : 0;
        mode = 2;
        values = std::move(dense);
        positions = scores ? bitmap : Bytes{};
    }
    Bytes out;
    put(out, 1, 1);
    put(out, rep, 1);
    put(out, mode, 1);
    put(out, scores, 1);
    put(out, flags, 1);
    put(out, 0, 3);
    put(out, x, 4);
    put(out, y, 4);
    put(out, w, 4);
    put(out, h, 4);
    put(out, count, 8);
    put(out, positions.size(), 4);
    put(out, values.size(), 4);
    append(out, positions);
    append(out, values);
    return out;
}

struct EncodedPage {
    uint32_t first, last, raw;
    Bytes stored;
};

EncodedPage encode_page(std::vector<std::pair<uint32_t, Bytes>> pending, int level) {
    check(!pending.empty(), "empty page");
    Bytes directory, payload;
    uint32_t previous = 0;
    for (size_t i = 0; i < pending.size(); ++i) {
        var(directory, i ? pending[i].first - previous : pending[i].first);
        var(directory, pending[i].second.size());
        previous = pending[i].first;
    }
    put(payload, directory.size(), 4);
    append(payload, directory);
    for (const auto &p : pending)
        append(payload, p.second);
    EncodedPage result{pending.front().first, pending.back().first, narrow(payload.size()), {}};
    magic(result.stored, "H10P");
    put(result.stored, 1, 1);
    put(result.stored, 1, 1);
    put(result.stored, 0, 2);
    put(result.stored, result.raw, 4);
    put(result.stored, pending.size(), 4);
    append(result.stored, compressed(payload, level));
    return result;
}
void canonical(Matrix &m, uint32_t a, uint32_t b, const Header &h, uint8_t u, uint32_t ri) {
    for (auto c : m.cells)
        check(c.x < h.bins(a, u, ri) && c.y < h.bins(b, u, ri) && (a != b || c.x <= c.y) &&
                  (m.scores ? c.value <= UINT32_MAX : c.value > 0),
              "invalid input cell");
    std::sort(m.cells.begin(), m.cells.end(),
              [](const Cell &a, const Cell &b) { return std::tie(a.y, a.x) < std::tie(b.y, b.x); });
    for (size_t i = 1; i < m.cells.size(); ++i)
        check(m.cells[i].x != m.cells[i - 1].x || m.cells[i].y != m.cells[i - 1].y,
              "duplicate input cell");
}
Matrix aggregate(const Matrix &source, uint32_t factor) {
    struct Sum {
        uint64_t count = 0;
        double score = 0;
    };
    std::map<std::pair<uint32_t, uint32_t>, Sum> sums;
    for (auto c : source.cells) {
        auto &s = sums[{c.y / factor, c.x / factor}];
        if (source.scores) {
            float v = floating(static_cast<uint32_t>(c.value));
            check(std::isfinite(v), "nonfinite score cannot be derived");
            s.score += double(v);
            check(std::isfinite(s.score), "derived score overflow");
        } else
            s.count = plus(s.count, c.value);
    }
    Matrix out;
    out.scores = source.scores;
    for (auto e : sums)
        out.cells.push_back(
            {e.first.second, e.first.first,
             source.scores ? bits(static_cast<float>(e.second.score)) : e.second.count});
    return out;
}
} // namespace
uint64_t Writer::position() const {
    auto p = ftello(file_);
    check(p >= 0, "cannot determine output position");
    return static_cast<uint64_t>(p);
}
uint64_t Writer::write(const Bytes &bytes) {
    auto p = position();
    check(std::fwrite(bytes.data(), 1, bytes.size(), file_) == bytes.size(), "write failed");
    return p;
}
void Writer::patch(uint64_t pos, const Bytes &bytes) {
    auto saved = position();
    check(fseeko(file_, pos, SEEK_SET) == 0, "seek failed");
    write(bytes);
    check(fseeko(file_, saved, SEEK_SET) == 0, "seek failed");
}
Writer::Writer(const std::string &output, Header header, const Options &options,
               std::shared_ptr<ThreadPool> pool)
    : output_(output), header_(std::move(header)), options_(options), pool_(std::move(pool)) {
    check(options_.blockBins > 0 && options_.blockBins <= 4096 && options_.pageBytes >= 1024 &&
              options_.pageBytes <= 16 * 1024 * 1024,
          "invalid block/page size");
    check(options_.level >= ZSTD_minCLevel() && options_.level <= ZSTD_maxCLevel(),
          "invalid Zstandard level");
    check(options_.threads > 0 && options_.threads <= 256, "invalid writer thread count");
    if (!pool_)
        pool_ = std::make_shared<ThreadPool>(options_.threads);
    check(!header_.chromosomes.empty(), "no chromosomes");
    for (auto &list : header_.resolutions) {
        std::sort(list.begin(), list.end(), [](auto a, auto b) { return a.bin < b.bin; });
        for (size_t i = 0; i < list.size(); ++i)
            check(list[i].bin > 0 && (!i || list[i].bin > list[i - 1].bin),
                  "invalid resolution list");
    }
    for (auto d : options_.derived) {
        auto target = header_.resolution(0, d.first), source = header_.resolution(0, d.second);
        check(d.first > d.second && d.first % d.second == 0,
              "derived target must be an integer multiple of finer source");
        auto &r = header_.resolutions[0][target];
        check(!r.mode, "duplicate derived target");
        r.mode = 1;
        r.source = source;
        r.aggregation = 1;
    }
    for (const auto &r : header_.resolutions[0])
        if (r.mode)
            check(!header_.resolutions[0][r.source].mode, "chained derivation is forbidden");
    std::set<std::string> names;
    for (uint32_t i = 0; i < header_.chromosomes.size(); ++i) {
        const auto &c = header_.chromosomes[i];
        check(c.length > 0 && c.length <= INT64_MAX && !c.name.empty() &&
                  names.insert(c.name).second,
              "invalid chromosome");
        uint64_t prev = 0;
        for (auto p : c.sites) {
            check(p > prev && p < c.length, "invalid fragment sites");
            prev = p;
        }
        for (uint8_t u = 0; u < 2; ++u)
            for (uint32_t r = 0; r < header_.resolutions[u].size(); ++r)
                check(header_.bins(i, u, r) <= UINT32_MAX, "too many bins");
    }
    names.clear();
    for (auto name : header_.norms)
        check(!name.empty() && name != "NONE" && names.insert(name).second,
              "invalid normalization dictionary");
    Bytes b;
    magic(b, "HIC\0");
    put(b, 10, 4);
    b.resize(88, 0);
    str(b, header_.genome);
    put(b, header_.attributes.size(), 4);
    for (const auto &a : header_.attributes) {
        str(b, a.first);
        str(b, a.second);
    }
    put(b, header_.chromosomes.size(), 4);
    for (const auto &c : header_.chromosomes) {
        str(b, c.name);
        put(b, c.length, 8);
    }
    for (const auto &list : header_.resolutions) {
        put(b, list.size(), 4);
        for (auto r : list) {
            put(b, r.bin, 4);
            put(b, r.mode, 1);
            put(b, r.aggregation, 1);
            put(b, 0, 2);
            put(b, r.source, 4);
        }
    }
    if (!header_.resolutions[1].empty())
        for (const auto &c : header_.chromosomes) {
            put(b, c.sites.size(), 4);
            for (auto site : c.sites)
                put(b, site, 8);
        }
    put(b, header_.norms.size(), 4);
    for (auto n : header_.norms)
        str(b, n);
    for (unsigned i = 0; i < 8; ++i)
        b[8 + i] = static_cast<uint8_t>(uint64_t(b.size()) >> (8 * i));
    // Stage beside the destination and publish only after every backpatch and
    // fclose succeeds. A failed conversion never clobbers an existing output.
    temporary_ = output + ".tmp.XXXXXX";
    std::vector<char> path(temporary_.begin(), temporary_.end());
    path.push_back(0);
    int fd = mkstemp(path.data());
    check(fd >= 0, "cannot create output temporary file");
    temporary_ = path.data();
    file_ = fdopen(fd, "w+b");
    if (!file_) {
        close(fd);
        std::remove(temporary_.c_str());
        throw std::runtime_error("V10: fdopen failed");
    }
    try {
        write(b);
    } catch (...) {
        std::fclose(file_);
        file_ = nullptr;
        std::remove(temporary_.c_str());
        throw;
    }
}
Writer::~Writer() {
    if (file_)
        std::fclose(file_);
    if (!temporary_.empty())
        std::remove(temporary_.c_str());
}
void Writer::matrix(uint32_t a, uint32_t b, const std::function<Matrix(uint8_t, uint32_t)> &load) {
    check(a <= b && b < header_.chromosomes.size(), "invalid matrix pair");
    for (auto e : entries_)
        check(e.a != a || e.b != b, "duplicate matrix pair");
    uint32_t n = narrow(header_.resolutions[0].size() + header_.resolutions[1].size());
    Bytes meta;
    magic(meta, "H10M");
    put(meta, 1, 4);
    put(meta, a, 4);
    put(meta, b, 4);
    put(meta, n, 4);
    put(meta, 0, 4);
    meta.resize(24 + uint64_t(n) * 76, 0);
    uint64_t metaPos = write(meta);
    entries_.push_back({a, b, metaPos, meta.size()});
    uint32_t ordinal = 0;
    for (uint8_t u = 0; u < 2; ++u)
        for (uint32_t ri = 0; ri < header_.resolutions[u].size(); ++ri, ++ordinal) {
            auto r = header_.resolutions[u][ri];
            Matrix m = load(u, ri);
            canonical(m, a, b, header_, u, ri);
            if (r.mode) {
                Matrix src = load(u, r.source);
                canonical(src, a, b, header_, u, r.source);
                Matrix derived = aggregate(src, r.bin / header_.resolutions[u][r.source].bin);
                check(derived.scores == m.scores && derived.cells.size() == m.cells.size(),
                      "requested derivation changes source matrix");
                for (size_t i = 0; i < m.cells.size(); ++i)
                    check(m.cells[i].x == derived.cells[i].x &&
                              m.cells[i].y == derived.cells[i].y &&
                              m.cells[i].value == derived.cells[i].value,
                          "requested derivation is not bitwise lossless");
            }
            const uint64_t occupied = m.cells.size();
            uint64_t countSum = 0;
            double scoreSum = 0;
            for (auto c : m.cells) {
                if (m.scores)
                    scoreSum += double(floating(static_cast<uint32_t>(c.value)));
                else
                    countSum = plus(countSum, c.value);
            }
            uint64_t columnBins = header_.bins(a, u, ri);
            uint64_t rowBins = header_.bins(b, u, ri);
            uint32_t blockBins = safe_block_bin_count(columnBins, rowBins, options_.blockBins);
            uint32_t columns = narrow(ceil_div(columnBins, blockBins));
            struct Page {
                uint32_t first, last, raw;
                uint64_t pos, len;
            };
            std::vector<Page> pages;
            uint64_t indexPos = 0, indexLen = 0, blockCount = 0;
            if (!r.mode && !m.cells.empty()) {
                auto block_number = [&](const Cell &c) {
                    return narrow(uint64_t(c.y / blockBins) * columns + c.x / blockBins);
                };
                // Reorder the existing matrix storage instead of duplicating every
                // cell into a map of block vectors.
                std::sort(m.cells.begin(), m.cells.end(), [&](const Cell &x, const Cell &y) {
                    return std::make_tuple(block_number(x), x.y, x.x) <
                           std::make_tuple(block_number(y), y.y, y.x);
                });
                std::vector<std::pair<uint32_t, Bytes>> pending;
                size_t pendingSize = 0;
                std::deque<std::future<EncodedPage>> queued;
                // Bound queued raw pages to about 64 MiB even with very large
                // pages or a very large requested thread count.
                const size_t byMemory = std::max<size_t>(
                    1, (64u * 1024u * 1024u) / options_.pageBytes);
                const size_t window = std::max<size_t>(
                    1, std::min<size_t>(options_.threads, byMemory));
                auto consume = [&]() {
                    EncodedPage encoded = queued.front().get();
                    queued.pop_front();
                    uint64_t pos = write(encoded.stored);
                    pages.push_back({encoded.first, encoded.last, encoded.raw, pos,
                                     encoded.stored.size()});
                };
                auto flush = [&]() {
                    if (pending.empty())
                        return;
                    auto page = std::move(pending);
                    queued.push_back(pool_->submit(
                        [page = std::move(page), level = options_.level]() mutable {
                            return encode_page(std::move(page), level);
                        }));
                    pending = {};
                    pendingSize = 0;
                    if (queued.size() >= window)
                        consume();
                };
                for (size_t begin = 0; begin < m.cells.size();) {
                    uint32_t number = block_number(m.cells[begin]);
                    size_t end = begin + 1;
                    while (end < m.cells.size() && block_number(m.cells[end]) == number)
                        ++end;
                    auto encoded = block(m.cells, begin, end, m.scores);
                    if (pendingSize && pendingSize + encoded.size() > options_.pageBytes)
                        flush();
                    pendingSize += encoded.size();
                    pending.emplace_back(number, std::move(encoded));
                    ++blockCount;
                    begin = end;
                }
                flush();
                while (!queued.empty())
                    consume();
                m.cells.clear();
                m.cells.shrink_to_fit();
                Bytes blob, checkpoints;
                constexpr uint32_t interval = 64;
                for (size_t i = 0; i < pages.size(); ++i) {
                    const auto &p = pages[i];
                    if (i % interval == 0) {
                        put(checkpoints, i, 4);
                        put(checkpoints, std::min<size_t>(interval, pages.size() - i), 4);
                        put(checkpoints, p.first, 4);
                        put(checkpoints, 0, 4);
                        put(checkpoints, p.pos, 8);
                        put(checkpoints, blob.size(), 8);
                    } else
                        var(blob, p.first - pages[i - 1].last - 1);
                    var(blob, p.last - p.first);
                    var(blob, p.len);
                    var(blob, p.raw);
                }
                Bytes idx;
                magic(idx, "H10I");
                put(idx, 1, 4);
                put(idx, pages.size(), 4);
                put(idx, interval, 4);
                put(idx, (pages.size() + interval - 1) / interval, 4);
                put(idx, 0, 4);
                put(idx, blob.size(), 8);
                append(idx, checkpoints);
                append(idx, blob);
                indexPos = write(idx);
                indexLen = idx.size();
            }
            Bytes desc;
            put(desc, u, 1);
            put(desc, r.mode, 1);
            put(desc, r.aggregation, 1);
            put(desc, m.scores, 1);
            put(desc, ri, 4);
            put(desc, r.bin, 4);
            put(desc, r.source, 4);
            put(desc, 0, 1);
            put(desc, 0, 3);
            put(desc, m.scores ? bits(scoreSum) : countSum, 8);
            put(desc, occupied, 8);
            put(desc, 0x7fc00000, 4);
            put(desc, 0x7fc00000, 4);
            put(desc, blockBins, 4);
            put(desc, columns, 4);
            put(desc, indexPos, 8);
            put(desc, indexLen, 8);
            put(desc, pages.size(), 4);
            put(desc, blockCount, 4);
            patch(metaPos + 24 + uint64_t(ordinal) * 76, desc);
        }
}
void Writer::finish(const std::vector<Vector> &vectors) {
    std::array<std::vector<const Vector *>, 3> groups;
    for (const auto &v : vectors) {
        check(v.kind <= 2, "unknown vector kind");
        groups[v.kind].push_back(&v);
    }
    for (uint8_t kind = 0; kind < 3; ++kind) {
        auto &group = groups[kind];
        if (group.empty())
            continue;
        auto key = [&](const Vector *v) {
            return std::make_tuple(kind == 1 ? 0 : v->norm, kind == 0 ? v->chr : 0, v->unit, v->ri);
        };
        std::sort(group.begin(), group.end(), [&](auto a, auto b) { return key(a) < key(b); });
        Bytes index;
        magic(index, kind == 0 ? "NVI0" : kind == 1 ? "EVI0" : "NEVI");
        put(index, 1, 4);
        put(index, group.size(), 4);
        put(index, 0, 4);
        for (size_t i = 0; i < group.size(); ++i) {
            const auto &v = *group[i];
            check(!i || key(group[i - 1]) < key(group[i]), "duplicate vector key");
            check(v.unit <= 1 && v.ri < header_.resolutions[v.unit].size() &&
                      (kind == 1 || v.norm < header_.norms.size()) &&
                      (kind != 0 || v.chr < header_.chromosomes.size()),
                  "invalid vector key");
            uint64_t required = kind == 0 ? header_.bins(v.chr, v.unit, v.ri) : 0;
            if (kind)
                for (uint32_t c = 0; c < header_.chromosomes.size(); ++c)
                    required = std::max(required, header_.bins(c, v.unit, v.ri));
            check(v.values.size() == required, "vector length does not match V10 bin count");
            Bytes chunks;
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
                            put(raw,
                                v.values[begin + j] ^ (t == 2 && j ? v.values[begin + j - 1] : 0),
                                4);
                    auto frame = compressed(raw, options_.level);
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
                auto pos = write(stored);
                put(chunks, begin, 8);
                put(chunks, n, 4);
                put(chunks, transform, 1);
                put(chunks, 1, 1);
                put(chunks, 0, 2);
                put(chunks, pos, 8);
                put(chunks, narrow(stored.size()), 4);
                put(chunks, uint64_t(n) * 4, 4);
            }
            Bytes entry;
            put(entry, 0, 4);
            if (kind != 1)
                put(entry, v.norm, 4);
            if (kind == 0)
                put(entry, v.chr, 4);
            put(entry, v.unit, 1);
            put(entry, 0, 3);
            put(entry, v.ri, 4);
            put(entry, header_.resolutions[v.unit][v.ri].bin, 4);
            put(entry, v.values.size(), 8);
            put(entry, nominal, 4);
            put(entry, chunks.size() / 32, 4);
            if (kind) {
                put(entry, v.scales.size(), 4);
                put(entry, 0, 4);
                for (auto s : v.scales) {
                    check(s.first < header_.chromosomes.size(), "invalid scale chromosome");
                    put(entry, s.first, 4);
                    put(entry, s.second, 4);
                }
            }
            append(entry, chunks);
            uint32_t len = narrow(entry.size());
            for (unsigned j = 0; j < 4; ++j)
                entry[j] = static_cast<uint8_t>(len >> (8 * j));
            append(index, entry);
        }
        uint64_t pos = write(index);
        Bytes loc;
        put(loc, pos, 8);
        put(loc, index.size(), 8);
        patch(32 + 16 * kind, loc);
    }
    std::sort(entries_.begin(), entries_.end(),
              [](auto a, auto b) { return std::tie(a.a, a.b) < std::tie(b.a, b.b); });
    Bytes footer;
    magic(footer, "H10F");
    put(footer, 1, 4);
    put(footer, 24 + uint64_t(entries_.size()) * 24, 8);
    put(footer, entries_.size(), 4);
    put(footer, 0, 4);
    for (auto e : entries_) {
        put(footer, e.a, 4);
        put(footer, e.b, 4);
        put(footer, e.pos, 8);
        put(footer, e.len, 8);
    }
    auto pos = write(footer);
    Bytes loc;
    put(loc, pos, 8);
    put(loc, footer.size(), 8);
    patch(16, loc);
    check(std::fflush(file_) == 0, "cannot flush output");
    check(fsync(fileno(file_)) == 0, "cannot sync output");
    auto f = file_;
    file_ = nullptr;
    check(std::fclose(f) == 0, "cannot close output");
    check(std::rename(temporary_.c_str(), output_.c_str()) == 0, "cannot publish output");
    temporary_.clear();
}
} // namespace hic10
