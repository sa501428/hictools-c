#include "writer.h"

#include "build.h"
#include "common/thread_pool.h"
#include "io.h"
#include "manifest.h"
#include "sort.h"
#include "vectors.h"
#include "v10/format.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <deque>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <sstream>
#include <tuple>
#include <unistd.h>
#include <zstd.h>

namespace hic10large {
namespace {

using hic10::Bytes;

struct BlockCell {
    uint32_t block = 0, x = 0, y = 0, reserved = 0;
    uint64_t count = 0;
};
static_assert(sizeof(BlockCell) == 24, "unexpected block-cell padding");

constexpr size_t BLOCK_RUN_HEADER_BYTES = 64;

struct BlockRun {
    std::string path;
    uint32_t chr1 = 0, chr2 = 0, resolution = 0;
    uint64_t records = 0, checksum = 0;
};

struct MatrixEntry {
    uint32_t a = 0, b = 0;
    uint64_t position = 0, length = 0;
};

class RemovePath {
  public:
    explicit RemovePath(std::string path) : path_(std::move(path)) {}
    ~RemovePath() { if (!path_.empty()) std::remove(path_.c_str()); }
    void release() { path_.clear(); }
  private:
    std::string path_;
};

bool block_less(const BlockCell &a, const BlockCell &b) {
    return std::tie(a.block, a.y, a.x) < std::tie(b.block, b.y, b.x);
}

std::array<unsigned char, BLOCK_RUN_HEADER_BYTES> block_run_header(const BlockRun &r) {
    std::array<unsigned char, BLOCK_RUN_HEADER_BYTES> b{};
    std::memcpy(b.data(), "H10Q", 4);
    put_u32(b.data() + 4, 1);
    put_u32(b.data() + 8, r.chr1);
    put_u32(b.data() + 12, r.chr2);
    put_u32(b.data() + 16, r.resolution);
    put_u32(b.data() + 20, sizeof(BlockCell));
    put_u64(b.data() + 24, r.records);
    put_u64(b.data() + 32, r.checksum);
    return b;
}

class BlockRunWriter {
  public:
    BlockRunWriter(const std::string &path, uint32_t a, uint32_t b, uint32_t resolution)
        : final_(path), temporary_(path + ".tmp-" + std::to_string(getpid())),
          info_{path, a, b, resolution, 0, 14695981039346656037ULL},
          output_(temporary_, "w+b") {
        auto header = block_run_header(info_);
        output_.write(header.data(), header.size());
        buffer_.reserve(65536);
    }
    ~BlockRunWriter() {
        if (!finished_) std::remove(temporary_.c_str());
    }
    void add(const BlockCell &cell) {
        require(cell.count, "zero block cell");
        buffer_.push_back(cell);
        info_.checksum = fnv1a(&cell, sizeof(cell), info_.checksum);
        ++info_.records;
        if (buffer_.size() == buffer_.capacity()) flush();
    }
    BlockRun finish() {
        require(!finished_, "block run writer finished twice");
        flush();
        auto header = block_run_header(info_);
        output_.seek(0);
        output_.write(header.data(), header.size());
        output_.flush();
        output_.close();
        atomic_rename(temporary_, final_);
        finished_ = true;
        return info_;
    }
  private:
    void flush() {
        if (buffer_.empty()) return;
        output_.write(buffer_.data(), buffer_.size() * sizeof(BlockCell));
        buffer_.clear();
    }
    std::string final_, temporary_;
    BlockRun info_;
    File output_;
    std::vector<BlockCell> buffer_;
    bool finished_ = false;
};

class BlockRunReader {
  public:
    explicit BlockRunReader(const BlockRun &run)
        : info_(run), input_(run.path, "rb"), buffer_(65536), remaining_(run.records) {
        std::array<unsigned char, BLOCK_RUN_HEADER_BYTES> b{};
        input_.read(b.data(), b.size());
        require(std::memcmp(b.data(), "H10Q", 4) == 0 && get_u32(b.data() + 4) == 1 &&
                    get_u32(b.data() + 20) == sizeof(BlockCell) &&
                    get_u32(b.data() + 8) == run.chr1 && get_u32(b.data() + 12) == run.chr2 &&
                    get_u32(b.data() + 16) == run.resolution &&
                    get_u64(b.data() + 24) == run.records && get_u64(b.data() + 32) == run.checksum,
                "invalid block run " + run.path);
    }
    bool next(BlockCell &cell) {
        if (at_ == available_) {
            if (!remaining_) {
                if (!verified_) {
                    require(checksum_ == info_.checksum,
                            "block run checksum mismatch " + info_.path);
                    verified_ = true;
                }
                return false;
            }
            available_ = static_cast<size_t>(std::min<uint64_t>(buffer_.size(), remaining_));
            input_.read(buffer_.data(), available_ * sizeof(BlockCell));
            checksum_ = fnv1a(buffer_.data(), available_ * sizeof(BlockCell), checksum_);
            remaining_ -= available_;
            at_ = 0;
        }
        cell = buffer_[at_++];
        return true;
    }
  private:
    BlockRun info_;
    File input_;
    std::vector<BlockCell> buffer_;
    size_t at_ = 0, available_ = 0;
    uint64_t remaining_ = 0;
    uint64_t checksum_ = 14695981039346656037ULL;
    bool verified_ = false;
};

void radix_sort_blocks(std::vector<BlockCell> &records) {
    if (records.size() < 2) return;
    std::vector<BlockCell> temporary(records.size());
    std::vector<BlockCell> *source = &records, *target = &temporary;
    for (unsigned pass = 0; pass < 12; ++pass) {
        std::array<size_t, 256> offsets{};
        auto byte = [&](const BlockCell &cell) -> unsigned {
            if (pass < 4) return (cell.x >> (8 * pass)) & 255;
            if (pass < 8) return (cell.y >> (8 * (pass - 4))) & 255;
            return (cell.block >> (8 * (pass - 8))) & 255;
        };
        for (const auto &record : *source) ++offsets[byte(record)];
        size_t total = 0;
        for (size_t &offset : offsets) {
            size_t count = offset;
            offset = total;
            total += count;
        }
        for (const auto &record : *source) (*target)[offsets[byte(record)]++] = record;
        std::swap(source, target);
    }
    if (source != &records) records.swap(*source);
}

uint32_t choose_block_bins(uint64_t column_bins, uint64_t row_bins, uint32_t resolution,
                           uint32_t requested, bool rotated) {
    uint32_t result = hic10::adaptive_block_bin_count(column_bins, row_bins, resolution,
                                                      requested, rotated);
    if (rotated) {
        while (true) {
            uint64_t columns = hic10::ceil_div(column_bins, result);
            uint64_t maximum = uint64_t(hic10::rotated_depth(row_bins - 1, result)) * columns +
                               (row_bins - 1) / result;
            if (maximum <= UINT32_MAX) break;
            require(result < UINT32_MAX, "cannot represent rotated V10 blocks");
            result = static_cast<uint32_t>(std::min<uint64_t>(UINT32_MAX, uint64_t(result) * 2));
        }
    }
    return result;
}

std::vector<BlockRun> create_block_runs(const RunInfo &cells, const std::string &directory,
                                        uint64_t memory_bytes, uint32_t block_bins,
                                        uint32_t columns, bool rotated) {
    uint64_t capacity = memory_bytes / (2 * sizeof(BlockCell));
    require(capacity >= 1024, "writer sort memory is too small");
    capacity = std::min<uint64_t>(capacity, std::numeric_limits<size_t>::max());
    std::vector<BlockCell> chunk;
    chunk.reserve(static_cast<size_t>(capacity));
    std::vector<BlockRun> runs;
    auto emit = [&]() {
        if (chunk.empty()) return;
        radix_sort_blocks(chunk);
        std::ostringstream name;
        name << "block-" << cells.chr1 << '-' << cells.chr2 << '-' << cells.resolution
             << '-' << runs.size() << ".h10q";
        BlockRunWriter writer(join_path(directory, name.str()), cells.chr1, cells.chr2,
                              cells.resolution);
        for (const auto &cell : chunk) writer.add(cell);
        runs.push_back(writer.finish());
        chunk.clear();
    };
    RunReader input(cells);
    CellRecord cell;
    while (input.next(cell)) {
        uint32_t number = rotated
            ? hic10::rotated_block_number(cell.x, cell.y, block_bins, columns)
            : hic10::narrow(uint64_t(cell.y / block_bins) * columns + cell.x / block_bins);
        chunk.push_back({number, cell.x, cell.y, 0, cell.count});
        if (chunk.size() == chunk.capacity()) emit();
    }
    emit();
    return runs;
}

BlockRun merge_block_group(const std::vector<BlockRun> &runs, const std::string &path) {
    require(!runs.empty(), "cannot merge empty block-run group");
    struct Head { BlockCell cell; size_t run; };
    auto later = [](const Head &a, const Head &b) {
        if (block_less(a.cell, b.cell)) return false;
        if (block_less(b.cell, a.cell)) return true;
        return a.run > b.run;
    };
    std::priority_queue<Head, std::vector<Head>, decltype(later)> heap(later);
    std::vector<std::unique_ptr<BlockRunReader>> readers;
    for (size_t i = 0; i < runs.size(); ++i) {
        readers.emplace_back(new BlockRunReader(runs[i]));
        BlockCell cell;
        if (readers.back()->next(cell)) heap.push({cell, i});
    }
    BlockRunWriter writer(path, runs[0].chr1, runs[0].chr2, runs[0].resolution);
    BlockCell previous{};
    bool have_previous = false;
    while (!heap.empty()) {
        Head head = heap.top();
        heap.pop();
        require(!have_previous || block_less(previous, head.cell),
                "duplicate or unordered cell during block merge");
        writer.add(head.cell);
        previous = head.cell;
        have_previous = true;
        BlockCell next;
        if (readers[head.run]->next(next)) heap.push({next, head.run});
    }
    return writer.finish();
}

BlockRun bounded_block_merge(std::vector<BlockRun> runs, const std::string &directory,
                             size_t fan_in) {
    require(!runs.empty() && fan_in >= 2, "invalid block merge");
    uint64_t pass = 0;
    while (runs.size() > fan_in) {
        std::vector<BlockRun> next;
        for (size_t begin = 0; begin < runs.size(); begin += fan_in) {
            size_t end = std::min(runs.size(), begin + fan_in);
            std::vector<BlockRun> group(runs.begin() + begin, runs.begin() + end);
            std::ostringstream name;
            name << "block-merge-" << pass << '-' << begin / fan_in << '-' << runs[0].chr1
                 << '-' << runs[0].chr2 << '-' << runs[0].resolution << ".h10q";
            next.push_back(merge_block_group(group, join_path(directory, name.str())));
            for (const auto &run : group) std::remove(run.path.c_str());
        }
        runs.swap(next);
        ++pass;
    }
    std::ostringstream name;
    name << "block-final-" << runs[0].chr1 << '-' << runs[0].chr2 << '-'
         << runs[0].resolution << ".h10q";
    BlockRun result = merge_block_group(runs, join_path(directory, name.str()));
    for (const auto &run : runs) std::remove(run.path.c_str());
    return result;
}

Bytes compress_bytes(const Bytes &input, int level) {
    Bytes output(ZSTD_compressBound(input.size()));
    size_t n = ZSTD_compress(output.data(), output.size(), input.data(), input.size(), level);
    require(!ZSTD_isError(n), "Zstandard block compression failed");
    output.resize(n);
    return output;
}

uint32_t varint_size(uint64_t value) {
    uint32_t result = 1;
    while (value >= 128) { value >>= 7; ++result; }
    return result;
}

struct BlockSpan {
    uint32_t number = 0;
    uint64_t first = 0, records = 0;
};

template <typename Callback>
void scan_block_span(const BlockRun &run, const BlockSpan &span, Callback callback) {
    require(span.records && span.first <= run.records &&
                span.records <= run.records - span.first,
            "invalid logical block span");
    File input(run.path, "rb");
    std::array<unsigned char, BLOCK_RUN_HEADER_BYTES> header{};
    input.read(header.data(), header.size());
    require(std::memcmp(header.data(), "H10Q", 4) == 0 && get_u32(header.data() + 4) == 1 &&
                get_u32(header.data() + 8) == run.chr1 &&
                get_u32(header.data() + 12) == run.chr2 &&
                get_u32(header.data() + 16) == run.resolution &&
                get_u64(header.data() + 24) == run.records,
            "logical block source changed");
    input.seek(BLOCK_RUN_HEADER_BYTES + span.first * sizeof(BlockCell));
    std::vector<BlockCell> buffer(65536);
    uint64_t remaining = span.records;
    while (remaining) {
        size_t n = static_cast<size_t>(std::min<uint64_t>(remaining, buffer.size()));
        input.read(buffer.data(), n * sizeof(BlockCell));
        for (size_t i = 0; i < n; ++i) {
            require(buffer[i].block == span.number, "logical block span crosses block number");
            callback(buffer[i]);
        }
        remaining -= n;
    }
}

class StreamingZstd {
  public:
    StreamingZstd(File &output, int level) : output_(output), compressed_(1 << 20) {
        stream_ = ZSTD_createCStream();
        require(stream_ != nullptr, "cannot create Zstandard stream");
        size_t result = ZSTD_initCStream(stream_, level);
        require(!ZSTD_isError(result), "cannot initialize Zstandard stream");
    }
    ~StreamingZstd() { if (stream_) ZSTD_freeCStream(stream_); }
    void write(const void *data, size_t size) {
        ZSTD_inBuffer input{data, size, 0};
        while (input.pos < input.size) {
            ZSTD_outBuffer output{compressed_.data(), compressed_.size(), 0};
            size_t result = ZSTD_compressStream2(stream_, &output, &input, ZSTD_e_continue);
            require(!ZSTD_isError(result), "Zstandard streaming compression failed");
            if (output.pos) output_.write(compressed_.data(), output.pos);
        }
    }
    void finish() {
        ZSTD_inBuffer input{nullptr, 0, 0};
        size_t remaining = 1;
        while (remaining) {
            ZSTD_outBuffer output{compressed_.data(), compressed_.size(), 0};
            remaining = ZSTD_compressStream2(stream_, &output, &input, ZSTD_e_end);
            require(!ZSTD_isError(remaining), "Zstandard stream finalization failed");
            if (output.pos) output_.write(compressed_.data(), output.pos);
        }
    }
  private:
    File &output_;
    ZSTD_CStream *stream_ = nullptr;
    std::vector<unsigned char> compressed_;
};

class RawEmitter {
  public:
    explicit RawEmitter(StreamingZstd &sink) : sink_(sink) { buffer_.reserve(1 << 20); }
    void byte(uint8_t value) {
        buffer_.push_back(value);
        if (buffer_.size() == buffer_.capacity()) flush();
    }
    void bytes(const void *data, size_t size) {
        flush();
        if (size) sink_.write(data, size);
    }
    void var(uint64_t value) {
        do {
            uint8_t word = value & 127;
            value >>= 7;
            byte(word | (value ? 128 : 0));
        } while (value);
    }
    void flush() {
        if (buffer_.empty()) return;
        sink_.write(buffer_.data(), buffer_.size());
        buffer_.clear();
    }
  private:
    StreamingZstd &sink_;
    std::vector<unsigned char> buffer_;
};

struct EncodedBlockFile {
    uint32_t number = 0;
    uint64_t bytes = 0;
    std::string path;
    EncodedBlockFile() = default;
    EncodedBlockFile(uint32_t n, uint64_t b, std::string p)
        : number(n), bytes(b), path(std::move(p)) {}
    EncodedBlockFile(const EncodedBlockFile &) = delete;
    EncodedBlockFile &operator=(const EncodedBlockFile &) = delete;
    EncodedBlockFile(EncodedBlockFile &&other) noexcept
        : number(other.number), bytes(other.bytes), path(std::move(other.path)) {
        other.path.clear();
    }
    EncodedBlockFile &operator=(EncodedBlockFile &&other) noexcept {
        if (this != &other) {
            if (!path.empty()) std::remove(path.c_str());
            number = other.number; bytes = other.bytes; path = std::move(other.path);
            other.path.clear();
        }
        return *this;
    }
    ~EncodedBlockFile() { if (!path.empty()) std::remove(path.c_str()); }
};

EncodedBlockFile encode_block_span(const BlockRun &run, const BlockSpan &span,
                                   int level) {
    uint32_t xmin = UINT32_MAX, ymin = UINT32_MAX, xmax = 0, ymax = 0;
    uint64_t common = 0, direct_bytes = 0;
    bool first = true, all_same = true;
    scan_block_span(run, span, [&](const BlockCell &cell) {
        xmin = std::min(xmin, cell.x); ymin = std::min(ymin, cell.y);
        xmax = std::max(xmax, cell.x); ymax = std::max(ymax, cell.y);
        if (first) common = cell.count;
        else all_same = all_same && cell.count == common;
        first = false;
        uint32_t width = varint_size(cell.count);
        require(direct_bytes <= UINT64_MAX - width, "logical block value size overflow");
        direct_bytes += width;
    });
    require(!first, "cannot encode empty logical block");
    const uint32_t width = hic10::narrow(uint64_t(xmax) - xmin + 1);
    const uint32_t height = hic10::narrow(uint64_t(ymax) - ymin + 1);
    const uint64_t area = uint64_t(width) * height;
    uint64_t sparse_bytes = 0, previous = 0, ordinal = 0;
    scan_block_span(run, span, [&](const BlockCell &cell) {
        uint64_t position = uint64_t(cell.y - ymin) * width + cell.x - xmin;
        require(!ordinal || position > previous, "unordered cell inside logical block");
        uint32_t bytes = varint_size(ordinal ? position - previous : position);
        require(sparse_bytes <= UINT64_MAX - bytes, "logical block position size overflow");
        sparse_bytes += bytes;
        previous = position;
        ++ordinal;
    });
    const uint64_t bitmap_bytes = area / 8 + (area % 8 != 0);
    const bool bitmap = bitmap_bytes < sparse_bytes && bitmap_bytes <= UINT32_MAX;
    const uint64_t position_bytes = bitmap ? bitmap_bytes : sparse_bytes;
    const uint64_t value_bytes = all_same ? varint_size(common) : direct_bytes;
    require(position_bytes <= UINT32_MAX && value_bytes <= UINT32_MAX &&
                position_bytes + value_bytes <= UINT32_MAX - 40,
            "logical block payload exceeds the V10 uint32 size limit");
    const uint32_t raw_bytes = static_cast<uint32_t>(40 + position_bytes + value_bytes);
    std::string path = run.path + ".encoded-" + std::to_string(span.number) + "-" +
                       std::to_string(getpid());
    RemovePath cleanup(path);
    File output(path, "w+b");
    Bytes stored_header;
    hic10::magic(stored_header, "H10B");
    hic10::put(stored_header, 1, 1); hic10::put(stored_header, 1, 1);
    hic10::put(stored_header, 0, 2); hic10::put(stored_header, raw_bytes, 4);
    hic10::put(stored_header, span.number, 4);
    output.write(stored_header.data(), stored_header.size());
    {
        StreamingZstd compressor(output, level);
        RawEmitter raw(compressor);
        Bytes payload_header;
        hic10::put(payload_header, 1, 1);
        hic10::put(payload_header, bitmap ? 1 : 0, 1);
        hic10::put(payload_header, all_same ? 0 : 2, 1);
        hic10::put(payload_header, 0, 1); // COUNT_UINT
        hic10::put(payload_header, bitmap ? 1 : 0, 1);
        hic10::put(payload_header, 0, 3);
        hic10::put(payload_header, xmin, 4); hic10::put(payload_header, ymin, 4);
        hic10::put(payload_header, width, 4); hic10::put(payload_header, height, 4);
        hic10::put(payload_header, span.records, 8);
        hic10::put(payload_header, position_bytes, 4); hic10::put(payload_header, value_bytes, 4);
        require(payload_header.size() == 40, "invalid streamed payload header");
        raw.bytes(payload_header.data(), payload_header.size());
        if (bitmap) {
            uint64_t byte_index = 0;
            uint8_t word = 0;
            scan_block_span(run, span, [&](const BlockCell &cell) {
                uint64_t position = uint64_t(cell.y - ymin) * width + cell.x - xmin;
                uint64_t target = position / 8;
                while (byte_index < target) { raw.byte(word); word = 0; ++byte_index; }
                word |= static_cast<uint8_t>(1u << (position % 8));
            });
            while (byte_index < bitmap_bytes) { raw.byte(word); word = 0; ++byte_index; }
        } else {
            uint64_t prior = 0, index = 0;
            scan_block_span(run, span, [&](const BlockCell &cell) {
                uint64_t position = uint64_t(cell.y - ymin) * width + cell.x - xmin;
                raw.var(index ? position - prior : position);
                prior = position;
                ++index;
            });
        }
        if (all_same) {
            raw.var(common);
        } else {
            scan_block_span(run, span, [&](const BlockCell &cell) { raw.var(cell.count); });
        }
        raw.flush();
        compressor.finish();
    }
    output.flush();
    uint64_t bytes = output.tell();
    require(bytes <= UINT32_MAX, "compressed logical block exceeds V10 uint32 size limit");
    output.close();
    cleanup.release();
    return EncodedBlockFile(span.number, bytes, path);
}

struct EncodedVectorChunk {
    uint64_t begin = 0;
    uint32_t count = 0;
    uint8_t transform = 0;
    Bytes stored;
};

EncodedVectorChunk encode_vector_chunk(uint64_t begin, const std::vector<uint32_t> &words,
                                       int level) {
    uint32_t count = hic10::narrow(words.size());
    Bytes best;
    uint8_t transform = 0;
    for (uint8_t candidate = 0; candidate < 3; ++candidate) {
        Bytes raw;
        raw.reserve(uint64_t(count) * 4);
        if (candidate == 1) {
            for (unsigned lane = 0; lane < 4; ++lane)
                for (uint32_t word : words) hic10::put(raw, word >> (8 * lane), 1);
        } else {
            for (uint32_t i = 0; i < count; ++i)
                hic10::put(raw, words[i] ^ (candidate == 2 && i ? words[i - 1] : 0), 4);
        }
        Bytes frame = compress_bytes(raw, level);
        if (best.empty() || frame.size() < best.size()) {
            best = std::move(frame);
            transform = candidate;
        }
    }
    Bytes stored;
    hic10::magic(stored, "H10V"); hic10::put(stored, 1, 1);
    hic10::put(stored, transform, 1); hic10::put(stored, 0, 2);
    hic10::put(stored, uint64_t(count) * 4, 4); hic10::put(stored, count, 4);
    hic10::append(stored, best);
    return {begin, count, transform, std::move(stored)};
}

template <typename Sink>
void write_vectors(Sink &output, const hic10::Header &header,
                   const VectorManifest &manifest, int level,
                   ThreadPool *compression_pool, int threads) {
    std::array<std::vector<const VectorInfo *>, 3> groups;
    for (const auto &vector : manifest.vectors) groups.at(vector.kind).push_back(&vector);
    for (uint8_t kind = 0; kind < 3; ++kind) {
        auto &group = groups[kind];
        if (group.empty()) continue;
        auto key = [&](const VectorInfo *v) {
            return std::make_tuple(kind == 1 ? 0 : v->norm, kind == 0 ? v->chr : 0, v->ri);
        };
        std::sort(group.begin(), group.end(), [&](auto a, auto b) { return key(a) < key(b); });
        Bytes index;
        hic10::magic(index, kind == 0 ? "NVI0" : kind == 1 ? "EVI0" : "NEVI");
        hic10::put(index, 1, 4); hic10::put(index, group.size(), 4); hic10::put(index, 0, 4);
        for (size_t i = 0; i < group.size(); ++i) {
            const VectorInfo &v = *group[i];
            require(!i || key(group[i - 1]) < key(group[i]), "duplicate vector manifest key");
            require(v.ri < header.resolutions[0].size() &&
                        header.resolutions[0][v.ri].bin == v.resolution,
                    "vector resolution index mismatch");
            require((kind == 1 || v.norm < header.norms.size()) &&
                        (kind != 0 || v.chr < header.chromosomes.size()),
                    "invalid vector manifest key");
            uint64_t required = kind == 0 ? header.bins(v.chr, 0, v.ri) : 0;
            if (kind)
                for (uint32_t chr = 0; chr < header.chromosomes.size(); ++chr)
                    required = std::max(required, header.bins(chr, 0, v.ri));
            require(v.words == required, "vector length does not match V10 geometry");
            VectorFileReader reader(v);
            Bytes chunks;
            constexpr uint32_t nominal = 65536;
            std::deque<std::future<EncodedVectorChunk>> pending;
            auto store = [&](EncodedVectorChunk encoded) {
                uint64_t position = output.write(encoded.stored);
                hic10::put(chunks, encoded.begin, 8); hic10::put(chunks, encoded.count, 4);
                hic10::put(chunks, encoded.transform, 1); hic10::put(chunks, 1, 1);
                hic10::put(chunks, 0, 2); hic10::put(chunks, position, 8);
                hic10::put(chunks, encoded.stored.size(), 4);
                hic10::put(chunks, uint64_t(encoded.count) * 4, 4);
            };
            auto consume = [&]() {
                store(pending.front().get());
                pending.pop_front();
            };
            for (uint64_t begin = 0; begin < v.words; begin += nominal) {
                uint32_t count = hic10::narrow(std::min<uint64_t>(nominal, v.words - begin));
                std::vector<uint32_t> words = reader.read(begin, count);
                if (compression_pool) {
                    pending.push_back(compression_pool->submit(
                        [level](uint64_t at, std::vector<uint32_t> values) {
                            return encode_vector_chunk(at, values, level);
                        }, begin, std::move(words)));
                    if (pending.size() >= static_cast<size_t>(threads)) consume();
                } else {
                    store(encode_vector_chunk(begin, words, level));
                }
            }
            while (!pending.empty()) consume();
            Bytes entry;
            hic10::put(entry, 0, 4);
            if (kind != 1) hic10::put(entry, v.norm, 4);
            if (kind == 0) hic10::put(entry, v.chr, 4);
            hic10::put(entry, 0, 1); hic10::put(entry, 0, 3);
            hic10::put(entry, v.ri, 4); hic10::put(entry, v.resolution, 4);
            hic10::put(entry, v.words, 8); hic10::put(entry, nominal, 4);
            hic10::put(entry, chunks.size() / 32, 4);
            if (kind) {
                hic10::put(entry, v.scales.size(), 4); hic10::put(entry, 0, 4);
                for (const auto &scale : v.scales) {
                    require(scale.first < header.chromosomes.size(), "vector scale chromosome outside table");
                    hic10::put(entry, scale.first, 4); hic10::put(entry, scale.second, 4);
                }
            }
            hic10::append(entry, chunks);
            uint32_t length = hic10::narrow(entry.size());
            for (unsigned lane = 0; lane < 4; ++lane)
                entry[lane] = static_cast<uint8_t>(length >> (8 * lane));
            hic10::append(index, entry);
        }
        uint64_t position = output.write(index);
        Bytes locator; hic10::put(locator, position, 8); hic10::put(locator, index.size(), 8);
        output.patch(32 + 16 * kind, locator);
    }
}

class Output {
  public:
    Output(const std::string &output, const Bytes &header)
        : final_(output), temporary_(output + ".large-tmp-" + std::to_string(getpid())),
          file_(temporary_, "w+b") {
        file_.write(header.data(), header.size());
    }
    ~Output() {
        if (!published_ && path_exists(temporary_)) std::remove(temporary_.c_str());
    }
    uint64_t position() const { return file_.tell(); }
    uint64_t write(const Bytes &bytes) {
        uint64_t at = position(); file_.write(bytes.data(), bytes.size()); return at;
    }
    uint64_t write(const void *data, size_t size) {
        uint64_t at = position(); file_.write(data, size); return at;
    }
    uint64_t copy(const std::string &path, uint64_t bytes) {
        require(file_size(path) == bytes, "encoded block sidecar length changed");
        File input(path, "rb");
        std::vector<unsigned char> buffer(8 * 1024 * 1024);
        uint64_t at = position(), remaining = bytes;
        while (remaining) {
            size_t n = static_cast<size_t>(std::min<uint64_t>(buffer.size(), remaining));
            input.read(buffer.data(), n);
            file_.write(buffer.data(), n);
            remaining -= n;
        }
        return at;
    }
    void patch(uint64_t at, const Bytes &bytes) {
        uint64_t saved = position(); file_.seek(at); file_.write(bytes.data(), bytes.size()); file_.seek(saved);
    }
    void publish() {
        file_.flush();
        require(fsync(fileno(file_.get())) == 0, "cannot sync V10 output");
        file_.close();
        atomic_rename(temporary_, final_);
        published_ = true;
    }
  private:
    std::string final_, temporary_;
    File file_;
    bool published_ = false;
};

Bytes make_header(hic10::Header &header) {
    Bytes b;
    hic10::magic(b, "HIC\0"); hic10::put(b, 10, 4); b.resize(88, 0);
    hic10::str(b, header.genome);
    hic10::put(b, header.attributes.size(), 4);
    for (const auto &a : header.attributes) { hic10::str(b, a.first); hic10::str(b, a.second); }
    hic10::put(b, header.chromosomes.size(), 4);
    for (const auto &c : header.chromosomes) { hic10::str(b, c.name); hic10::put(b, c.length, 8); }
    for (const auto &list : header.resolutions) {
        hic10::put(b, list.size(), 4);
        for (const auto &r : list) {
            hic10::put(b, r.bin, 4); hic10::put(b, r.mode, 1); hic10::put(b, r.aggregation, 1);
            hic10::put(b, 0, 2); hic10::put(b, r.source, 4);
        }
    }
    hic10::put(b, header.norms.size(), 4);
    for (const auto &name : header.norms) hic10::str(b, name);
    for (unsigned i = 0; i < 8; ++i) b[8 + i] = static_cast<uint8_t>(uint64_t(b.size()) >> (8 * i));
    return b;
}

std::pair<uint64_t, uint64_t> run_statistics(const RunInfo &run,
                                             uint64_t columns, uint64_t rows,
                                             bool cis) {
    RunReader reader(run);
    CellRecord cell;
    uint64_t occupied = 0, sum = 0;
    while (reader.next(cell)) {
        require(cell.x < columns && cell.y < rows && (!cis || cell.x <= cell.y),
                "cell outside V10 matrix geometry");
        require(cell.count <= UINT64_MAX - sum, "matrix count sum exceeds uint64");
        sum += cell.count;
        ++occupied;
    }
    require(occupied == run.records, "cell run record count changed");
    return {occupied, sum};
}

} // namespace

void write_v10(const std::string &stage_path, const std::string &build_path,
               const std::string &output_path, const WriteOptions &options) {
    StageManifest stage = read_stage_manifest(stage_path);
    BuildManifest build = read_build_manifest(build_path);
    require(build.source_fingerprint == stage.source_fingerprint,
            "stage and build manifests refer to different sources");
    require(options.memory_bytes >= 1024 * 1024, "writer memory budget is too small");
    require(options.block_bins && options.block_bins <= 4096, "invalid block bin minimum");
    require(options.threads > 0, "writer thread count must be positive");
    require(options.compression_level >= ZSTD_minCLevel() &&
                options.compression_level <= ZSTD_maxCLevel(), "invalid Zstandard level");
    std::string temporary_directory = options.temporary_directory.empty()
        ? join_path(build_path.substr(0, build_path.find_last_of('/')), "block-work")
        : options.temporary_directory;
    if (temporary_directory == "block-work") temporary_directory = "./block-work";
    make_directory(temporary_directory);

    hic10::Header header;
    header.genome = options.genome;
    header.attributes = {{"software", "hictools-c hic_v10_large"}, {"hicFileScaling", "1.0"}};
    for (const auto &chromosome : stage.chromosomes)
        header.chromosomes.push_back({chromosome.name, chromosome.length, {}});
    for (uint32_t resolution : build.resolutions) header.resolutions[0].push_back({resolution});
    VectorManifest vectors;
    if (!options.vector_manifest.empty()) {
        vectors = read_vector_manifest(options.vector_manifest);
        require(vectors.source_fingerprint == stage.source_fingerprint,
                "vector manifest refers to a different source");
        header.norms = vectors.norms;
    }
    auto derive = [&](uint32_t target_bin, uint32_t source_bin) {
        uint32_t target = header.resolution(0, target_bin), source = header.resolution(0, source_bin);
        require(target_bin > source_bin && target_bin % source_bin == 0,
                "derived target must be an integer multiple of source");
        require(!hic10::required_materialized_resolution(target_bin),
                "500 kb must remain materialized");
        auto &r = header.resolutions[0][target];
        require(!r.mode || r.source == source, "conflicting derived source");
        r.mode = 1; r.source = source; r.aggregation = 1;
    };
    for (const auto &d : options.derived) derive(d.first, d.second);
    for (const auto &r : header.resolutions[0])
        if (uint32_t source = hic10::required_derived_source(r.bin)) derive(r.bin, source);
    for (const auto &r : header.resolutions[0])
        if (r.mode) require(!header.resolutions[0][r.source].mode, "chained derivation is forbidden");
    require(hic10::required_bp_resolution_policy(header.resolutions[0]),
            "mandatory V10 resolution policy is not satisfied");

    Output output(output_path, make_header(header));
    std::unique_ptr<ThreadPool> compression_pool(
        options.threads > 1 ? new ThreadPool(static_cast<size_t>(options.threads)) : nullptr);
    std::vector<MatrixEntry> matrices;
    for (const PairInfo &pair : stage.pairs) {
        uint32_t a = pair.chr1, b = pair.chr2;
        const uint32_t n = hic10::narrow(header.resolutions[0].size());
        Bytes meta;
        hic10::magic(meta, "H10M"); hic10::put(meta, 1, 4);
        hic10::put(meta, a, 4); hic10::put(meta, b, 4); hic10::put(meta, n, 4); hic10::put(meta, 0, 4);
        meta.resize(24 + uint64_t(n) * 76, 0);
        uint64_t meta_position = output.write(meta);
        matrices.push_back({a, b, meta_position, meta.size()});
        std::cerr << "Writing V10 matrix " << stage.chromosomes[a].name << " x "
                  << stage.chromosomes[b].name << "\n";
        for (uint32_t ri = 0; ri < header.resolutions[0].size(); ++ri) {
            const auto &resolution = header.resolutions[0][ri];
            auto found = build.cells.find(std::make_tuple(a, b, resolution.bin));
            require(found != build.cells.end(), "build manifest is missing a matrix resolution");
            const RunInfo &cells = found->second;
            uint64_t column_bins = header.bins(a, 0, ri), row_bins = header.bins(b, 0, ri);
            bool rotated = a == b;
            auto statistics = run_statistics(cells, column_bins, row_bins, rotated);
            uint32_t block_bins = choose_block_bins(column_bins, row_bins, resolution.bin,
                                                    options.block_bins, rotated);
            uint32_t columns = hic10::narrow(hic10::ceil_div(column_bins, block_bins));
            uint32_t index_count = 0;
            uint64_t index_position = 0, index_length = 0;
            if (!resolution.mode && cells.records) {
                auto runs = create_block_runs(cells, temporary_directory, options.memory_bytes,
                                              block_bins, columns, rotated);
                BlockRun blocks = bounded_block_merge(std::move(runs), temporary_directory,
                                                      options.merge_fan_in);
                RemovePath remove_blocks(blocks.path);
                const std::string index_path = blocks.path + ".index-tmp-" +
                                               std::to_string(getpid());
                RemovePath remove_index(index_path);
                {
                    File index_file(index_path, "w+b");
                    {
                        BlockRunReader reader(blocks);
                        BlockCell cell;
                        std::deque<std::future<EncodedBlockFile>> pending;
                        BlockSpan span;
                        uint64_t record_index = 0;
                        auto store = [&](EncodedBlockFile encoded) {
                            uint64_t position = output.copy(encoded.path, encoded.bytes);
                            require(index_count < UINT32_MAX,
                                    "too many V10 blocks in one resolution");
                            ++index_count;
                            std::array<unsigned char, 16> entry{};
                            put_u32(entry.data(), encoded.number);
                            put_u32(entry.data() + 4, hic10::narrow(encoded.bytes));
                            put_u64(entry.data() + 8, position);
                            index_file.write(entry.data(), entry.size());
                        };
                        auto consume = [&]() {
                            store(pending.front().get());
                            pending.pop_front();
                        };
                        auto submit = [&]() {
                            if (!span.records) return;
                            if (compression_pool) {
                                const int level = options.compression_level;
                                pending.push_back(compression_pool->submit(
                                    [level](BlockRun source, BlockSpan block) {
                                        return encode_block_span(source, block, level);
                                    }, blocks, span));
                                if (pending.size() >= static_cast<size_t>(options.threads))
                                    consume();
                            } else {
                                store(encode_block_span(blocks, span,
                                                        options.compression_level));
                            }
                        };
                        while (reader.next(cell)) {
                            if (!span.records || cell.block != span.number) {
                                submit();
                                span = {cell.block, record_index, 0};
                            }
                            ++span.records;
                            ++record_index;
                        }
                        submit();
                        while (!pending.empty()) consume();
                    }
                    index_file.flush();
                    index_file.seek(0);
                    Bytes header;
                    hic10::magic(header, "H10I"); hic10::put(header, 2, 4);
                    index_length = uint64_t(24) + uint64_t(index_count) * 16;
                    hic10::put(header, index_length, 8);
                    hic10::put(header, index_count, 4); hic10::put(header, 0, 4);
                    index_position = output.write(header);
                    std::vector<unsigned char> copy(8 * 1024 * 1024);
                    uint64_t remaining = uint64_t(index_count) * 16;
                    while (remaining) {
                        size_t n = static_cast<size_t>(std::min<uint64_t>(copy.size(), remaining));
                        index_file.read(copy.data(), n);
                        output.write(copy.data(), n);
                        remaining -= n;
                    }
                    index_file.close();
                }
            }
            Bytes descriptor;
            hic10::put(descriptor, 0, 1); hic10::put(descriptor, resolution.mode, 1);
            hic10::put(descriptor, resolution.aggregation, 1); hic10::put(descriptor, 0, 1);
            hic10::put(descriptor, ri, 4); hic10::put(descriptor, resolution.bin, 4);
            hic10::put(descriptor, resolution.source, 4); hic10::put(descriptor, rotated, 1);
            hic10::put(descriptor, 0, 3); hic10::put(descriptor, statistics.second, 8);
            hic10::put(descriptor, statistics.first, 8);
            hic10::put(descriptor, 0x7fc00000, 4); hic10::put(descriptor, 0x7fc00000, 4);
            hic10::put(descriptor, block_bins, 4); hic10::put(descriptor, columns, 4);
            hic10::put(descriptor, index_position, 8); hic10::put(descriptor, index_length, 8);
            hic10::put(descriptor, index_count, 4); hic10::put(descriptor, 0, 4);
            output.patch(meta_position + 24 + uint64_t(ri) * 76, descriptor);
        }
    }
    if (!options.vector_manifest.empty())
        write_vectors(output, header, vectors, options.compression_level,
                      compression_pool.get(), options.threads);
    std::sort(matrices.begin(), matrices.end(), [](const auto &x, const auto &y) {
        return std::tie(x.a, x.b) < std::tie(y.a, y.b);
    });
    Bytes footer;
    hic10::magic(footer, "H10F"); hic10::put(footer, 1, 4);
    hic10::put(footer, uint64_t(24) + uint64_t(matrices.size()) * 24, 8);
    hic10::put(footer, matrices.size(), 4); hic10::put(footer, 0, 4);
    for (const auto &matrix : matrices) {
        hic10::put(footer, matrix.a, 4); hic10::put(footer, matrix.b, 4);
        hic10::put(footer, matrix.position, 8); hic10::put(footer, matrix.length, 8);
    }
    uint64_t footer_position = output.write(footer);
    Bytes locator; hic10::put(locator, footer_position, 8); hic10::put(locator, footer.size(), 8);
    output.patch(16, locator);
    output.publish();
    std::cerr << "Large-data V10 matrix assembly complete: " << output_path << "\n";
}

} // namespace hic10large
