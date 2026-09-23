#include "writer.h"

#include "build.h"
#include "common/thread_pool.h"
#include "io.h"
#include "manifest.h"
#include "sort.h"
#include "vectors.h"
#include "v10/format.h"
#include "v10/reader.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
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

class MatrixStatistics {
  public:
    MatrixStatistics(uint64_t columns, uint64_t rows, bool cis)
        : columns_(columns), rows_(rows), cis_(cis) {}
    void add(const CellRecord &cell) {
        require(cell.x < columns_ && cell.y < rows_ && (!cis_ || cell.x <= cell.y),
                "cell outside V10 matrix geometry");
        require(cell.count <= UINT64_MAX - sum_, "matrix count sum exceeds uint64");
        sum_ += cell.count;
        require(occupied_ < UINT64_MAX, "matrix occupancy exceeds uint64");
        ++occupied_;
    }
    std::pair<uint64_t, uint64_t> result() const { return {occupied_, sum_}; }
  private:
    uint64_t columns_ = 0, rows_ = 0, occupied_ = 0, sum_ = 0;
    bool cis_ = false;
};

struct BlockBuild {
    std::vector<BlockRun> runs;
    std::pair<uint64_t, uint64_t> statistics;
};

class BlockRunBuilder {
  public:
    BlockRunBuilder(RunInfo metadata, std::string directory, uint64_t memory_bytes,
                    uint32_t block_bins, uint32_t block_columns,
                    uint64_t coordinate_columns, uint64_t coordinate_rows, bool rotated)
        : metadata_(std::move(metadata)), directory_(std::move(directory)),
          block_bins_(block_bins), block_columns_(block_columns), rotated_(rotated),
          statistics_(coordinate_columns, coordinate_rows, rotated) {
        uint64_t capacity = memory_bytes / (2 * sizeof(BlockCell));
        require(capacity >= 1024, "writer sort memory is too small");
        capacity = std::min<uint64_t>(capacity, std::numeric_limits<size_t>::max());
        chunk_.reserve(static_cast<size_t>(capacity));
    }
    void add(const CellRecord &cell) {
        statistics_.add(cell);
        uint32_t number = rotated_
            ? hic10::rotated_block_number(cell.x, cell.y, block_bins_, block_columns_)
            : hic10::narrow(uint64_t(cell.y / block_bins_) * block_columns_ +
                            cell.x / block_bins_);
        chunk_.push_back({number, cell.x, cell.y, 0, cell.count});
        if (chunk_.size() == chunk_.capacity()) emit();
    }
    BlockBuild finish() {
        emit();
        return {std::move(runs_), statistics_.result()};
    }
  private:
    void emit() {
        if (chunk_.empty()) return;
        radix_sort_blocks(chunk_);
        std::ostringstream name;
        name << "block-" << metadata_.chr1 << '-' << metadata_.chr2 << '-'
             << metadata_.resolution << '-' << runs_.size() << ".h10q";
        BlockRunWriter writer(join_path(directory_, name.str()), metadata_.chr1,
                              metadata_.chr2, metadata_.resolution);
        for (const auto &cell : chunk_) writer.add(cell);
        runs_.push_back(writer.finish());
        chunk_.clear();
    }
    RunInfo metadata_;
    std::string directory_;
    uint32_t block_bins_ = 0, block_columns_ = 0;
    bool rotated_ = false;
    MatrixStatistics statistics_;
    std::vector<BlockCell> chunk_;
    std::vector<BlockRun> runs_;
};

BlockBuild create_block_runs(const RunInfo &cells, const std::string &directory,
                             uint64_t memory_bytes, uint32_t block_bins,
                             uint32_t block_columns, uint64_t coordinate_columns,
                             uint64_t coordinate_rows, bool rotated) {
    BlockRunBuilder builder(cells, directory, memory_bytes, block_bins, block_columns,
                            coordinate_columns, coordinate_rows, rotated);
    RunReader input(cells);
    CellRecord cell;
    while (input.next(cell)) builder.add(cell);
    BlockBuild result = builder.finish();
    require(result.statistics.first == cells.records, "cell run record count changed");
    return result;
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
    if (runs.size() == 1) return std::move(runs.front());
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
    Bytes data;
    EncodedBlockFile() = default;
    EncodedBlockFile(uint32_t n, uint64_t b, std::string p)
        : number(n), bytes(b), path(std::move(p)) {}
    EncodedBlockFile(uint32_t n, Bytes b)
        : number(n), bytes(b.size()), data(std::move(b)) {}
    EncodedBlockFile(const EncodedBlockFile &) = delete;
    EncodedBlockFile &operator=(const EncodedBlockFile &) = delete;
    EncodedBlockFile(EncodedBlockFile &&other) noexcept
        : number(other.number), bytes(other.bytes), path(std::move(other.path)),
          data(std::move(other.data)) {
        other.path.clear();
    }
    EncodedBlockFile &operator=(EncodedBlockFile &&other) noexcept {
        if (this != &other) {
            if (!path.empty()) std::remove(path.c_str());
            number = other.number; bytes = other.bytes; path = std::move(other.path);
            data = std::move(other.data);
            other.path.clear();
        }
        return *this;
    }
    ~EncodedBlockFile() { if (!path.empty()) std::remove(path.c_str()); }
};

EncodedBlockFile encode_block_span(const BlockRun &run, const BlockSpan &span,
                                   int level) {
    // Most logical blocks are small. Reading each one once and reusing it for
    // the sizing/representation/emission passes avoids millions of open/seek
    // cycles on deep maps. Exceptionally large blocks retain the fully
    // streaming path, so memory remains bounded independently of occupancy.
    constexpr uint64_t CACHE_BYTES = 8ULL * 1024 * 1024;
    std::vector<BlockCell> cached;
    if (span.records <= CACHE_BYTES / sizeof(BlockCell)) {
        cached.reserve(static_cast<size_t>(span.records));
        scan_block_span(run, span, [&](const BlockCell &cell) { cached.push_back(cell); });
    }
    auto scan = [&](auto callback) {
        if (!cached.empty()) {
            for (const BlockCell &cell : cached) callback(cell);
        } else {
            scan_block_span(run, span, callback);
        }
    };
    uint32_t xmin = UINT32_MAX, ymin = UINT32_MAX, xmax = 0, ymax = 0;
    uint64_t common = 0, direct_bytes = 0;
    bool first = true, all_same = true;
    scan([&](const BlockCell &cell) {
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
    scan([&](const BlockCell &cell) {
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

    if (!cached.empty()) {
        Bytes raw = payload_header;
        raw.reserve(raw_bytes);
        if (bitmap) {
            raw.resize(40 + static_cast<size_t>(bitmap_bytes), 0);
            for (const BlockCell &cell : cached) {
                uint64_t position = uint64_t(cell.y - ymin) * width + cell.x - xmin;
                raw[40 + position / 8] |= static_cast<uint8_t>(1u << (position % 8));
            }
        } else {
            uint64_t prior = 0;
            for (size_t i = 0; i < cached.size(); ++i) {
                uint64_t position = uint64_t(cached[i].y - ymin) * width +
                                    cached[i].x - xmin;
                hic10::var(raw, i ? position - prior : position);
                prior = position;
            }
        }
        if (all_same) hic10::var(raw, common);
        else for (const BlockCell &cell : cached) hic10::var(raw, cell.count);
        require(raw.size() == raw_bytes, "cached logical block size mismatch");
        const size_t bound = ZSTD_compressBound(raw.size());
        Bytes stored(16 + bound);
        std::memcpy(stored.data(), "H10B", 4);
        stored[4] = 1; stored[5] = 1;
        put_u32(stored.data() + 8, raw_bytes);
        put_u32(stored.data() + 12, span.number);
        size_t compressed = ZSTD_compress(stored.data() + 16, bound,
                                          raw.data(), raw.size(), level);
        require(!ZSTD_isError(compressed), "Zstandard block compression failed");
        stored.resize(16 + compressed);
        require(stored.size() <= UINT32_MAX,
                "compressed logical block exceeds V10 uint32 size limit");
        return EncodedBlockFile(span.number, std::move(stored));
    }

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
        raw.bytes(payload_header.data(), payload_header.size());
        if (bitmap) {
            uint64_t byte_index = 0;
            uint8_t word = 0;
            scan([&](const BlockCell &cell) {
                uint64_t position = uint64_t(cell.y - ymin) * width + cell.x - xmin;
                uint64_t target = position / 8;
                while (byte_index < target) { raw.byte(word); word = 0; ++byte_index; }
                word |= static_cast<uint8_t>(1u << (position % 8));
            });
            while (byte_index < bitmap_bytes) { raw.byte(word); word = 0; ++byte_index; }
        } else {
            uint64_t prior = 0, index = 0;
            scan([&](const BlockCell &cell) {
                uint64_t position = uint64_t(cell.y - ymin) * width + cell.x - xmin;
                raw.var(index ? position - prior : position);
                prior = position;
                ++index;
            });
        }
        if (all_same) {
            raw.var(common);
        } else {
            scan([&](const BlockCell &cell) { raw.var(cell.count); });
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

Bytes make_header(const hic10::Header &header) {
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

struct WriterState {
    hic10::Header header;
    VectorManifest vectors;
};

WriterState prepare_writer_state(const StageManifest &stage, const BuildManifest &build,
                                 const WriteOptions &options) {
    WriterState state;
    auto &header = state.header;
    header.genome = options.genome;
    header.attributes = {{"software", "hictools-c hic_v10_large"},
                         {"hicFileScaling", "1.0"}};
    for (const auto &chromosome : stage.chromosomes)
        header.chromosomes.push_back({chromosome.name, chromosome.length, {}});
    for (uint32_t resolution : build.resolutions)
        header.resolutions[0].push_back({resolution});
    if (!options.vector_manifest.empty()) {
        state.vectors = read_vector_manifest(options.vector_manifest);
        require(state.vectors.source_fingerprint == stage.source_fingerprint,
                "vector manifest refers to a different source");
        header.norms = state.vectors.norms;
    }
    auto derive = [&](uint32_t target_bin, uint32_t source_bin) {
        uint32_t target = header.resolution(0, target_bin);
        uint32_t source = header.resolution(0, source_bin);
        require(target_bin > source_bin && target_bin % source_bin == 0,
                "derived target must be an integer multiple of source");
        require(!hic10::required_materialized_resolution(target_bin),
                "500 kb must remain materialized");
        auto &resolution = header.resolutions[0][target];
        require(!resolution.mode || resolution.source == source,
                "conflicting derived source");
        resolution.mode = 1;
        resolution.source = source;
        resolution.aggregation = 1;
    };
    for (const auto &d : options.derived) derive(d.first, d.second);
    for (const auto &resolution : header.resolutions[0])
        if (uint32_t source = hic10::required_derived_source(resolution.bin))
            derive(resolution.bin, source);
    for (const auto &resolution : header.resolutions[0])
        if (resolution.mode)
            require(!header.resolutions[0][resolution.source].mode,
                    "chained derivation is forbidden");
    require(hic10::required_bp_resolution_policy(header.resolutions[0]),
            "mandatory V10 resolution policy is not satisfied");
    return state;
}

std::pair<uint64_t, uint64_t> derived_statistics(const RunInfo &source,
                                                 uint32_t target_resolution,
                                                 uint64_t columns, uint64_t rows,
                                                 bool cis) {
    MatrixStatistics statistics(columns, rows, cis);
    RollupAccumulator rollup(source.resolution, target_resolution);
    RunReader reader(source);
    CellRecord cell;
    while (reader.next(cell))
        rollup.accept(cell, [&](const CellRecord &rolled) { statistics.add(rolled); });
    rollup.finish([&](const CellRecord &rolled) { statistics.add(rolled); });
    return statistics.result();
}

struct WrittenBlocks {
    uint32_t count = 0;
    uint64_t position = 0, length = 0;
};

WrittenBlocks write_block_build(Output &output, BlockBuild build,
                                const std::string &temporary_directory,
                                const WriteOptions &options, ThreadPool *compression_pool) {
    WrittenBlocks result;
    if (!build.statistics.first) return result;
    BlockRun blocks = bounded_block_merge(std::move(build.runs), temporary_directory,
                                          options.merge_fan_in);
    RemovePath remove_blocks(blocks.path);
    const std::string index_path = blocks.path + ".index-tmp-" +
                                   std::to_string(getpid());
    RemovePath remove_index(index_path);
    File index_file(index_path, "w+b");
    {
        BlockRunReader reader(blocks);
        BlockCell cell;
        std::deque<std::future<EncodedBlockFile>> pending;
        BlockSpan span;
        uint64_t record_index = 0;
        auto store = [&](EncodedBlockFile encoded) {
            uint64_t position = encoded.data.empty()
                ? output.copy(encoded.path, encoded.bytes)
                : output.write(encoded.data);
            require(result.count < UINT32_MAX, "too many V10 blocks in one resolution");
            ++result.count;
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
                if (pending.size() >= static_cast<size_t>(options.threads)) consume();
            } else {
                store(encode_block_span(blocks, span, options.compression_level));
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
    result.length = uint64_t(24) + uint64_t(result.count) * 16;
    hic10::put(header, result.length, 8);
    hic10::put(header, result.count, 4); hic10::put(header, 0, 4);
    result.position = output.write(header);
    std::vector<unsigned char> copy(8 * 1024 * 1024);
    uint64_t remaining = uint64_t(result.count) * 16;
    while (remaining) {
        size_t n = static_cast<size_t>(std::min<uint64_t>(copy.size(), remaining));
        index_file.read(copy.data(), n);
        output.write(copy.data(), n);
        remaining -= n;
    }
    index_file.close();
    return result;
}

std::string pair_part_path(const std::string &directory, const PairInfo &pair) {
    std::ostringstream name;
    name << "pair-" << std::setfill('0') << std::setw(5) << pair.chr1 << '-'
         << std::setw(5) << pair.chr2 << ".v10.hic";
    return join_path(directory, name.str());
}

void validate_fragment_header(const hic10::Header &header,
                              const StageManifest &stage,
                              const BuildManifest &build) {
    require(header.chromosomes.size() == stage.chromosomes.size(),
            "pair fragment chromosome table differs from stage manifest");
    for (size_t i = 0; i < stage.chromosomes.size(); ++i)
        require(header.chromosomes[i].name == stage.chromosomes[i].name &&
                    header.chromosomes[i].length == stage.chromosomes[i].length,
                "pair fragment chromosome table differs from stage manifest");
    require(header.resolutions[1].empty() &&
                header.resolutions[0].size() == build.resolutions.size(),
            "pair fragment resolutions differ from build manifest");
    for (size_t i = 0; i < build.resolutions.size(); ++i)
        require(header.resolutions[0][i].bin == build.resolutions[i],
                "pair fragment resolutions differ from build manifest");
    require(header.norms.empty(),
            "pair fragments must omit vectors; supply --vectors to merge-pairs");
}

Bytes canonical_fragment_header(hic10::Reader &reader) {
    Bytes header = reader.read_bytes(0, reader.header_length());
    require(header.size() >= 32, "truncated pair fragment header");
    std::fill(header.begin() + 16, header.begin() + 32, 0);
    return header;
}

void copy_relocated_matrix(const std::string &path, hic10::Reader &reader,
                           Output &output, uint64_t old_begin, uint64_t old_end,
                           uint64_t new_begin) {
    require(old_begin <= old_end && old_end <= reader.file_size(),
            "invalid pair fragment matrix section");
    std::vector<uint64_t> fields = reader.matrix_relocation_fields();
    for (uint64_t field : fields)
        require(field >= old_begin && field <= old_end - 8,
                "pair fragment relocation field lies outside matrix section");

    File input(path, "rb");
    input.seek(old_begin);
    std::vector<unsigned char> buffer(8 * 1024 * 1024);
    uint64_t position = old_begin;
    size_t field_index = 0;
    while (position < old_end) {
        uint64_t length = std::min<uint64_t>(buffer.size(), old_end - position);
        // A buffer must not end in the middle of a relocation field. One buffer
        // usually consumes many fields, so the one that can straddle the end is
        // the last field starting before it, not the first unconsumed one.
        {
            const auto first = fields.begin() + static_cast<std::ptrdiff_t>(field_index);
            const auto bound = std::lower_bound(first, fields.end(), position + length);
            if (bound != first && *(bound - 1) + 8 > position + length)
                length = *(bound - 1) - position;
        }
        require(length, "cannot align pair fragment relocation buffer");
        input.read(buffer.data(), static_cast<size_t>(length));
        const uint64_t end = position + length;
        while (field_index < fields.size() && fields[field_index] + 8 <= end) {
            const uint64_t field = fields[field_index++];
            require(field >= position, "unordered pair fragment relocation fields");
            unsigned char *word = buffer.data() + (field - position);
            const uint64_t old_value = get_u64(word);
            require(old_value >= old_begin && old_value < old_end,
                    "pair fragment pointer lies outside matrix section");
            require(old_value - old_begin <= UINT64_MAX - new_begin,
                    "relocated V10 pointer overflow");
            put_u64(word, new_begin + (old_value - old_begin));
        }
        output.write(buffer.data(), static_cast<size_t>(length));
        position = end;
    }
    require(field_index == fields.size(), "unprocessed pair fragment relocation fields");
}

void write_footer(Output &output, std::vector<MatrixEntry> matrices) {
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
    require(options.merge_fan_in >= 2, "writer merge fan-in must be at least two");
    require(options.resolution_batch > 0, "resolution batch must be positive");
    require(options.pair_index == std::numeric_limits<size_t>::max() ||
                options.pair_index < stage.pairs.size(),
            "writer pair index is outside stage manifest");
    require(options.compression_level >= ZSTD_minCLevel() &&
                options.compression_level <= ZSTD_maxCLevel(), "invalid Zstandard level");
    std::string temporary_directory = options.temporary_directory.empty()
        ? join_path(build_path.substr(0, build_path.find_last_of('/')), "block-work")
        : options.temporary_directory;
    if (temporary_directory == "block-work") temporary_directory = "./block-work";
    make_directory(temporary_directory);

    WriterState state = prepare_writer_state(stage, build, options);
    hic10::Header &header = state.header;

    Output output(output_path, make_header(header));
    std::unique_ptr<ThreadPool> compression_pool(
        options.threads > 1 ? new ThreadPool(static_cast<size_t>(options.threads)) : nullptr);
    std::vector<MatrixEntry> matrices;
    for (size_t pair_index = 0; pair_index < stage.pairs.size(); ++pair_index) {
        if (options.pair_index != std::numeric_limits<size_t>::max() &&
            pair_index != options.pair_index)
            continue;
        const PairInfo &pair = stage.pairs[pair_index];
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
        const bool rotated = a == b;
        auto publish = [&](uint32_t ri, std::pair<uint64_t, uint64_t> statistics,
                           BlockBuild block_build) {
            const auto &resolution = header.resolutions[0][ri];
            uint64_t column_bins = header.bins(a, 0, ri), row_bins = header.bins(b, 0, ri);
            uint32_t block_bins = choose_block_bins(column_bins, row_bins, resolution.bin,
                                                    options.block_bins, rotated);
            uint32_t columns = hic10::narrow(hic10::ceil_div(column_bins, block_bins));
            WrittenBlocks written;
            if (!resolution.mode)
                written = write_block_build(output, std::move(block_build),
                                            temporary_directory, options,
                                            compression_pool.get());
            Bytes descriptor;
            hic10::put(descriptor, 0, 1); hic10::put(descriptor, resolution.mode, 1);
            hic10::put(descriptor, resolution.aggregation, 1); hic10::put(descriptor, 0, 1);
            hic10::put(descriptor, ri, 4); hic10::put(descriptor, resolution.bin, 4);
            hic10::put(descriptor, resolution.source, 4); hic10::put(descriptor, rotated, 1);
            hic10::put(descriptor, 0, 3); hic10::put(descriptor, statistics.second, 8);
            hic10::put(descriptor, statistics.first, 8);
            hic10::put(descriptor, 0x7fc00000, 4); hic10::put(descriptor, 0x7fc00000, 4);
            hic10::put(descriptor, block_bins, 4); hic10::put(descriptor, columns, 4);
            hic10::put(descriptor, written.position, 8); hic10::put(descriptor, written.length, 8);
            hic10::put(descriptor, written.count, 4); hic10::put(descriptor, 0, 4);
            output.patch(meta_position + 24 + uint64_t(ri) * 76, descriptor);
        };

        if (!build.root_only) {
            for (uint32_t ri = 0; ri < header.resolutions[0].size(); ++ri) {
                const auto &resolution = header.resolutions[0][ri];
                uint64_t column_bins = header.bins(a, 0, ri), row_bins = header.bins(b, 0, ri);
                if (resolution.mode) {
                    const auto &source = header.resolutions[0][resolution.source];
                    auto found = build.cells.find(std::make_tuple(a, b, source.bin));
                    require(found != build.cells.end(),
                            "build manifest is missing a derived source resolution");
                    publish(ri, derived_statistics(found->second, resolution.bin,
                                                    column_bins, row_bins, rotated), {});
                } else {
                    auto found = build.cells.find(std::make_tuple(a, b, resolution.bin));
                    require(found != build.cells.end(),
                            "build manifest is missing a matrix resolution");
                    uint32_t block_bins = choose_block_bins(column_bins, row_bins,
                                                            resolution.bin,
                                                            options.block_bins, rotated);
                    uint32_t columns = hic10::narrow(
                        hic10::ceil_div(column_bins, block_bins));
                    BlockBuild prepared = create_block_runs(
                        found->second, temporary_directory, options.memory_bytes,
                        block_bins, columns, column_bins, row_bins, rotated);
                    auto statistics = prepared.statistics;
                    publish(ri, statistics, std::move(prepared));
                }
            }
            continue;
        }

        auto root_found = build.cells.find(
            std::make_tuple(a, b, build.resolutions.front()));
        require(root_found != build.cells.end(),
                "root-only build is missing a chromosome pair");
        const RunInfo &root = root_found->second;
        std::vector<uint32_t> materialized, derived;
        for (uint32_t ri = 0; ri < header.resolutions[0].size(); ++ri)
            (header.resolutions[0][ri].mode ? derived : materialized).push_back(ri);

        struct BlockPipeline {
            uint32_t ri = 0;
            std::unique_ptr<RollupAccumulator> rollup;
            std::unique_ptr<BlockRunBuilder> builder;
        };
        struct StatisticPipeline {
            uint32_t ri = 0;
            std::unique_ptr<RollupAccumulator> rollup;
            std::unique_ptr<MatrixStatistics> statistics;
        };
        for (size_t begin = 0; begin < materialized.size(); begin += options.resolution_batch) {
            const size_t end = std::min(materialized.size(), begin + options.resolution_batch);
            const uint64_t per_resolution_memory = options.memory_bytes / (end - begin);
            std::vector<BlockPipeline> blocks;
            for (size_t at = begin; at < end; ++at) {
                uint32_t ri = materialized[at];
                const auto &resolution = header.resolutions[0][ri];
                uint64_t column_bins = header.bins(a, 0, ri), row_bins = header.bins(b, 0, ri);
                uint32_t block_bins = choose_block_bins(column_bins, row_bins,
                                                        resolution.bin,
                                                        options.block_bins, rotated);
                uint32_t columns = hic10::narrow(hic10::ceil_div(column_bins, block_bins));
                RunInfo metadata = root;
                metadata.resolution = resolution.bin;
                metadata.records = metadata.checksum = 0;
                BlockPipeline pipeline;
                pipeline.ri = ri;
                if (resolution.bin != root.resolution)
                    pipeline.rollup.reset(new RollupAccumulator(root.resolution,
                                                                resolution.bin));
                pipeline.builder.reset(new BlockRunBuilder(
                    metadata, temporary_directory, per_resolution_memory,
                    block_bins, columns, column_bins, row_bins, rotated));
                blocks.push_back(std::move(pipeline));
            }
            std::vector<StatisticPipeline> statistics;
            if (begin == 0) {
                for (uint32_t ri : derived) {
                    const auto &resolution = header.resolutions[0][ri];
                    StatisticPipeline pipeline;
                    pipeline.ri = ri;
                    pipeline.rollup.reset(new RollupAccumulator(root.resolution,
                                                                resolution.bin));
                    pipeline.statistics.reset(new MatrixStatistics(
                        header.bins(a, 0, ri), header.bins(b, 0, ri), rotated));
                    statistics.push_back(std::move(pipeline));
                }
            }
            RunReader input(root);
            CellRecord cell;
            while (input.next(cell)) {
                for (auto &pipeline : blocks) {
                    if (pipeline.rollup)
                        pipeline.rollup->accept(cell, [&](const CellRecord &rolled) {
                            pipeline.builder->add(rolled);
                        });
                    else
                        pipeline.builder->add(cell);
                }
                for (auto &pipeline : statistics)
                    pipeline.rollup->accept(cell, [&](const CellRecord &rolled) {
                        pipeline.statistics->add(rolled);
                    });
            }
            for (auto &pipeline : blocks) {
                if (pipeline.rollup)
                    pipeline.rollup->finish([&](const CellRecord &rolled) {
                        pipeline.builder->add(rolled);
                    });
                BlockBuild prepared = pipeline.builder->finish();
                auto result = prepared.statistics;
                publish(pipeline.ri, result, std::move(prepared));
            }
            for (auto &pipeline : statistics) {
                pipeline.rollup->finish([&](const CellRecord &rolled) {
                    pipeline.statistics->add(rolled);
                });
                publish(pipeline.ri, pipeline.statistics->result(), {});
            }
        }
    }
    if (!options.vector_manifest.empty())
        write_vectors(output, header, state.vectors, options.compression_level,
                      compression_pool.get(), options.threads);
    write_footer(output, std::move(matrices));
    output.publish();
    std::cerr << "Large-data V10 matrix assembly complete: " << output_path << "\n";
}

void print_write_tasks(const std::string &stage_path, const std::string &build_path,
                       const std::string &parts_directory, const std::string &output) {
    StageManifest stage = read_stage_manifest(stage_path);
    BuildManifest build = read_build_manifest(build_path, false);
    require(build.source_fingerprint == stage.source_fingerprint,
            "stage and build manifests refer to different sources");
    for (size_t pair = 0; pair < stage.pairs.size(); ++pair)
        std::cout << "write-pair\t" << pair << '\t'
                  << pair_part_path(parts_directory, stage.pairs[pair]) << "\n";
    std::cout << "merge-pairs\t0\t" << output << "\n";
}

void write_pair_v10(const std::string &stage_path, const std::string &build_path,
                    const std::string &parts_directory, size_t pair_index,
                    const WriteOptions &options) {
    require(options.vector_manifest.empty(),
            "write-pair does not duplicate normalization vectors; pass them to merge-pairs");
    StageManifest stage = read_stage_manifest(stage_path);
    BuildManifest build = read_build_manifest(build_path, false);
    require(build.source_fingerprint == stage.source_fingerprint,
            "stage and build manifests refer to different sources");
    require(pair_index < stage.pairs.size(), "writer pair index is outside stage manifest");
    make_directory(parts_directory);
    const PairInfo &pair = stage.pairs[pair_index];
    const std::string output = pair_part_path(parts_directory, pair);
    WriteOptions pair_options = options;
    pair_options.pair_index = pair_index;
    if (path_exists(output)) {
        hic10::Reader reader(output);
        WriterState expected = prepare_writer_state(stage, build, pair_options);
        Bytes expected_header = make_header(expected.header);
        Bytes actual_header = canonical_fragment_header(reader);
        require(actual_header == expected_header,
                "existing pair fragment was produced with different writer options "
                "(actual header " + std::to_string(actual_header.size()) + " bytes/" +
                hex64(fnv1a(actual_header.data(), actual_header.size())) +
                ", expected " + std::to_string(expected_header.size()) + " bytes/" +
                hex64(fnv1a(expected_header.data(), expected_header.size())) + ")");
        require(reader.matrices().size() == 1 &&
                    reader.matrices()[0].chr1 == pair.chr1 &&
                    reader.matrices()[0].chr2 == pair.chr2,
                "existing pair fragment contains the wrong matrix");
        for (const auto &locator : reader.vector_indexes())
            require(!locator.position && !locator.length,
                    "existing pair fragment contains normalization vectors");
        require(reader.footer().position + reader.footer().length == reader.file_size(),
                "existing pair fragment has unexpected trailing data");
        reader.matrix_relocation_fields();
        std::cerr << "Pair writer task " << pair_index << " already complete\n";
        return;
    }
    write_v10(stage_path, build_path, output, pair_options);
}

void merge_pair_v10(const std::string &stage_path, const std::string &build_path,
                    const std::string &parts_directory, const std::string &output_path,
                    const WriteOptions &options) {
    StageManifest stage = read_stage_manifest(stage_path);
    BuildManifest build = read_build_manifest(build_path, false);
    require(build.source_fingerprint == stage.source_fingerprint,
            "stage and build manifests refer to different sources");
    require(!stage.pairs.empty(), "stage manifest has no chromosome pairs");
    require(options.threads > 0, "writer thread count must be positive");
    require(options.compression_level >= ZSTD_minCLevel() &&
                options.compression_level <= ZSTD_maxCLevel(), "invalid Zstandard level");
    require(options.derived.empty(),
            "merge-pairs inherits resolution derivations from pair fragments");

    for (const PairInfo &pair : stage.pairs)
        require(output_path != pair_part_path(parts_directory, pair),
                "merge output must not overwrite a pair fragment");

    const std::string first_path = pair_part_path(parts_directory, stage.pairs.front());
    hic10::Reader first(first_path);
    validate_fragment_header(first.header(), stage, build);
    Bytes fragment_header = canonical_fragment_header(first);
    hic10::Header final_header = first.header();
    VectorManifest vectors;
    if (!options.vector_manifest.empty()) {
        vectors = read_vector_manifest(options.vector_manifest);
        require(vectors.source_fingerprint == stage.source_fingerprint,
                "vector manifest refers to a different source");
        final_header.norms = vectors.norms;
    }

    Output output(output_path, make_header(final_header));
    std::vector<MatrixEntry> matrices;
    matrices.reserve(stage.pairs.size());
    for (size_t pair_index = 0; pair_index < stage.pairs.size(); ++pair_index) {
        const PairInfo &expected = stage.pairs[pair_index];
        const std::string path = pair_part_path(parts_directory, expected);
        hic10::Reader fragment(path);
        validate_fragment_header(fragment.header(), stage, build);
        require(fragment.header_length() == fragment_header.size() &&
                    canonical_fragment_header(fragment) == fragment_header,
                "pair fragments have different headers");
        require(fragment.matrices().size() == 1 &&
                    fragment.matrices()[0].chr1 == expected.chr1 &&
                    fragment.matrices()[0].chr2 == expected.chr2,
                "pair fragment contains the wrong matrix: " + path);
        for (const auto &locator : fragment.vector_indexes())
            require(!locator.position && !locator.length,
                    "pair fragment contains normalization vectors: " + path);
        require(fragment.footer().position + fragment.footer().length == fragment.file_size(),
                "pair fragment has unexpected trailing data: " + path);
        const uint64_t old_begin = fragment.header_length();
        const uint64_t old_end = fragment.footer().position;
        const uint64_t new_begin = output.position();
        hic10::MatrixKey key{expected.chr1, expected.chr2};
        hic10::FileLocator matrix = fragment.matrix_location(key);
        require(matrix.position >= old_begin && matrix.position <= old_end &&
                    matrix.length <= old_end - matrix.position,
                "pair fragment matrix descriptor lies outside matrix section");
        copy_relocated_matrix(path, fragment, output, old_begin, old_end, new_begin);
        matrices.push_back({expected.chr1, expected.chr2,
                            new_begin + (matrix.position - old_begin), matrix.length});
        std::cerr << "Merged V10 matrix " << stage.chromosomes[expected.chr1].name
                  << " x " << stage.chromosomes[expected.chr2].name << "\n";
    }

    std::unique_ptr<ThreadPool> compression_pool(
        options.threads > 1 ? new ThreadPool(static_cast<size_t>(options.threads)) : nullptr);
    if (!options.vector_manifest.empty())
        write_vectors(output, final_header, vectors, options.compression_level,
                      compression_pool.get(), options.threads);
    write_footer(output, std::move(matrices));
    output.publish();
    std::cerr << "Merged pair fragments into V10 file: " << output_path << "\n";
}

} // namespace hic10large
