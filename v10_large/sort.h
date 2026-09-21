#pragma once

#include "io.h"
#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace hic10large {

struct CellRecord {
    uint32_t x = 0;
    uint32_t y = 0;
    uint64_t count = 0;
};
static_assert(sizeof(CellRecord) == 16, "cell records must remain compact");

constexpr size_t RUN_HEADER_BYTES = 64;
constexpr uint32_t RUN_VERSION = 1;

struct RunInfo {
    std::string path;
    uint32_t chr1 = 0, chr2 = 0, resolution = 0;
    uint64_t records = 0, checksum = 0;
};

std::array<unsigned char, RUN_HEADER_BYTES> encode_run_header(const RunInfo &run);
RunInfo decode_run_header(const unsigned char *bytes, size_t size, const std::string &path);
RunInfo inspect_run(const std::string &path, bool checksum);

void radix_sort_cells(std::vector<CellRecord> &records);
void aggregate_sorted_cells(std::vector<CellRecord> &records);

class RunWriter {
  public:
    RunWriter(const std::string &path, uint32_t chr1, uint32_t chr2, uint32_t resolution);
    ~RunWriter();
    void add(const CellRecord &record);
    RunInfo finish();
  private:
    void flush();
    std::string final_, temporary_;
    RunInfo info_;
    File output_;
    std::vector<CellRecord> buffer_;
    bool finished_ = false;
};

class RunReader {
  public:
    explicit RunReader(const RunInfo &run, size_t buffer_records = 65536);
    bool next(CellRecord &record);
  private:
    RunInfo info_;
    File input_;
    std::vector<CellRecord> buffer_;
    size_t at_ = 0, available_ = 0;
    uint64_t remaining_ = 0;
    uint64_t checksum_ = 14695981039346656037ULL;
    bool verified_ = false;
};

namespace detail {
// Sparse per-row accumulator shared by durable rollups and the fused V10
// writer. Coordinates can span hundreds of millions of 1 bp bins, while one
// active coarse row usually touches only a tiny fraction of those columns.
class SparseRow {
  public:
    void add(uint32_t key, uint64_t value) {
        if (slots_.empty()) rehash(16);
        if ((used_ + 1) * 10 >= slots_.size() * 7) rehash(slots_.size() * 2);
        size_t at = locate(key);
        if (slots_[at].key == EMPTY) {
            slots_[at] = {key, value};
            touched_.push_back(key);
            ++used_;
        } else {
            require(value <= UINT64_MAX - slots_[at].value,
                    "cell count exceeds uint64 during streaming rollup");
            slots_[at].value += value;
        }
    }

    template <class Emit> void flush(Emit &&emit) {
        std::sort(touched_.begin(), touched_.end());
        std::vector<size_t> occupied;
        occupied.reserve(touched_.size());
        for (uint32_t key : touched_) {
            size_t at = locate(key);
            require(slots_[at].key == key && slots_[at].value,
                    "invalid sparse rollup row");
            emit(key, slots_[at].value);
            occupied.push_back(at);
        }
        // Open-addressed probe chains must remain intact until all lookups for
        // this row are complete.
        for (size_t at : occupied) slots_[at] = {};
        touched_.clear();
        used_ = 0;
    }

  private:
    static constexpr uint32_t EMPTY = UINT32_MAX;
    struct Slot {
        uint32_t key = EMPTY;
        uint64_t value = 0;
    };
    size_t locate(uint32_t key) const {
        size_t at = (uint64_t(key) * 11400714819323198485ULL) & (slots_.size() - 1);
        while (slots_[at].key != EMPTY && slots_[at].key != key)
            at = (at + 1) & (slots_.size() - 1);
        return at;
    }
    void rehash(size_t size) {
        std::vector<Slot> old = std::move(slots_);
        slots_.assign(size, {});
        for (const Slot &slot : old) {
            if (slot.key == EMPTY) continue;
            size_t at = locate(slot.key);
            slots_[at] = slot;
        }
    }
    std::vector<Slot> slots_;
    std::vector<uint32_t> touched_;
    size_t used_ = 0;
};
} // namespace detail

class RollupAccumulator {
  public:
    RollupAccumulator(uint32_t source_resolution, uint32_t output_resolution)
        : factor_(output_resolution / source_resolution) {
        require(source_resolution && output_resolution > source_resolution &&
                    output_resolution % source_resolution == 0,
                "invalid rollup resolution");
    }
    template <class Emit> void accept(const CellRecord &record, Emit &&emit) {
        const uint32_t y = record.y / factor_;
        const uint32_t x = record.x / factor_;
        if (!have_row_ || y != active_row_) {
            require(!have_row_ || y > active_row_,
                    "source rows are not ordered during rollup");
            flush(std::forward<Emit>(emit));
            active_row_ = y;
            have_row_ = true;
        }
        row_.add(x, record.count);
    }
    template <class Emit> void finish(Emit &&emit) {
        flush(std::forward<Emit>(emit));
        have_row_ = false;
    }
  private:
    template <class Emit> void flush(Emit &&emit) {
        if (!have_row_) return;
        row_.flush([&](uint32_t x, uint64_t count) {
            emit(CellRecord{x, active_row_, count});
        });
    }
    uint32_t factor_ = 0, active_row_ = 0;
    bool have_row_ = false;
    detail::SparseRow row_;
};

template <class Emit>
void stream_rollup_cells(const RunInfo &source, uint32_t output_resolution, Emit &&emit) {
    require(output_resolution > source.resolution &&
                output_resolution % source.resolution == 0,
            "invalid rollup resolution");
    RollupAccumulator accumulator(source.resolution, output_resolution);
    RunReader reader(source);
    CellRecord record;
    while (reader.next(record)) accumulator.accept(record, emit);
    accumulator.finish(emit);
}

RunInfo merge_runs(const std::vector<RunInfo> &runs, const std::string &output,
                   const std::function<void(const CellRecord &)> &tap = {});
RunInfo bounded_merge(std::vector<RunInfo> runs, const std::string &directory,
                      const std::string &prefix, size_t fan_in,
                      const std::function<void(const CellRecord &)> &tap = {},
                      bool remove_inputs = false);

std::vector<RunInfo> sort_staged_shard(const std::string &shard_path,
                                       const std::string &run_directory,
                                       uint32_t output_resolution,
                                       uint64_t memory_bytes,
                                       const std::string &prefix);

RunInfo rollup_cell_file(const RunInfo &source,
                         const std::string &output_directory,
                         uint32_t output_resolution,
                         const std::string &prefix);

} // namespace hic10large
