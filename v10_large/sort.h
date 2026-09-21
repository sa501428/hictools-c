#pragma once

#include "io.h"
#include <cstdint>
#include <functional>
#include <string>
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
