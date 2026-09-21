#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace hic10large {

struct WriteOptions {
    std::string genome;
    std::string temporary_directory;
    std::string vector_manifest;
    uint64_t memory_bytes = 8ULL * 1024 * 1024 * 1024;
    uint32_t block_bins = 256;
    int compression_level = 6;
    int threads = 4;
    size_t merge_fan_in = 128;
    size_t resolution_batch = 4;
    size_t pair_index = std::numeric_limits<size_t>::max();
    std::vector<std::pair<uint32_t, uint32_t>> derived;
};

void write_v10(const std::string &stage_manifest, const std::string &build_manifest,
               const std::string &output, const WriteOptions &options);

void print_write_tasks(const std::string &stage_manifest,
                       const std::string &build_manifest,
                       const std::string &parts_directory,
                       const std::string &output);
void write_pair_v10(const std::string &stage_manifest,
                    const std::string &build_manifest,
                    const std::string &parts_directory, size_t pair_index,
                    const WriteOptions &options);
void merge_pair_v10(const std::string &stage_manifest,
                    const std::string &build_manifest,
                    const std::string &parts_directory,
                    const std::string &output, const WriteOptions &options);

} // namespace hic10large
