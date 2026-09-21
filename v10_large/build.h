#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <vector>
#include "sort.h"

namespace hic10large {

struct BuildOptions {
    std::vector<uint32_t> resolutions;
    uint64_t memory_bytes = 8ULL * 1024 * 1024 * 1024;
    size_t merge_fan_in = 128;
};

void build_cells(const std::string &stage_manifest, const std::string &build_directory,
                 const BuildOptions &options);
void map_root_shard(const std::string &stage_manifest, const std::string &build_directory,
                    size_t shard_index, const BuildOptions &options);
void reduce_root_group(const std::string &stage_manifest, const std::string &build_directory,
                       size_t pair_index, size_t group_index, const BuildOptions &options);
void build_pair_cells(const std::string &stage_manifest, const std::string &build_directory,
                      size_t pair_index, const BuildOptions &options);
void finalize_build(const std::string &stage_manifest, const std::string &build_directory,
                    const BuildOptions &options);
void print_build_tasks(const std::string &stage_manifest, const std::string &build_directory,
                       const BuildOptions &options);

struct BuildManifest {
    std::string stage_manifest;
    uint64_t source_fingerprint = 0;
    std::vector<uint32_t> resolutions;
    std::map<std::tuple<uint32_t, uint32_t, uint32_t>, RunInfo> cells;
};

BuildManifest read_build_manifest(const std::string &path, bool inspect_files = true);

} // namespace hic10large
