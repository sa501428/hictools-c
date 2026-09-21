#pragma once

#include "scale.h"
#include <string>

namespace hic10large {

struct NormalizeOptions {
    bool vc = true;
    bool vc_sqrt = true;
    bool scale = true;
    uint64_t memory_bytes = 8ULL * 1024 * 1024 * 1024;
    std::string temporary_directory;
    ScaleOptions scale_options;
};

void normalize_cells(const std::string &stage_manifest, const std::string &build_manifest,
                     const std::string &output_directory, const NormalizeOptions &options);
void normalize_chromosome(const std::string &stage_manifest, const std::string &build_manifest,
                          const std::string &output_directory, uint32_t chromosome,
                          const NormalizeOptions &options);
void expected_resolution(const std::string &stage_manifest, const std::string &build_manifest,
                         const std::string &output_directory, uint32_t resolution_index,
                         const NormalizeOptions &options);
void finalize_vectors(const std::string &stage_manifest, const std::string &build_manifest,
                      const std::string &output_directory, const NormalizeOptions &options);
void print_normalize_tasks(const std::string &stage_manifest, const std::string &build_manifest,
                           const std::string &output_directory, const NormalizeOptions &options);

} // namespace hic10large
