#pragma once

#include "sort.h"
#include <cstdint>
#include <string>
#include <vector>

namespace hic10large {

struct ScaleOptions {
    double tolerance = 1.0e-4;
    double row_sum_tolerance = 5.0e-2;
    double delta = 5.0e-2;
    int max_iterations = 2000;
    int threads = 1;
    uint64_t sort_memory_bytes = 8ULL * 1024 * 1024 * 1024;
    size_t merge_fan_in = 128;
    std::string temporary_directory;
};

struct ScaleResult {
    bool success = false;
    int iterations = 0;
    std::string reason;
    std::vector<float> values;
};

ScaleResult scale_cis(const RunInfo &cells, uint32_t bins, const ScaleOptions &options);

} // namespace hic10large
