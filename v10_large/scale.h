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
    // Optional non-owning inputs. scale_cis consumes them synchronously.
    const std::vector<double> *raw_vc = nullptr;
    const std::vector<float> *warm_norm = nullptr;
    uint32_t warm_resolution = 0;
    double warm_excluded_fraction = 0;
    double warm_vc_exponent = 0.5;
};

struct ScaleResult {
    bool success = false;
    int iterations = 0;
    int cutoff = 1;
    int percentile_cutoff = 1;
    int zscore_cutoff = 1;
    double excluded_fraction = 0;
    bool used_warm_start = false;
    std::string reason;
    std::vector<float> values;
};

ScaleResult scale_cis(const RunInfo &cells, uint32_t bins, const ScaleOptions &options);

} // namespace hic10large
