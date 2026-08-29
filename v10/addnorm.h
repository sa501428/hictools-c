#pragma once
#include <string>
namespace hic10 {
struct AddNormOptions {
    bool vc = true, vc_sqrt = true, scale = true;
    int threads = 4;
    double tolerance = 1.0e-4;
    int max_iterations = 2000;
    int minimum_scale_resolution = 0;
    int compression_level = 3;
};
// Atomically replaces the V10 file with an equivalent file containing fresh
// EVI0, NVI0, and NEVI indexes. Matrix pages and metadata are copied unchanged.
void add_norm_v10(const std::string &path, const AddNormOptions &options);
} // namespace hic10
