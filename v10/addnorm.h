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
// Atomically adds missing type/resolution bundles. Existing dictionary IDs,
// vectors, expected values, and matrix block bytes are preserved.
void add_norm_v10(const std::string &path, const AddNormOptions &options);
} // namespace hic10
