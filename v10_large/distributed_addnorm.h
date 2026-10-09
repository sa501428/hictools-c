#pragma once
#include "normalize.h"
#include <vector>
namespace hic10large {
struct AddnormPlanOptions {
    uint32_t rows_per_task = 1048576;
    uint32_t chunks_per_task = 64;
    uint32_t norms_per_task = 1;
    int level = 3;
};
void plan_addnorm(const std::string &source, const std::vector<std::string> &texts,
                  const std::string &directory, const AddnormPlanOptions &options);
// Phase names: cells, expected, reduce, compress-norm, compress-expected.
// The ordinal is the phase's zero-based task index; the Slurm launcher maps
// compact array indices to these immutable ordinals.
void run_addnorm_task(const std::string &plan, const std::string &phase, uint64_t ordinal,
                      const NormalizeOptions &options);
void print_addnorm_tasks(const std::string &plan, const std::string &phase, bool count_only,
                         bool dependencies = false);
void assemble_addnorm(const std::string &plan, uint32_t threads);
} // namespace hic10large
