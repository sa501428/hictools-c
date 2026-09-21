#pragma once

#include <cstdint>
#include <string>

namespace hic10large {

struct StageOptions {
    uint64_t chunk_records = 250000000;
    uint64_t progress_records = 100000000;
};

void stage_hbs(const std::string &input, const std::string &output_directory,
               const StageOptions &options);
void inspect_stage(const std::string &manifest_path, bool verify_shards);

} // namespace hic10large
