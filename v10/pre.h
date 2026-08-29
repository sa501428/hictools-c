#pragma once
#include "common/pair_parser.h"
#include "format.h"
namespace hic10 {
struct PreOptions {
    std::vector<uint32_t> resolutions;
    std::string genome, tmpDir = "/tmp";
    InputFormat format = InputFormat::AUTO;
    uint32_t readAhead = 0; // 0 uses the writer thread count
    int mapq = 0;
    bool intra = false, nearDiagonal = false;
};
void pre(const std::string &input, const std::string &output, const std::string &genomeSpec,
         const Options &options, const PreOptions &preOptions);
} // namespace hic10
