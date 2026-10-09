#pragma once
#include "normalize.h"
#include "sort.h"
#include "v10/reader.h"
namespace hic10large {
RunInfo extract_addnorm_cis_run(hic10::Reader &reader, uint32_t chr, uint8_t unit, uint32_t ri,
                                const std::string &directory, const NormalizeOptions &options);
void add_text_norm_file(const std::string &path, const std::string &text,
                        const NormalizeOptions &options, int level);
void add_norm_file(const std::string &path, const NormalizeOptions &options, int level);
}
