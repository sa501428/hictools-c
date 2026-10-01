#pragma once
#include "normalize.h"
namespace hic10large {
void add_text_norm_file(const std::string &path, const std::string &text,
                        const NormalizeOptions &options, int level);
void add_norm_file(const std::string &path, const NormalizeOptions &options, int level);
}
