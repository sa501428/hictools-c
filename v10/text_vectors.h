#pragma once
#include "reader.h"

namespace hic10 {
using ExpectedCellConsumer = std::function<void(uint32_t, uint32_t, long double)>;
using ExpectedCellStream = std::function<void(Reader &, uint32_t, uint8_t, uint32_t,
                                             const ExpectedCellConsumer &)>;
// HIC_NORM_VECTORS 1 text import: exact supplied divisors plus computed NEVI.
// The default stream uses full chromosome matrices. The large executable
// supplies a disk-backed derived-resolution stream through the callback.
void add_text_vectors_v10(const std::string &path, const std::string &text_path,
                         int level, const std::string &temporary_directory = "/tmp",
                         const ExpectedCellStream &stream = {});
} // namespace hic10
