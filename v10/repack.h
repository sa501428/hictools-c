#pragma once
// Atomic addnorm needs to expand the variable header when it introduces the
// first normalization names. This helper copies the matrix section byte for
// byte and relocates only its absolute file-position fields.
#include "reader.h"
#include <cstdio>

namespace hic10 {
// Writes a freshly serialized header followed by the unchanged matrix section.
// The returned footer has its matrix-record positions relocated but is not yet
// written; the caller writes vector data/indexes before appending it.
Bytes repack_matrix_prefix(FILE *output, Reader &reader, const Header &header);
} // namespace hic10
