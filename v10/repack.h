#pragma once
// Atomic addnorm needs to expand the variable header when it introduces the
// normalization names. This helper copies referenced matrix storage byte for
// byte and relocates only its absolute file-position fields.
#include "reader.h"
#include <cstdio>

namespace hic10 {
// Writes a freshly serialized header followed by all referenced matrix storage.
// Copies block bytes exactly and relocates metadata/index pointers; works with
// interleaved matrix/vector layouts and omits obsolete unreferenced data.
// The returned footer has its matrix-record positions relocated but is not yet
// written; the caller writes vector data/indexes before appending it.
Bytes repack_matrix_prefix(FILE *output, Reader &reader, const Header &header);
} // namespace hic10
