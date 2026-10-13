#pragma once
// Atomic addnorm needs to expand the variable header when it introduces the
// normalization names. This helper copies referenced matrix storage byte for
// byte and relocates only its absolute file-position fields.
#include "reader.h"
#include <cstdio>

class ThreadPool;

namespace hic10 {
// Writes a freshly serialized header followed by all referenced matrix storage.
// Copies block bytes exactly and relocates metadata/index pointers; works with
// interleaved matrix/vector layouts and omits obsolete unreferenced data.
// The returned footer has its matrix-record positions relocated but is not yet
// written; the caller writes vector data/indexes before appending it.
// A supplied pool copies matrix ranges concurrently to disjoint output offsets.
Bytes repack_matrix_prefix(FILE *output, Reader &reader, const Header &header,
                           ThreadPool *pool = nullptr, Reader *extra = nullptr);
// Large-file repair: preserve the source body at a constant offset, stream
// pointer updates, append optional new matrices, validate, then rename atomically.
// No per-block relocation map or storage sorting; memory is independent of the
// total stored block count. Unreferenced source bytes are intentionally retained.
void append_matrix_atomic(Reader &reader, const Header &header, const std::string &path,
                          uint32_t threads, Reader *extra,
                          const std::function<void(const std::string &)> &validate);
} // namespace hic10
