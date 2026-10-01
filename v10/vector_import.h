#pragma once
#include <string>
namespace hic10 {
// Import BP vectors produced by the large normalizer or another algorithm using
// its documented vectors.manifest/H10W sidecar format. Existing entries win.
void add_vector_manifest_v10(const std::string &path, const std::string &manifest, int level);
}
