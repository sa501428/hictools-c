#pragma once
#include "reader.h"

namespace hic10 {
// Transactional, additive vector writer shared by the regular and large tools.
// Existing keys win; normalization IDs are never renumbered. New algorithms
// supply exact float32 words through Vector::values or a chunk loader.
struct StoredVectorChunk {
    uint64_t begin = 0, offset = 0, checksum = 0;
    uint32_t count = 0, stored = 0;
    uint8_t transform = 0;
    std::string path;
};
class VectorOutput {
  public:
    VectorOutput(Reader &reader, const Header &header, const std::string &path, int level,
                 uint32_t threads = 1);
    ~VectorOutput();
    bool contains(uint8_t kind, uint32_t norm, uint32_t chr, uint8_t unit, uint32_t ri) const;
    void add(const Vector &vector);
    void add_precompressed(const Vector &vector, uint64_t count,
                           const std::vector<StoredVectorChunk> &chunks);
    void finish(const Header &header,
                const std::function<void(const std::string &)> &before_replace = {});

  private:
    std::unique_ptr<class VectorOutputImpl> impl_;
};
} // namespace hic10
