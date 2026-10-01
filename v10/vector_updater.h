#pragma once
#include "reader.h"

namespace hic10 {
// Transactional, additive vector writer shared by the regular and large tools.
// Existing keys win; normalization IDs are never renumbered. New algorithms
// supply exact float32 words through Vector::values or a chunk loader.
class VectorOutput {
  public:
    VectorOutput(Reader &reader, const Header &header, const std::string &path, int level);
    ~VectorOutput();
    bool contains(uint8_t kind, uint32_t norm, uint32_t chr, uint8_t unit, uint32_t ri) const;
    void add(const Vector &vector);
    void finish(const Header &header);
  private:
    std::unique_ptr<class VectorOutputImpl> impl_;
};
} // namespace hic10
