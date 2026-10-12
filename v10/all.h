#pragma once
#include "reader.h"

namespace hic10 {
// Scaled overview convention uses ordinary V10 header attributes and matrices.
struct AllOptions {
    uint32_t bins = 500;
    uint32_t source = 2500000;
};
bool is_all(const std::string &name);
class AllBuilder {
  public:
    AllBuilder(Header &header, const AllOptions &options = {});
    uint32_t index() const { return index_; }
    uint32_t source() const { return source_; }
    void add(uint32_t a, uint32_t b, uint64_t x, uint64_t y, uint32_t input_bin,
             uint64_t value, bool scores = false);
    void collect(Reader &reader);
    Matrix matrix(uint8_t unit, uint32_t ri) const;
    void write(Writer &writer) const;
  private:
    Header header_;
    uint32_t index_, source_, finest_;
    uint64_t group_;
    std::vector<uint64_t> offsets_;
    std::map<std::pair<uint32_t, uint32_t>, uint64_t> counts_;
    std::map<std::pair<uint32_t, uint32_t>, double> scores_;
    bool scored_ = false;
};
// Existing All is preserved. An absent All is appended without renumbering IDs.
void add_all(const std::string &input, const std::string &output,
             const AllOptions &options = {}, int level = 6, uint32_t threads = 4);
} // namespace hic10
