#pragma once
// Local V10 reader with full-matrix and block-streaming interfaces. Intentionally
// separate from both straw and the V9 hic_addnorm reader.
#include "format.h"
#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace hic10 {
struct FileLocator {
    uint64_t position = 0, length = 0;
};
struct MatrixKey {
    uint32_t chr1, chr2;
    bool operator<(const MatrixKey &o) const {
        return std::tie(chr1, chr2) < std::tie(o.chr1, o.chr2);
    }
};

struct VectorEntry {
    uint8_t kind = 0, unit = 0;
    uint32_t norm = 0, chr = 0, ri = 0;
    Bytes bytes;
    size_t descriptors = 0;
    uint32_t chunks = 0;
};

class Reader {
  public:
    explicit Reader(const std::string &path);
    ~Reader();
    Reader(const Reader &) = delete;
    Reader &operator=(const Reader &) = delete;

    const Header &header() const { return header_; }
    const std::vector<MatrixKey> &matrices() const { return matrix_keys_; }
    FileLocator matrix_location(MatrixKey key) const;
    uint64_t header_length() const { return header_length_; }
    uint64_t file_size() const { return file_size_; }
    const std::array<FileLocator, 3> &vector_indexes() const { return vector_indexes_; }
    const FileLocator &footer() const { return footer_; }
    Bytes read_bytes(uint64_t position, uint64_t length);

    // First byte occupied by an active vector chunk or vector index. Everything
    // before this point contains the header and matrix section and can be copied
    // when rebuilding the vector section in place.
    uint64_t vector_data_start();

    // Absolute u64 fields inside the matrix section that point elsewhere in
    // that section. Optionally returns the referenced matrix storage intervals.
    std::vector<uint64_t> matrix_relocation_fields(std::vector<FileLocator> *storage = nullptr);
    std::vector<VectorEntry> vector_entries();
    // Decode one materialized block at a time, without retaining a full matrix.
    void stream_materialized(uint32_t chr1, uint32_t chr2, uint8_t unit, uint32_t ri,
                             const std::function<void(const Cell &)> &emit);

    // Returns canonical cells (x <= y for cis). Derived resolutions are summed
    // exactly from their declared materialized source before conversion to float.
    Matrix matrix(uint32_t chr1, uint32_t chr2, uint8_t unit, uint32_t resolution_index);

  private:
    struct Zoom {
        uint8_t unit = 0, mode = 0, aggregation = 0, type = 0, grid = 0;
        uint32_t resolution = 0, bin = 0, source = UINT32_MAX, block_bins = 0, columns = 0;
        uint64_t occupied = 0;
        FileLocator block_index;
        uint64_t block_index_position_field = 0;
        uint32_t blocks = 0;
    };
    struct MatrixMeta {
        std::vector<Zoom> zooms;
    };

    FILE *file_ = nullptr;
    uint64_t file_size_ = 0, header_length_ = 0;
    Header header_;
    FileLocator footer_;
    std::array<FileLocator, 3> vector_indexes_{};
    std::map<MatrixKey, FileLocator> matrix_locations_;
    std::map<MatrixKey, MatrixMeta> matrix_metadata_;
    std::vector<MatrixKey> matrix_keys_;

    const MatrixMeta &metadata(MatrixKey key);
    Matrix materialized(MatrixKey key, const Zoom &zoom,
                        const std::function<void(const Cell &)> &emit = {});
};
} // namespace hic10
