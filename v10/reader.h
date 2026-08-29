#pragma once
// Local, full-matrix V10 reader used by hic_v10 addnorm. This is intentionally
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

class Reader {
  public:
    explicit Reader(const std::string &path);
    ~Reader();
    Reader(const Reader &) = delete;
    Reader &operator=(const Reader &) = delete;

    const Header &header() const { return header_; }
    const std::vector<MatrixKey> &matrices() const { return matrix_keys_; }
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
    // that section. A header repacker shifts these fields by the header delta.
    std::vector<uint64_t> matrix_relocation_fields();

    // Returns canonical cells (x <= y for cis). Derived resolutions are summed
    // exactly from their declared materialized source before conversion to float.
    Matrix matrix(uint32_t chr1, uint32_t chr2, uint8_t unit, uint32_t resolution_index);

  private:
    struct Zoom {
        uint8_t unit = 0, mode = 0, aggregation = 0, type = 0, grid = 0;
        uint32_t resolution = 0, bin = 0, source = UINT32_MAX, block_bins = 0, columns = 0;
        uint64_t occupied = 0;
        FileLocator page_index;
        uint64_t page_index_position_field = 0;
        uint32_t pages = 0, blocks = 0;
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
    Matrix materialized(MatrixKey key, const Zoom &zoom);
};
} // namespace hic10
