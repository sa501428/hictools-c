#pragma once
// HiC V9 binary format writer.
//
// Usage pattern:
//   HicWriter w("out.hic");
//   w.write_header(genome, resolutions, attributes);
//   // For each chr-pair matrix (in sorted order):
//   w.begin_matrix(chr1_idx, chr2_idx);
//   for each resolution:
//       w.write_zoom(zoom_data);  // writes block index placeholder + blocks
//   w.end_matrix();              // patches block index positions
//   // After all matrices:
//   w.write_footer(expected_values, norm_vectors, norm_expected_values);
//   // Patches the master index position in the header.

#include "hic_file_def.h"
#include "expected_value.h"
#include "block_writer.h"
#include "genome.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdio>
#include <cstdint>
#include <functional>
#include <memory>

// One worker-owned file containing compressed V9 blocks. The ordered writer
// streams ranges from it and deletion happens automatically after the matrix is
// appended (or during exception unwinding).
struct StagedBlockFile {
    std::string path;
    explicit StagedBlockFile(std::string p) : path(std::move(p)) {}
    ~StagedBlockFile();
};

struct CompressedBlock {
    int block_number;
    std::vector<uint8_t> bytes;
    uint64_t staged_offset = 0;
    uint32_t staged_size = 0;
};

// Per-zoom-level data ready to write. Compression is completed by preprocessing
// workers so the single output writer only performs ordered I/O and back-patching.
struct ZoomWriteData {
    int  bin_size;
    int  block_bin_count;
    int  block_column_count;
    float sum_counts;
    std::shared_ptr<StagedBlockFile> staged_file;
    // Blocks sorted by block_number
    std::vector<CompressedBlock> blocks;
};

class HicWriter {
public:
    explicit HicWriter(const std::string& path);
    ~HicWriter();

    // ---- Header ----
    // resolutions: list of BP bin sizes, coarsest first.
    // attributes:  key-value pairs stored in header (e.g. "software" → "hic_pre v1.0").
    void write_header(
        const Genome& genome,
        const std::vector<int>& resolutions,
        const std::unordered_map<std::string, std::string>& attributes
    );

    // ---- Matrix body ----
    // Write the metadata section for a chromosome pair then its block data.
    // Call for each (chr1, chr2) pair in sorted order (chr1 <= chr2).
    void write_matrix(
        int chr1_idx, int chr2_idx,
        const std::vector<ZoomWriteData>& zooms
    );

    // ---- Footer ----
    // Call after all matrices have been written.
    // expected_vals: one per resolution (from ExpectedValueCalculation).
    // norm_vecs: normalization vectors (VC, VC_SQRT, SCALE) keyed by (norm_type, chr_idx, bin_size).
    // norm_expected: normalized expected values.
    void write_footer(
        const std::vector<ExpectedValueCalculation*>& expected_vals,
        const std::vector<NormVector>& norm_vecs,
        const std::vector<NormExpectedValueVector>& norm_expected
    );

    int64_t current_position() const;

private:
    FILE*   file_ = nullptr;
    std::string path_;

    // Positions needing back-patching
    int64_t master_index_pos_position_ = 0;   // where we wrote the placeholder
    int64_t norm_vi_pos_position_      = 0;   // normVectorIndexPosition placeholder
    int64_t norm_vi_len_position_      = 0;   // normVectorIndexLength placeholder

    // Master index: maps chr-pair key → (position, size)
    std::vector<MatrixIndexEntry> master_index_;
};
