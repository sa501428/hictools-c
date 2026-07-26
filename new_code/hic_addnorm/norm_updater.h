#pragma once
// hic_addnorm: reads an existing V9 .hic file and appends/updates
// normalization vectors in-place.
//
// Strategy:
//   1. Parse the .hic header → genome, resolutions, master index
//   2. For each resolution, for each chromosome (intra-only for SCALE):
//      a. Read the raw contact matrix from the .hic blocks
//      b. Compute VC (row sums)
//      c. Compute VC_SQRT (sqrt of VC)
//      d. Compute SCALE (run balance() algorithm)
//   3. Compute normalized expected value vectors
//   4. Rewrite the footer section of the .hic file with the new norms
//      (preserves all existing data in the body)

#include "../common/hic_file_def.h"
#include "../common/genome.h"
#include "scale_norm.h"
#include <string>
#include <vector>
#include <unordered_map>

struct AddNormOptions {
    bool   build_vc        = true;
    bool   build_vc_sqrt   = true;
    bool   build_scale     = true;
    int    num_threads     = 4;
    double scale_tolerance = 1.0e-4;
    int    scale_max_iter  = 2000;
    // Minimum resolution to compute SCALE for (below this resolution,
    // data is too sparse; skip SCALE to save time)
    int    min_scale_res   = 0; // 0 = compute for all
};

// Main entry point
void add_norm(const std::string& hic_path, const AddNormOptions& opts);

// ---- Lightweight .hic reader (just enough for addnorm) ----
// No dependency on straw or libcurl; reads local files only.

struct HicHeader {
    int                     version;
    int64_t                 footer_position;
    std::string             genome_id;
    int64_t                 norm_vi_position;          // value stored in header
    int64_t                 norm_vi_length;            // value stored in header
    int64_t                 norm_vi_position_offset;   // file offset of norm_vi_position field
    int64_t                 norm_vi_length_offset;     // file offset of norm_vi_length field
    std::vector<Chromosome> chromosomes;
    std::vector<int>        bp_resolutions;
    // Attributes
    std::unordered_map<std::string, std::string> attributes;
};

struct HicZoomMeta {
    std::string unit;
    int         bin_size;
    float       sum_counts;
    int         block_bin_count;
    int         block_column_count;
    // Block index for random access
    std::vector<BlockIndexEntry> block_index;
    int64_t     block_index_file_pos;
};

struct HicMatrixMeta {
    int chr1_idx;
    int chr2_idx;
    std::vector<HicZoomMeta> zooms;
    int64_t                  file_position; // position of this matrix metadata
    int32_t                  meta_size;
};

// Parse the header from a .hic file (read-only).
HicHeader read_hic_header(FILE* f);

// Parse master index + matrix metadata from footer.
std::vector<HicMatrixMeta> read_matrix_meta(FILE* f, const HicHeader& header);

// Read one block (decompress from zlib) → contact records.
std::vector<ContactRecord> read_block(FILE* f, const BlockIndexEntry& entry);

// Read all intra-chromosomal contacts for (chr, resolution) as sparse arrays.
// Returns (row[], col[], val[]) for the upper triangle, and n_bins (chr_len/res+1).
struct SparseMatrix {
    std::vector<uint32_t> row;
    std::vector<uint32_t> col;
    std::vector<float>    val;
    uint32_t              n_bins;
};

SparseMatrix read_intra_contacts(
    FILE* f,
    const HicMatrixMeta& meta,
    int bin_size,
    uint32_t expected_n_bins
);
