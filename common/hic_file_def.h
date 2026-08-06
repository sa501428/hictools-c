#pragma once
// HiC V9 file format constants and core data structures.
// All multi-byte values are little-endian.

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

// ---- Format constants ----

static constexpr int HIC_VERSION = 9;
static constexpr int BLOCK_CAPACITY = 1000;   // max blocks in memory before spilling
static constexpr int V9_DEPTH_BASE  = 2;      // base for intra-chr block depth calculation

// Default BP resolutions (coarsest first, matching Java defaults)
static const int DEFAULT_BP_RESOLUTIONS[] = {
    2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000
};
static constexpr int N_DEFAULT_BP_RESOLUTIONS = 10;

// Block column count sizing thresholds (from MatrixPP.java)
static constexpr int INTRA_CUTOFF = 500;   // bin size below which nColumns is reduced for intra
static constexpr int INTER_CUTOFF = 5000;  // same for inter

// ---- Data structures ----

struct Chromosome {
    std::string name;
    int32_t     index;   // 0-based index in chromosome array
    int64_t     length;  // bp
};

// A single contact record (binned)
struct ContactRecord {
    int32_t bin1;
    int32_t bin2;
    float   count;
};

// Master index entry: maps "chr1_chr2" key → file position
struct MatrixIndexEntry {
    std::string key;
    int64_t     position;
    int32_t     size;     // size of the matrix METADATA section (not blocks)
};

// Block index entry within a zoom level
struct BlockIndexEntry {
    int32_t block_number;
    int64_t position;
    int32_t size_bytes;
};

// Metadata for a single zoom level within a matrix
struct ZoomMeta {
    int32_t bin_size;
    int32_t block_bin_count;    // blockSize: dimension along diagonal (intra) or side (inter)
    int32_t block_column_count; // number of block columns
    float   sum_counts;
    std::vector<BlockIndexEntry> block_index;
    int64_t block_index_file_pos; // file position where block index starts (for back-patching)
};

// Expected value vector for one resolution
struct ExpectedValueVector {
    std::string unit;     // "BP"
    int32_t     bin_size;
    std::vector<float> values;                        // indexed by genomic distance in bins
    std::unordered_map<int32_t, float> chr_scale;    // chr_index → scale factor
};

// Normalization type strings (as stored in file)
static constexpr const char* NORM_VC       = "VC";
static constexpr const char* NORM_VC_SQRT  = "VC_SQRT";
static constexpr const char* NORM_SCALE    = "SCALE";

// Normalization vector index entry
struct NormVectorIndexEntry {
    std::string norm_type;   // e.g. "VC", "SCALE"
    int32_t     chr_idx;
    std::string unit;        // "BP"
    int32_t     bin_size;
    int64_t     position;    // file position of the norm array
    int64_t     n_bytes;
};

// A normalization vector (one per chromosome per norm type per resolution)
struct NormVector {
    std::string norm_type;
    int32_t     chr_idx;
    int32_t     bin_size;
    std::vector<float> values;  // NaN means unmappable bin
};

// Normalized expected value entry
struct NormExpectedValueVector {
    std::string norm_type;
    std::string unit;
    int32_t     bin_size;
    std::vector<float> values;
    std::unordered_map<int32_t, float> chr_scale;
};
