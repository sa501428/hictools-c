#pragma once
// Block number calculation (V9 format) and zlib block compression.
//
// V9 introduces:
//   - Intra-chromosomal: rotated block structure along the diagonal.
//   - Inter-chromosomal: standard row-major grid.
//
// Block binary format (after decompression):
//   int32 nRecords
//   int32 binColumnOffset
//   int32 binRowOffset
//   byte  useFloatContact  (1 = float, 0 = short)
//   byte  useIntXPos       (1 = int,   0 = short)
//   byte  useIntYPos       (1 = int,   0 = short)
//   byte  matrixRepresentation (1 = list-of-rows)
//   [list-of-rows data]

#include "hic_file_def.h"
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <cmath>
#include <string>
#include <cstdio>

// --------------------------------------------------------------------------
// V9 block number calculation
// --------------------------------------------------------------------------

// V9 depth for intra-chromosomal blocks.
// depth = floor( log2(1 + |binX - binY| / sqrt(2) / blockBinCount) )
// This mirrors the Java V9Depth class.
inline int v9_depth(int bin_x, int bin_y, int block_bin_count, int depth_base = V9_DEPTH_BASE) {
    int dist = std::abs(bin_x - bin_y);
    if (dist == 0) return 0;
    if (depth_base <= 1) depth_base = 2;
    double d = std::log(1.0 + dist / (std::sqrt(2.0) * block_bin_count))
             / std::log((double)depth_base);
    return static_cast<int>(d);
}

// Intra-chromosomal block number (V9)
inline int intra_block_number(int bin_x, int bin_y, int block_bin_count,
                              int block_column_count, int depth_base = V9_DEPTH_BASE) {
    int depth = v9_depth(bin_x, bin_y, block_bin_count, depth_base);
    int pos_along_diag = (bin_x + bin_y) / 2 / block_bin_count;
    return depth * block_column_count + pos_along_diag;
}

// Inter-chromosomal block number (V9 / standard row-major)
inline int inter_block_number(int bin_x, int bin_y, int block_bin_count, int block_column_count) {
    int block_col = bin_x / block_bin_count;
    int block_row = bin_y / block_bin_count;
    return block_column_count * block_row + block_col;
}

// Compute blockBinCount and blockColumnCount for a chromosome pair.
// Mirrors MatrixPP.getNumColumnsFromNumBins() / MatrixZoomDataPP constructor.
struct BlockParams {
    int block_bin_count;    // dimension of block in bins
    int block_column_count; // number of block columns in grid
};

BlockParams compute_block_params(
    int64_t chr1_len, int64_t chr2_len,
    int bin_size, bool is_intra,
    int block_capacity = BLOCK_CAPACITY);

// --------------------------------------------------------------------------
// In-memory block accumulator
// --------------------------------------------------------------------------

// A single block's contact data, collected during the pass over pairs.
struct BlockData {
    int block_number;
    // Map from (bin_x, bin_y) to count using a flat key: (bin_x << 20) | bin_y
    // or just a vector of contacts sorted later.
    std::unordered_map<int64_t, float> contacts;  // key = bin_x | ((int64_t)bin_y << 32)

    BlockData() = default;
    explicit BlockData(int num) : block_number(num) {}

    void add(int bin_x, int bin_y, float count) {
        int64_t key = (int64_t)bin_x | ((int64_t)bin_y << 32);
        contacts[key] += count;
    }

    int num_records() const { return static_cast<int>(contacts.size()); }

    // Merge another block into this one
    void merge(const BlockData& other) {
        for (auto& [k, v] : other.contacts) contacts[k] += v;
    }
};

// --------------------------------------------------------------------------
// Block compression and serialization (V9 "list of rows" format)
// --------------------------------------------------------------------------

// Compress a BlockData into a zlib-deflated byte buffer.
// Returns compressed bytes ready to be written to the .hic file.
// Also fills in the index entry (size is set after compression).
std::vector<uint8_t> compress_block(const BlockData& block);

// Write compressed block bytes directly to a FILE*, return number of bytes written.
size_t write_compressed_block(FILE* f, const BlockData& block);

// --------------------------------------------------------------------------
// Temp-file-backed block storage
// --------------------------------------------------------------------------
// During preprocessing, blocks for a single (chr1, chr2, resolution) combo
// may not fit in memory.  We spill to per-chr-pair temp files.

// Write a list of BlockData to a temp file.
// Format per block: int32 block_num, int32 n_contacts,
//                   (int32 bin_x, int32 bin_y, float count) × n
void write_blocks_to_tempfile(FILE* f, const std::vector<BlockData>& blocks);

// Read blocks from a temp file (appending to existing map).
// Returns number of blocks read.
int merge_blocks_from_tempfile(FILE* f, std::unordered_map<int, BlockData>& block_map);
