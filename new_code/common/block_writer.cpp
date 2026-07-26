#include "block_writer.h"
#include "little_endian.h"
#include <zlib.h>
#include <algorithm>
#include <map>
#include <cstring>
#include <stdexcept>
#include <cassert>

// --------------------------------------------------------------------------
// BlockParams computation  (mirrors MatrixPP / MatrixZoomDataPP in Java)
// --------------------------------------------------------------------------

BlockParams compute_block_params(
    int64_t chr1_len, int64_t chr2_len,
    int bin_size, bool is_intra,
    int block_capacity)
{
    // Maximum sqrt for int overflow safety
    const int MAX_SQRT = (int)std::sqrt((double)INT32_MAX);

    int64_t len = std::max(chr1_len, chr2_len);
    int n_bins  = (int)(len / bin_size + 1);

    int cutoff = is_intra ? INTRA_CUTOFF : INTER_CUTOFF;
    int n_cols;
    if (bin_size < cutoff) {
        int64_t num   = (int64_t)n_bins * bin_size;
        int64_t denom = (int64_t)block_capacity * cutoff;
        n_cols = (int)(num / denom) + 1;
    } else {
        n_cols = n_bins / block_capacity + 1;
    }
    n_cols = std::min(n_cols, MAX_SQRT - 1);

    int block_bin_count = n_bins / n_cols + 1;

    BlockParams bp;
    bp.block_column_count = n_cols;
    bp.block_bin_count    = block_bin_count;
    return bp;
}

// --------------------------------------------------------------------------
// Block serialization (V9 "list of rows" format)
// --------------------------------------------------------------------------

// Build the uncompressed byte buffer for a block.
static std::vector<uint8_t> serialize_block(const BlockData& block) {
    // Collect contacts, determine offsets and data type flags
    // Sort by (bin_y, bin_x) for row-major list-of-rows encoding
    struct Record { int32_t bx, by; float v; };
    std::vector<Record> recs;
    recs.reserve(block.contacts.size());
    for (auto& [key, val] : block.contacts) {
        int32_t bx = static_cast<int32_t>(key & 0xFFFFFFFF);
        int32_t by = static_cast<int32_t>((key >> 32) & 0xFFFFFFFF);
        recs.push_back({bx, by, val});
    }
    // Sort rows first, then columns
    std::sort(recs.begin(), recs.end(), [](const Record& a, const Record& b){
        return (a.by != b.by) ? (a.by < b.by) : (a.bx < b.bx);
    });

    if (recs.empty()) {
        // Empty block: write a minimal valid header
        std::vector<uint8_t> buf;
        write_int32(buf, 0);      // nRecords
        write_int32(buf, 0);      // binColumnOffset
        write_int32(buf, 0);      // binRowOffset
        write_byte(buf, 1);       // useFloatContact
        write_byte(buf, 0);       // useIntXPos
        write_byte(buf, 0);       // useIntYPos
        write_byte(buf, 1);       // matrixRepresentation = list of rows
        write_int16(buf, 0);      // rowCount = 0
        return buf;
    }

    // Compute min offsets
    int32_t bin_x_offset = recs[0].bx;
    int32_t bin_y_offset = recs[0].by;
    int32_t bin_x_max    = recs[0].bx;
    int32_t bin_y_max    = recs[0].by;
    for (auto& r : recs) {
        bin_x_offset = std::min(bin_x_offset, r.bx);
        bin_y_offset = std::min(bin_y_offset, r.by);
        bin_x_max    = std::max(bin_x_max, r.bx);
        bin_y_max    = std::max(bin_y_max, r.by);
    }

    // Determine data type flags
    // useFloatContact: true unless all values are integers < Short.MAX_VALUE
    bool all_integer = true;
    float max_count  = 0.0f;
    for (auto& r : recs) {
        if (std::floor(r.v) != r.v) all_integer = false;
        max_count = std::max(max_count, r.v);
    }
    bool use_short_contact = all_integer && (max_count < 32767.0f);
    bool use_short_x = (bin_x_max - bin_x_offset) < 32767;
    bool use_short_y = (bin_y_max - bin_y_offset) < 32767;

    // Build list-of-rows: group contacts by bin_y
    std::map<int32_t, std::vector<std::pair<int32_t,float>>> rows;
    for (auto& r : recs) {
        rows[r.by].emplace_back(r.bx, r.v);
    }

    // Compute sizes for list-of-rows
    int value_size = use_short_contact ? 2 : 4;
    // Each row: rowNumber (2 or 4) + recordCount (2 or 4) + records*(binX(2 or 4) + value)
    size_t lor_size = 0;
    for (auto& [ry, row_recs] : rows) {
        size_t row_num_sz  = use_short_y ? 2 : 4;
        size_t rec_cnt_sz  = use_short_x ? 2 : 4;
        size_t per_record  = (use_short_x ? 2 : 4) + value_size;
        lor_size += row_num_sz + rec_cnt_sz + row_recs.size() * per_record;
    }

    std::vector<uint8_t> buf;
    buf.reserve(16 + lor_size);

    write_int32(buf, (int32_t)recs.size());
    write_int32(buf, bin_x_offset);
    write_int32(buf, bin_y_offset);
    write_byte(buf, use_short_contact ? 0 : 1);   // useFloatContact
    write_byte(buf, use_short_x       ? 0 : 1);   // useIntXPos
    write_byte(buf, use_short_y       ? 0 : 1);   // useIntYPos
    write_byte(buf, 1);                             // matrixRepresentation = list of rows

    // rowCount
    if (use_short_y) write_int16(buf, (int16_t)rows.size());
    else             write_int32(buf, (int32_t)rows.size());

    for (auto& [ry, row_recs] : rows) {
        int32_t rel_y = ry - bin_y_offset;
        int32_t n_rec = (int32_t)row_recs.size();

        // rowNumber
        if (use_short_y) write_int16(buf, (int16_t)rel_y);
        else             write_int32(buf, rel_y);

        // recordCount
        if (use_short_x) write_int16(buf, (int16_t)n_rec);
        else             write_int32(buf, n_rec);

        for (auto& [bx, val] : row_recs) {
            int32_t rel_x = bx - bin_x_offset;
            // binColumn
            if (use_short_x) write_int16(buf, (int16_t)rel_x);
            else             write_int32(buf, rel_x);
            // value
            if (use_short_contact) write_int16(buf, (int16_t)val);
            else                   write_float(buf, val);
        }
    }
    return buf;
}

// --------------------------------------------------------------------------
// zlib compression
// --------------------------------------------------------------------------

static std::vector<uint8_t> zlib_compress(const std::vector<uint8_t>& data) {
    uLong bound = compressBound(data.size());
    std::vector<uint8_t> out(bound);
    uLong out_len = bound;
    int rc = compress2(out.data(), &out_len, data.data(), data.size(), Z_DEFAULT_COMPRESSION);
    if (rc != Z_OK) throw std::runtime_error("zlib compress2 failed: " + std::to_string(rc));
    out.resize(out_len);
    return out;
}

std::vector<uint8_t> compress_block(const BlockData& block) {
    auto raw = serialize_block(block);
    return zlib_compress(raw);
}

size_t write_compressed_block(FILE* f, const BlockData& block) {
    auto compressed = compress_block(block);
    size_t n = fwrite(compressed.data(), 1, compressed.size(), f);
    if (n != compressed.size()) throw std::runtime_error("Error writing compressed block");
    return n;
}

// --------------------------------------------------------------------------
// Temp file I/O
// --------------------------------------------------------------------------

void write_blocks_to_tempfile(FILE* f, const std::vector<BlockData>& blocks) {
    for (auto& block : blocks) {
        int32_t n = (int32_t)block.contacts.size();
        fwrite_int32(f, block.block_number);
        fwrite_int32(f, n);
        for (auto& [key, val] : block.contacts) {
            int32_t bx = static_cast<int32_t>(key & 0xFFFFFFFF);
            int32_t by = static_cast<int32_t>((key >> 32) & 0xFFFFFFFF);
            fwrite_int32(f, bx);
            fwrite_int32(f, by);
            fwrite_float(f, val);
        }
    }
}

int merge_blocks_from_tempfile(FILE* f, std::unordered_map<int, BlockData>& block_map) {
    int count = 0;
    while (!feof(f)) {
        // Attempt to read block_num + n_contacts
        uint8_t b4[4];
        if (fread(b4, 1, 4, f) != 4) break; // EOF
        int32_t block_num = (int32_t)((uint32_t)b4[0] | ((uint32_t)b4[1]<<8) | ((uint32_t)b4[2]<<16) | ((uint32_t)b4[3]<<24));
        if (fread(b4, 1, 4, f) != 4) break;
        int32_t n = (int32_t)((uint32_t)b4[0] | ((uint32_t)b4[1]<<8) | ((uint32_t)b4[2]<<16) | ((uint32_t)b4[3]<<24));

        auto& bd = block_map.try_emplace(block_num, block_num).first->second;
        for (int i = 0; i < n; i++) {
            int32_t bx = fread_int32(f);
            int32_t by = fread_int32(f);
            float   v  = fread_float(f);
            bd.add(bx, by, v);
        }
        count++;
    }
    return count;
}
