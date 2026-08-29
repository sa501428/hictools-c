#include "hic_writer.h"
#include "little_endian.h"
#include <algorithm>
#include <stdexcept>
#include <cstring>

StagedBlockFile::~StagedBlockFile() {
    if (!path.empty()) std::remove(path.c_str());
}

namespace {
uint64_t block_size(const CompressedBlock& block) {
    return block.staged_size ? block.staged_size : block.bytes.size();
}
void copy_staged(FILE* source, FILE* destination, uint64_t offset, uint64_t length,
                 std::vector<uint8_t>& buffer) {
    if (fseeko(source, static_cast<off_t>(offset), SEEK_SET) != 0)
        throw std::runtime_error("Cannot seek staged V9 block section");
    while (length) {
        size_t wanted = static_cast<size_t>(std::min<uint64_t>(buffer.size(), length));
        if (std::fread(buffer.data(), 1, wanted, source) != wanted)
            throw std::runtime_error("Cannot read staged V9 block section");
        if (std::fwrite(buffer.data(), 1, wanted, destination) != wanted)
            throw std::runtime_error("Cannot append staged V9 block section");
        length -= wanted;
    }
}
} // namespace

// ============================================================
//  Constructor / Destructor
// ============================================================

HicWriter::HicWriter(const std::string& path) : path_(path) {
    file_ = fopen(path.c_str(), "wb");
    if (!file_) throw std::runtime_error("Cannot open output file: " + path);
}

HicWriter::~HicWriter() {
    if (file_) { fclose(file_); file_ = nullptr; }
}

int64_t HicWriter::current_position() const {
    return (int64_t)ftell(file_);
}

// ============================================================
//  Header
// ============================================================

void HicWriter::write_header(
    const Genome& genome,
    const std::vector<int>& resolutions,
    const std::unordered_map<std::string, std::string>& attributes)
{
    FILE* f = file_;

    // Magic: "HIC\0"
    fwrite("HIC", 1, 3, f);
    fwrite_byte(f, 0);

    // Version = 9
    fwrite_int32(f, HIC_VERSION);

    // Placeholder for master index position (back-patched in write_footer)
    master_index_pos_position_ = current_position();
    fwrite_int64(f, 0LL);

    // Genome ID: we use an empty string here; caller can pass via attributes
    // Actually the Java code writes genomeId directly. We'll write a placeholder.
    // The genome ID is the first attribute "genome" or we can derive from filename.
    // For now, write an empty string (user can specify via --genome flag stored in attributes).
    std::string genome_id;
    auto git = attributes.find("genomeID");
    if (git != attributes.end()) genome_id = git->second;
    fwrite_string(f, genome_id);

    // normVectorIndexPosition placeholder
    norm_vi_pos_position_ = current_position();
    fwrite_int64(f, 0LL);

    // normVectorIndexLength placeholder
    norm_vi_len_position_ = current_position();
    fwrite_int64(f, 0LL);

    // Attribute dictionary (exclude genomeID, it's written separately above)
    // Standard attributes expected by readers:
    //   "software", "statistics", "graphs", "nviIndex", "nviLength", "hicFileScaling", "v9BaseDepth"
    std::vector<std::pair<std::string,std::string>> attrs;
    for (auto& [k, v] : attributes) {
        if (k == "genomeID") continue;
        attrs.push_back({k, v});
    }
    // Always include software tag
    bool has_software = false;
    for (auto& a : attrs) if (a.first == "software") { has_software = true; break; }
    if (!has_software) attrs.push_back({"software", "hic_pre C++ v1.0"});
    // Always include hicFileScaling
    bool has_scaling = false;
    for (auto& a : attrs) if (a.first == "hicFileScaling") { has_scaling = true; break; }
    if (!has_scaling) attrs.push_back({"hicFileScaling", "1.0"});

    fwrite_int32(f, (int32_t)attrs.size());
    for (auto& [k, v] : attrs) {
        fwrite_string(f, k);
        fwrite_string(f, v);
    }

    // Chromosome list (all chromosomes including "All" at index 0)
    const auto& chroms = genome.chromosomes();
    fwrite_int32(f, (int32_t)chroms.size());
    for (auto& c : chroms) {
        fwrite_string(f, c.name);
        fwrite_int64(f, c.length);
    }

    // BP resolutions (coarsest first, then finer)
    fwrite_int32(f, (int32_t)resolutions.size());
    for (int r : resolutions) fwrite_int32(f, r);

    // Fragment resolutions = 0 (BP only)
    fwrite_int32(f, 0);
}

// ============================================================
//  Matrix body
// ============================================================

void HicWriter::write_matrix(
    int chr1_idx, int chr2_idx,
    const std::vector<ZoomWriteData>& zooms)
{
    FILE* f = file_;

    int64_t matrix_start = current_position();

    fwrite_int32(f, chr1_idx);
    fwrite_int32(f, chr2_idx);

    // Count non-empty zoom levels
    int n_res = (int)zooms.size();
    fwrite_int32(f, n_res);

    // For each zoom level, write the metadata + block index placeholder,
    // then remember where the block index is so we can back-patch it.

    // We store: (position of block index in file, reference to zoom blocks)
    struct ZoomInfo {
        int64_t block_index_pos; // file position of first block index entry
        int     n_blocks;
        const ZoomWriteData* zoom;
    };
    std::vector<ZoomInfo> zoom_infos;
    zoom_infos.reserve(zooms.size());

    for (size_t zoom_index = 0; zoom_index < zooms.size(); ++zoom_index) {
        const auto& z = zooms[zoom_index];
        // unit = "BP"
        fwrite_string(f, "BP");
        fwrite_int32(f, static_cast<int32_t>(zoom_index));
        fwrite_float(f, z.sum_counts);
        fwrite_float(f, 0.0f);    // occupiedCellCount (unused)
        fwrite_float(f, 0.0f);    // stdDev (unused)
        fwrite_float(f, 0.0f);    // percent95 (unused)
        fwrite_int32(f, z.bin_size);
        fwrite_int32(f, z.block_bin_count);
        fwrite_int32(f, z.block_column_count);

        // Block count
        fwrite_int32(f, (int32_t)z.blocks.size());

        // Block index: write placeholders (blockNumber=num, position=0, size=0)
        int64_t idx_pos = current_position();
        for (const auto& block : z.blocks) {
            fwrite_int32(f, block.block_number);
            fwrite_int64(f, 0LL);  // position placeholder
            fwrite_int32(f, 0);    // size placeholder
        }
        zoom_infos.push_back({idx_pos, (int)z.blocks.size(), &z});
    }

    // Record matrix metadata size (up to here)
    int64_t matrix_meta_end = current_position();
    int64_t meta_size = matrix_meta_end - matrix_start;

    // Build master index key: "chr1Idx_chr2Idx"
    std::string key = std::to_string(chr1_idx) + "_" + std::to_string(chr2_idx);
    master_index_.push_back({key, matrix_start, (int32_t)meta_size});

    // Now write the actual compressed block data for each zoom level,
    // back-patching the block index entries.
    for (size_t zi = 0; zi < zooms.size(); zi++) {
        auto& info = zoom_infos[zi];
        auto& z    = zooms[zi];

        int64_t idx_ptr = info.block_index_pos;
        std::unique_ptr<FILE, decltype(&fclose)> staged(nullptr, &fclose);
        if (info.zoom->staged_file) {
            staged.reset(fopen(info.zoom->staged_file->path.c_str(), "rb"));
            if (!staged)
                throw std::runtime_error("Cannot open staged V9 block section: " +
                                         info.zoom->staged_file->path);
        }
        std::vector<uint8_t> copy_buffer(staged ? 1024 * 1024 : 0);
        for (const auto& block : z.blocks) {
            int64_t block_start = current_position();
            uint64_t block_length = block_size(block);
            if (block_length > INT32_MAX)
                throw std::runtime_error("Compressed V9 block exceeds signed 32-bit index size");
            size_t block_bytes = static_cast<size_t>(block_length);
            if (staged) {
                copy_staged(staged.get(), f, block.staged_offset, block_bytes, copy_buffer);
            } else if (fwrite(block.bytes.data(), 1, block_bytes, f) != block_bytes) {
                throw std::runtime_error("Error writing compressed block");
            }

            // Back-patch this block's index entry
            // Entry layout: int32 blockNum + int64 position + int32 size = 16 bytes
            long saved = ftell(f);
            fseek(f, (long)idx_ptr, SEEK_SET);
            fwrite_int32(f, block.block_number);
            fwrite_int64(f, block_start);
            fwrite_int32(f, (int32_t)block_bytes);
            fseek(f, saved, SEEK_SET);

            idx_ptr += 16; // advance to next block index entry
        }
    }
}

// ============================================================
//  Footer
// ============================================================

void HicWriter::write_footer(
    const std::vector<ExpectedValueCalculation*>& expected_vals,
    const std::vector<NormVector>& norm_vecs,
    const std::vector<NormExpectedValueVector>& norm_expected)
{
    FILE* f = file_;

    // The master index position is here
    int64_t footer_start = current_position();

    // Back-patch the master index position in the header
    patch_int64(f, master_index_pos_position_, footer_start);

    // ---- nBytesV5 placeholder (long) ----
    // We'll compute how many bytes the "V5" portion takes (master index + expected values)
    // then back-patch it.
    int64_t n_bytes_v5_pos = current_position();
    fwrite_int64(f, 0LL);  // placeholder

    // ---- Master index ----
    fwrite_int32(f, (int32_t)master_index_.size());
    for (auto& entry : master_index_) {
        fwrite_string(f, entry.key);
        fwrite_int64(f, entry.position);
        fwrite_int32(f, entry.size);
    }

    // ---- Expected value vectors (non-normalized) ----
    // Only include resolutions that have data.
    std::vector<ExpectedValueCalculation*> non_empty_ev;
    for (auto* ev : expected_vals) {
        if (ev && ev->has_data()) {
            const_cast<ExpectedValueCalculation*>(ev)->compute_density();
            non_empty_ev.push_back(ev);
        }
    }

    fwrite_int32(f, (int32_t)non_empty_ev.size());
    for (auto* ev : non_empty_ev) {
        fwrite_string(f, "BP");
        fwrite_int32(f, ev->bin_size());
        int64_t n_vals = ev->n_values();
        fwrite_int64(f, n_vals);
        for (float v : ev->density()) fwrite_float(f, v);
        // chrScaleFactors
        auto& sf = ev->chr_scale_factors();
        fwrite_int32(f, (int32_t)sf.size());
        for (auto& [ci, factor] : sf) {
            fwrite_int32(f, ci);
            fwrite_float(f, factor);
        }
    }

    // Compute nBytesV5 = bytes from start of footer to here (before norm expected vecs)
    int64_t v5_end = current_position();
    int64_t n_bytes_v5 = v5_end - (n_bytes_v5_pos + 8); // exclude the nBytesV5 field itself
    patch_int64(f, n_bytes_v5_pos, n_bytes_v5);

    // ---- Normalized expected value vectors ----
    fwrite_int32(f, (int32_t)norm_expected.size());
    for (auto& nev : norm_expected) {
        fwrite_string(f, nev.norm_type);
        fwrite_string(f, nev.unit);
        fwrite_int32(f, nev.bin_size);
        fwrite_int64(f, (int64_t)nev.values.size());
        for (float v : nev.values) fwrite_float(f, v);
        fwrite_int32(f, (int32_t)nev.chr_scale.size());
        for (auto& [ci, sf] : nev.chr_scale) {
            fwrite_int32(f, ci);
            fwrite_float(f, sf);
        }
    }

    // ---- Normalization vector index ----
    // We write all norm vectors contiguously after the index.
    // First write the index with placeholder positions, then write the data,
    // then back-patch positions.

    int64_t norm_vi_start = current_position();

    // Back-patch normVectorIndexPosition
    patch_int64(f, norm_vi_pos_position_, norm_vi_start);

    fwrite_int32(f, (int32_t)norm_vecs.size());

    // Write index entries with placeholder positions
    struct NVIEntry { int64_t pos_placeholder; };
    std::vector<int64_t> nv_pos_placeholders;
    nv_pos_placeholders.reserve(norm_vecs.size());

    for (auto& nv : norm_vecs) {
        fwrite_string(f, nv.norm_type);
        fwrite_int32(f, nv.chr_idx);
        fwrite_string(f, "BP");
        fwrite_int32(f, nv.bin_size);
        nv_pos_placeholders.push_back(current_position());
        fwrite_int64(f, 0LL);  // position placeholder
        int64_t n_bytes = (int64_t)nv.values.size() * 4 + 8; // nValues(8) + values*4
        fwrite_int64(f, n_bytes);
    }

    int64_t norm_vi_end = current_position();
    int64_t norm_vi_len = norm_vi_end - norm_vi_start;
    patch_int64(f, norm_vi_len_position_, norm_vi_len);

    // ---- Normalization vector arrays ----
    for (size_t i = 0; i < norm_vecs.size(); i++) {
        auto& nv = norm_vecs[i];
        int64_t nv_start = current_position();

        // Back-patch position in index
        patch_int64(f, nv_pos_placeholders[i], nv_start);

        fwrite_int64(f, (int64_t)nv.values.size());
        for (float v : nv.values) fwrite_float(f, v);
    }
}
