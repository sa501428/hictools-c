#include "norm_updater.h"
#include "../common/little_endian.h"
#include "../common/expected_value.h"
#include <zlib.h>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <cstdio>
#include <unistd.h>  // ftruncate(), fileno()

// ============================================================
//  Lightweight .hic reader
// ============================================================

HicHeader read_hic_header(FILE* f) {
    HicHeader h;
    fseek(f, 0, SEEK_SET);

    // Magic
    char magic[4];
    if (fread(magic, 1, 4, f) != 4) throw std::runtime_error("Cannot read magic");
    if (magic[0] != 'H' || magic[1] != 'I' || magic[2] != 'C')
        throw std::runtime_error("Not a HiC file (bad magic)");

    h.version = fread_int32(f);
    if (h.version != 9) throw std::runtime_error("Only HiC V9 supported");

    h.footer_position = fread_int64(f);
    h.genome_id       = fread_string(f);
    h.norm_vi_position_offset = (int64_t)ftell(f);
    h.norm_vi_position = fread_int64(f);
    h.norm_vi_length_offset = (int64_t)ftell(f);
    h.norm_vi_length   = fread_int64(f);

    // Attributes
    int n_attr = fread_int32(f);
    for (int i = 0; i < n_attr; i++) {
        std::string k = fread_string(f);
        std::string v = fread_string(f);
        h.attributes[k] = v;
    }

    // Chromosomes
    int n_chrs = fread_int32(f);
    h.chromosomes.resize(n_chrs);
    for (int i = 0; i < n_chrs; i++) {
        h.chromosomes[i].name   = fread_string(f);
        h.chromosomes[i].length = fread_int64(f);
        h.chromosomes[i].index  = i;
    }

    // BP resolutions
    int n_res = fread_int32(f);
    h.bp_resolutions.resize(n_res);
    for (int i = 0; i < n_res; i++) h.bp_resolutions[i] = fread_int32(f);

    // Fragment resolutions (skip)
    int n_frag = fread_int32(f);
    for (int i = 0; i < n_frag; i++) fread_int32(f);

    return h;
}

std::vector<HicMatrixMeta> read_matrix_meta(FILE* f, const HicHeader& header) {
    // Seek to footer
    fseek(f, (long)header.footer_position, SEEK_SET);

    int64_t n_bytes_v5 = fread_int64(f);

    // Master index
    int n_entries = fread_int32(f);
    std::unordered_map<std::string, std::pair<int64_t,int32_t>> master_index;
    for (int i = 0; i < n_entries; i++) {
        std::string key  = fread_string(f);
        int64_t     pos  = fread_int64(f);
        int32_t     size = fread_int32(f);
        master_index[key] = {pos, size};
    }

    // Read each matrix metadata
    std::vector<HicMatrixMeta> matrices;
    for (auto& [key, ps] : master_index) {
        auto [pos, sz] = ps;
        fseek(f, (long)pos, SEEK_SET);

        HicMatrixMeta meta;
        meta.chr1_idx     = fread_int32(f);
        meta.chr2_idx     = fread_int32(f);
        meta.file_position = pos;
        meta.meta_size    = sz;

        int n_res = fread_int32(f);
        meta.zooms.resize(n_res);
        for (int ri = 0; ri < n_res; ri++) {
            HicZoomMeta& z = meta.zooms[ri];
            z.unit        = fread_string(f);
            /*int res_idx =*/ fread_int32(f);
            z.sum_counts  = fread_float(f);
            /*float occ =*/ fread_float(f);
            /*float std =*/ fread_float(f);
            /*float p95 =*/ fread_float(f);
            z.bin_size            = fread_int32(f);
            z.block_bin_count     = fread_int32(f);
            z.block_column_count  = fread_int32(f);

            int n_blocks = fread_int32(f);
            z.block_index_file_pos = (int64_t)ftell(f);
            z.block_index.resize(n_blocks);
            for (int bi = 0; bi < n_blocks; bi++) {
                z.block_index[bi].block_number = fread_int32(f);
                z.block_index[bi].position     = fread_int64(f);
                z.block_index[bi].size_bytes   = fread_int32(f);
            }
        }
        matrices.push_back(std::move(meta));
    }
    return matrices;
}

// ============================================================
//  Block decompression
// ============================================================

std::vector<ContactRecord> read_block(FILE* f, const BlockIndexEntry& entry) {
    if (entry.size_bytes == 0) return {};

    // Read compressed bytes
    std::vector<uint8_t> compressed(entry.size_bytes);
    fseek(f, (long)entry.position, SEEK_SET);
    if (fread(compressed.data(), 1, entry.size_bytes, f) != (size_t)entry.size_bytes)
        throw std::runtime_error("Failed to read block");

    // Decompress with zlib
    // Try to estimate uncompressed size
    uLong uncomp_size = (uLong)entry.size_bytes * 10;
    std::vector<uint8_t> uncomp(uncomp_size);
    int rc;
    while (true) {
        rc = uncompress(uncomp.data(), &uncomp_size, compressed.data(), (uLong)entry.size_bytes);
        if (rc == Z_BUF_ERROR) {
            uncomp_size *= 2;
            uncomp.resize(uncomp_size);
        } else break;
    }
    if (rc != Z_OK) throw std::runtime_error("zlib decompression failed: " + std::to_string(rc));
    uncomp.resize(uncomp_size);

    // Parse block data
    const uint8_t* p = uncomp.data();
    const uint8_t* end = p + uncomp.size();
    auto read4 = [&]() -> int32_t {
        if (p + 4 > end) throw std::runtime_error("Block parse overflow");
        int32_t v = (int32_t)((uint32_t)p[0] | ((uint32_t)p[1]<<8) | ((uint32_t)p[2]<<16) | ((uint32_t)p[3]<<24));
        p += 4; return v;
    };
    auto read2 = [&]() -> int16_t {
        if (p + 2 > end) throw std::runtime_error("Block parse overflow");
        int16_t v = (int16_t)((uint16_t)p[0] | ((uint16_t)p[1]<<8));
        p += 2; return v;
    };
    auto read1 = [&]() -> uint8_t { return *p++; };
    auto readf = [&]() -> float {
        if (p + 4 > end) throw std::runtime_error("Block parse overflow");
        uint32_t bits = (uint32_t)p[0] | ((uint32_t)p[1]<<8) | ((uint32_t)p[2]<<16) | ((uint32_t)p[3]<<24);
        float v; std::memcpy(&v, &bits, 4); p += 4; return v;
    };

    int32_t n_records       = read4();
    int32_t bin_x_offset    = read4();
    int32_t bin_y_offset    = read4();
    uint8_t use_float       = read1();
    uint8_t use_int_x       = read1();
    uint8_t use_int_y       = read1();
    uint8_t matrix_rep      = read1();

    std::vector<ContactRecord> records;
    records.reserve(n_records);

    auto read_x = [&]() -> int32_t { return use_int_x ? read4() : (int32_t)read2(); };
    auto read_y = [&]() -> int32_t { return use_int_y ? read4() : (int32_t)read2(); };
    auto read_v = [&]() -> float   { return use_float ? readf() : (float)read2(); };

    if (matrix_rep == 1) {
        // List of rows
        int32_t row_count = read_y();
        for (int32_t r = 0; r < row_count; r++) {
            int32_t row_num  = read_y();
            int32_t rec_cnt  = read_x();
            for (int32_t c = 0; c < rec_cnt; c++) {
                int32_t bin_col = read_x();
                float   val     = read_v();
                ContactRecord cr;
                cr.bin1  = bin_x_offset + bin_col;
                cr.bin2  = bin_y_offset + row_num;
                cr.count = val;
                records.push_back(cr);
            }
        }
    } else if (matrix_rep == 2) {
        // Dense
        int32_t n_dense = read4();
        int16_t w       = read2();
        for (int32_t idx = 0; idx < n_dense; idx++) {
            float val = read_v();
            bool missing = use_float ? std::isnan(val) : (val == -32768.0f);
            if (!missing && val != 0.0f) {
                int32_t row = idx / w;
                int32_t col = idx - row * w;
                ContactRecord cr;
                cr.bin1  = bin_x_offset + col;
                cr.bin2  = bin_y_offset + row;
                cr.count = val;
                records.push_back(cr);
            }
        }
    }
    return records;
}

// ============================================================
//  Read all intra-chromosomal contacts for (chr, bin_size)
// ============================================================

SparseMatrix read_intra_contacts(
    FILE* f, const HicMatrixMeta& meta, int bin_size, uint32_t expected_n_bins)
{
    // Find the right zoom level
    const HicZoomMeta* zoom = nullptr;
    for (auto& z : meta.zooms) {
        if (z.bin_size == bin_size && z.unit == "BP") { zoom = &z; break; }
    }
    if (!zoom) throw std::runtime_error("Resolution not found in matrix");

    SparseMatrix sm;
    sm.n_bins = expected_n_bins;

    for (auto& entry : zoom->block_index) {
        auto recs = read_block(f, entry);
        for (auto& cr : recs) {
            if (cr.count <= 0) continue;
            // Ensure upper triangle
            int32_t b1 = cr.bin1, b2 = cr.bin2;
            if (b1 > b2) std::swap(b1, b2);
            sm.row.push_back((uint32_t)b1);
            sm.col.push_back((uint32_t)b2);
            sm.val.push_back(cr.count);
        }
    }
    return sm;
}

// ============================================================
//  VC and VC_SQRT computation
// ============================================================

static std::vector<float> compute_vc(const SparseMatrix& sm) {
    if (sm.n_bins == 0) return {};
    std::vector<float> row_sums(sm.n_bins, 0.0f);
    for (size_t i = 0; i < sm.row.size(); i++) {
        row_sums[sm.row[i]] += sm.val[i];
        if (sm.row[i] != sm.col[i])
            row_sums[sm.col[i]] += sm.val[i];
    }
    return row_sums;
}

static std::vector<float> compute_vc_sqrt(const std::vector<float>& vc) {
    std::vector<float> vs(vc.size());
    for (size_t i = 0; i < vc.size(); i++) {
        vs[i] = std::isnan(vc[i]) ? std::numeric_limits<float>::quiet_NaN() : std::sqrt(vc[i]);
    }
    return vs;
}

// Juicer rescales every saved normalization vector so the sum of the
// normalized matrix equals the sum of the raw matrix.
static void fix_by_sum_factor(const SparseMatrix& sm, std::vector<float>& norm) {
    double raw_sum = 0.0;
    double normalized_sum = 0.0;
    for (size_t p = 0; p < sm.row.size(); ++p) {
        const uint32_t x = sm.row[p], y = sm.col[p];
        const float nx = norm[x], ny = norm[y];
        if (!(nx > 0.0f) || !(ny > 0.0f) ||
            !std::isfinite(nx) || !std::isfinite(ny)) {
            continue;
        }
        const double multiplier = (x == y) ? 1.0 : 2.0;
        raw_sum += multiplier * sm.val[p];
        normalized_sum += multiplier * sm.val[p] / ((double)nx * ny);
    }
    if (raw_sum <= 0.0 || normalized_sum <= 0.0) return;
    const float factor = (float)std::sqrt(normalized_sum / raw_sum);
    for (float& v : norm) if (std::isfinite(v) && v > 0.0f) v *= factor;
}

// ============================================================
//  Compute normalized expected value
// ============================================================

static std::vector<float> compute_norm_expected(
    const SparseMatrix& sm,
    const std::vector<float>& norm_vec,
    int n_bins, int bin_size,
    const Genome& genome, int chr_idx)
{
    // Apply norm vector to contacts, compute distance decay
    std::vector<double> actual(n_bins, 0.0);
    std::vector<double> possible(n_bins, 0.0);
    double chr_len_bins = (double)genome.at(chr_idx).length / bin_size;

    for (long p = 0; p < (long)sm.row.size(); p++) {
        int32_t b1 = sm.row[p], b2 = sm.col[p];
        float n1 = (b1 < (int32_t)norm_vec.size()) ? norm_vec[b1] : std::numeric_limits<float>::quiet_NaN();
        float n2 = (b2 < (int32_t)norm_vec.size()) ? norm_vec[b2] : std::numeric_limits<float>::quiet_NaN();
        if (std::isnan(n1) || std::isnan(n2) || n1 == 0 || n2 == 0) continue;
        int32_t dist = b2 - b1;
        if (dist >= n_bins) continue;
        double normalized = sm.val[p] / (n1 * n2);
        actual[dist] += normalized;
        if (b1 != b2) actual[dist] += normalized; // mirror
    }

    // possibleDistances
    for (int64_t i = 0; i < (int64_t)chr_len_bins && i < n_bins; i++) {
        possible[i] += (chr_len_bins - i);
    }

    std::vector<float> expected(n_bins, 0.0f);
    for (int i = 0; i < n_bins; i++) {
        expected[i] = (possible[i] > 0) ? (float)(actual[i] / possible[i]) : 0.0f;
    }
    return expected;
}

// ============================================================
//  Footer rewriting
// ============================================================

// Read the existing "V5" footer bytes (master index + expected values)
// These need to be preserved exactly.
static std::vector<uint8_t> read_v5_footer(FILE* f, const HicHeader& header) {
    fseek(f, (long)header.footer_position, SEEK_SET);
    int64_t n_bytes_v5 = fread_int64(f);
    std::vector<uint8_t> data((size_t)n_bytes_v5);
    if (fread(data.data(), 1, data.size(), f) != data.size())
        throw std::runtime_error("Cannot read V5 footer");
    return data;
}

// ============================================================
//  Main add_norm function
// ============================================================

void add_norm(const std::string& hic_path, const AddNormOptions& opts) {
    FILE* f = fopen(hic_path.c_str(), "r+b");
    if (!f) throw std::runtime_error("Cannot open hic file: " + hic_path);

    fprintf(stderr, "Reading hic header...\n");
    HicHeader header = read_hic_header(f);
    fprintf(stderr, "  version=%d, genome=%s, %zu chromosomes, %zu resolutions\n",
            header.version, header.genome_id.c_str(),
            header.chromosomes.size(), header.bp_resolutions.size());

    fprintf(stderr, "Reading matrix metadata...\n");
    std::vector<HicMatrixMeta> matrices = read_matrix_meta(f, header);

    // Build a lookup: (chr1, chr2) -> matrix meta
    std::unordered_map<int64_t, const HicMatrixMeta*> meta_map;
    for (auto& m : matrices) {
        int64_t key = (int64_t)m.chr1_idx * 100000 + m.chr2_idx;
        meta_map[key] = &m;
    }

    // Rebuild genome from header chromosomes
    std::vector<std::pair<std::string,int64_t>> chr_list;
    for (auto& c : header.chromosomes) {
        if (c.name != "ALL" && c.name != "All") chr_list.push_back({c.name, c.length});
    }
    Genome genome = Genome::from_list(chr_list);

    // Storage for all norm vectors and norm expected values
    std::vector<NormVector>             all_norm_vecs;
    std::vector<NormExpectedValueVector> all_norm_expected;

    // Process each chromosome at each resolution
    for (int res : header.bp_resolutions) {
        fprintf(stderr, "\nProcessing resolution %d bp\n", res);

        // Expected value accumulator for each norm type at this resolution
        ExpectedValueCalculation ev_vc(genome, res);
        ExpectedValueCalculation ev_vc_sqrt(genome, res);
        ExpectedValueCalculation ev_scale(genome, res);

        // Process each chromosome (intra-chromosomal only for SCALE)
        auto chrs = genome.chromosomes_without_all();
        for (auto* chr : chrs) {
            int ci = chr->index;
            int64_t key = (int64_t)ci * 100000 + ci;
            auto mit = meta_map.find(key);
            if (mit == meta_map.end()) continue;
            const HicMatrixMeta& meta = *mit->second;

            fprintf(stderr, "  chr %s ...", chr->name.c_str());
            fflush(stderr);

            // Read sparse contacts
            SparseMatrix sm;
            try {
                uint32_t expected_n_bins = (uint32_t)(chr->length / res + 1);
                sm = read_intra_contacts(f, meta, res, expected_n_bins);
            } catch (std::exception& e) {
                fprintf(stderr, " (skipped: %s)\n", e.what());
                continue;
            }
            if (sm.row.empty()) { fprintf(stderr, " (empty)\n"); continue; }

            int n_bins = (int)sm.n_bins;

            // ---- VC ----
            std::vector<float> vc_vec;
            if (opts.build_vc || opts.build_vc_sqrt || opts.build_scale) {
                vc_vec = compute_vc(sm);
            }
            if (opts.build_vc) {
                fix_by_sum_factor(sm, vc_vec);
                // Accumulate expected
                for (long p = 0; p < (long)sm.row.size(); p++) {
                    int32_t b1 = sm.row[p], b2 = sm.col[p];
                    float n1 = vc_vec[b1], n2 = vc_vec[b2];
                    if (std::isnan(n1) || std::isnan(n2) || n1 == 0 || n2 == 0) continue;
                    float norm_v = sm.val[p] / (n1 * n2);
                    ev_vc.add_distance(ci, b1, b2, norm_v);
                }
                all_norm_vecs.push_back({NORM_VC, ci, res, vc_vec});
                fprintf(stderr, " VC");
            }

            // ---- VC_SQRT ----
            if (opts.build_vc_sqrt && !vc_vec.empty()) {
                // VC_SQRT starts from unscaled row sums, not the saved/scaled VC.
                auto raw_vc = compute_vc(sm);
                auto vcs_vec = compute_vc_sqrt(raw_vc);
                fix_by_sum_factor(sm, vcs_vec);
                // Accumulate expected
                for (long p = 0; p < (long)sm.row.size(); p++) {
                    int32_t b1 = sm.row[p], b2 = sm.col[p];
                    float n1 = vcs_vec[b1], n2 = vcs_vec[b2];
                    if (std::isnan(n1) || std::isnan(n2) || n1 == 0 || n2 == 0) continue;
                    float norm_v = sm.val[p] / (n1 * n2);
                    ev_vc_sqrt.add_distance(ci, b1, b2, norm_v);
                }
                all_norm_vecs.push_back({NORM_VC_SQRT, ci, res, vcs_vec});
                fprintf(stderr, " VC_SQRT");
            }

            // ---- SCALE ----
            if (opts.build_scale && (opts.min_scale_res == 0 || res >= opts.min_scale_res)) {
                ScaleParams sp;
                sp.tolerance       = opts.scale_tolerance;
                sp.total_max_iter  = opts.scale_max_iter;
                sp.num_threads     = opts.num_threads;

                std::vector<double> scale_b(n_bins, std::numeric_limits<double>::quiet_NaN());
                scale_balance((long)sm.row.size(), sm.row, sm.col, sm.val,
                              (uint32_t)n_bins, scale_b, sp);
                pp_norm_vector((long)sm.row.size(), sm.row, sm.col, sm.val,
                               (uint32_t)n_bins, scale_b, opts.num_threads);

                // Convert to float
                std::vector<float> scale_f(n_bins);
                for (int i = 0; i < n_bins; i++) {
                    scale_f[i] = std::isnan(scale_b[i])
                                 ? std::numeric_limits<float>::quiet_NaN()
                                 : (float)scale_b[i];
                }
                fix_by_sum_factor(sm, scale_f);
                bool has_finite_scale = false;
                for (float v : scale_f) {
                    if (std::isfinite(v) && v > 0.0f) {
                        has_finite_scale = true;
                        break;
                    }
                }
                if (!has_finite_scale) {
                    fprintf(stderr, " SCALE_FAILED");
                    fprintf(stderr, "\n");
                    continue;
                }
                // Accumulate expected
                for (long p = 0; p < (long)sm.row.size(); p++) {
                    int32_t b1 = sm.row[p], b2 = sm.col[p];
                    float n1 = (b1 < n_bins) ? scale_f[b1] : std::numeric_limits<float>::quiet_NaN();
                    float n2 = (b2 < n_bins) ? scale_f[b2] : std::numeric_limits<float>::quiet_NaN();
                    if (std::isnan(n1) || std::isnan(n2) || n1 == 0 || n2 == 0) continue;
                    float norm_v = sm.val[p] / (n1 * n2);
                    ev_scale.add_distance(ci, b1, b2, norm_v);
                }
                all_norm_vecs.push_back({NORM_SCALE, ci, res, scale_f});
                fprintf(stderr, " SCALE");
            }
            fprintf(stderr, "\n");
        }

        // Compute expected value vectors for this resolution
        if (opts.build_vc && ev_vc.has_data()) {
            ev_vc.compute_density();
            NormExpectedValueVector nev;
            nev.norm_type = NORM_VC;
            nev.unit      = "BP";
            nev.bin_size  = res;
            nev.values    = ev_vc.density();
            nev.chr_scale = ev_vc.chr_scale_factors();
            all_norm_expected.push_back(std::move(nev));
        }
        if (opts.build_vc_sqrt && ev_vc_sqrt.has_data()) {
            ev_vc_sqrt.compute_density();
            NormExpectedValueVector nev;
            nev.norm_type = NORM_VC_SQRT;
            nev.unit      = "BP";
            nev.bin_size  = res;
            nev.values    = ev_vc_sqrt.density();
            nev.chr_scale = ev_vc_sqrt.chr_scale_factors();
            all_norm_expected.push_back(std::move(nev));
        }
        if (opts.build_scale && ev_scale.has_data()) {
            ev_scale.compute_density();
            NormExpectedValueVector nev;
            nev.norm_type = NORM_SCALE;
            nev.unit      = "BP";
            nev.bin_size  = res;
            nev.values    = ev_scale.density();
            nev.chr_scale = ev_scale.chr_scale_factors();
            all_norm_expected.push_back(std::move(nev));
        }
    }

    // ============================================================
    //  Rewrite footer: preserve V5 section, replace norm section
    // ============================================================

    fprintf(stderr, "\nRewriting footer with %zu norm vectors...\n",
            all_norm_vecs.size());

    // Read the V5 footer bytes (master index + expected values)
    std::vector<uint8_t> v5_data = read_v5_footer(f, header);

    // Truncate file at footer position (discard old footer)
    // On most systems: truncate to footer_position + 8 (nBytesV5 field) + v5_data
    int64_t new_footer_start = header.footer_position;
    fseek(f, 0, SEEK_END);
    long file_end = ftell(f);

    // We'll rewrite everything from footer_position onward
    fseek(f, (long)new_footer_start, SEEK_SET);

    // nBytesV5
    fwrite_int64(f, (int64_t)v5_data.size());
    fwrite(v5_data.data(), 1, v5_data.size(), f);

    // ---- Normalized expected value vectors ----
    fwrite_int32(f, (int32_t)all_norm_expected.size());
    for (auto& nev : all_norm_expected) {
        fwrite_string(f, nev.norm_type);
        fwrite_string(f, nev.unit);
        fwrite_int32(f, nev.bin_size);
        fwrite_int64(f, (int64_t)nev.values.size());
        for (float v : nev.values) fwrite_float(f, v);
        fwrite_int32(f, (int32_t)nev.chr_scale.size());
        for (auto& [ci, sf] : nev.chr_scale) { fwrite_int32(f, ci); fwrite_float(f, sf); }
    }

    // ---- Normalization vector index ----
    int64_t nvi_pos = (int64_t)ftell(f);

    fwrite_int32(f, (int32_t)all_norm_vecs.size());

    // Write index entries with placeholder positions
    std::vector<int64_t> nv_placeholders;
    nv_placeholders.reserve(all_norm_vecs.size());
    for (auto& nv : all_norm_vecs) {
        fwrite_string(f, nv.norm_type);
        fwrite_int32(f, nv.chr_idx);
        fwrite_string(f, "BP");
        fwrite_int32(f, nv.bin_size);
        nv_placeholders.push_back((int64_t)ftell(f));
        fwrite_int64(f, 0LL);
        int64_t nbytes = (int64_t)nv.values.size() * 4 + 8;
        fwrite_int64(f, nbytes);
    }

    int64_t nvi_end = (int64_t)ftell(f);
    int64_t nvi_len = nvi_end - nvi_pos;

    // ---- Normalization vector arrays ----
    for (size_t i = 0; i < all_norm_vecs.size(); i++) {
        auto& nv = all_norm_vecs[i];
        int64_t pos = (int64_t)ftell(f);
        // Back-patch position
        patch_int64(f, nv_placeholders[i], pos);
        fwrite_int64(f, (int64_t)nv.values.size());
        for (float v : nv.values) fwrite_float(f, v);
    }

    // Truncate file here (remove any old data beyond new footer)
    long new_end = ftell(f);
    fflush(f);
    ftruncate(fileno(f), new_end);

    // ---- Back-patch header: normVectorIndexPosition and normVectorIndexLength ----
    patch_int64(f, header.norm_vi_position_offset, nvi_pos);
    patch_int64(f, header.norm_vi_length_offset,   nvi_len);

    fclose(f);
    fprintf(stderr, "Done. Wrote %zu norm vectors.\n", all_norm_vecs.size());
}
