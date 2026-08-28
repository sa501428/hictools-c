#include "preprocessor.h"
#include "../common/thread_pool.h"
#include "../common/little_endian.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <cassert>
#include <deque>
#include <future>
#include <sstream>
#include <cerrno>
#include <unistd.h>  // close(), ftruncate()

// ============================================================
//  Constructor
// ============================================================

Preprocessor::Preprocessor(const Genome& genome, const PreprocessorOptions& opts)
    : genome_(genome), opts_(opts) {
    if (opts_.resolutions.empty()) {
        opts_.resolutions.assign(DEFAULT_BP_RESOLUTIONS,
                                  DEFAULT_BP_RESOLUTIONS + N_DEFAULT_BP_RESOLUTIONS);
    }
    // One expected value calculator per resolution
    for (int r : opts_.resolutions) {
        ev_calcs_.push_back(std::make_unique<ExpectedValueCalculation>(genome_, r));
    }
}

Preprocessor::~Preprocessor() {
    cleanup_temp_files();
}

// ============================================================
//  Main entry point
// ============================================================

void Preprocessor::run(const std::string& input_path, const std::string& output_path) {
  try {
    fprintf(stderr, "hic_pre: distributing contacts to temp files...\n");
    distribute_pairs(input_path);

    fprintf(stderr, "hic_pre: %zu chromosome pairs found\n", chr_pair_order_.size());

    // Open writer
    writer_ = std::make_unique<HicWriter>(output_path);

    auto attrs = opts_.attributes;
    if (attrs.find("software") == attrs.end())
        attrs["software"] = "hic_pre C++ v1.0";
    if (attrs.find("hicFileScaling") == attrs.end())
        attrs["hicFileScaling"] = "1.0";
    if (attrs.find("v9BaseDepth") == attrs.end() && opts_.v9_depth_base != 2)
        attrs["v9BaseDepth"] = std::to_string(opts_.v9_depth_base);

    fprintf(stderr, "hic_pre: writing header...\n");
    writer_->write_header(genome_, opts_.resolutions, attrs);

    // Process the whole-genome (All × All) matrix first
    fprintf(stderr, "hic_pre: computing whole-genome matrix...\n");
    write_whole_genome_matrix();

    // Process chromosome pairs in order using thread pool
    fprintf(stderr, "hic_pre: writing chromosome pair matrices...\n");

    // Prepare pairs in parallel and write them in deterministic input order.
    // Keep a bounded rolling window so slow pairs do not create a barrier after
    // every N items and completed compressed data cannot grow without bound.
    int n_threads = std::max(1, opts_.num_threads);
    ThreadPool pool(n_threads);

    // Collect pairs with data
    std::vector<std::pair<int,int>> pairs_with_data;
    for (int key : chr_pair_order_) {
        int c1 = key / genome_.size();
        int c2 = key % genome_.size();
        pairs_with_data.push_back({c1, c2});
    }

    using PairResult = std::pair<std::pair<int,int>, std::vector<ZoomWriteData>>;
    std::deque<std::future<PairResult>> pending;
    size_t next_pair = 0;
    const size_t window_size = std::max<size_t>(1, (size_t)n_threads * 2);

    auto submit_pair = [&](size_t pair_index) {
        int c1 = pairs_with_data[pair_index].first;
        int c2 = pairs_with_data[pair_index].second;
        int key = chr_pair_index(c1, c2);
        auto& stream = chr_pair_streams_.at(key);
        pending.push_back(pool.submit([this, c1, c2, &stream]() -> PairResult {
            auto zooms = process_chr_pair(c1, c2, stream);
            return {{c1, c2}, std::move(zooms)};
        }));
    };

    while (next_pair < pairs_with_data.size() && pending.size() < window_size) {
        submit_pair(next_pair++);
    }
    while (!pending.empty()) {
        PairResult result = pending.front().get();
        pending.pop_front();

        int c1 = result.first.first;
        int c2 = result.first.second;
        fprintf(stderr, "  writing chr %s x chr %s\n",
                genome_.at(c1).name.c_str(), genome_.at(c2).name.c_str());
        writer_->write_matrix(c1, c2, result.second);

        if (next_pair < pairs_with_data.size()) submit_pair(next_pair++);
    }

    cleanup_temp_files();

    // Write footer
    fprintf(stderr, "hic_pre: writing footer...\n");
    std::vector<ExpectedValueCalculation*> ev_ptrs;
    for (auto& ev : ev_calcs_) ev_ptrs.push_back(ev.get());

    writer_->write_footer(ev_ptrs, {}, {});
    fprintf(stderr, "hic_pre: done.\n");
  } catch (...) {
    cleanup_temp_files();
    throw;
  }
}

// ============================================================
//  Phase 1: distribute pairs
// ============================================================

void Preprocessor::flush_stream(ChrPairStream& s) {
    if (s.buf.empty()) return;
    if (!s.file) {
        throw std::runtime_error("Internal error: chromosome-pair temp stream is closed");
    }
    size_t written = fwrite(s.buf.data(), 1, s.buf.size(), s.file);
    if (written != s.buf.size()) {
        throw std::runtime_error("Cannot write temp file " + s.path + ": "
                                 + std::strerror(errno));
    }
    s.buf.clear();
}

void Preprocessor::close_active_stream() {
    if (active_stream_key_ < 0) return;

    auto it = chr_pair_streams_.find(active_stream_key_);
    if (it == chr_pair_streams_.end()) {
        throw std::runtime_error("Internal error: active chromosome-pair stream is missing");
    }

    ChrPairStream& s = it->second;
    flush_stream(s);
    if (fflush(s.file) != 0) {
        throw std::runtime_error("Cannot flush temp file " + s.path + ": "
                                 + std::strerror(errno));
    }
    if (fclose(s.file) != 0) {
        s.file = nullptr;
        throw std::runtime_error("Cannot close temp file " + s.path + ": "
                                 + std::strerror(errno));
    }
    s.file = nullptr;
    active_stream_key_ = -1;
}

void Preprocessor::cleanup_temp_files() noexcept {
    for (auto& [key, stream] : chr_pair_streams_) {
        if (stream.file) {
            fclose(stream.file);
            stream.file = nullptr;
        }
        if (!stream.path.empty()) remove(stream.path.c_str());
    }
    chr_pair_streams_.clear();
    chr_pair_order_.clear();
    whole_genome_blocks_.clear();
    whole_genome_sum_counts_ = 0.0f;
    active_stream_key_ = -1;
}

void Preprocessor::write_contact(int c1, int c2, int pos1, int pos2, float score) {
    // Ensure upper triangle: c1 <= c2 (and pos1 <= pos2 for intra)
    if (c1 > c2 || (c1 == c2 && pos1 > pos2)) {
        std::swap(c1, c2);
        std::swap(pos1, pos2);
    }

    int key = chr_pair_index(c1, c2);
    if (key != active_stream_key_) {
        close_active_stream();

        // Java hic_tools pre requires each chromosome combination to appear in
        // one contiguous input block. Enforce the same contract here.
        if (chr_pair_streams_.find(key) != chr_pair_streams_.end()) {
            throw std::runtime_error(
                "Chromosome combination " + genome_.at(c1).name + "_"
                + genome_.at(c2).name
                + " appears in multiple blocks; input must be grouped by chromosome pair");
        }

        ChrPairStream s;
        std::string tmp_path = opts_.tmp_dir + "/hic_pre_c"
                             + std::to_string(c1) + "_" + std::to_string(c2) + "_XXXXXX";
        std::vector<char> tmpl(tmp_path.begin(), tmp_path.end());
        tmpl.push_back('\0');
        int fd = mkstemp(tmpl.data());
        if (fd < 0) {
            throw std::runtime_error("Cannot create temp file " + tmp_path + ": "
                                     + std::strerror(errno));
        }
        s.path = tmpl.data();
        s.file = fdopen(fd, "wb");
        if (!s.file) {
            int saved_errno = errno;
            close(fd);
            remove(s.path.c_str());
            throw std::runtime_error("Cannot open temp file " + s.path + ": "
                                     + std::strerror(saved_errno));
        }
        chr_pair_streams_.emplace(key, std::move(s));
        chr_pair_order_.push_back(key);
        active_stream_key_ = key;
    }

    ChrPairStream& s = chr_pair_streams_.at(key);

    // Write (pos1, pos2, score) = 12 bytes to buffer
    uint8_t rec[12];
    // pos1
    rec[0] = (uint8_t)(pos1 & 0xFF); rec[1] = (uint8_t)((pos1>>8) & 0xFF);
    rec[2] = (uint8_t)((pos1>>16) & 0xFF); rec[3] = (uint8_t)((pos1>>24) & 0xFF);
    // pos2
    rec[4] = (uint8_t)(pos2 & 0xFF); rec[5] = (uint8_t)((pos2>>8) & 0xFF);
    rec[6] = (uint8_t)((pos2>>16) & 0xFF); rec[7] = (uint8_t)((pos2>>24) & 0xFF);
    // score
    uint32_t sbits; std::memcpy(&sbits, &score, 4);
    rec[8]  = (uint8_t)(sbits & 0xFF); rec[9]  = (uint8_t)((sbits>>8) & 0xFF);
    rec[10] = (uint8_t)((sbits>>16) & 0xFF); rec[11] = (uint8_t)((sbits>>24) & 0xFF);

    s.buf.insert(s.buf.end(), rec, rec + 12);
    s.n_contacts++;

    if (s.buf.size() >= ChrPairStream::BUF_SIZE) flush_stream(s);
}

bool Preprocessor::should_skip(int chr1, int chr2, int pos1, int pos2,
                                int mapq1, int mapq2) const {
    if (opts_.intra_only && chr1 != chr2) return true;
    if (opts_.near_diagonal_only && chr1 == chr2 && std::abs(pos1 - pos2) > 10000000) return true;
    if (std::min(mapq1, mapq2) < opts_.mapq_threshold) return true;
    return false;
}

void Preprocessor::distribute_pairs(const std::string& input_path) {
    auto iter = open_pair_iterator(input_path, genome_, opts_.input_format);
    if (auto source = iter->source_resolution()) {
        for (auto r : opts_.resolutions)
            if (r <= 0 || static_cast<uint32_t>(r) % source)
                throw std::runtime_error("HBS: output resolutions must be multiples of the source resolution");
    }
    AlignmentPair pair;
    int64_t n_total = 0, n_skipped = 0;

    const int64_t genome_length = genome_.at(0).length;
    const int genome_bin_size = std::max<int64_t>(1, genome_length / 500);
    auto genome_bp = compute_block_params(
        genome_length, genome_length, genome_bin_size, false, BLOCK_CAPACITY);

    while (iter->next(pair)) {
        n_total++;
        if (should_skip(pair.chr1, pair.chr2, pair.pos1, pair.pos2,
                        pair.mapq1, pair.mapq2)) {
            n_skipped++;
            continue;
        }
        pair.pos1 = (int)std::min<int64_t>(pair.pos1, genome_.at(pair.chr1).length);
        pair.pos2 = (int)std::min<int64_t>(pair.pos2, genome_.at(pair.chr2).length);

        float score = pair.score * (float)opts_.scaling_factor;

        write_contact(pair.chr1, pair.chr2, pair.pos1, pair.pos2, score);

        int64_t gw1 = genome_.genome_position(pair.chr1, pair.pos1);
        int64_t gw2 = genome_.genome_position(pair.chr2, pair.pos2);
        if (gw1 >= 0 && gw2 >= 0) {
            if (gw1 > gw2) std::swap(gw1, gw2);
            int bin1 = (int)(gw1 / genome_bin_size);
            int bin2 = (int)(gw2 / genome_bin_size);
            int bnum = inter_block_number(
                bin1, bin2, genome_bp.block_bin_count, genome_bp.block_column_count);
            whole_genome_blocks_.try_emplace(bnum, bnum).first->second.add(
                bin1, bin2, score);
            whole_genome_sum_counts_ += score;
        }

        // Accumulate expected values (intra only)
        if (pair.chr1 == pair.chr2) {
            int c = (pair.chr1 <= pair.chr2) ? pair.chr1 : pair.chr2;
            int b1 = pair.pos1;
            int b2 = pair.pos2;
            if (b1 > b2) std::swap(b1, b2);
            for (size_t ri = 0; ri < opts_.resolutions.size(); ri++) {
                int r = opts_.resolutions[ri];
                ev_calcs_[ri]->add_distance(c, b1 / r, b2 / r, score);
            }
        }
    }
    iter->close();

    close_active_stream();

    fprintf(stderr, "  total pairs read: %lld, skipped: %lld\n",
            (long long)n_total, (long long)n_skipped);
}

// ============================================================
//  Whole-genome (All × All) matrix
// ============================================================

void Preprocessor::write_whole_genome_matrix() {
    // The whole-genome matrix uses a single coarse resolution binned at
    // genome-wide coordinates (like the Java code).
    const Chromosome& all_chr = genome_.at(0); // "All"
    int64_t total_len = all_chr.length;

    // Same binSize calculation as HiCFileBuilder.getInitialGenomeWideMatrixPP
    int bin_size = (int)(total_len / 500);
    if (bin_size == 0) bin_size = 1;
    auto bp = compute_block_params(total_len, total_len, bin_size, false, BLOCK_CAPACITY);

    // Build ZoomWriteData
    ZoomWriteData z;
    z.bin_size          = bin_size;
    z.block_bin_count   = bp.block_bin_count;
    z.block_column_count = bp.block_column_count;
    z.sum_counts        = whole_genome_sum_counts_;

    z.blocks.reserve(whole_genome_blocks_.size());
    for (auto& [bnum, bd] : whole_genome_blocks_) {
        z.blocks.push_back({bnum, compress_block(bd)});
    }
    std::sort(z.blocks.begin(), z.blocks.end(),
              [](const auto& a, const auto& b){ return a.block_number < b.block_number; });
    whole_genome_blocks_.clear();

    writer_->write_matrix(0, 0, {z});
}

// ============================================================
//  Phase 2: process one chr-pair temp file
// ============================================================

ZoomWriteData Preprocessor::bin_one_resolution(
    int chr1_idx, int chr2_idx,
    int bin_size, int block_bin_count, int block_column_count,
    FILE* contacts_file, bool is_intra)
{
    // Block accumulator (in-memory blocks, spill to disk when full)
    std::unordered_map<int, BlockData> block_map;
    block_map.reserve(BLOCK_CAPACITY);

    // Spill file for when block_map exceeds capacity
    // (using a per-resolution temp file)
    std::string spill_path;
    FILE* spill_file = nullptr;
    std::vector<std::string> spill_paths;

    auto maybe_spill = [&]() {
        if ((int)block_map.size() > BLOCK_CAPACITY) {
            // Sort by block number and write to spill
            if (!spill_file) {
                spill_path = opts_.tmp_dir + "/hic_pre_spill_XXXXXX";
                std::vector<char> tmpl(spill_path.begin(), spill_path.end());
                tmpl.push_back('\0');
                int fd = mkstemp(tmpl.data());
                if (fd < 0) throw std::runtime_error("mkstemp failed");
                close(fd);
                spill_path = tmpl.data();
                spill_file = fopen(spill_path.c_str(), "wb");
                if (!spill_file) throw std::runtime_error("Cannot open spill file");
                spill_paths.push_back(spill_path);
            }
            std::vector<BlockData> to_spill;
            to_spill.reserve(block_map.size());
            for (auto& [k, v] : block_map) to_spill.push_back(std::move(v));
            block_map.clear();
            std::sort(to_spill.begin(), to_spill.end(),
                      [](const BlockData& a, const BlockData& b){ return a.block_number < b.block_number; });
            write_blocks_to_tempfile(spill_file, to_spill);
            fflush(spill_file);
        }
    };

    // Read contacts from temp file (pos1, pos2, score) = 12 bytes each
    rewind(contacts_file);
    uint8_t rec[12];
    float sum_counts = 0.0f;

    while (fread(rec, 1, 12, contacts_file) == 12) {
        int32_t p1 = (int32_t)((uint32_t)rec[0] | ((uint32_t)rec[1]<<8)
                              | ((uint32_t)rec[2]<<16) | ((uint32_t)rec[3]<<24));
        int32_t p2 = (int32_t)((uint32_t)rec[4] | ((uint32_t)rec[5]<<8)
                              | ((uint32_t)rec[6]<<16) | ((uint32_t)rec[7]<<24));
        uint32_t sbits = (uint32_t)rec[8] | ((uint32_t)rec[9]<<8)
                       | ((uint32_t)rec[10]<<16) | ((uint32_t)rec[11]<<24);
        float score;
        std::memcpy(&score, &sbits, 4);

        int bin1 = p1 / bin_size;
        int bin2 = p2 / bin_size;

        // For intra, ensure bin1 <= bin2
        if (is_intra && bin1 > bin2) std::swap(bin1, bin2);

        int bnum;
        if (is_intra) {
            bnum = intra_block_number(bin1, bin2, block_bin_count,
                                      block_column_count, opts_.v9_depth_base);
        } else {
            bnum = inter_block_number(bin1, bin2, block_bin_count, block_column_count);
        }

        block_map.try_emplace(bnum, bnum).first->second.add(bin1, bin2, score);
        sum_counts += score;
        // For intra off-diagonal, count twice (upper+lower triangle)
        if (is_intra && bin1 != bin2) sum_counts += score;

        maybe_spill();
    }

    if (spill_file) { fclose(spill_file); spill_file = nullptr; }

    // Merge spill files back into block_map
    for (auto& sp : spill_paths) {
        FILE* sf = fopen(sp.c_str(), "rb");
        if (sf) {
            merge_blocks_from_tempfile(sf, block_map);
            fclose(sf);
            remove(sp.c_str());
        }
    }

    // Build ZoomWriteData
    ZoomWriteData z;
    z.bin_size           = bin_size;
    z.block_bin_count    = block_bin_count;
    z.block_column_count = block_column_count;
    z.sum_counts         = sum_counts;

    z.blocks.reserve(block_map.size());
    for (auto& [bnum, bd] : block_map) {
        z.blocks.push_back({bnum, compress_block(bd)});
    }
    std::sort(z.blocks.begin(), z.blocks.end(),
              [](const auto& a, const auto& b){ return a.block_number < b.block_number; });

    return z;
}

std::vector<ZoomWriteData> Preprocessor::process_chr_pair(
    int chr1_idx, int chr2_idx, ChrPairStream& stream)
{
    bool is_intra = (chr1_idx == chr2_idx);
    const Chromosome& chr1 = genome_.at(chr1_idx);
    const Chromosome& chr2 = genome_.at(chr2_idx);

    std::unique_ptr<FILE, decltype(&fclose)> cf(fopen(stream.path.c_str(), "rb"), &fclose);
    if (!cf) {
        throw std::runtime_error("Cannot open temp file " + stream.path + ": "
                                 + std::strerror(errno));
    }

    std::vector<ZoomWriteData> zooms;
    zooms.reserve(opts_.resolutions.size());

    for (size_t ri = 0; ri < opts_.resolutions.size(); ri++) {
        int bin_size = opts_.resolutions[ri];
        auto bp = compute_block_params(chr1.length, chr2.length, bin_size, is_intra);
        auto z = bin_one_resolution(
            chr1_idx, chr2_idx,
            bin_size, bp.block_bin_count, bp.block_column_count,
            cf.get(), is_intra
        );
        zooms.push_back(std::move(z));
    }
    return zooms;
}
