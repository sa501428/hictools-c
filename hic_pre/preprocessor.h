#pragma once
// Core preprocessing pipeline: contacts pairs → HiC V9 file.
//
// Memory model:
//   Phase 1 (distribute): Stream input pairs into per-chr-pair temp files.
//                          Each temp file stores raw (int32 pos1, int32 pos2, float score).
//                          Input must group each normalized chromosome pair into one
//                          contiguous block, matching Java hic_tools pre. Only the
//                          current chromosome-pair temp file is open for writing.
//
//   Phase 2 (bin+write):  For each chr-pair temp file (one at a time):
//                          - Read all contacts
//                          - For each resolution, accumulate into BlockData hash map
//                          - When in-memory blocks exceed BLOCK_CAPACITY, spill to a second
//                            temp file per resolution (block-level temp files)
//                          - After all contacts read, merge temp blocks and write compressed
//                            blocks to the output .hic file
//
//   Phase 3 (footer):     Write expected values + empty norm vector placeholders.
//
// Parallelism:
//   - Phase 2 bins and compresses different chromosome pairs in parallel using
//     a bounded rolling work queue.
//   - Ordered output I/O overlaps with preparation of later chromosome pairs.

#include "../common/hic_file_def.h"
#include "../common/genome.h"
#include "../common/pair_parser.h"
#include "../common/expected_value.h"
#include "../common/hic_writer.h"
#include "../common/block_writer.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>

struct PreprocessorOptions {
    std::vector<int> resolutions;    // BP resolutions (coarsest first)
    int  mapq_threshold      = 0;    // min(mapq1, mapq2) threshold
    int  count_threshold     = 0;    // min contact count to include in block
    bool intra_only          = false;// only intra-chromosomal contacts
    bool near_diagonal_only  = false;// only contacts within 10 Mb of diagonal
    int  num_threads         = 4;    // number of worker threads
    std::string tmp_dir      = "/tmp"; // directory for temp files
    int  v9_depth_base       = 2;    // V9 block depth base (default 2)
    double scaling_factor    = 1.0;  // global scaling factor applied to scores
    InputFormat input_format = InputFormat::AUTO;
    std::unordered_map<std::string, std::string> attributes; // header attributes
};

class Preprocessor {
public:
    Preprocessor(const Genome& genome, const PreprocessorOptions& opts);
    ~Preprocessor();

    // Run the full pipeline: input_path → output_path.
    void run(const std::string& input_path, const std::string& output_path);

private:
    const Genome&        genome_;
    PreprocessorOptions  opts_;

    // ---- Phase 1: distribute contacts ----

    // Metadata for one buffered output stream per chr-pair index. Only the stream
    // identified by active_stream_key_ is open during phase 1.
    // chr-pair index = chr1_idx * n_chroms + chr2_idx (always chr1 <= chr2)
    struct ChrPairStream {
        FILE*               file    = nullptr;
        std::string         path;
        int64_t             n_contacts = 0;
        // Write buffer
        std::vector<uint8_t> buf;
        static constexpr size_t BUF_SIZE = 256 * 1024; // 256 KB per stream
    };

    std::unordered_map<int, ChrPairStream> chr_pair_streams_;
    int                                     active_stream_key_ = -1;

    // Expected value calculators: one per resolution
    std::vector<std::unique_ptr<ExpectedValueCalculation>> ev_calcs_;

    // Whole-genome blocks are accumulated during phase 1 to avoid a second
    // complete scan of the input.
    std::unordered_map<int, BlockData> whole_genome_blocks_;
    float                              whole_genome_sum_counts_ = 0.0f;

    // Master writer (single-threaded ordered output)
    std::unique_ptr<HicWriter>  writer_;

    // Ordered list of chr-pair indices that have data (for deterministic output order)
    std::vector<int>  chr_pair_order_;

    // ---- Helpers ----

    int chr_pair_index(int c1, int c2) const {
        return c1 * genome_.size() + c2;
    }

    // Phase 1: read all input pairs, write to temp files
    void distribute_pairs(const std::string& input_path);

    // Flush a stream's write buffer to disk
    void flush_stream(ChrPairStream& s);

    // Flush and close the active phase-1 chromosome-pair stream
    void close_active_stream();

    // Close and remove every phase-1 temp file; safe to call during unwinding
    void cleanup_temp_files() noexcept;

    // Write a contact record to the appropriate temp file
    void write_contact(int c1, int c2, int pos1, int pos2, float score);

    // Phase 2: process one chr-pair temp file → write matrix to .hic
    // Returns ZoomWriteData for all resolutions (to be written by write_matrix_to_hic).
    std::vector<ZoomWriteData> process_chr_pair(
        int chr1_idx, int chr2_idx, ChrPairStream& stream
    );

    // Inner: read contacts from temp file, bin at one resolution, compress blocks.
    ZoomWriteData bin_one_resolution(
        int chr1_idx, int chr2_idx,
        int bin_size, int block_bin_count, int block_column_count,
        FILE* contacts_file, bool is_intra
    );

    // Phase 3: finalize
    void write_whole_genome_matrix();

    // Check whether a contact should be skipped
    bool should_skip(int chr1, int chr2, int pos1, int pos2, int mapq1, int mapq2) const;
};
