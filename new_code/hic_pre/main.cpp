// hic_pre: Convert aligned read pairs to HiC V9 format.
//
// Usage:
//   hic_pre [options] <input_pairs> <output.hic> <chrom.sizes>
//
// Input formats auto-detected from extension / header:
//   .mnd / merged_nodups   (Juicer format)
//   .pairs                 (4DN format)
//   short format           (chr1 pos1 chr2 pos2)
//   All formats support .gz compression.
//
// Options:
//   -r <resolutions>  Comma-separated BP resolutions (default: 2500000,...,1000)
//   -q <mapq>         MAPQ threshold (default: 0)
//   -t <threads>      Number of threads (default: 4)
//   -T <tmpdir>       Temp directory (default: /tmp)
//   -g <genome>       Genome ID string stored in header (e.g. hg38)
//   --intra           Only include intra-chromosomal contacts
//   --near-diagonal   Only include contacts within 10 Mb of diagonal
//   -d <depth>        V9 depth base (default: 2)

#include "preprocessor.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>

static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [options] <input_pairs> <output.hic> <chrom.sizes>\n"
        "\n"
        "Options:\n"
        "  -r <res,...>   Comma-separated BP resolutions\n"
        "                 Default: 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000\n"
        "  -q <mapq>      Min MAPQ threshold (default 0)\n"
        "  -t <threads>   Number of worker threads (default 4)\n"
        "  -T <tmpdir>    Temp file directory (default /tmp)\n"
        "  -g <genome>    Genome ID (e.g. hg38); stored in file header\n"
        "  --intra        Only keep intra-chromosomal contacts\n"
        "  --near-diag    Only keep contacts within 10 Mb of diagonal\n"
        "  -d <depth>     V9 block depth base (default 2, range 1-10)\n"
        "  -f <format>    Input format: auto|mnd|short|pairs (default auto)\n"
        "  -h             Show this help\n",
        prog);
}

static std::vector<int> parse_resolutions(const std::string& s) {
    std::vector<int> res;
    std::istringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) res.push_back(std::stoi(tok));
    }
    // Sort descending (coarsest first) to match Java convention
    std::sort(res.begin(), res.end(), std::greater<int>());
    return res;
}

int main(int argc, char* argv[]) {
    if (argc < 4) { usage(argv[0]); return 1; }

    PreprocessorOptions opts;
    std::string genome_id;
    std::string fmt_str = "auto";

    // Parse options
    int i = 1;
    for (; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") { usage(argv[0]); return 0; }
        else if (arg == "-r" && i+1 < argc) { opts.resolutions = parse_resolutions(argv[++i]); }
        else if (arg == "-q" && i+1 < argc) { opts.mapq_threshold = std::atoi(argv[++i]); }
        else if (arg == "-t" && i+1 < argc) { opts.num_threads = std::atoi(argv[++i]); }
        else if (arg == "-T" && i+1 < argc) { opts.tmp_dir = argv[++i]; }
        else if (arg == "-g" && i+1 < argc) { genome_id = argv[++i]; }
        else if (arg == "-d" && i+1 < argc) { opts.v9_depth_base = std::atoi(argv[++i]); }
        else if (arg == "-f" && i+1 < argc) { fmt_str = argv[++i]; }
        else if (arg == "--intra")         { opts.intra_only = true; }
        else if (arg == "--near-diag")     { opts.near_diagonal_only = true; }
        else break;
    }

    if (argc - i < 3) {
        fprintf(stderr, "Error: need <input_pairs> <output.hic> <chrom.sizes>\n");
        usage(argv[0]);
        return 1;
    }

    std::string input_path  = argv[i];
    std::string output_path = argv[i+1];
    std::string chrom_sizes = argv[i+2];

    // Input format
    if      (fmt_str == "mnd")   opts.input_format = InputFormat::MND;
    else if (fmt_str == "short") opts.input_format = InputFormat::SHORT;
    else if (fmt_str == "pairs") opts.input_format = InputFormat::PAIRS;
    else                          opts.input_format = InputFormat::AUTO;

    // Genome ID attribute
    if (!genome_id.empty()) opts.attributes["genomeID"] = genome_id;

    // Load genome
    Genome genome;
    try {
        genome = Genome::from_chrom_sizes(chrom_sizes);
    } catch (std::exception& e) {
        fprintf(stderr, "Error loading chrom.sizes: %s\n", e.what());
        return 1;
    }

    fprintf(stderr, "Loaded %d chromosomes from %s\n",
            genome.size() - 1, chrom_sizes.c_str());

    if (opts.resolutions.empty()) {
        opts.resolutions.assign(DEFAULT_BP_RESOLUTIONS,
                                 DEFAULT_BP_RESOLUTIONS + N_DEFAULT_BP_RESOLUTIONS);
    }

    fprintf(stderr, "Resolutions:");
    for (int r : opts.resolutions) fprintf(stderr, " %d", r);
    fprintf(stderr, "\n");
    fprintf(stderr, "Threads: %d\n", opts.num_threads);
    fprintf(stderr, "Temp dir: %s\n", opts.tmp_dir.c_str());

    try {
        Preprocessor pre(genome, opts);
        pre.run(input_path, output_path);
    } catch (std::exception& e) {
        fprintf(stderr, "Error: %s\n", e.what());
        return 1;
    }

    return 0;
}
