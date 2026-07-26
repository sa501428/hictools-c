// hic_addnorm: Add normalization vectors to a HiC V9 file.
//
// Computes and appends:
//   VC     - vanilla coverage (row sums)
//   VC_SQRT - square root of VC
//   SCALE  - matrix balancing (intra-chromosomal only)
//
// Normalized expected value vectors are also computed and stored.
//
// Usage:
//   hic_addnorm [options] <hic_file>
//
// The file is modified IN-PLACE. Make a backup first if needed.

#include "norm_updater.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <stdexcept>

static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [options] <hic_file>\n"
        "\n"
        "Adds normalization vectors (VC, VC_SQRT, SCALE) to an existing V9 .hic file.\n"
        "File is modified in-place; make a backup first if needed.\n"
        "\n"
        "Options:\n"
        "  --no-vc        Skip VC normalization\n"
        "  --no-vc-sqrt   Skip VC_SQRT normalization\n"
        "  --no-scale     Skip SCALE normalization\n"
        "  -t <threads>   Number of threads for SCALE (default 4)\n"
        "  --tol <val>    SCALE convergence tolerance (default 1e-4)\n"
        "  --iter <n>     SCALE max total iterations (default 2000)\n"
        "  --min-res <bp> Minimum resolution to run SCALE (default all)\n"
        "  -h             Show this help\n",
        prog);
}

int main(int argc, char* argv[]) {
    if (argc < 2) { usage(argv[0]); return 1; }

    AddNormOptions opts;
    std::string hic_path;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if      (arg == "-h" || arg == "--help")   { usage(argv[0]); return 0; }
        else if (arg == "--no-vc")                  opts.build_vc       = false;
        else if (arg == "--no-vc-sqrt")             opts.build_vc_sqrt  = false;
        else if (arg == "--no-scale")               opts.build_scale    = false;
        else if (arg == "-t" && i+1 < argc)         opts.num_threads    = std::atoi(argv[++i]);
        else if (arg == "--tol" && i+1 < argc)      opts.scale_tolerance = std::atof(argv[++i]);
        else if (arg == "--iter" && i+1 < argc)     opts.scale_max_iter = std::atoi(argv[++i]);
        else if (arg == "--min-res" && i+1 < argc)  opts.min_scale_res  = std::atoi(argv[++i]);
        else if (arg[0] != '-') {
            if (!hic_path.empty()) {
                fprintf(stderr, "Error: multiple positional arguments\n");
                usage(argv[0]);
                return 1;
            }
            hic_path = arg;
        } else {
            fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            usage(argv[0]);
            return 1;
        }
    }

    if (hic_path.empty()) {
        fprintf(stderr, "Error: no hic file specified\n");
        usage(argv[0]);
        return 1;
    }

    fprintf(stderr, "hic_addnorm: processing %s\n", hic_path.c_str());
    fprintf(stderr, "  VC=%s  VC_SQRT=%s  SCALE=%s  threads=%d\n",
            opts.build_vc ? "yes" : "no",
            opts.build_vc_sqrt ? "yes" : "no",
            opts.build_scale ? "yes" : "no",
            opts.num_threads);

    try {
        add_norm(hic_path, opts);
    } catch (std::exception& e) {
        fprintf(stderr, "Error: %s\n", e.what());
        return 1;
    }

    return 0;
}
