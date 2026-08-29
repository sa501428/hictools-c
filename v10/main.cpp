#include "addnorm.h"
#include "pre.h"
#include <climits>
#include <cmath>
#include <iostream>
#include <sstream>

static void usage() {
    std::cout
        << "Usage:\n"
           "  hic_v10 pre [options] <pairs> <output.hic> <chrom.sizes|genome-id>\n"
           "  hic_v10 convert [options] <input.v9.hic> <output.v10.hic>\n"
           "  hic_v10 addnorm [options] <file.v10.hic>\n\n"
           "Common options:\n"
           "  --level N          Zstandard compression level (default 3)\n"
           "\nPre/convert writer options:\n"
           "  --block-bins N     Minimum rectangular block width (default 256; max 4096)\n"
           "  --page-bytes N     Target uncompressed page bytes (default 131072)\n"
           "  --derive T:S       Derive target BP resolution T from source S; repeatable\n"
           "                     Both must exist; exact equivalence is verified\n"
           "  --scores           Force SCORE_FLOAT32, even for integral values\n\n"
           "Pre options (same parsers and MAPQ filtering as hic_pre):\n"
           "  -r N,N,...         BP resolutions (default: existing V9 resolution set)\n"
           "  -q N               Minimum MAPQ (default 0)\n"
           "  -f FORMAT          auto|pairs|short|mnd|bin|bn|hbs\n"
           "  -g GENOME          Genome ID stored in header\n"
           "  -T DIR             Pair-spool directory (default /tmp)\n"
           "  --intra            Retain cis contacts only\n"
           "  --near-diag        Discard cis contacts beyond 10 Mb\n"
           "\nAddnorm options (V10 file is replaced atomically in place):\n"
           "  --no-vc            Skip VC normalization\n"
           "  --no-vc-sqrt       Skip VC_SQRT normalization\n"
           "  --no-scale         Skip SCALE normalization\n"
           "  -t N               SCALE threads (default 4)\n"
           "  --tol X            SCALE convergence tolerance (default 1e-4)\n"
           "  --iter N           SCALE maximum total iterations (default 2000)\n"
           "  --min-res N        Minimum BP resolution for SCALE (default all)\n"
           "  -h, --help         Show this help\n\n"
           "Raw and normalized expected vectors are rebuilt at every stored or derived\n"
           "resolution. Default hic_pre and hic_addnorm remain V9-only.\n";
}
static uint32_t number(const std::string &s) {
    hic10::check(!s.empty() && s.find_first_not_of("0123456789") == std::string::npos,
                 "invalid unsigned integer: " + s);
    return hic10::narrow(std::stoull(s));
}
int main(int argc, char **argv) {
    try {
        if (argc < 2) {
            usage();
            return 1;
        }
        std::string command = argv[1];
        if (command == "-h" || command == "--help") {
            usage();
            return 0;
        }
        hic10::check(command == "pre" || command == "convert" || command == "addnorm",
                     "subcommand must be pre, convert, or addnorm");
        hic10::Options opts;
        hic10::PreOptions pre;
        hic10::AddNormOptions addnorm;
        std::vector<std::string> args;
        for (int i = 2; i < argc; ++i) {
            std::string arg = argv[i];
            auto value = [&]() {
                hic10::check(i + 1 < argc, "missing value for " + arg);
                return std::string(argv[++i]);
            };
            if (arg == "-h" || arg == "--help") {
                usage();
                return 0;
            }
            if (command != "addnorm" && arg == "--scores")
                opts.scores = true;
            else if (arg == "--level") {
                auto s = value();
                size_t n = 0;
                opts.level = std::stoi(s, &n);
                hic10::check(n == s.size(), "invalid compression level");
                addnorm.compression_level = opts.level;
            } else if (command != "addnorm" && arg == "--block-bins")
                opts.blockBins = number(value());
            else if (command != "addnorm" && arg == "--page-bytes")
                opts.pageBytes = number(value());
            else if (command != "addnorm" && arg == "--derive") {
                auto s = value();
                auto colon = s.find(':');
                hic10::check(colon != std::string::npos, "--derive needs target:source");
                opts.derived.emplace_back(number(s.substr(0, colon)), number(s.substr(colon + 1)));
            } else if (command == "pre" && arg == "-r") {
                std::istringstream in(value());
                std::string s;
                while (std::getline(in, s, ','))
                    pre.resolutions.push_back(number(s));
                hic10::check(!pre.resolutions.empty(), "empty resolution list");
            } else if (command == "pre" && arg == "-q") {
                auto n = number(value());
                hic10::check(n <= INT_MAX, "MAPQ too large");
                pre.mapq = n;
            } else if (command == "pre" && arg == "-T")
                pre.tmpDir = value();
            else if (command == "pre" && arg == "-g")
                pre.genome = value();
            else if (command == "pre" && arg == "--intra")
                pre.intra = true;
            else if (command == "pre" && arg == "--near-diag")
                pre.nearDiagonal = true;
            else if (command == "addnorm" && arg == "--no-vc")
                addnorm.vc = false;
            else if (command == "addnorm" && arg == "--no-vc-sqrt")
                addnorm.vc_sqrt = false;
            else if (command == "addnorm" && arg == "--no-scale")
                addnorm.scale = false;
            else if (command == "addnorm" && arg == "-t") {
                auto n = number(value());
                hic10::check(n > 0 && n <= INT_MAX, "thread count must be positive");
                addnorm.threads = static_cast<int>(n);
            } else if (command == "addnorm" && arg == "--tol") {
                auto s = value();
                size_t n = 0;
                addnorm.tolerance = std::stod(s, &n);
                hic10::check(n == s.size() && addnorm.tolerance > 0 &&
                                 std::isfinite(addnorm.tolerance),
                             "tolerance must be positive");
            } else if (command == "addnorm" && arg == "--iter") {
                auto n = number(value());
                hic10::check(n > 0 && n <= INT_MAX, "iteration count must be positive");
                addnorm.max_iterations = static_cast<int>(n);
            } else if (command == "addnorm" && arg == "--min-res") {
                auto n = number(value());
                hic10::check(n <= INT_MAX, "minimum resolution too large");
                addnorm.minimum_scale_resolution = static_cast<int>(n);
            } else if (command == "pre" && arg == "-f") {
                auto f = value();
                const std::map<std::string, InputFormat> formats = {
                    {"auto", InputFormat::AUTO},   {"pairs", InputFormat::PAIRS},
                    {"short", InputFormat::SHORT}, {"mnd", InputFormat::MND},
                    {"bin", InputFormat::BIN},     {"bn", InputFormat::BN},
                    {"hbs", InputFormat::HBS}};
                hic10::check(formats.count(f), "unknown input format " + f);
                pre.format = formats.at(f);
            } else {
                hic10::check(arg.empty() || arg[0] != '-', "unknown option " + arg);
                args.push_back(arg);
            }
        }
        size_t wanted = command == "pre" ? 3 : command == "convert" ? 2 : 1;
        hic10::check(args.size() == wanted,
                     "incorrect positional arguments (see --help)");
        if (command == "pre")
            hic10::pre(args[0], args[1], args[2], opts, pre);
        else if (command == "convert")
            hic10::convert(args[0], args[1], opts);
        else
            hic10::add_norm_v10(args[0], addnorm);
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
