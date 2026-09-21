#include "io.h"
#include "build.h"
#include "stage.h"
#include "writer.h"
#include "normalize.h"

#include <climits>
#include <iostream>
#include <string>
#include <vector>

using namespace hic10large;

namespace {
void usage() {
    std::cout
        << "Usage:\n"
           "  hic_v10_large stage [--chunk-records N] [--progress-records N]\n"
           "      <input.hbs.gz> <work-dir>\n"
           "  hic_v10_large inspect [--verify] <work-dir/stage.manifest>\n\n"
           "  hic_v10_large build-cells -r N,N,... [--root-only] [--memory SIZE] [--fan-in N]\n"
           "      <stage.manifest> <build-dir>\n\n"
           "  hic_v10_large plan -r N,N,... [--root-only] [--fan-in N] <stage.manifest> <build-dir>\n"
           "  hic_v10_large map-root -r N,N,... [--memory SIZE] [--fan-in N] <stage.manifest> <build-dir> <shard-id>\n"
           "  hic_v10_large reduce-root -r N,N,... [--fan-in N] <stage.manifest> <build-dir> <pair-id> <group-id>\n"
           "  hic_v10_large build-pair -r N,N,... [--root-only] [--memory SIZE] <stage.manifest> <build-dir> <pair-id>\n"
           "  hic_v10_large finalize-build -r N,N,... [--root-only] <stage.manifest> <build-dir>\n\n"
           "  hic_v10_large write [--genome NAME] [--memory SIZE] [--tmp DIR] [-t N]\n"
           "      [--vectors vectors.manifest]\n"
           "      [--derive TARGET:SOURCE] <stage.manifest> <build.manifest> <output.hic>\n\n"
           "  hic_v10_large normalize [--memory SIZE] [--tmp DIR] [-t N]\n"
           "      [--no-vc] [--no-vc-sqrt] [--no-scale]\n"
           "      <stage.manifest> <build.manifest> <vectors-dir>\n\n"
           "  hic_v10_large plan-normalize [normalization options]\n"
           "      <stage.manifest> <build.manifest> <vectors-dir>\n"
           "  hic_v10_large normalize-chr [normalization options]\n"
           "      <stage.manifest> <build.manifest> <vectors-dir> <chromosome-id>\n"
           "  hic_v10_large expected-res [normalization options]\n"
           "      <stage.manifest> <build.manifest> <vectors-dir> <resolution-index>\n"
           "  hic_v10_large finalize-vectors [normalization options]\n"
           "      <stage.manifest> <build.manifest> <vectors-dir>\n\n"
           "The large-data builder is separate from hic_v10 pre. Stage validates the\n"
           "entire gzip/HBS stream, checks chromosome-pair contiguity, and writes\n"
           "record-aligned shards plus a durable manifest.\n";
}

uint64_t number(const std::string &text) {
    require(!text.empty() && text.find_first_not_of("0123456789") == std::string::npos,
            "invalid unsigned integer: " + text);
    return std::stoull(text);
}

uint32_t number32(const std::string &text, const std::string &what) {
    uint64_t value = number(text);
    require(value <= UINT32_MAX, what + " is outside uint32 range");
    return static_cast<uint32_t>(value);
}

std::vector<uint32_t> resolution_list(const std::string &text) {
    std::vector<uint32_t> result;
    size_t begin = 0;
    while (begin <= text.size()) {
        size_t end = text.find(',', begin);
        std::string word = text.substr(begin, end == std::string::npos ? end : end - begin);
        uint64_t value = number(word);
        require(value && value <= UINT32_MAX, "resolution outside uint32 range");
        result.push_back(static_cast<uint32_t>(value));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return result;
}

std::pair<uint32_t, uint32_t> derivation(const std::string &text) {
    size_t colon = text.find(':');
    require(colon != std::string::npos, "--derive requires TARGET:SOURCE");
    uint64_t target = number(text.substr(0, colon));
    uint64_t source = number(text.substr(colon + 1));
    require(target <= UINT32_MAX && source <= UINT32_MAX, "derived resolution is too large");
    return {static_cast<uint32_t>(target), static_cast<uint32_t>(source)};
}
} // namespace

int main(int argc, char **argv) {
    try {
        if (argc < 2 || std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
            usage();
            return argc < 2 ? 1 : 0;
        }
        std::string command = argv[1];
        if (command == "stage") {
            StageOptions options;
            std::vector<std::string> args;
            for (int i = 2; i < argc; ++i) {
                std::string arg = argv[i];
                if (arg == "--chunk-records") {
                    require(i + 1 < argc, "missing value for --chunk-records");
                    options.chunk_records = number(argv[++i]);
                } else if (arg == "--progress-records") {
                    require(i + 1 < argc, "missing value for --progress-records");
                    options.progress_records = number(argv[++i]);
                } else {
                    require(arg.empty() || arg[0] != '-', "unknown option " + arg);
                    args.push_back(arg);
                }
            }
            require(args.size() == 2, "stage needs input.hbs.gz and work-dir");
            stage_hbs(args[0], args[1], options);
        } else if (command == "inspect") {
            bool verify = false;
            std::vector<std::string> args;
            for (int i = 2; i < argc; ++i) {
                std::string arg = argv[i];
                if (arg == "--verify") verify = true;
                else {
                    require(arg.empty() || arg[0] != '-', "unknown option " + arg);
                    args.push_back(arg);
                }
            }
            require(args.size() == 1, "inspect needs stage.manifest");
            inspect_stage(args[0], verify);
        } else if (command == "build-cells" || command == "plan" || command == "map-root" ||
                   command == "reduce-root" ||
                   command == "build-pair" || command == "finalize-build") {
            BuildOptions options;
            std::vector<std::string> args;
            for (int i = 2; i < argc; ++i) {
                std::string arg = argv[i];
                if (arg == "-r") {
                    require(i + 1 < argc, "missing value for -r");
                    options.resolutions = resolution_list(argv[++i]);
                } else if (arg == "--memory") {
                    require(i + 1 < argc, "missing value for --memory");
                    options.memory_bytes = parse_size(argv[++i]);
                } else if (arg == "--fan-in") {
                    require(i + 1 < argc, "missing value for --fan-in");
                    options.merge_fan_in = static_cast<size_t>(number(argv[++i]));
                } else if (arg == "--root-only") {
                    options.root_only = true;
                } else {
                    require(arg.empty() || arg[0] != '-', "unknown option " + arg);
                    args.push_back(arg);
                }
            }
            size_t wanted = command == "reduce-root" ? 4 :
                            command == "map-root" || command == "build-pair" ? 3 : 2;
            require(args.size() == wanted, command + " received the wrong number of arguments");
            if (command == "build-cells") build_cells(args[0], args[1], options);
            else if (command == "plan") print_build_tasks(args[0], args[1], options);
            else if (command == "map-root")
                map_root_shard(args[0], args[1], static_cast<size_t>(number(args[2])), options);
            else if (command == "reduce-root")
                reduce_root_group(args[0], args[1], static_cast<size_t>(number(args[2])),
                                  static_cast<size_t>(number(args[3])), options);
            else if (command == "build-pair")
                build_pair_cells(args[0], args[1], static_cast<size_t>(number(args[2])), options);
            else finalize_build(args[0], args[1], options);
        } else if (command == "write") {
            WriteOptions options;
            std::vector<std::string> args;
            for (int i = 2; i < argc; ++i) {
                std::string arg = argv[i];
                auto value = [&]() {
                    require(i + 1 < argc, "missing value for " + arg);
                    return std::string(argv[++i]);
                };
                if (arg == "--genome") options.genome = value();
                else if (arg == "--memory") options.memory_bytes = parse_size(value());
                else if (arg == "--tmp") options.temporary_directory = value();
                else if (arg == "--vectors") options.vector_manifest = value();
                else if (arg == "--block-bins") options.block_bins = static_cast<uint32_t>(number(value()));
                else if (arg == "--level") options.compression_level = std::stoi(value());
                else if (arg == "-t") {
                    uint64_t threads = number(value());
                    require(threads && threads <= INT_MAX, "writer thread count is invalid");
                    options.threads = static_cast<int>(threads);
                }
                else if (arg == "--fan-in") options.merge_fan_in = static_cast<size_t>(number(value()));
                else if (arg == "--derive") options.derived.push_back(derivation(value()));
                else {
                    require(arg.empty() || arg[0] != '-', "unknown option " + arg);
                    args.push_back(arg);
                }
            }
            require(args.size() == 3, "write needs stage.manifest, build.manifest, and output.hic");
            write_v10(args[0], args[1], args[2], options);
        } else if (command == "normalize" || command == "plan-normalize" ||
                   command == "normalize-chr" || command == "expected-res" ||
                   command == "finalize-vectors") {
            NormalizeOptions options;
            std::vector<std::string> args;
            for (int i = 2; i < argc; ++i) {
                std::string arg = argv[i];
                auto value = [&]() {
                    require(i + 1 < argc, "missing value for " + arg);
                    return std::string(argv[++i]);
                };
                if (arg == "--memory") options.memory_bytes = parse_size(value());
                else if (arg == "--tmp") options.temporary_directory = value();
                else if (arg == "-t") options.scale_options.threads = static_cast<int>(number(value()));
                else if (arg == "--tol") options.scale_options.tolerance = std::stod(value());
                else if (arg == "--iter") options.scale_options.max_iterations = static_cast<int>(number(value()));
                else if (arg == "--fan-in") options.scale_options.merge_fan_in = static_cast<size_t>(number(value()));
                else if (arg == "--no-vc") options.vc = false;
                else if (arg == "--no-vc-sqrt") options.vc_sqrt = false;
                else if (arg == "--no-scale") options.scale = false;
                else {
                    require(arg.empty() || arg[0] != '-', "unknown option " + arg);
                    args.push_back(arg);
                }
            }
            size_t wanted = command == "normalize-chr" || command == "expected-res" ? 4 : 3;
            require(args.size() == wanted, command + " received the wrong number of arguments");
            if (command == "normalize") normalize_cells(args[0], args[1], args[2], options);
            else if (command == "plan-normalize")
                print_normalize_tasks(args[0], args[1], args[2], options);
            else if (command == "normalize-chr")
                normalize_chromosome(args[0], args[1], args[2],
                                     number32(args[3], "chromosome ID"), options);
            else if (command == "expected-res")
                expected_resolution(args[0], args[1], args[2],
                                    number32(args[3], "resolution index"), options);
            else finalize_vectors(args[0], args[1], args[2], options);
        } else {
            fail("unknown command " + command);
        }
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
