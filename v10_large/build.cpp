#include "build.h"

#include "io.h"
#include "manifest.h"
#include "sort.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <unistd.h>

namespace hic10large {
namespace {

std::string pair_prefix(uint32_t a, uint32_t b) {
    std::ostringstream out;
    out << "pair-" << std::setfill('0') << std::setw(5) << a << '-' << std::setw(5) << b;
    return out.str();
}

std::string map_manifest_path(const std::string &directory, size_t shard) {
    std::ostringstream out;
    out << "map-" << std::setfill('0') << std::setw(8) << shard << ".manifest";
    return join_path(join_path(directory, "maps"), out.str());
}

std::string pair_manifest_path(const std::string &directory, uint32_t a, uint32_t b) {
    return join_path(join_path(directory, "pairs"), pair_prefix(a, b) + ".manifest");
}

std::string reduce_manifest_path(const std::string &directory, size_t pair, size_t group,
                                 size_t fan_in) {
    std::ostringstream out;
    out << "reduce-f" << std::setfill('0') << std::setw(5) << fan_in << "-p"
        << std::setw(5) << pair << "-g"
        << std::setw(5) << group << ".manifest";
    return join_path(join_path(directory, "reductions"), out.str());
}

size_t reduction_group_count(const PairInfo &pair, const BuildOptions &options) {
    return pair.shards.size() > options.merge_fan_in
        ? (pair.shards.size() + options.merge_fan_in - 1) / options.merge_fan_in : 0;
}

uint32_t parent_resolution(const std::vector<uint32_t> &completed, uint32_t target) {
    uint32_t parent = 0;
    for (uint32_t candidate : completed)
        if (target % candidate == 0)
            parent = std::max(parent, candidate);
    require(parent, "resolution " + std::to_string(target) +
                    " has no completed divisor; include the HBS source resolution");
    return parent;
}

std::vector<uint32_t> checked_resolutions(const StageManifest &stage,
                                          const BuildOptions &options) {
    require(!options.resolutions.empty(), "no build resolutions");
    auto resolutions = options.resolutions;
    std::sort(resolutions.begin(), resolutions.end());
    require(std::adjacent_find(resolutions.begin(), resolutions.end()) == resolutions.end(),
            "duplicate build resolution");
    require(resolutions.front() == stage.source_resolution,
            "the finest build resolution must equal the HBS source resolution");
    for (uint32_t resolution : resolutions)
        require(resolution && resolution % stage.source_resolution == 0,
                "build resolution must be a multiple of HBS resolution");
    return resolutions;
}

void write_map_manifest(const std::string &path, uint64_t fingerprint, size_t shard,
                        const std::vector<RunInfo> &runs) {
    std::string temporary = path + ".tmp-" + std::to_string(getpid());
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    require(bool(out), "cannot create map manifest");
    out << "HIC_V10_LARGE_MAP 1\nsource " << hex64(fingerprint) << "\nshard " << shard
        << "\nruns " << runs.size() << "\n";
    for (const auto &run : runs)
        out << "run " << run.chr1 << ' ' << run.chr2 << ' ' << run.resolution << ' '
            << run.records << ' ' << hex64(run.checksum) << ' ' << std::quoted(run.path) << "\n";
    out << "end\n";
    out.flush(); require(bool(out), "cannot write map manifest"); out.close();
    atomic_rename(temporary, path);
}

std::vector<RunInfo> read_map_manifest(const std::string &path, uint64_t fingerprint,
                                       size_t shard) {
    std::ifstream in(path, std::ios::binary);
    require(bool(in), "missing map task result " + path);
    std::string magic, tag, hash;
    unsigned version = 0;
    in >> magic >> version;
    require(magic == "HIC_V10_LARGE_MAP" && version == 1, "invalid map manifest");
    in >> tag >> hash;
    require(tag == "source" && unhex64(hash) == fingerprint, "map manifest source mismatch");
    size_t actual_shard = 0, count = 0;
    in >> tag;
    require(tag == "shard", "invalid map manifest");
    in >> actual_shard >> tag >> count;
    require(actual_shard == shard && tag == "runs", "map manifest task mismatch");
    std::vector<RunInfo> runs;
    for (size_t i = 0; i < count; ++i) {
        RunInfo run;
        std::string checksum;
        in >> tag >> run.chr1 >> run.chr2 >> run.resolution >> run.records >> checksum
           >> std::quoted(run.path);
        require(tag == "run", "invalid map run entry");
        run.checksum = unhex64(checksum);
        RunInfo actual = inspect_run(run.path, false);
        require(actual.chr1 == run.chr1 && actual.chr2 == run.chr2 &&
                    actual.resolution == run.resolution && actual.records == run.records &&
                    actual.checksum == run.checksum,
                "map run changed after task completion");
        runs.push_back(std::move(run));
    }
    in >> tag;
    require(tag == "end", "truncated map manifest");
    return runs;
}

void write_reduce_manifest(const std::string &path, uint64_t fingerprint, size_t pair,
                           size_t group, const RunInfo &run) {
    std::string temporary = path + ".tmp-" + std::to_string(getpid());
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    require(bool(out), "cannot create reduction manifest");
    out << "HIC_V10_LARGE_REDUCE 1\nsource " << hex64(fingerprint)
        << "\npair " << pair << "\ngroup " << group << "\nrun " << run.chr1 << ' '
        << run.chr2 << ' ' << run.resolution << ' ' << run.records << ' '
        << hex64(run.checksum) << ' ' << std::quoted(run.path) << "\nend\n";
    out.flush(); require(bool(out), "cannot write reduction manifest"); out.close();
    atomic_rename(temporary, path);
}

RunInfo read_reduce_manifest(const std::string &path, uint64_t fingerprint,
                             size_t pair, size_t group) {
    std::ifstream in(path, std::ios::binary);
    require(bool(in), "missing reduction task result " + path);
    std::string magic, tag, hash, checksum;
    unsigned version = 0;
    size_t actual_pair = 0, actual_group = 0;
    RunInfo run;
    in >> magic >> version;
    require(magic == "HIC_V10_LARGE_REDUCE" && version == 1,
            "invalid reduction manifest " + path);
    in >> tag >> hash;
    require(tag == "source" && unhex64(hash) == fingerprint,
            "reduction manifest source mismatch");
    in >> tag >> actual_pair;
    require(tag == "pair" && actual_pair == pair, "reduction manifest pair mismatch");
    in >> tag >> actual_group;
    require(tag == "group" && actual_group == group, "reduction manifest group mismatch");
    in >> tag >> run.chr1 >> run.chr2 >> run.resolution >> run.records >> checksum
       >> std::quoted(run.path);
    require(tag == "run", "invalid reduction run entry");
    run.checksum = unhex64(checksum);
    in >> tag;
    require(tag == "end", "truncated reduction manifest " + path);
    RunInfo actual = inspect_run(run.path, false);
    require(actual.chr1 == run.chr1 && actual.chr2 == run.chr2 &&
                actual.resolution == run.resolution && actual.records == run.records &&
                actual.checksum == run.checksum,
            "reduction run changed after task completion");
    return run;
}

void write_pair_manifest(const std::string &path, uint64_t fingerprint,
                         const std::vector<uint32_t> &resolutions,
                         const std::map<uint32_t, RunInfo> &cells) {
    std::string temporary = path + ".tmp-" + std::to_string(getpid());
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    require(bool(out), "cannot create pair manifest");
    out << "HIC_V10_LARGE_PAIR 1\nsource " << hex64(fingerprint) << "\nresolutions "
        << resolutions.size();
    for (auto resolution : resolutions) out << ' ' << resolution;
    out << "\n";
    for (const auto &entry : cells) {
        const RunInfo &run = entry.second;
        out << "cell " << run.chr1 << ' ' << run.chr2 << ' ' << run.resolution << ' '
            << run.records << ' ' << hex64(run.checksum) << ' ' << std::quoted(run.path) << "\n";
    }
    out << "end\n";
    out.flush(); require(bool(out), "cannot write pair manifest"); out.close();
    atomic_rename(temporary, path);
}

std::vector<RunInfo> read_pair_manifest(const std::string &path, uint64_t fingerprint,
                                        const std::vector<uint32_t> &resolutions) {
    std::ifstream in(path, std::ios::binary);
    require(bool(in), "missing pair task result " + path);
    std::string magic, tag, hash;
    unsigned version = 0;
    in >> magic >> version >> tag >> hash;
    require(magic == "HIC_V10_LARGE_PAIR" && version == 1 && tag == "source" &&
                unhex64(hash) == fingerprint,
            "invalid pair manifest");
    size_t count = 0;
    in >> tag >> count;
    require(tag == "resolutions" && count == resolutions.size(), "pair resolution list mismatch");
    for (uint32_t resolution : resolutions) {
        uint32_t actual = 0; in >> actual;
        require(actual == resolution, "pair resolution list mismatch");
    }
    std::vector<RunInfo> result;
    while (in >> tag && tag == "cell") {
        RunInfo run;
        std::string checksum;
        in >> run.chr1 >> run.chr2 >> run.resolution >> run.records >> checksum >> std::quoted(run.path);
        run.checksum = unhex64(checksum);
        RunInfo actual = inspect_run(run.path, false);
        require(actual.records == run.records && actual.checksum == run.checksum &&
                    actual.resolution == run.resolution,
                "pair cell artifact changed after completion");
        result.push_back(std::move(run));
    }
    require(tag == "end" && result.size() == resolutions.size(), "truncated pair manifest");
    return result;
}

void cleanup_pair_inputs(const StageManifest &stage, const std::string &directory,
                         size_t pair_index, const BuildOptions &options) {
    const PairInfo &pair = stage.pairs.at(pair_index);
    for (size_t shard : pair.shards) {
        std::string path = map_manifest_path(directory, shard);
        if (!path_exists(path)) continue;
        auto mapped = read_map_manifest(path, stage.source_fingerprint, shard);
        for (const auto &run : mapped) std::remove(run.path.c_str());
        std::remove(path.c_str());
    }
    for (size_t group = 0; group < reduction_group_count(pair, options); ++group) {
        std::string path = reduce_manifest_path(directory, pair_index, group,
                                                options.merge_fan_in);
        if (!path_exists(path)) continue;
        RunInfo reduced = read_reduce_manifest(path, stage.source_fingerprint,
                                               pair_index, group);
        std::remove(reduced.path.c_str());
        std::remove(path.c_str());
    }
}

} // namespace

void build_cells(const std::string &stage_path, const std::string &directory,
                 const BuildOptions &options) {
    StageManifest stage = read_stage_manifest(stage_path);
    require(options.merge_fan_in >= 2, "merge fan-in must be at least two");
    checked_resolutions(stage, options);
    make_directory(directory);
    make_directory(join_path(directory, "runs"));
    make_directory(join_path(directory, "cells"));
    make_directory(join_path(directory, "maps"));
    make_directory(join_path(directory, "pairs"));
    make_directory(join_path(directory, "reductions"));
    std::set<std::pair<uint32_t, uint32_t>> completed_pairs;
    for (const auto &pair : stage.pairs)
        if (path_exists(pair_manifest_path(directory, pair.chr1, pair.chr2)))
            completed_pairs.insert({pair.chr1, pair.chr2});
    for (size_t shard = 0; shard < stage.shards.size(); ++shard)
        if (!completed_pairs.count({stage.shards[shard].chr1, stage.shards[shard].chr2}) &&
            !path_exists(map_manifest_path(directory, shard)))
            map_root_shard(stage_path, directory, shard, options);
    for (size_t pair = 0; pair < stage.pairs.size(); ++pair) {
        if (completed_pairs.count({stage.pairs[pair].chr1, stage.pairs[pair].chr2})) continue;
        for (size_t group = 0; group < reduction_group_count(stage.pairs[pair], options); ++group)
            if (!path_exists(reduce_manifest_path(directory, pair, group, options.merge_fan_in)))
                reduce_root_group(stage_path, directory, pair, group, options);
    }
    for (size_t pair = 0; pair < stage.pairs.size(); ++pair)
        if (!path_exists(pair_manifest_path(directory, stage.pairs[pair].chr1, stage.pairs[pair].chr2)))
            build_pair_cells(stage_path, directory, pair, options);
    finalize_build(stage_path, directory, options);
}

void map_root_shard(const std::string &stage_path, const std::string &directory,
                    size_t shard_index, const BuildOptions &options) {
    StageManifest stage = read_stage_manifest(stage_path);
    auto resolutions = checked_resolutions(stage, options);
    require(shard_index < stage.shards.size(), "map shard index outside manifest");
    make_directory(directory); make_directory(join_path(directory, "runs"));
    make_directory(join_path(directory, "maps"));
    std::string result_path = map_manifest_path(directory, shard_index);
    if (path_exists(result_path)) {
        read_map_manifest(result_path, stage.source_fingerprint, shard_index);
        std::cerr << "Map shard " << shard_index << " already complete\n";
        return;
    }
    std::string stage_directory = stage_path.substr(0, stage_path.find_last_of('/'));
    if (stage_directory == stage_path) stage_directory = ".";
    const ShardInfo &shard = stage.shards[shard_index];
    if (path_exists(pair_manifest_path(directory, shard.chr1, shard.chr2))) {
        std::cerr << "Map shard " << shard_index << " is already consumed by a completed pair\n";
        return;
    }
    std::string prefix = pair_prefix(shard.chr1, shard.chr2) + "-root-s" +
                         std::to_string(shard.part);
    auto runs = sort_staged_shard(join_path(stage_directory, shard.path),
                                  join_path(directory, "runs"), resolutions.front(),
                                  options.memory_bytes, prefix);
    // Map-side combine: publish one sorted run per staged shard even when a
    // small RAM budget forced several radix-sort spills inside the task.
    if (runs.size() > 1) {
        RunInfo combined = bounded_merge(std::move(runs), join_path(directory, "runs"),
                                         prefix + "-map", options.merge_fan_in, {}, true);
        runs = {std::move(combined)};
    }
    write_map_manifest(result_path, stage.source_fingerprint, shard_index, runs);
}

void reduce_root_group(const std::string &stage_path, const std::string &directory,
                       size_t pair_index, size_t group_index, const BuildOptions &options) {
    StageManifest stage = read_stage_manifest(stage_path);
    auto resolutions = checked_resolutions(stage, options);
    require(options.merge_fan_in >= 2, "merge fan-in must be at least two");
    require(pair_index < stage.pairs.size(), "reduction pair index outside manifest");
    const PairInfo &pair = stage.pairs[pair_index];
    const size_t groups = reduction_group_count(pair, options);
    require(group_index < groups, "reduction group index outside task graph");
    make_directory(directory); make_directory(join_path(directory, "runs"));
    make_directory(join_path(directory, "reductions"));
    const std::string result_path = reduce_manifest_path(directory, pair_index, group_index,
                                                         options.merge_fan_in);
    if (path_exists(result_path)) {
        read_reduce_manifest(result_path, stage.source_fingerprint, pair_index, group_index);
        std::cerr << "Reduction task " << pair_index << ':' << group_index
                  << " already complete\n";
        return;
    }
    if (path_exists(pair_manifest_path(directory, pair.chr1, pair.chr2))) {
        std::cerr << "Reduction task " << pair_index << ':' << group_index
                  << " is superseded by a completed pair\n";
        return;
    }
    const size_t begin = group_index * options.merge_fan_in;
    const size_t end = std::min(pair.shards.size(), begin + options.merge_fan_in);
    std::vector<RunInfo> inputs;
    inputs.reserve(end - begin);
    for (size_t at = begin; at < end; ++at) {
        const size_t shard = pair.shards[at];
        auto runs = read_map_manifest(map_manifest_path(directory, shard),
                                      stage.source_fingerprint, shard);
        inputs.insert(inputs.end(), runs.begin(), runs.end());
    }
    std::string prefix = pair_prefix(pair.chr1, pair.chr2) + "-root-reduce-g" +
                         std::to_string(group_index) + "-f" +
                         std::to_string(options.merge_fan_in);
    RunInfo reduced = bounded_merge(std::move(inputs), join_path(directory, "runs"), prefix,
                                    options.merge_fan_in);
    require(reduced.resolution == resolutions.front(), "root reduction resolution mismatch");
    write_reduce_manifest(result_path, stage.source_fingerprint, pair_index, group_index, reduced);
}

void build_pair_cells(const std::string &stage_path, const std::string &directory,
                      size_t pair_index, const BuildOptions &options) {
    StageManifest stage = read_stage_manifest(stage_path);
    auto resolutions = checked_resolutions(stage, options);
    require(pair_index < stage.pairs.size(), "pair task index outside manifest");
    make_directory(directory); make_directory(join_path(directory, "runs"));
    make_directory(join_path(directory, "cells")); make_directory(join_path(directory, "pairs"));
    make_directory(join_path(directory, "reductions"));
    const PairInfo &pair = stage.pairs[pair_index];
    std::string result_path = pair_manifest_path(directory, pair.chr1, pair.chr2);
    if (path_exists(result_path)) {
        read_pair_manifest(result_path, stage.source_fingerprint, resolutions);
        cleanup_pair_inputs(stage, directory, pair_index, options);
        std::cerr << "Pair task " << pair_index << " already complete\n";
        return;
    }
    const std::string prefix = pair_prefix(pair.chr1, pair.chr2);
    std::map<uint32_t, RunInfo> cells;
    std::vector<uint32_t> completed;
    std::cerr << "Building cells for " << stage.chromosomes[pair.chr1].name << " x "
              << stage.chromosomes[pair.chr2].name << " (" << pair.records << " records)\n";
    for (uint32_t resolution : resolutions) {
        std::string resolution_prefix = prefix + "-r" + std::to_string(resolution);
        RunInfo merged;
        if (resolution == resolutions.front()) {
            std::vector<RunInfo> runs;
            const size_t groups = reduction_group_count(pair, options);
            if (groups) {
                for (size_t group = 0; group < groups; ++group)
                    runs.push_back(read_reduce_manifest(
                        reduce_manifest_path(directory, pair_index, group, options.merge_fan_in),
                        stage.source_fingerprint, pair_index, group));
            } else {
                for (size_t shard_index : pair.shards) {
                    auto mapped = read_map_manifest(map_manifest_path(directory, shard_index),
                                                    stage.source_fingerprint, shard_index);
                    runs.insert(runs.end(), mapped.begin(), mapped.end());
                }
            }
            merged = bounded_merge(std::move(runs), join_path(directory, "cells"),
                                   resolution_prefix + "-final", options.merge_fan_in);
        } else {
            uint32_t parent = parent_resolution(completed, resolution);
            // The parent is already row-sorted. Aggregate one coarsened row at
            // a time, avoiding another external sort or raw-spool scan.
            merged = rollup_cell_file(cells.at(parent), join_path(directory, "cells"),
                                      resolution, resolution_prefix + "-final");
        }
        cells.emplace(resolution, merged);
        completed.push_back(resolution);
        std::cerr << "  " << resolution << " bp: " << merged.records << " occupied cells\n";
    }
    write_pair_manifest(result_path, stage.source_fingerprint, resolutions, cells);
    // Once the durable pair manifest exists, root map runs are superseded by
    // the canonical finest-resolution cell file and can be reclaimed safely.
    cleanup_pair_inputs(stage, directory, pair_index, options);
}

void finalize_build(const std::string &stage_path, const std::string &directory,
                    const BuildOptions &options) {
    StageManifest stage = read_stage_manifest(stage_path);
    auto resolutions = checked_resolutions(stage, options);
    std::string temporary = join_path(directory, "build.manifest.tmp-" + std::to_string(getpid()));
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    require(bool(out), "cannot create build manifest");
    out << "HIC_V10_LARGE_BUILD 1\nstage " << std::quoted(stage_path) << ' '
        << hex64(stage.source_fingerprint) << "\nresolutions " << resolutions.size();
    for (auto resolution : resolutions) out << ' ' << resolution;
    out << "\n";
    for (const PairInfo &pair : stage.pairs) {
        auto cells = read_pair_manifest(pair_manifest_path(directory, pair.chr1, pair.chr2),
                                        stage.source_fingerprint, resolutions);
        for (const RunInfo &run : cells)
            out << "cell " << run.chr1 << ' ' << run.chr2 << ' ' << run.resolution << ' '
                << run.records << ' ' << hex64(run.checksum) << ' ' << std::quoted(run.path) << "\n";
    }
    out << "end\n"; out.flush(); require(bool(out), "cannot write build manifest"); out.close();
    atomic_rename(temporary, join_path(directory, "build.manifest"));
    for (size_t pair = 0; pair < stage.pairs.size(); ++pair)
        cleanup_pair_inputs(stage, directory, pair, options);
}

void print_build_tasks(const std::string &stage_path, const std::string &directory,
                       const BuildOptions &options) {
    StageManifest stage = read_stage_manifest(stage_path);
    checked_resolutions(stage, options);
    for (size_t shard = 0; shard < stage.shards.size(); ++shard)
        std::cout << "map-root\t" << shard << '\t' << map_manifest_path(directory, shard) << "\n";
    for (size_t pair = 0; pair < stage.pairs.size(); ++pair)
        for (size_t group = 0; group < reduction_group_count(stage.pairs[pair], options); ++group)
            std::cout << "reduce-root\t" << pair << ':' << group << '\t'
                      << reduce_manifest_path(directory, pair, group, options.merge_fan_in) << "\n";
    for (size_t pair = 0; pair < stage.pairs.size(); ++pair)
        std::cout << "build-pair\t" << pair << '\t'
                  << pair_manifest_path(directory, stage.pairs[pair].chr1, stage.pairs[pair].chr2) << "\n";
    std::cout << "finalize-build\t0\t" << join_path(directory, "build.manifest") << "\n";
}

BuildManifest read_build_manifest(const std::string &path, bool inspect_files) {
    std::ifstream in(path, std::ios::binary);
    require(bool(in), "cannot open build manifest " + path);
    std::string magic;
    unsigned version = 0;
    in >> magic >> version;
    require(magic == "HIC_V10_LARGE_BUILD" && version == 1, "invalid build manifest");
    BuildManifest result;
    std::string tag;
    size_t expected_resolutions = 0;
    std::string directory = path.substr(0, path.find_last_of('/'));
    if (directory == path) directory = ".";
    while (in >> tag) {
        if (tag == "stage") {
            std::string fingerprint;
            in >> std::quoted(result.stage_manifest) >> fingerprint;
            result.source_fingerprint = unhex64(fingerprint);
        } else if (tag == "resolutions") {
            in >> expected_resolutions;
            result.resolutions.resize(expected_resolutions);
            for (auto &resolution : result.resolutions) in >> resolution;
        } else if (tag == "cell") {
            uint32_t a = 0, b = 0, resolution = 0;
            uint64_t records = 0;
            std::string checksum, cell_path;
            in >> a >> b >> resolution >> records >> checksum >> std::quoted(cell_path);
            RunInfo info{cell_path, a, b, resolution, records, unhex64(checksum)};
            if (inspect_files) {
                RunInfo actual = inspect_run(cell_path, false);
                require(actual.chr1 == a && actual.chr2 == b &&
                            actual.resolution == resolution && actual.records == records &&
                            actual.checksum == info.checksum,
                        "cell file does not match build manifest: " + cell_path);
            }
            require(result.cells.emplace(std::make_tuple(a, b, resolution), info).second,
                    "duplicate cell entry in build manifest");
        } else if (tag == "end") {
            break;
        } else {
            fail("unknown build manifest field " + tag);
        }
        require(bool(in), "truncated build manifest " + path);
    }
    require(!result.stage_manifest.empty() && result.resolutions.size() == expected_resolutions &&
                !result.cells.empty(),
            "incomplete build manifest " + path);
    return result;
}

} // namespace hic10large
