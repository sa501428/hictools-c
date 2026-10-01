#include "addnorm.h"
#include "build.h"
#include "manifest.h"
#include "v10/reader.h"
#include "v10/vector_import.h"
#include "v10/text_vectors.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <unistd.h>

namespace hic10large {
namespace {
struct Workspace {
    std::string path;
    explicit Workspace(const std::string &parent) {
        make_directory(parent);
        std::string pattern = join_path(parent, "hic-v10-addnorm-XXXXXX");
        std::vector<char> name(pattern.begin(), pattern.end()); name.push_back(0);
        require(mkdtemp(name.data()) != nullptr, "cannot create addnorm workspace");
        path = name.data();
    }
    ~Workspace() { std::error_code error; std::filesystem::remove_all(path, error); }
};
RunInfo extract_cis_run(hic10::Reader &reader, uint32_t chr, uint8_t unit, uint32_t ri,
                        const std::string &directory, const NormalizeOptions &options) {
    auto r = reader.header().resolutions[unit][ri];
    std::cerr << "Extracting " << reader.header().chromosomes[chr].name << " at " << r.bin << (unit ? " FRAG for addnorm\n" : " bp for addnorm\n");
    std::vector<CellRecord> buffer;
    const uint64_t capacity = options.memory_bytes / (2 * sizeof(CellRecord));
    std::vector<RunInfo> runs;
    auto flush = [&]() {
        if (buffer.empty()) return;
        radix_sort_cells(buffer);
        auto name = "u" + std::to_string(unit) + "-c" + std::to_string(chr) + "-r" + std::to_string(ri) + "-part" + std::to_string(runs.size());
        RunWriter output(join_path(directory, name + ".h10r"), chr, chr, r.bin);
        for (auto cell : buffer) output.add(cell);
        runs.push_back(output.finish()); buffer.clear();
    };
    reader.stream_materialized(chr, chr, unit, ri, [&](const hic10::Cell &cell) {
        buffer.push_back({cell.x, cell.y, cell.value});
        if (buffer.size() >= capacity) flush();
    });
    flush();
    if (runs.empty()) return {};
    auto merged = bounded_merge(std::move(runs), directory,
        "u" + std::to_string(unit) + "-c" + std::to_string(chr) + "-r" + std::to_string(ri),
        options.scale_options.merge_fan_in, {}, true);
    return merged;
}
} // namespace
void add_text_norm_file(const std::string &path, const std::string &text,
                        const NormalizeOptions &options, int level) {
    require(options.memory_bytes >= 32768, "addnorm sort memory must be at least 32 KiB");
    require(options.scale_options.merge_fan_in >= 2, "merge fan-in must be at least 2");
    Workspace work(options.temporary_directory.empty() ? "/tmp" : options.temporary_directory);
    using Key = std::tuple<uint32_t, uint8_t, uint32_t>;
    std::map<Key, RunInfo> cache;
    hic10::ExpectedCellStream stream = [&](hic10::Reader &reader, uint32_t chr,
                                          uint8_t unit, uint32_t ri,
                                          const hic10::ExpectedCellConsumer &emit) {
        auto resolution = reader.header().resolutions[unit][ri];
        if (!resolution.mode) {
            reader.stream_materialized(chr, chr, unit, ri, [&](const hic10::Cell &cell) {
                emit(cell.x, cell.y, static_cast<long double>(cell.value));
            });
            return;
        }
        auto target = cache.find({chr, unit, ri});
        if (target == cache.end()) {
            auto source = cache.find({chr, unit, resolution.source});
            if (source == cache.end())
                source = cache.emplace(Key{chr, unit, resolution.source},
                    extract_cis_run(reader, chr, unit, resolution.source, work.path, options)).first;
            RunInfo run;
            if (!source->second.path.empty())
                run = rollup_cell_file(source->second, work.path, resolution.bin,
                    "u" + std::to_string(unit) + "-c" + std::to_string(chr) + "-r" + std::to_string(ri));
            target = cache.emplace(Key{chr, unit, ri}, std::move(run)).first;
        }
        if (target->second.path.empty()) return;
        RunReader cells(target->second);
        CellRecord cell;
        while (cells.next(cell)) emit(cell.x, cell.y, static_cast<long double>(cell.count));
    };
    hic10::add_text_vectors_v10(path, text, level, work.path, stream);
}
void add_norm_file(const std::string &path, const NormalizeOptions &options, int level) {
    require(options.memory_bytes >= 32768, "addnorm sort memory must be at least 32 KiB");
    require(options.scale_options.merge_fan_in >= 2, "merge fan-in must be at least 2");
    hic10::Reader reader(path);
    const auto &header = reader.header();
    require(header.resolutions[1].empty(), "large addnorm supports BP count files; use hic_v10 addnorm for FRAG");
    require(!header.resolutions[0].empty(), "no BP resolutions to normalize");
    auto existing = reader.vector_entries();
    bool work_needed = false;
    for (uint32_t ri = 0; ri < header.resolutions[0].size(); ++ri) {
        bool raw_present = false;
        for (const auto &e : existing)
            if (e.kind == 1 && e.unit == 0 && e.ri == ri) raw_present = true;
        work_needed |= !raw_present;
        for (const auto &item : std::vector<std::pair<std::string, bool>>{
                 {"VC", options.vc}, {"VC_SQRT", options.vc_sqrt}, {"SCALE", options.scale}}) {
            if (!item.second) continue;
            auto found = std::find(header.norms.begin(), header.norms.end(), item.first);
            bool present = false;
            if (found != header.norms.end())
                for (const auto &e : existing)
                    if (e.kind != 1 && e.unit == 0 && e.ri == ri &&
                        e.norm == uint32_t(found - header.norms.begin())) present = true;
            work_needed |= !present;
        }
    }
    if (!work_needed) {
        std::cerr << "Requested V10 normalization bundles already exist: " << path << '\n';
        return;
    }
    Workspace work(options.temporary_directory.empty() ? "/tmp" : options.temporary_directory);
    StageManifest stage;
    stage.source = std::filesystem::absolute(path).string();
    stage.source_bytes = reader.file_size();
    auto header_bytes = reader.read_bytes(0, reader.header_length());
    stage.source_fingerprint = fnv1a(header_bytes.data(), header_bytes.size());
    stage.source_resolution = header.resolutions[0].front().bin;
    for (const auto &chr : header.chromosomes) stage.chromosomes.push_back({chr.name, chr.length});
    const auto stage_path = join_path(work.path, "stage.manifest");
    write_stage_manifest(stage, stage_path);
    const auto build_path = join_path(work.path, "build.manifest");
    std::ofstream build(build_path);
    require(bool(build), "cannot create addnorm build manifest");
    build << "HIC_V10_LARGE_BUILD 1\nstage " << std::quoted(stage_path) << ' '
          << hex64(stage.source_fingerprint) << "\nresolutions " << header.resolutions[0].size();
    for (auto r : header.resolutions[0]) build << ' ' << r.bin;
    build << "\nstorage all\n";
    bool have_cells = false;
    for (auto key : reader.matrices()) {
        if (key.chr1 != key.chr2) continue;
        const auto &name = header.chromosomes[key.chr1].name;
        if (name == "ALL" || name == "All" || name == "all") continue;
        std::map<uint32_t, RunInfo> runs_by_index;
        for (uint32_t ri = 0; ri < header.resolutions[0].size(); ++ri) {
            auto r = header.resolutions[0][ri];
            if (r.mode) continue;
            auto merged = extract_cis_run(reader, key.chr1, 0, ri, work.path, options);
            if (!merged.path.empty()) runs_by_index.emplace(ri, std::move(merged));
        }
        for (uint32_t ri = 0; ri < header.resolutions[0].size(); ++ri) {
            auto r = header.resolutions[0][ri];
            if (!r.mode) continue;
            auto source = runs_by_index.find(r.source);
            if (source == runs_by_index.end()) continue;
            runs_by_index.emplace(ri, rollup_cell_file(source->second, work.path, r.bin,
                "c" + std::to_string(key.chr1) + "-r" + std::to_string(ri)));
        }
        for (auto &entry : runs_by_index) {
            const auto &run = entry.second;
            build << "cell " << run.chr1 << ' ' << run.chr2 << ' ' << run.resolution << ' '
                  << run.records << ' ' << hex64(run.checksum) << ' ' << std::quoted(run.path) << '\n';
            have_cells = true;
        }
    }
    build << "end\n"; build.flush(); require(bool(build), "cannot write addnorm cell manifest"); build.close();
    require(have_cells, "no cis counts to normalize");
    NormalizeOptions local = options;
    local.temporary_directory = join_path(work.path, "normalize-work");
    const auto vectors = join_path(work.path, "vectors");
    normalize_cells(stage_path, build_path, vectors, local);
    hic10::add_vector_manifest_v10(path, join_path(vectors, "vectors.manifest"), level);
    std::cerr << "V10 normalization added: " << path << '\n';
}
}
