#include "normalize.h"

#include "build.h"
#include "io.h"
#include "manifest.h"
#include "sort.h"
#include "vectors.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <tuple>

namespace hic10large {
namespace {

using U128 = unsigned __int128;
constexpr uint32_t NAN_WORD = 0x7fc00000;

U128 checked_add(U128 a, uint64_t b) {
    U128 maximum = ~U128(0);
    require(a <= maximum - b, "uint128 contact total overflow");
    return a + b;
}

U128 checked_add(U128 a, U128 b) {
    U128 maximum = ~U128(0);
    require(a <= maximum - b, "uint128 contact total overflow");
    return a + b;
}

long double as_long_double(U128 value) {
    constexpr U128 low_mask = (U128(1) << 64) - 1;
    uint64_t low = static_cast<uint64_t>(value & low_mask);
    uint64_t high = static_cast<uint64_t>(value >> 64);
    return std::ldexp(static_cast<long double>(high), 64) + low;
}

std::string vector_name(uint8_t kind, uint32_t norm, uint32_t chr, uint32_t resolution) {
    std::ostringstream out;
    out << "vector-k" << unsigned(kind) << "-n" << norm << "-c" << chr
        << "-r" << resolution << ".h10w";
    return out.str();
}

std::string chromosome_manifest_path(const std::string &directory, uint32_t chromosome) {
    std::ostringstream out;
    out << "chr-" << std::setfill('0') << std::setw(5) << chromosome << ".manifest";
    return join_path(directory, out.str());
}

std::string expected_manifest_path(const std::string &directory, uint32_t resolution_index) {
    std::ostringstream out;
    out << "expected-r" << std::setfill('0') << std::setw(8) << resolution_index << ".manifest";
    return join_path(directory, out.str());
}

VectorManifest configured_manifest(uint64_t source_fingerprint, const NormalizeOptions &options) {
    VectorManifest result;
    result.source_fingerprint = source_fingerprint;
    if (options.vc) result.norms.push_back("VC");
    if (options.vc_sqrt) result.norms.push_back("VC_SQRT");
    if (options.scale) result.norms.push_back("SCALE");
    return result;
}

uint32_t norm_id(const VectorManifest &manifest, const std::string &name) {
    auto found = std::find(manifest.norms.begin(), manifest.norms.end(), name);
    return found == manifest.norms.end() ? NO_ID : static_cast<uint32_t>(found - manifest.norms.begin());
}

void validate_subset(const VectorManifest &subset, const VectorManifest &configured,
                     const std::string &path) {
    require(subset.source_fingerprint == configured.source_fingerprint,
            "normalization task source mismatch: " + path);
    require(subset.norms == configured.norms,
            "normalization task options mismatch: " + path);
}

uint32_t chromosome_bins(uint64_t length, uint32_t resolution) {
    uint64_t bins = length / resolution + (length % resolution != 0);
    require(bins <= UINT32_MAX, "chromosome has more than uint32 bins at " +
                               std::to_string(resolution) + " bp");
    return static_cast<uint32_t>(bins);
}

bool valid_norm(const std::vector<float> &values) {
    for (float value : values)
        if (value > 0 && std::isfinite(value)) return true;
    return false;
}

float rescale_norm(const RunInfo &cells, std::vector<float> &norm) {
    RunReader reader(cells);
    CellRecord cell;
    long double raw = 0, normalized = 0;
    while (reader.next(cell)) {
        float a = norm[cell.x], b = norm[cell.y];
        if (!(a > 0) || !(b > 0) || !std::isfinite(a) || !std::isfinite(b)) continue;
        long double multiple = cell.x == cell.y ? 1 : 2;
        raw += multiple * cell.count;
        normalized += multiple * cell.count / (static_cast<long double>(a) * b);
    }
    if (!(raw > 0) || !(normalized > 0)) return 0;
    float factor = static_cast<float>(std::sqrt(normalized / raw));
    for (float &value : norm)
        if (value > 0 && std::isfinite(value)) value *= factor;
    return factor;
}

std::vector<uint32_t> length_order(const std::map<uint32_t, uint32_t> &bins,
                                   const std::map<uint32_t, long double> &observed) {
    std::vector<uint32_t> lengths;
    for (const auto &entry : bins)
        if (observed.count(entry.first) && observed.at(entry.first) > 0)
            lengths.push_back(entry.second);
    std::sort(lengths.begin(), lengths.end());
    return lengths;
}

struct Possible {
    std::vector<uint32_t> lengths;
    std::vector<uint64_t> suffix;
    explicit Possible(std::vector<uint32_t> input) : lengths(std::move(input)), suffix(lengths.size() + 1) {
        for (size_t i = lengths.size(); i-- > 0;) suffix[i] = suffix[i + 1] + lengths[i];
    }
    long double operator()(uint64_t distance) const {
        size_t first = std::upper_bound(lengths.begin(), lengths.end(), distance) - lengths.begin();
        return static_cast<long double>(suffix[first]) -
               static_cast<long double>(lengths.size() - first) * distance;
    }
    uint64_t support() const { return lengths.empty() ? 0 : lengths.back(); }
};

struct RawExpected {
    explicit RawExpected(uint32_t size) : actual(size) {}
    std::vector<U128> actual;
    std::map<uint32_t, U128> observed;
    std::map<uint32_t, uint32_t> bins;
    void add(uint32_t chr, const RunInfo &cells, uint32_t chromosome_bins) {
        bins[chr] = chromosome_bins;
        U128 total = 0;
        RunReader reader(cells);
        CellRecord cell;
        while (reader.next(cell)) {
            uint32_t distance = cell.y - cell.x;
            actual[distance] = checked_add(actual[distance], cell.count);
            total = checked_add(total, cell.count);
        }
        if (total) observed[chr] = checked_add(observed[chr], total);
    }
    bool has_data() const { return !observed.empty(); }
    VectorInfo finish(uint32_t ri, uint32_t resolution, const std::string &directory) {
        std::map<uint32_t, long double> observed_float;
        for (const auto &entry : observed) observed_float[entry.first] = as_long_double(entry.second);
        Possible possible(length_order(bins, observed_float));
        VectorInfo info;
        info.kind = 1; info.norm = NO_ID; info.chr = NO_ID;
        info.ri = ri; info.resolution = resolution;
        const uint64_t support = possible.support();
        U128 numerator = support ? actual[0] : 0;
        long double denominator = support ? possible(0) : 0;
        uint64_t low = 0, high = 0;
        std::vector<uint32_t> words;
        words.reserve(actual.size());
        for (uint64_t distance = 0; distance < support; ++distance) {
            if (numerator < 400) {
                while (numerator < 400 && high + 1 < support) {
                    ++high;
                    numerator = checked_add(numerator, actual[high]);
                    denominator += possible(high);
                }
            } else {
                while (high > low) {
                    U128 endpoints = checked_add(actual[low], actual[high]);
                    require(endpoints <= numerator, "invalid raw expected smoothing window");
                    if (numerator - endpoints < 400) break;
                    numerator -= endpoints;
                    denominator -= possible(low) + possible(high);
                    ++low; --high;
                }
            }
            float value = denominator > 0 ? static_cast<float>(as_long_double(numerator) / denominator) : 0;
            words.push_back(float_bits(value));
            if (high + 2 < support) {
                numerator = checked_add(numerator, actual[high + 1]);
                numerator = checked_add(numerator, actual[high + 2]);
                denominator += possible(high + 1) + possible(high + 2);
                high += 2;
            } else if (high + 1 < support) {
                ++high; numerator = checked_add(numerator, actual[high]);
                denominator += possible(high);
            }
        }
        words.resize(actual.size(), NAN_WORD);
        for (const auto &entry : observed) {
            long double expected = 0;
            uint32_t n = bins[entry.first];
            for (uint32_t d = 0; d < n && d < words.size(); ++d)
                expected += static_cast<long double>(n - d) * bits_float(words[d]);
            if (expected > 0 && entry.second)
                info.scales[entry.first] = float_bits(static_cast<float>(expected / as_long_double(entry.second)));
        }
        VectorFileWriter output(join_path(directory, vector_name(1, NO_ID, NO_ID, resolution)), info);
        output.add(words);
        return output.finish();
    }
};

struct NormalizedExpected {
    explicit NormalizedExpected(uint32_t size) : actual(size, 0) {}
    std::vector<long double> actual;
    std::map<uint32_t, long double> observed;
    std::map<uint32_t, uint32_t> bins;
    void add(uint32_t chr, const RunInfo &cells, uint32_t chromosome_bins,
             const std::vector<float> &norm) {
        bins[chr] = chromosome_bins;
        long double total = 0;
        RunReader reader(cells);
        CellRecord cell;
        while (reader.next(cell)) {
            float a = norm[cell.x], b = norm[cell.y];
            if (!(a > 0) || !(b > 0) || !std::isfinite(a) || !std::isfinite(b)) continue;
            long double value = static_cast<long double>(cell.count) /
                                (static_cast<long double>(a) * b);
            actual[cell.y - cell.x] += value;
            total += value;
        }
        if (total > 0) observed[chr] += total;
    }
    bool has_data() const { return !observed.empty(); }
    VectorInfo finish(uint32_t norm, uint32_t ri, uint32_t resolution,
                      const std::string &directory) {
        Possible possible(length_order(bins, observed));
        const uint64_t support = possible.support();
        long double numerator = support ? actual[0] : 0;
        long double denominator = support ? possible(0) : 0;
        uint64_t low = 0, high = 0;
        std::vector<uint32_t> words;
        words.reserve(actual.size());
        for (uint64_t distance = 0; distance < support; ++distance) {
            if (numerator < 400) {
                while (numerator < 400 && high + 1 < support) {
                    ++high; numerator += actual[high]; denominator += possible(high);
                }
            } else {
                while (high > low && numerator - actual[low] - actual[high] >= 400) {
                    numerator -= actual[low] + actual[high];
                    denominator -= possible(low) + possible(high);
                    ++low; --high;
                }
            }
            words.push_back(float_bits(denominator > 0 ? static_cast<float>(numerator / denominator) : 0));
            if (high + 2 < support) {
                numerator += actual[high + 1] + actual[high + 2];
                denominator += possible(high + 1) + possible(high + 2); high += 2;
            } else if (high + 1 < support) {
                ++high; numerator += actual[high]; denominator += possible(high);
            }
        }
        words.resize(actual.size(), NAN_WORD);
        VectorInfo info;
        info.kind = 2; info.norm = norm; info.chr = NO_ID;
        info.ri = ri; info.resolution = resolution;
        for (const auto &entry : observed) {
            long double expected = 0;
            uint32_t n = bins[entry.first];
            for (uint32_t d = 0; d < n && d < words.size(); ++d)
                expected += static_cast<long double>(n - d) * bits_float(words[d]);
            if (expected > 0 && entry.second > 0)
                info.scales[entry.first] = float_bits(static_cast<float>(expected / entry.second));
        }
        VectorFileWriter output(join_path(directory, vector_name(2, norm, NO_ID, resolution)), info);
        output.add(words);
        return output.finish();
    }
};

VectorInfo write_norm(uint32_t norm, uint32_t chr, uint32_t ri, uint32_t resolution,
                      const std::vector<float> &values, const std::string &directory) {
    VectorInfo info;
    info.kind = 0; info.norm = norm; info.chr = chr;
    info.ri = ri; info.resolution = resolution;
    VectorFileWriter output(join_path(directory, vector_name(0, norm, chr, resolution)), info);
    for (float value : values) output.add(float_bits(value));
    return output.finish();
}

} // namespace

void normalize_chromosome(const std::string &stage_path, const std::string &build_path,
                          const std::string &directory, uint32_t chr,
                          const NormalizeOptions &options) {
    StageManifest stage = read_stage_manifest(stage_path);
    BuildManifest build = read_build_manifest(build_path);
    require(stage.source_fingerprint == build.source_fingerprint,
            "stage and build manifests refer to different inputs");
    require(chr < stage.chromosomes.size(), "chromosome task index outside manifest");
    require(options.memory_bytes >= 1024 * sizeof(CellRecord),
            "normalization memory budget is too small");
    make_directory(directory);
    const std::string result_path = chromosome_manifest_path(directory, chr);
    VectorManifest manifest = configured_manifest(stage.source_fingerprint, options);
    if (path_exists(result_path)) {
        validate_subset(read_vector_manifest(result_path), manifest, result_path);
        std::cerr << "Chromosome normalization task " << chr << " already complete\n";
        return;
    }
    std::string temporary_root = options.temporary_directory.empty()
        ? join_path(directory, "scale-work") : options.temporary_directory;
    std::string temporary = join_path(temporary_root, "chr-" + std::to_string(chr));
    make_directory(temporary_root);
    make_directory(temporary);
    const uint32_t vc_id = norm_id(manifest, "VC");
    const uint32_t vcs_id = norm_id(manifest, "VC_SQRT");
    const uint32_t scale_id = norm_id(manifest, "SCALE");
    bool scale_active = options.scale;

    // SCALE failure propagates toward finer bins, so resolutions are processed
    // from largest bin size to smallest.
    for (size_t reverse = build.resolutions.size(); reverse-- > 0;) {
        const uint32_t ri = static_cast<uint32_t>(reverse);
        const uint32_t resolution = build.resolutions[ri];
        auto found = build.cells.find(std::make_tuple(chr, chr, resolution));
        if (found == build.cells.end()) continue;
        const RunInfo &cells = found->second;
        const uint32_t bins = chromosome_bins(stage.chromosomes[chr].length, resolution);
        std::cerr << "Normalizing " << stage.chromosomes[chr].name << " at "
                  << resolution << " bp\n";
        std::vector<U128> coverage(bins);
        RunReader reader(cells);
        CellRecord cell;
        while (reader.next(cell)) {
            require(cell.x < bins && cell.y < bins, "cis cell lies outside chromosome bins");
            coverage[cell.x] = checked_add(coverage[cell.x], cell.count);
            if (cell.x != cell.y) coverage[cell.y] = checked_add(coverage[cell.y], cell.count);
        }
        if (options.vc) {
            std::vector<float> norm(bins);
            for (uint32_t i = 0; i < bins; ++i)
                norm[i] = static_cast<float>(as_long_double(coverage[i]));
            rescale_norm(cells, norm);
            if (valid_norm(norm))
                manifest.vectors.push_back(write_norm(vc_id, chr, ri, resolution, norm, directory));
        }
        if (options.vc_sqrt) {
            std::vector<float> norm(bins);
            for (uint32_t i = 0; i < bins; ++i)
                norm[i] = std::sqrt(static_cast<float>(as_long_double(coverage[i])));
            rescale_norm(cells, norm);
            if (valid_norm(norm))
                manifest.vectors.push_back(write_norm(vcs_id, chr, ri, resolution, norm, directory));
        }
        // Release the 128-bit exact coverage array before constructing SCALE's
        // CSR and iterative working vectors.
        std::vector<U128>().swap(coverage);
        if (scale_active) {
            ScaleOptions scale_options = options.scale_options;
            scale_options.temporary_directory = temporary;
            scale_options.sort_memory_bytes = options.memory_bytes;
            ScaleResult scale = scale_cis(cells, bins, scale_options);
            if (scale.success) {
                rescale_norm(cells, scale.values);
                manifest.vectors.push_back(write_norm(scale_id, chr, ri, resolution,
                                                      scale.values, directory));
                std::cerr << "  SCALE converged in " << scale.iterations << " iterations\n";
            } else {
                scale_active = false;
                std::cerr << "  SCALE failed: " << scale.reason
                          << "; skipping finer resolutions for this chromosome\n";
            }
        }
    }
    write_vector_manifest(manifest, result_path);
}

void expected_resolution(const std::string &stage_path, const std::string &build_path,
                         const std::string &directory, uint32_t ri,
                         const NormalizeOptions &options) {
    StageManifest stage = read_stage_manifest(stage_path);
    BuildManifest build = read_build_manifest(build_path);
    require(stage.source_fingerprint == build.source_fingerprint,
            "stage and build manifests refer to different inputs");
    require(ri < build.resolutions.size(), "expected-value resolution index outside manifest");
    make_directory(directory);
    VectorManifest manifest = configured_manifest(stage.source_fingerprint, options);
    const std::string result_path = expected_manifest_path(directory, ri);
    if (path_exists(result_path)) {
        validate_subset(read_vector_manifest(result_path), manifest, result_path);
        std::cerr << "Expected-value task " << ri << " already complete\n";
        return;
    }
    std::map<std::tuple<uint32_t, uint32_t>, VectorInfo> norms;
    for (uint32_t chr = 0; chr < stage.chromosomes.size(); ++chr) {
        const std::string path = chromosome_manifest_path(directory, chr);
        VectorManifest subset = read_vector_manifest(path);
        validate_subset(subset, manifest, path);
        for (const VectorInfo &vector : subset.vectors) {
            require(vector.kind == 0 && vector.chr == chr,
                    "unexpected entry in chromosome vector manifest " + path);
            if (vector.ri == ri)
                require(norms.emplace(std::make_tuple(vector.norm, chr), vector).second,
                        "duplicate chromosome normalization vector");
        }
    }
    const uint32_t resolution = build.resolutions[ri];
    uint32_t maximum_bins = 0;
    for (const auto &chromosome : stage.chromosomes)
        maximum_bins = std::max(maximum_bins, chromosome_bins(chromosome.length, resolution));
    RawExpected raw(maximum_bins);
    std::vector<std::unique_ptr<NormalizedExpected>> expected;
    expected.reserve(manifest.norms.size());
    for (size_t norm = 0; norm < manifest.norms.size(); ++norm)
        expected.emplace_back(new NormalizedExpected(maximum_bins));
    for (uint32_t chr = 0; chr < stage.chromosomes.size(); ++chr) {
        auto cell_it = build.cells.find(std::make_tuple(chr, chr, resolution));
        if (cell_it == build.cells.end()) continue;
        const uint32_t bins = chromosome_bins(stage.chromosomes[chr].length, resolution);
        raw.add(chr, cell_it->second, bins);
        for (uint32_t norm = 0; norm < manifest.norms.size(); ++norm) {
            auto vector_it = norms.find(std::make_tuple(norm, chr));
            if (vector_it == norms.end()) continue;
            require(vector_it->second.words == bins, "normalization vector length mismatch");
            VectorFileReader input(vector_it->second);
            std::vector<uint32_t> words = input.read(0, bins);
            std::vector<float> values(bins);
            std::transform(words.begin(), words.end(), values.begin(), bits_float);
            expected[norm]->add(chr, cell_it->second, bins, values);
        }
    }
    if (raw.has_data()) manifest.vectors.push_back(raw.finish(ri, resolution, directory));
    for (uint32_t norm = 0; norm < manifest.norms.size(); ++norm)
        if (expected[norm]->has_data())
            manifest.vectors.push_back(expected[norm]->finish(norm, ri, resolution, directory));
    write_vector_manifest(manifest, result_path);
}

void finalize_vectors(const std::string &stage_path, const std::string &build_path,
                      const std::string &directory, const NormalizeOptions &options) {
    StageManifest stage = read_stage_manifest(stage_path);
    BuildManifest build = read_build_manifest(build_path);
    require(stage.source_fingerprint == build.source_fingerprint,
            "stage and build manifests refer to different inputs");
    VectorManifest manifest = configured_manifest(stage.source_fingerprint, options);
    std::set<std::tuple<uint8_t, uint32_t, uint32_t, uint32_t>> keys;
    auto append = [&](const std::string &path) {
        VectorManifest subset = read_vector_manifest(path);
        validate_subset(subset, manifest, path);
        for (VectorInfo &vector : subset.vectors) {
            auto key = std::make_tuple(vector.kind, vector.norm, vector.chr, vector.ri);
            require(keys.insert(key).second, "duplicate normalization artifact while finalizing");
            manifest.vectors.push_back(std::move(vector));
        }
    };
    for (uint32_t chr = 0; chr < stage.chromosomes.size(); ++chr)
        append(chromosome_manifest_path(directory, chr));
    for (uint32_t ri = 0; ri < build.resolutions.size(); ++ri)
        append(expected_manifest_path(directory, ri));
    std::sort(manifest.vectors.begin(), manifest.vectors.end(), [](const VectorInfo &a, const VectorInfo &b) {
        return std::tie(a.kind, a.norm, a.chr, a.ri) < std::tie(b.kind, b.norm, b.chr, b.ri);
    });
    const std::string result = join_path(directory, "vectors.manifest");
    write_vector_manifest(manifest, result);
    std::cerr << "Normalization sidecars complete: " << result << "\n";
}

void normalize_cells(const std::string &stage_path, const std::string &build_path,
                     const std::string &directory, const NormalizeOptions &options) {
    StageManifest stage = read_stage_manifest(stage_path);
    BuildManifest build = read_build_manifest(build_path);
    require(stage.source_fingerprint == build.source_fingerprint,
            "stage and build manifests refer to different inputs");
    for (uint32_t chr = 0; chr < stage.chromosomes.size(); ++chr)
        normalize_chromosome(stage_path, build_path, directory, chr, options);
    for (uint32_t ri = 0; ri < build.resolutions.size(); ++ri)
        expected_resolution(stage_path, build_path, directory, ri, options);
    finalize_vectors(stage_path, build_path, directory, options);
}

void print_normalize_tasks(const std::string &stage_path, const std::string &build_path,
                           const std::string &directory, const NormalizeOptions &) {
    StageManifest stage = read_stage_manifest(stage_path);
    BuildManifest build = read_build_manifest(build_path, false);
    require(stage.source_fingerprint == build.source_fingerprint,
            "stage and build manifests refer to different inputs");
    for (uint32_t chr = 0; chr < stage.chromosomes.size(); ++chr)
        std::cout << "normalize-chr\t" << chr << '\t'
                  << chromosome_manifest_path(directory, chr) << "\n";
    for (uint32_t ri = 0; ri < build.resolutions.size(); ++ri)
        std::cout << "expected-res\t" << ri << '\t'
                  << expected_manifest_path(directory, ri) << "\n";
    std::cout << "finalize-vectors\t0\t" << join_path(directory, "vectors.manifest") << "\n";
}

} // namespace hic10large
