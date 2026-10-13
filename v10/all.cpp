#include "all.h"
#include "repack.h"
#include <iostream>
#include <cctype>
#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <set>

namespace hic10 {
bool is_all(const std::string &name) {
    return name.size() == 3 && std::tolower(static_cast<unsigned char>(name[0])) == 'a' &&
        std::tolower(static_cast<unsigned char>(name[1])) == 'l' &&
        std::tolower(static_cast<unsigned char>(name[2])) == 'l';
}
namespace {
uint32_t choose_source(const Header &h, uint64_t preference, uint64_t physical_scale,
                       const std::vector<uint32_t> &allowed = {}) {
    auto eligible = [&](const Resolution &r) {
        return !r.mode && (allowed.empty() || std::find(allowed.begin(), allowed.end(), r.bin) != allowed.end());
    };
    uint32_t below = 0, above = 0, first = 0;
    preference = std::min(preference, physical_scale);
    for (const auto &r : h.resolutions[0]) if (eligible(r)) {
        if (!first) first = r.bin;
        if (r.bin <= preference) below = r.bin;
        if (!above && r.bin >= preference && r.bin <= physical_scale) above = r.bin;
    }
    check(first > 0, "All has no usable source resolution");
    if (!below) return above ? above : first;
    using Wide = unsigned __int128;
    return above && Wide(above) * below < Wide(preference) * preference ? above : below;
}
}
AllBuilder::AllBuilder(Header &h, const AllOptions &o) {
    check(o.bins > 0 && o.bins <= 4096 && o.source > 0, "invalid All overview options");
    check(!h.resolutions[0].empty(), "All requires a BP resolution");
    finest_ = h.resolutions[0].front().bin;
    preference_ = o.source;
    index_ = narrow(h.chromosomes.size());
    for (uint32_t i = 0; i < h.chromosomes.size(); ++i) if (is_all(h.chromosomes[i].name)) {
        check(o.rebuild && index_ == h.chromosomes.size(), "All already exists or is duplicated");
        index_ = i;
    }
    uint64_t total = 0, maximum_length = 0;
    offsets_.resize(h.chromosomes.size());
    for (uint32_t i = 0; i < h.chromosomes.size(); ++i) {
        const auto &c = h.chromosomes[i];
        if (i == index_) continue;
        offsets_[i] = total;
        total = plus(total, c.length);
        maximum_length = std::max(maximum_length, c.length);
    }
    check(total > 0, "All requires real chromosomes");
    // V9 uses unpadded genome-wide kilobase coordinates. Raise the scale only
    // for unusual genomes where kb would enlarge expected vectors or overflow
    // desktop boundary arrays. Ordinary mammalian maps retain exactly scale 1000.
    coordinate_scale_ = std::max<uint64_t>(1000, std::max(ceil_div(total, maximum_length),
                                                        ceil_div(total, uint64_t(INT32_MAX))));
    span_ = std::max<uint64_t>(1, total / coordinate_scale_);
    // Match V9's floor(All length / 500) overview width, rounded upward only
    // when the advertised finest resolution cannot represent that interval.
    uint64_t nominal = std::max<uint64_t>(1, span_ / o.bins);
    stride_ = std::max<uint64_t>(1, ceil_div(nominal, finest_));
    display_resolution_ = stride_ * finest_;
    if (ceil_div(span_, display_resolution_) > uint64_t(o.bins) + 1) {
        stride_ = ceil_div(ceil_div(span_, uint64_t(o.bins) + 1), finest_);
        display_resolution_ = stride_ * finest_;
    }
    bins_ = narrow(ceil_div(span_, display_resolution_));
    check(display_resolution_ <= UINT64_MAX / coordinate_scale_, "All physical scale overflow");
    source_ = choose_source(h, o.source, display_resolution_ * coordinate_scale_, o.sources);
    std::string lengths;
    for (uint32_t i = 0; i < h.chromosomes.size(); ++i) if (i != index_) {
        if (!lengths.empty()) lengths += ',';
        lengths += std::to_string(ceil_div(h.chromosomes[i].length, source_));
    }
    if (index_ == h.chromosomes.size()) h.chromosomes.push_back({"All", span_, {}});
    else { h.chromosomes[index_].length = span_; h.chromosomes[index_].sites.clear(); }
    h.attributes.erase(std::remove_if(h.attributes.begin(), h.attributes.end(), [](const auto &a) {
        return a.first.compare(0, 7, "hicAll.") == 0;
    }), h.attributes.end());
    h.attributes.insert(h.attributes.end(), {
        {"hicAll.version", "2"}, {"hicAll.sourceResolution", std::to_string(source_)},
        {"hicAll.groupFactor", "1"},
        {"hicAll.finestResolution", std::to_string(finest_)},
        {"hicAll.binCount", std::to_string(bins_)}, {"hicAll.chromosomeBins", lengths},
        {"hicAll.coordinateScale", std::to_string(coordinate_scale_)},
        {"hicAll.storageStride", std::to_string(stride_)}});
    header_ = h;
}
void AllBuilder::add(uint32_t a, uint32_t b, uint64_t x, uint64_t y,
                     uint32_t bin, uint64_t value, bool score) {
    check(a < offsets_.size() && b < offsets_.size() && a != index_ && b != index_ && bin > 0, "invalid All source chromosome/bin");
    auto coordinate = [&](uint32_t chr, uint64_t p) {
        check(p <= UINT64_MAX / bin, "All source coordinate overflow");
        uint64_t bp = p * bin;
        check(bp < header_.chromosomes[chr].length, "All source coordinate out of bounds");
        uint64_t pseudo = std::min(span_ - 1, plus(offsets_[chr], bp) / coordinate_scale_);
        return narrow(pseudo / display_resolution_);
    };
    uint32_t i = coordinate(a, x), j = coordinate(b, y);
    if (i > j) std::swap(i, j);
    auto key = std::make_pair(i, j);
    if (score) {
        if (!scored_) {
            for (auto c : counts_) scores_[c.first] = static_cast<double>(c.second);
            counts_.clear();
            scored_ = true;
        }
        scores_[key] += floating(narrow(value));
    } else if (scored_) scores_[key] += static_cast<double>(value);
    else counts_[key] = plus(counts_[key], value);
}
void AllBuilder::collect(Reader &reader) {
    for (auto key : reader.matrices()) {
        if (key.chr1 == index_ || key.chr2 == index_) continue;
        auto populated = reader.populated_bp_resolutions(key);
        if (populated.empty()) continue;
        uint32_t bin = std::find(populated.begin(), populated.end(), source_) != populated.end()
            ? source_ : choose_source(reader.header(), preference_, display_resolution_ * coordinate_scale_, populated);
        uint32_t ri = reader.header().resolution(0, bin);
        if (bin != source_) std::cerr << "All repair: pair " << key.chr1 << '_' << key.chr2
            << " uses populated fallback " << bin << " bp\n";
        // stream_materialized emits canonical raw cells and preserves u64 counts.
        // The type is a matrix property; Reader exposes it through this small
        // descriptor read rather than materializing the source matrix.
        auto loc = reader.matrix_location(key);
        auto desc = reader.read_bytes(loc.position + 24 + uint64_t(ri) * 76, 4);
        bool scored = desc.at(3) == 1;
        reader.stream_materialized(key.chr1, key.chr2, 0, ri, [&](const Cell &c) {
            add(key.chr1, key.chr2, c.x, c.y, bin, c.value, scored);
        }, true);
    }
}
Matrix AllBuilder::matrix(uint8_t unit, uint32_t ri) const {
    Matrix result;
    if (unit) return result; // No restriction-fragment interpretation for All.
    uint32_t bin = header_.resolutions[0].at(ri).bin;
    auto key = [&](std::pair<uint32_t, uint32_t> c) {
        return std::make_pair(narrow(uint64_t(c.first) * display_resolution_ / bin),
                              narrow(uint64_t(c.second) * display_resolution_ / bin));
    };
    result.scores = scored_;
    if (scored_) {
        std::map<std::pair<uint32_t, uint32_t>, double> sums;
        for (auto c : scores_) sums[key(c.first)] += c.second;
        for (auto c : sums) result.cells.push_back({c.first.first, c.first.second,
                                                  bits(static_cast<float>(c.second))});
    } else {
        std::map<std::pair<uint32_t, uint32_t>, uint64_t> sums;
        for (auto c : counts_) sums[key(c.first)] = plus(sums[key(c.first)], c.second);
        for (auto c : sums) if (c.second)
            result.cells.push_back({c.first.first, c.first.second, c.second});
    }
    return result;
}
void AllBuilder::write(Writer &writer) const {
    writer.matrix(index_, index_, [&](uint8_t u, uint32_t ri) { return matrix(u, ri); });
}
void add_all(const std::string &input, const std::string &output,
             const AllOptions &options, int level, uint32_t threads) {
    check(options.bins > 0 && options.bins <= 4096 && options.source > 0, "invalid All overview options");
    check(threads > 0 && threads <= 256, "All repair copy threads must be 1..256");
    Reader reader(input);
    for (const auto &c : reader.header().chromosomes) if (is_all(c.name) && !options.rebuild) {
        // Refuse a misleading successful repair of an existing broken All.
        uint32_t chr = narrow(&c - reader.header().chromosomes.data());
        auto loc = reader.matrix_location({chr, chr});
        check(loc.length > 0, "All header exists without an All/All matrix");
        if (input != output) {
            append_matrix_atomic(reader, reader.header(), output, threads, nullptr,
                [](const std::string &path) {
                    Reader validated(path); validated.validate_matrix_storage();
                    validated.for_each_vector_entry([](uint64_t, VectorEntry) {});
                });
        }
        return;
    }
    if (options.rebuild) reader.for_each_vector_entry([&](uint64_t, VectorEntry v) {
        if (v.kind == 0) check(!is_all(reader.header().chromosomes[v.chr].name),
            "cannot rebuild an All with normalization vectors; remove those vectors first");
    });
    Header header = reader.header();
    AllOptions actual = options;
    std::set<uint32_t> populated;
    for (auto key : reader.matrices()) {
        if (is_all(reader.header().chromosomes[key.chr1].name) || is_all(reader.header().chromosomes[key.chr2].name)) continue;
        for (auto bin : reader.populated_bp_resolutions(key)) populated.insert(bin);
    }
    actual.sources.assign(populated.begin(), populated.end());
    AllBuilder all(header, actual);
    std::cerr << "All repair: aggregating stored " << all.source() << " bp matrices\n";
    all.collect(reader);
    std::string pattern = output + ".all-matrix.XXXXXX";
    std::vector<char> name(pattern.begin(), pattern.end()); name.push_back(0);
    int fd = mkstemp(name.data()); check(fd >= 0, "cannot create All matrix temporary");
    close(fd);
    struct Cleanup { std::string path; ~Cleanup() { std::remove(path.c_str()); } } cleanup{name.data()};
    Options writer_options; writer_options.level = level; writer_options.threads = threads;
    writer_options.verifyDerived = false;
    // Preserve the original global storage policy, including explicit overrides.
    for (const auto &r : header.resolutions[0]) {
        if (r.mode) writer_options.derived.emplace_back(r.bin, header.resolutions[0][r.source].bin);
        else writer_options.materialized.push_back(r.bin);
    }
    { Writer writer(cleanup.path, header, writer_options); all.write(writer); writer.finish({}); }
    Reader extra(cleanup.path);
    append_matrix_atomic(reader, header, output, threads, &extra, [&](const std::string &path) {
        Reader validated(path);
        check(validated.header().chromosomes.size() == header.chromosomes.size(), "All repair header mismatch");
        validated.validate_matrix_storage();
        validated.for_each_vector_entry([](uint64_t, VectorEntry) {});
        for (uint32_t ri = 0; ri < header.resolutions[0].size(); ++ri)
            validated.matrix(all.index(), all.index(), 0, ri);
    });
}
} // namespace hic10
