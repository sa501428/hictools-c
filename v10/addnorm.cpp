#include "addnorm.h"
#include "hic_addnorm/scale_norm.h"
#include "reader.h"
#include "vector_updater.h"
#include "expected_vector.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>

namespace hic10 {
namespace {
struct Sparse {
    std::vector<uint32_t> row, col;
    std::vector<double> value;
    uint32_t bins = 0;
};
Sparse sparse_matrix(const Matrix &matrix, uint32_t bins) {
    Sparse s;
    s.bins = bins;
    s.row.reserve(matrix.cells.size());
    s.col.reserve(matrix.cells.size());
    s.value.reserve(matrix.cells.size());
    for (auto cell : matrix.cells) {
        double value = matrix.scores ? floating(static_cast<uint32_t>(cell.value))
                                     : static_cast<double>(cell.value);
        if (!std::isfinite(value) || value <= 0)
            continue;
        check(cell.x < bins && cell.y < bins, "normalization cell outside chromosome");
        s.row.push_back(cell.x);
        s.col.push_back(cell.y);
        s.value.push_back(value);
    }
    return s;
}
std::vector<float> raw_vc(const Sparse &s) {
    std::vector<double> sums(s.bins, 0);
    for (size_t i = 0; i < s.row.size(); ++i) {
        sums[s.row[i]] += s.value[i];
        if (s.row[i] != s.col[i])
            sums[s.col[i]] += s.value[i];
    }
    std::vector<float> result(s.bins);
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = static_cast<float>(sums[i]);
    return result;
}
void fix_sum(const Sparse &s, const std::vector<std::vector<float>*> &norms) {
    std::vector<double> raw(norms.size(), 0.0), normalized(norms.size(), 0.0);
    for (size_t i = 0; i < s.row.size(); ++i) {
        uint32_t x = s.row[i], y = s.col[i];
        double multiple = x == y ? 1 : 2;
        for (size_t n = 0; n < norms.size(); ++n) {
            float nx = (*norms[n])[x], ny = (*norms[n])[y];
            if (!(nx > 0) || !(ny > 0) || !std::isfinite(nx) || !std::isfinite(ny))
                continue;
            raw[n] += multiple * s.value[i];
            normalized[n] += multiple * s.value[i] / (double(nx) * ny);
        }
    }
    for (size_t n = 0; n < norms.size(); ++n) {
        if (raw[n] <= 0 || normalized[n] <= 0)
            continue;
        float factor = static_cast<float>(std::sqrt(normalized[n] / raw[n]));
        for (float &v : *norms[n])
            if (v > 0 && std::isfinite(v))
                v *= factor;
    }
}
struct ExpectedAccumulator {
    std::vector<double> actual;
    std::map<uint32_t, double> observed;
    std::map<uint32_t, uint32_t> chromosome_bins;
    explicit ExpectedAccumulator(uint32_t size) : actual(size, 0) {}
    static void add_many(
        uint32_t chr, const Sparse &s,
        const std::vector<std::pair<ExpectedAccumulator*, const std::vector<float>*>>& work) {
        std::vector<double> totals(work.size(), 0.0);
        for (const auto& item : work) item.first->chromosome_bins[chr] = s.bins;
        for (size_t i = 0; i < s.row.size(); ++i) {
            for (size_t n = 0; n < work.size(); ++n) {
                double value = s.value[i];
                if (work[n].second) {
                    float a = (*work[n].second)[s.row[i]];
                    float b = (*work[n].second)[s.col[i]];
                    if (!(a > 0) || !(b > 0) || !std::isfinite(a) || !std::isfinite(b))
                        continue;
                    value /= double(a) * b;
                }
                work[n].first->actual[s.col[i] - s.row[i]] += value;
                totals[n] += value;
            }
        }
        for (size_t n = 0; n < work.size(); ++n)
            if (totals[n] > 0) work[n].first->observed[chr] += totals[n];
    }
    Vector finish(uint8_t kind, uint32_t norm, uint8_t unit, uint32_t ri, bool smooth) {
        return finish_expected(actual, observed, chromosome_bins, kind, norm, unit, ri, smooth);
    }
};
uint32_t norm_id(const Header &h, const std::string &name) {
    auto it = std::find(h.norms.begin(), h.norms.end(), name);
    check(it != h.norms.end(), "normalization dictionary lacks requested type " + name);
    return static_cast<uint32_t>(it - h.norms.begin());
}
bool real_chromosome(const Chromosome &c) {
    return c.name != "ALL" && c.name != "All" && c.name != "all";
}
} // namespace
void add_norm_v10(const std::string &path, const AddNormOptions &options) {
    Reader reader(path);
    Header h = reader.header();
    bool build_scale = false;
    if (options.scale)
        for (uint8_t unit = 0; unit < 2; ++unit)
            for (const auto &resolution : h.resolutions[unit])
                if (unit != 0 || options.minimum_scale_resolution == 0 ||
                    resolution.bin >= uint32_t(options.minimum_scale_resolution))
                    build_scale = true;
    if (options.vc && std::find(h.norms.begin(), h.norms.end(), "VC") == h.norms.end())
        h.norms.push_back("VC");
    if (options.vc_sqrt && std::find(h.norms.begin(), h.norms.end(), "VC_SQRT") == h.norms.end())
        h.norms.push_back("VC_SQRT");
    if (build_scale && std::find(h.norms.begin(), h.norms.end(), "SCALE") == h.norms.end())
        h.norms.push_back("SCALE");
    uint32_t vc = options.vc ? norm_id(h, "VC") : UINT32_MAX;
    uint32_t vcs = options.vc_sqrt ? norm_id(h, "VC_SQRT") : UINT32_MAX;
    uint32_t scale = build_scale ? norm_id(h, "SCALE") : UINT32_MAX;
    auto existing = reader.vector_entries();
    auto present = [&](uint32_t norm, uint8_t unit, uint32_t ri) {
        for (const auto &e : existing)
            if (e.kind != 1 && e.norm == norm && e.unit == unit && e.ri == ri) return true;
        return false;
    };
    bool work = false;
    for (uint8_t unit = 0; unit < 2; ++unit)
        for (uint32_t ri = 0; ri < h.resolutions[unit].size(); ++ri) {
            bool raw_present = false;
            for (const auto &e : existing)
                if (e.kind == 1 && e.unit == unit && e.ri == ri) raw_present = true;
            bool scale_requested = build_scale && (unit || !options.minimum_scale_resolution ||
                h.resolutions[unit][ri].bin >= uint32_t(options.minimum_scale_resolution));
            work |= !raw_present || (options.vc && !present(vc, unit, ri)) ||
                (options.vc_sqrt && !present(vcs, unit, ri)) ||
                (scale_requested && !present(scale, unit, ri));
        }
    if (!work) {
        std::fprintf(stderr, "Requested V10 normalization bundles already exist: %s\n", path.c_str());
        return;
    }
    VectorOutput output(reader, h, path, options.compression_level);
    for (uint8_t unit = 0; unit < 2; ++unit)
        for (uint32_t ri = 0; ri < h.resolutions[unit].size(); ++ri) {
            uint32_t bin = h.resolutions[unit][ri].bin;
            bool derived = h.resolutions[unit][ri].mode != 0;
            std::fprintf(stderr, "\n%s resolution %u %s%s\n", unit ? "FRAG" : "BP", bin,
                         derived ? "(derived on the fly) " : "", "");
            uint32_t maximum = 0;
            for (uint32_t chr = 0; chr < h.chromosomes.size(); ++chr)
                maximum = std::max(maximum, narrow(h.bins(chr, unit, ri)));
            check(maximum != 0, "cannot create an empty expected vector");
            const bool do_vc = options.vc && !present(vc, unit, ri);
            const bool do_vcs = options.vc_sqrt && !present(vcs, unit, ri);
            const bool do_raw = !output.contains(1, 0, 0, unit, ri);
            const bool do_scale =
                build_scale && !present(scale, unit, ri) && (unit != 0 || options.minimum_scale_resolution == 0 ||
                                bin >= uint32_t(options.minimum_scale_resolution));

            if (!do_vc && !do_vcs && !do_scale && !do_raw) continue;

            ExpectedAccumulator raw(do_raw ? maximum : 0);
            ExpectedAccumulator expected_vc(do_vc ? maximum : 0), expected_vcs(do_vcs ? maximum : 0),
                                expected_scale(do_scale ? maximum : 0);
            uint32_t written_vc = 0, written_vcs = 0, written_scale = 0;

            auto save_norm = [&](uint32_t chr, uint32_t id,
                                 const std::vector<float>& norm, uint32_t& written) {
                bool valid = false;
                for (float value : norm)
                    if (value > 0 && std::isfinite(value)) {
                        valid = true;
                        break;
                    }
                if (!valid)
                    return;
                Vector v;
                v.kind = 0;
                v.norm = id;
                v.chr = chr;
                v.unit = unit;
                v.ri = ri;
                v.values.reserve(norm.size());
                for (float value : norm)
                    v.values.push_back(bits(value));
                output.add(v);
                ++written;
            };

            // Materializing a derived matrix can itself be expensive. Keep the
            // chromosome outermost so raw expected, VC, VC_SQRT, and SCALE all
            // consume the same sparse matrix and the same raw coverage vector.
            for (uint32_t chr = 0; chr < h.chromosomes.size(); ++chr) {
                if (!real_chromosome(h.chromosomes[chr]))
                    continue;
                auto s = sparse_matrix(reader.matrix(chr, chr, unit, ri),
                                       narrow(h.bins(chr, unit, ri)));
                std::vector<float> coverage, vc_norm, vcs_norm, scale_norm;
                if (!s.row.empty() && (do_vc || do_vcs))
                    coverage = raw_vc(s);

                if (do_vc) vc_norm = coverage;
                if (do_vcs) {
                    vcs_norm = coverage;
                    for (float& value : vcs_norm) value = std::sqrt(value);
                }
                if (!s.row.empty() && do_scale) {
                    ScaleParams p;
                    p.tolerance = options.tolerance;
                    p.total_max_iter = options.max_iterations;
                    p.num_threads = options.threads;
                    std::vector<float> values;
                    values.reserve(s.value.size());
                    std::vector<double> scale_vc(s.bins, 0.0);
                    for (size_t i = 0; i < s.value.size(); ++i) {
                        double value = s.value[i];
                        values.push_back(static_cast<float>(value));
                        const double stored = values.back();
                        scale_vc[s.row[i]] += stored;
                        if (s.row[i] != s.col[i]) scale_vc[s.col[i]] += stored;
                    }
                    std::vector<double> b(s.bins,
                                          std::numeric_limits<double>::quiet_NaN());
                    scale_balance(s.row.size(), s.row, s.col, values, s.bins, b, p,
                                  &scale_vc);
                    pp_norm_vector(s.row.size(), s.row, s.col, values, s.bins, b,
                                   options.threads);
                    scale_norm.resize(b.size());
                    for (size_t i = 0; i < b.size(); ++i)
                        scale_norm[i] = static_cast<float>(b[i]);
                }

                std::vector<std::vector<float>*> norms;
                if (!vc_norm.empty()) norms.push_back(&vc_norm);
                if (!vcs_norm.empty()) norms.push_back(&vcs_norm);
                if (!scale_norm.empty()) norms.push_back(&scale_norm);
                fix_sum(s, norms);

                std::vector<std::pair<ExpectedAccumulator*, const std::vector<float>*>>
                    expected_work;
                if (do_raw) expected_work.push_back({&raw, nullptr});
                if (!vc_norm.empty()) expected_work.push_back({&expected_vc, &vc_norm});
                if (!vcs_norm.empty()) expected_work.push_back({&expected_vcs, &vcs_norm});
                if (!scale_norm.empty())
                    expected_work.push_back({&expected_scale, &scale_norm});
                ExpectedAccumulator::add_many(chr, s, expected_work);

                if (!vc_norm.empty()) save_norm(chr, vc, vc_norm, written_vc);
                if (!vcs_norm.empty()) save_norm(chr, vcs, vcs_norm, written_vcs);
                if (!scale_norm.empty()) save_norm(chr, scale, scale_norm, written_scale);
            }

            if (do_raw) output.add(raw.finish(1, 0, unit, ri, true));
            if (written_vc)
                output.add(expected_vc.finish(2, vc, unit, ri, true));
            if (written_vcs)
                output.add(expected_vcs.finish(2, vcs, unit, ri, true));
            if (written_scale)
                output.add(expected_scale.finish(2, scale, unit, ri, true));
            if (do_vc)
                std::fprintf(stderr, "  VC (%u chromosomes)\n", written_vc);
            if (do_vcs)
                std::fprintf(stderr, "  VC_SQRT (%u chromosomes)\n", written_vcs);
            if (do_scale)
                std::fprintf(stderr, "  SCALE (%u chromosomes)\n", written_scale);
        }
    output.finish(h);
    std::fprintf(stderr, "\nV10 normalization complete: %s\n", path.c_str());
}
} // namespace hic10
