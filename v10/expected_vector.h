#pragma once
#include "format.h"

namespace hic10 {
// Shared full-distance expected-vector smoothing and chromosome-scale logic.
template<class Real>
Vector finish_expected(const std::vector<Real> &actual,
                       const std::map<uint32_t, Real> &observed,
                       const std::map<uint32_t, uint32_t> &chromosome_bins,
                       uint8_t kind, uint32_t norm, uint8_t unit, uint32_t ri, bool smooth) {
    std::vector<uint32_t> lengths;
    for (auto e : chromosome_bins)
        if (observed.count(e.first))
            lengths.push_back(e.second);
    std::sort(lengths.begin(), lengths.end());
    std::vector<uint64_t> suffix(lengths.size() + 1, 0);
    for (size_t i = lengths.size(); i-- > 0;)
        suffix[i] = suffix[i + 1] + lengths[i];
    auto possible = [&](uint64_t distance) {
        size_t first = std::upper_bound(lengths.begin(), lengths.end(), distance) -
                       lengths.begin();
        return Real(suffix[first]) - Real(lengths.size() - first) * distance;
    };
    Vector v;
    v.kind = kind;
    v.norm = norm;
    v.unit = unit;
    v.ri = ri;
    v.values.reserve(actual.size());
    uint64_t support = lengths.empty() ? 0 : lengths.back();
    if (smooth && support) {
        const Real minimum = 400;
        Real numerator = actual[0], denominator = possible(0);
        uint64_t lo = 0, hi = 0;
        for (uint64_t i = 0; i < support; ++i) {
            if (numerator < minimum) {
                while (numerator < minimum && hi + 1 < support) {
                    ++hi;
                    numerator += actual[hi];
                    denominator += possible(hi);
                }
            } else
                while (hi > lo && numerator - actual[lo] - actual[hi] >= minimum) {
                    numerator -= actual[lo] + actual[hi];
                    denominator -= possible(lo) + possible(hi);
                    ++lo;
                    --hi;
                }
            v.values.push_back(
                bits(denominator > 0 ? static_cast<float>(numerator / denominator) : 0));
            if (hi + 2 < support) {
                numerator += actual[hi + 1] + actual[hi + 2];
                denominator += possible(hi + 1) + possible(hi + 2);
                hi += 2;
            } else if (hi + 1 < support) {
                ++hi;
                numerator += actual[hi];
                denominator += possible(hi);
            }
        }
    } else
        for (uint64_t i = 0; i < support; ++i) {
            Real denominator = possible(i);
            v.values.push_back(
                bits(denominator > 0 ? static_cast<float>(actual[i] / denominator) : 0));
        }
    v.values.resize(actual.size(), 0x7fc00000);
    for (auto e : observed) {
        Real total = 0;
        uint32_t n = chromosome_bins.at(e.first);
        for (uint32_t d = 0; d < n && d < v.values.size(); ++d)
            total += Real(n - d) * floating(v.values[d]);
        if (total > 0 && e.second > 0)
            v.scales[e.first] = bits(static_cast<float>(total / e.second));
    }
    return v;

}
} // namespace hic10
