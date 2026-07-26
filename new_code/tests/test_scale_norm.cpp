#include "hic_addnorm/scale_norm.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

int main() {
    constexpr uint32_t k = 101;
    std::vector<uint32_t> row;
    std::vector<uint32_t> col;
    std::vector<float> val;

    auto add = [&](uint32_t r, uint32_t c, float value) {
        row.push_back(r);
        col.push_back(c);
        val.push_back(value);
    };

    // Row 0 is nonempty but has only one nonzero. The other 100 rows form a
    // diagonally supported chain and have at least two nonzeros. With 101
    // nonempty rows, the former C++ 1st-percentile initialization selected a
    // cutoff of 2 and masked row 0 before attempting to balance.
    add(0, 0, 2.0f);
    for (uint32_t i = 1; i < k; i++) {
        add(i, i, 2.0f);
        if (i + 1 < k) add(i, i + 1, 0.25f);
    }

    ScaleParams params;
    params.num_threads = 2;
    std::vector<double> b;
    scale_balance((long)row.size(), row, col, val, k, b, params);

    if (b.size() != k || !std::isfinite(b[0]) || b[0] <= 0.0) {
        std::cerr << "The nonempty sparse row was incorrectly excluded\n";
        return 1;
    }

    std::vector<double> ab(k, 0.0);
    for (size_t p = 0; p < row.size(); p++) {
        ab[row[p]] += val[p] * b[col[p]];
        if (row[p] != col[p]) ab[col[p]] += val[p] * b[row[p]];
    }

    double max_error = 0.0;
    for (uint32_t i = 0; i < k; i++) {
        if (!std::isfinite(b[i])) continue;
        max_error = std::max(max_error, std::fabs(b[i] * ab[i] - 1.0));
    }
    if (max_error > params.row_sum_tolerance) {
        std::cerr << "Accepted SCALE vector has row-sum error "
                  << max_error << '\n';
        return 1;
    }

    return 0;
}
