#pragma once
// SCALE normalization (matrix balancing) for intra-chromosomal contact matrices.
// Adapted from SCALE-normalize-main/FINITO/finito.c and d_thMul.c / d_ppNormVector.c.
//
// The algorithm:
//   Given a symmetric sparse matrix (upper triangle only) with m non-zero entries:
//     rows i[], cols j[], values x[]
//   Find vector b[] such that (D * A * D) ≈ J  (Jacobian / doubly stochastic)
//   where D = diag(b).
//
// The first attempt includes every non-empty row and starts from sqrt(VC).
// If convergence or the balanced row-sum check fails, progressively stricter
// low-coverage cutoffs are tried. A converged vector with isolated offending
// rows gets one targeted masking retry at the same cutoff. Excluded rows are
// masked as NaN.
// Returns the normalization vector b[] of length k (n_bins).

#include <cstdint>
#include <vector>
#include <cmath>

struct ScaleParams {
    double tolerance      = 1.0e-4;  // convergence tolerance
    double row_sum_tolerance = 5.0e-2; // max |b[i] * (A*b)[i] - 1|
    double delta          = 5.0e-2;  // convergence rate threshold
    int    max_iter       = 1000;    // max iterations per perc level
    int    total_max_iter = 2000;    // absolute max iterations
    int    num_threads    = 1;       // threads for matrix-vector multiply
    int    diag_width     = -1;      // diagonal width (-1 = full matrix)
};

// Run the SCALE (balance) algorithm.
// Input:
//   m: number of non-zero entries (upper triangle including diagonal)
//   row[m], col[m]: 0-based row/col indices
//   val[m]: contact values
//   k: number of bins (size of normalization vector)
// Output:
//   norm_vec[k]: normalization vector (NaN for unmappable bins)
// Returns final iteration count.
int scale_balance(
    long m,
    const std::vector<uint32_t>& row,
    const std::vector<uint32_t>& col,
    const std::vector<float>&    val,
    uint32_t k,
    std::vector<double>& norm_vec,
    const ScaleParams& params = ScaleParams{}
);

// Post-process: normalize b[] so that sqrt(b^T A b / 1^T A 1) = 1.
// This is ppNormVector from d_ppNormVector.c.
void pp_norm_vector(
    long m,
    const std::vector<uint32_t>& row,
    const std::vector<uint32_t>& col,
    const std::vector<float>&    val,
    uint32_t k,
    std::vector<double>& b,
    int num_threads = 1
);
