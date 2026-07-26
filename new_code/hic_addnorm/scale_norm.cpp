// SCALE matrix-balancing normalization.
// Direct port of SCALE-normalize-main/FINITO/finito.c, d_thMul.c, d_ppNormVector.c.

#include "scale_norm.h"
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <vector>
#include <thread>
#include <stdexcept>

// --------------------------------------------------------------------------
// Threaded symmetric upper-triangle matrix × vector multiply
// result[i] = sum_j (A_ij * v[j])  for upper-triangle entries
// (diagonal entries contribute once to result[i], off-diagonal contribute twice)
// --------------------------------------------------------------------------

static void utmv_mul_chunk(
    const uint32_t* row, const uint32_t* col, const float* val,
    long m, const double* v, double* res)
{
    for (long p = 0; p < m; p++) {
        res[row[p]] += (double)val[p] * v[col[p]];
        if (row[p] != col[p])
            res[col[p]] += (double)val[p] * v[row[p]];
    }
}

static void utmv_mul(
    const std::vector<uint32_t>& row,
    const std::vector<uint32_t>& col,
    const std::vector<float>&    val,
    long m, const double* v, uint32_t k,
    double* res, int threads,
    std::vector<std::vector<double>>& workspace)
{
    // Partition work among threads
    if (threads <= 1) {
        std::fill(res, res + k, 0.0);
        utmv_mul_chunk(row.data(), col.data(), val.data(), m, v, res);
        return;
    }

    // Ensure workspace has enough buffers
    workspace.resize(threads, std::vector<double>(k, 0.0));
    for (auto& ws : workspace) std::fill(ws.begin(), ws.end(), 0.0);

    long chunk = m / threads;
    std::vector<std::thread> ths;
    ths.reserve(threads);

    for (int t = 0; t < threads; t++) {
        long start = t * chunk;
        long end   = (t == threads - 1) ? m : start + chunk;
        double* ws = workspace[t].data();
        ths.emplace_back([&row, &col, &val, v, ws, start, end](){
            utmv_mul_chunk(
                row.data() + start, col.data() + start, val.data() + start,
                end - start, v, ws);
        });
    }
    for (auto& th : ths) th.join();

    // Reduce
    std::fill(res, res + k, 0.0);
    for (int t = 0; t < threads; t++) {
        for (uint32_t p = 0; p < k; p++) res[p] += workspace[t][p];
    }
}

// --------------------------------------------------------------------------
// Post-process: scale b so that (b^T A b) / (1^T A 1) = 1
// From d_ppNormVector.c
// --------------------------------------------------------------------------

void pp_norm_vector(
    long m,
    const std::vector<uint32_t>& row,
    const std::vector<uint32_t>& col,
    const std::vector<float>&    val,
    uint32_t k,
    std::vector<double>& b,
    int num_threads)
{
    std::vector<std::vector<double>> ws;
    std::vector<double> one(k), u(k, 0.0), v_vec(k, 0.0);

    // Replace zero b[p] with NaN
    for (uint32_t p = 0; p < k; p++) if (b[p] == 0.0) b[p] = std::numeric_limits<double>::quiet_NaN();
    for (uint32_t p = 0; p < k; p++) one[p] = std::isnan(b[p]) ? 0.0 : 1.0;

    // u = A * one
    utmv_mul(row, col, val, m, one.data(), k, u.data(), num_threads, ws);
    double s1 = 0.0;
    for (uint32_t p = 0; p < k; p++) s1 += u[p] * one[p];

    // v = A * b (treating NaN as 0)
    for (uint32_t p = 0; p < k; p++) v_vec[p] = std::isnan(b[p]) ? 0.0 : b[p];
    utmv_mul(row, col, val, m, v_vec.data(), k, u.data(), num_threads, ws);
    double s2 = 0.0;
    for (uint32_t p = 0; p < k; p++) s2 += u[p] * v_vec[p];

    if (s1 == 0.0 || s2 == 0.0) return;
    double s = std::sqrt(s2 / s1);
    for (uint32_t p = 0; p < k; p++) {
        if (!std::isnan(b[p])) b[p] = s / b[p];
    }
}

// --------------------------------------------------------------------------
// Main SCALE balance algorithm — direct port of finito.c
// --------------------------------------------------------------------------

static int cmpint(const void* a, const void* b) {
    if (*(int*)a < *(int*)b) return -1;
    if (*(int*)a > *(int*)b) return 1;
    return 0;
}

int scale_balance(
    long m,
    const std::vector<uint32_t>& ii_vec,
    const std::vector<uint32_t>& jj_vec,
    const std::vector<float>&    xx_vec,
    uint32_t k,
    std::vector<double>& norm_vec,
    const ScaleParams& params)
{
    const double tol      = params.tolerance;
    const double del      = params.delta;
    int threads           = params.num_threads;
    int maxiter           = params.max_iter;
    int total_max         = params.total_max_iter;
    int width             = params.diag_width;
    double perc_in        = params.percentile;

    const uint32_t* ii = ii_vec.data();
    const uint32_t* jj = jj_vec.data();
    const float*    xx = xx_vec.data();

    norm_vec.assign(k, std::numeric_limits<double>::quiet_NaN());
    double* b = norm_vec.data();

    // Allocate working arrays
    std::vector<double> current(k), row_s(k), col_s(k), dr(k), dc(k), s_v(k), one_v(k);
    std::vector<int>    bad(k, 0), bad0(k, 0);
    std::vector<double> b_conv(k), b0(k);
    std::vector<int>    bad_conv(k, 0);
    std::vector<std::vector<double>> ws; // for threaded mul

    // Count non-zeros per row
    std::vector<int> nz(k, 0);
    for (long p = 0; p < m; p++) {
        nz[ii[p]]++;
        if (ii[p] != jj[p]) nz[jj[p]]++;
    }

    // Sort nz values to find percentile cutoffs
    uint32_t n0 = 0;
    std::vector<int> nz0;
    nz0.reserve(k);
    for (uint32_t p = 0; p < k; p++) if (nz[p] > 0) nz0.push_back(nz[p]);
    n0 = (uint32_t)nz0.size();
    std::sort(nz0.begin(), nz0.end());

    int junk = (int)(n0 * 0.2);
    int bound = (n0 > 0) ? nz0[junk] : 0;
    junk = (int)(n0 * perc_in);
    int low   = (n0 > 0) ? nz0[junk] : 0;
    int low0  = low;

    // Diagonal-only filtering (optional)
    std::vector<float> diag_v(k, 0.0f);
    if (width > 0) {
        for (long p = 0; p < m; p++) {
            uint32_t dist = (jj[p] >= ii[p]) ? (jj[p] - ii[p]) : (ii[p] - jj[p]);
            if (dist <= (uint32_t)width) {
                diag_v[ii[p]] += xx[p];
                if (ii[p] != jj[p]) diag_v[jj[p]] += xx[p];
            }
        }
        std::vector<float> diag_sorted;
        diag_sorted.reserve(k);
        for (uint32_t p = 0; p < k; p++) if (diag_v[p] > 0) diag_sorted.push_back(diag_v[p]);
        std::sort(diag_sorted.begin(), diag_sorted.end());
        if (!diag_sorted.empty()) {
            double med   = diag_sorted[diag_sorted.size() / 2];
            int lowMed   = (int)(med * 0.01);
            for (uint32_t p = 0; p < k; p++) if (diag_v[p] < (float)lowMed) bad0[p] = 1;
        }
    }

    // Initial bad rows (too few non-zeros)
    for (uint32_t p = 0; p < k; p++) if (nz[p] < low) bad[p] = 1;
    for (uint32_t p = 0; p < k; p++) bad[p] |= bad0[p];
    for (uint32_t p = 0; p < k; p++) one_v[p] = 1.0 - bad[p];
    for (uint32_t p = 0; p < k; p++) dr[p] = 1.0 - bad[p];
    for (uint32_t p = 0; p < k; p++) dc[p] = 1.0 - bad[p];

    for (uint32_t p = 0; p < k; p++) current[p] = std::sqrt(dr[p] * dc[p]);

    // Initial row sums
    utmv_mul(ii_vec, jj_vec, xx_vec, m, dc.data(), k, row_s.data(), threads, ws);
    for (uint32_t p = 0; p < k; p++) row_s[p] *= dr[p];

    double ber     = 10.0 * (1.0 + tol);
    int    iter    = 0;
    int    all_iter = 0;
    bool   conv    = false, div_flag = false;
    double low_conv = 1000.0, low_div = 0.0;
    double ber_conv = tol;
    bool   yes     = true;
    const double erez = 1.0e-4;

    std::vector<double> report(total_max + 10, 0.0);

    while (ber > tol && iter < maxiter && all_iter < total_max) {
        iter++;
        all_iter++;

        // dr = dr / row_s (for non-bad rows)
        for (uint32_t p = 0; p < k; p++) dr[p] = (bad[p] == 0 && row_s[p] != 0) ? dr[p] / row_s[p] : 0.0;

        // col = A * dr
        utmv_mul(ii_vec, jj_vec, xx_vec, m, dr.data(), k, col_s.data(), threads, ws);
        for (uint32_t p = 0; p < k; p++) col_s[p] *= dc[p];
        for (uint32_t p = 0; p < k; p++) dc[p] = (bad[p] == 0 && col_s[p] != 0) ? dc[p] / col_s[p] : 0.0;

        // row = A * dc
        utmv_mul(ii_vec, jj_vec, xx_vec, m, dc.data(), k, row_s.data(), threads, ws);
        for (uint32_t p = 0; p < k; p++) row_s[p] *= dr[p];

        // b = sqrt(dr * dc)
        for (uint32_t p = 0; p < k; p++) b[p] = std::sqrt(dr[p] * dc[p]);

        // Error
        ber = 0.0;
        int numBad = 0;
        for (uint32_t p = 0; p < k; p++) {
            if (bad[p]) continue;
            double t1 = std::fabs((b[p] - current[p]) / (b[p] + current[p] + 1e-300));
            if (t1 > ber) ber = t1;
            if (t1 > tol) numBad++;
        }
        if (all_iter - 1 < (int)report.size()) report[all_iter - 1] = ber;
        for (uint32_t p = 0; p < k; p++) b0[p] = current[p];
        for (uint32_t p = 0; p < k; p++) current[p] = b[p];

        if (ber < tol) {
            // Converged
            ber_conv  = ber;
            low_conv  = low;
            yes       = true;
            if (low <= low0) break;
            conv = true;
            for (uint32_t p = 0; p < k; p++) b_conv[p] = b[p];
            for (uint32_t p = 0; p < k; p++) bad_conv[p] = bad[p];

            if (div_flag) {
                if (low_conv - low_div <= 1) break;
                low = (int)((low_conv + low_div) / 2);
            } else {
                low = (int)(low_conv / 2);
            }
            // Reset
            for (uint32_t p = 0; p < k; p++) { bad[p] = 0; one_v[p] = 1.0; }
            for (uint32_t p = 0; p < k; p++) if (nz[p] < low) bad[p] = 1;
            for (uint32_t p = 0; p < k; p++) { bad[p] |= bad0[p]; one_v[p] = 1 - bad[p]; }
            iter = 0; ber = 10.0;
            for (uint32_t p = 0; p < k; p++) dr[p] = dc[p] = 1.0 - bad[p];
            utmv_mul(ii_vec, jj_vec, xx_vec, m, dc.data(), k, row_s.data(), threads, ws);
            for (uint32_t p = 0; p < k; p++) row_s[p] *= dr[p];
            continue;
        }

        if (iter <= 5) continue;

        // Check convergence rate
        int prev = all_iter - 6;
        if (prev >= 0 && (report[all_iter - 1] * (1.0 + del) < report[prev])) continue;

        // Diverged
        div_flag  = true;
        low_div   = low;

        if (conv) {
            if (low_conv - low_div <= 1) {
                for (uint32_t p = 0; p < k; p++) b[p] = b_conv[p];
                for (uint32_t p = 0; p < k; p++) bad[p] = bad_conv[p];
                ber = ber_conv;
                break;
            } else if (((double)numBad) / n0 < erez && yes) {
                // Erez's trick: mark slow-converging rows as bad
                for (uint32_t p = 0; p < k; p++) {
                    if (bad[p]) continue;
                    double t1 = std::fabs((b[p] - b0[p]) / (b[p] + b0[p] + 1e-300));
                    if (t1 > tol) { bad[p] = 1; one_v[p] = 0; }
                }
                yes = false;
                goto next_iter;
            } else {
                low = (int)((low_div + low_conv) / 2);
                yes = true;
            }
        } else if (((double)numBad) / n0 < erez && yes) {
            for (uint32_t p = 0; p < k; p++) {
                if (bad[p]) continue;
                double t1 = std::fabs((b[p] - b0[p]) / (b[p] + b0[p] + 1e-300));
                if (t1 > tol) { bad[p] = 1; one_v[p] = 0; }
            }
            yes = false;
            goto next_iter;
        } else {
            low = 2 * low;
            yes = true;
        }

        for (uint32_t p = 0; p < k; p++) bad[p] = 0;
        for (uint32_t p = 0; p < k; p++) if (nz[p] < low) bad[p] = 1;
        for (uint32_t p = 0; p < k; p++) { bad[p] |= bad0[p]; one_v[p] = 1 - bad[p]; }

    next_iter:
        ber = 10.0; iter = 0;
        for (uint32_t p = 0; p < k; p++) dr[p] = dc[p] = 1.0 - bad[p];
        utmv_mul(ii_vec, jj_vec, xx_vec, m, dc.data(), k, row_s.data(), threads, ws);
        for (uint32_t p = 0; p < k; p++) row_s[p] *= dr[p];

        if (low > bound) break;
        if (all_iter >= total_max) break;
    }

    // Final: mark bad rows as NaN
    for (uint32_t p = 0; p < k; p++) if (bad[p]) b[p] = std::numeric_limits<double>::quiet_NaN();

    // If >20% of rows excluded, set all to NaN (convergence failure)
    int final_junk = 0;
    for (int v : nz0) if (v < (int)low_conv) final_junk++;
    double final_perc = (n0 > 0) ? (double)final_junk / n0 : 0.0;
    if (final_perc > 0.2) {
        for (uint32_t p = 0; p < k; p++) b[p] = std::numeric_limits<double>::quiet_NaN();
    }

    return all_iter;
}
