#include "expected_value.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

ExpectedValueCalculation::ExpectedValueCalculation(const Genome& genome, int bin_size)
    : genome_(genome), bin_size_(bin_size) {
    // Compute max bins needed: max chromosome length / bin_size + 1
    int64_t max_len = 0;
    for (auto* c : genome.chromosomes_without_all()) {
        max_len = std::max(max_len, c->length);
    }
    n_bins_ = max_len / bin_size + 1;
    actual_distances_.assign(n_bins_, 0.0);
}

void ExpectedValueCalculation::add_distance(int chr_idx, int bin1, int bin2, double weight) {
    if (weight <= 0.0) return;
    std::lock_guard<std::mutex> lock(mutex_);
    // Accumulate chromosome total
    chr_counts_[chr_idx] += weight;
    // Accumulate distance bin
    int64_t dist = std::abs((int64_t)bin1 - bin2);
    if (dist < n_bins_) {
        actual_distances_[dist] += weight;
    }
}

// Port of ExpectedValueCalculation.computeDensity() from Java.
// Uses adaptive window smoothing to compute density_avg.
void ExpectedValueCalculation::compute_density() {
    if (computed_) return;
    computed_ = true;

    // possibleDistances[d] = sum over all chromosomes of (nChrBins - d) for d < nChrBins
    // Only chromosomes that have observed contacts contribute.
    std::vector<double> possible_distances(n_bins_, 0.0);
    int64_t max_num_bins = 0;

    for (auto* chr : genome_.chromosomes_without_all()) {
        if (chr_counts_.find(chr->index) == chr_counts_.end()) continue;
        int64_t n_chr_bins = chr->length / bin_size_;
        max_num_bins = std::max(max_num_bins, n_chr_bins);
        for (int64_t i = 0; i < n_chr_bins && i < n_bins_; i++) {
            possible_distances[i] += (double)(n_chr_bins - i);
        }
    }

    if (max_num_bins == 0) {
        density_.assign(1, 0.0f);
        return;
    }

    density_.resize(max_num_bins, 0.0f);

    // Adaptive window smoothing (from Java computeDensity):
    // Window expands until it contains at least shotNoiseMinimum counts,
    // then shrinks as possible while staying above threshold.
    const double shot_noise_minimum = 400.0;

    double num_sum = actual_distances_[0];
    double den_sum = possible_distances[0];
    int64_t bound1 = 0;
    int64_t bound2 = 0;

    for (int64_t ii = 0; ii < max_num_bins; ii++) {
        if (num_sum < shot_noise_minimum) {
            // Expand window to the right until we have enough counts
            while (num_sum < shot_noise_minimum && bound2 < max_num_bins) {
                bound2++;
                if (bound2 < n_bins_) {
                    num_sum += actual_distances_[bound2];
                    den_sum += possible_distances[bound2];
                }
            }
        } else if (num_sum >= shot_noise_minimum && bound2 - bound1 > 0) {
            // Shrink window from both ends while maintaining threshold
            while (bound2 - bound1 > 0
                   && bound2 < n_bins_
                   && bound1 < n_bins_
                   && num_sum - actual_distances_[bound1] - actual_distances_[bound2] >= shot_noise_minimum) {
                num_sum -= actual_distances_[bound1] + actual_distances_[bound2];
                den_sum -= possible_distances[bound1] + possible_distances[bound2];
                bound1++;
                bound2--;
            }
        }

        density_[ii] = (den_sum > 0.0) ? (float)(num_sum / den_sum) : 0.0f;

        // Advance window center by 1 (bump bound2 up by 2 to keep centered)
        if (bound2 + 2 < max_num_bins) {
            if (bound2 + 1 < n_bins_) { num_sum += actual_distances_[bound2 + 1]; den_sum += possible_distances[bound2 + 1]; }
            if (bound2 + 2 < n_bins_) { num_sum += actual_distances_[bound2 + 2]; den_sum += possible_distances[bound2 + 2]; }
            bound2 += 2;
        } else if (bound2 + 1 < max_num_bins) {
            if (bound2 + 1 < n_bins_) { num_sum += actual_distances_[bound2 + 1]; den_sum += possible_distances[bound2 + 1]; }
            bound2++;
        }
    }

    // Compute per-chromosome scale factors:
    // scale = expected_total / observed_total
    // where expected_total = sum over d of (nChrBins - d) * density_[d]
    for (auto& [chr_idx, observed] : chr_counts_) {
        const Chromosome& chr = genome_.at(chr_idx);
        int64_t n_chr_bins = chr.length / bin_size_;
        double expected = 0.0;
        for (int64_t n = 0; n < n_chr_bins && n < (int64_t)density_.size(); n++) {
            expected += (double)(n_chr_bins - n) * density_[n];
        }
        if (observed > 0.0 && expected > 0.0) {
            chr_scale_[chr_idx] = (float)(expected / observed);
        }
    }
}

float ExpectedValueCalculation::expected_at(int64_t dist_bins) const {
    if (dist_bins < 0 || dist_bins >= (int64_t)density_.size()) return 0.0f;
    return density_[dist_bins];
}
