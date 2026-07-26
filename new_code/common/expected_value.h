#pragma once
// Expected value (distance-decay) computation for intra-chromosomal contacts.
//
// For each resolution, we track:
//   actualDistances[d] = total contact count at genomic distance d bins
//   chromosomeCounts[chr_idx] = total intra contact count for that chr
//
// After all pairs are processed, computeDensity() smooths actualDistances
// to produce the expected value vector (one value per distance bin).
// chrScaleFactors are then computed so that the integral of expected * chr_length
// matches the observed total for each chromosome.
//
// This mirrors ExpectedValueCalculation.java.

#include "genome.h"
#include "hic_file_def.h"
#include <vector>
#include <unordered_map>
#include <mutex>
#include <cstdint>

class ExpectedValueCalculation {
public:
    ExpectedValueCalculation(const Genome& genome, int bin_size);

    // Thread-safe: add a contact at (bin1, bin2) for chromosome chr_idx.
    // Call only for intra-chromosomal contacts (chr1 == chr2).
    void add_distance(int chr_idx, int bin1, int bin2, double weight = 1.0);

    // Compute the density (expected value) vector from accumulated distances.
    // Must be called after all add_distance() calls.
    void compute_density();

    // The expected value at distance d bins (after compute_density()).
    float expected_at(int64_t dist_bins) const;

    int  bin_size()   const { return bin_size_; }
    int64_t n_values() const { return static_cast<int64_t>(density_.size()); }
    const std::vector<float>& density() const { return density_; }
    const std::unordered_map<int32_t, float>& chr_scale_factors() const { return chr_scale_; }
    // Data availability must be testable before compute_density(); both pre and
    // addnorm use this predicate to decide whether to serialize an expected vector.
    bool has_data() const { return !chr_counts_.empty(); }

private:
    const Genome&                      genome_;
    int                                bin_size_;
    int64_t                            n_bins_;         // max chromosome / bin_size + 1
    std::vector<double>                actual_distances_; // indexed by distance in bins
    std::unordered_map<int, double>    chr_counts_;     // chr_idx → total count
    std::vector<float>                 density_;        // output expected values
    std::unordered_map<int32_t, float> chr_scale_;      // chr_idx → scale factor
    bool                               computed_ = false;
    std::mutex                         mutex_;
};
