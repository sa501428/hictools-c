#include "scale.h"

#include "common/thread_pool.h"
#include "io.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <future>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace hic10large {
namespace {

struct CsrEntry {
    uint32_t column = 0;
    float value = 0;
};
static_assert(sizeof(CsrEntry) == 8, "CSR entries must remain compact");
constexpr size_t CSR_HEADER_BYTES = 64;

class Mapping {
  public:
    Mapping(const std::string &path, size_t expected) : size_(expected) {
        fd_ = open(path.c_str(), O_RDONLY);
        require(fd_ >= 0, "cannot open CSR sidecar " + path);
        struct stat s{};
        require(fstat(fd_, &s) == 0 && s.st_size >= 0 && uint64_t(s.st_size) == expected,
                "CSR sidecar size mismatch " + path);
        data_ = mmap(nullptr, expected, PROT_READ, MAP_SHARED, fd_, 0);
        require(data_ != MAP_FAILED, "cannot map CSR sidecar " + path);
    }
    ~Mapping() {
        if (data_ != MAP_FAILED) munmap(data_, size_);
        if (fd_ >= 0) close(fd_);
    }
    const unsigned char *bytes() const { return static_cast<const unsigned char *>(data_); }
  private:
    int fd_ = -1;
    void *data_ = MAP_FAILED;
    size_t size_ = 0;
};

struct CsrFiles {
    std::string index, data;
    uint32_t bins = 0;
    uint64_t entries = 0;
};

std::string scale_prefix(const RunInfo &cells) {
    return "scale-" + std::to_string(cells.chr1) + "-" + std::to_string(cells.resolution);
}

void emit_transpose_chunk(std::vector<CellRecord> &chunk, std::vector<RunInfo> &runs,
                          const RunInfo &source, const ScaleOptions &options) {
    if (chunk.empty()) return;
    radix_sort_cells(chunk);
    aggregate_sorted_cells(chunk);
    std::string path = join_path(options.temporary_directory,
        scale_prefix(source) + "-edge-" + std::to_string(runs.size()) + ".h10r");
    RunWriter writer(path, source.chr1, source.chr2, source.resolution);
    for (const auto &record : chunk) writer.add(record);
    runs.push_back(writer.finish());
    chunk.clear();
}

CsrFiles build_csr(const RunInfo &cells, uint32_t bins, const ScaleOptions &options) {
    uint64_t capacity = options.sort_memory_bytes / (2 * sizeof(CellRecord));
    require(capacity >= 1024, "SCALE sort memory is too small");
    capacity = std::min<uint64_t>(capacity, std::numeric_limits<size_t>::max());
    std::vector<CellRecord> chunk;
    chunk.reserve(static_cast<size_t>(capacity));
    std::vector<RunInfo> runs;
    // The source is already ordered by (row=y,column=x).  Copy that half
    // directly and externally sort only the transposed off-diagonal half.
    // This avoids sorting roughly half of the directed CSR entries.
    const std::string direct_path = join_path(options.temporary_directory,
        scale_prefix(cells) + "-direct.h10r");
    RunWriter direct(direct_path, cells.chr1, cells.chr2, cells.resolution);
    RunReader input(cells);
    CellRecord cell;
    while (input.next(cell)) {
        require(cell.x < bins && cell.y < bins && cell.x <= cell.y,
                "SCALE cell outside cis matrix");
        direct.add(cell);
        if (cell.x != cell.y) {
            chunk.push_back({cell.y, cell.x, cell.count});
            if (chunk.size() == chunk.capacity())
                emit_transpose_chunk(chunk, runs, cells, options);
        }
    }
    emit_transpose_chunk(chunk, runs, cells, options);
    runs.push_back(direct.finish());
    RunInfo directed = bounded_merge(std::move(runs), options.temporary_directory,
                                     scale_prefix(cells) + "-directed",
                                     options.merge_fan_in, {}, true);

    CsrFiles result;
    result.bins = bins;
    result.entries = directed.records;
    result.index = join_path(options.temporary_directory, scale_prefix(cells) + ".csr-index");
    result.data = join_path(options.temporary_directory, scale_prefix(cells) + ".csr-data");
    File index(result.index + ".tmp", "w+b"), data(result.data + ".tmp", "w+b");
    std::array<unsigned char, CSR_HEADER_BYTES> header{};
    std::memcpy(header.data(), "H10C", 4);
    put_u32(header.data() + 4, 1);
    put_u32(header.data() + 8, bins);
    put_u32(header.data() + 12, sizeof(CsrEntry));
    put_u64(header.data() + 16, directed.records);
    index.write(header.data(), header.size());
    data.write(header.data(), header.size());
    std::vector<uint64_t> offsets;
    offsets.reserve(65536);
    std::vector<CsrEntry> entries;
    entries.reserve(65536);
    auto flush_offsets = [&]() {
        if (!offsets.empty()) { index.write(offsets.data(), offsets.size() * 8); offsets.clear(); }
    };
    auto flush_entries = [&]() {
        if (!entries.empty()) { data.write(entries.data(), entries.size() * sizeof(CsrEntry)); entries.clear(); }
    };
    uint32_t next_row = 0;
    uint64_t emitted = 0;
    RunReader reader(directed);
    while (reader.next(cell)) {
        require(cell.y < bins && cell.x < bins, "directed SCALE edge outside matrix");
        while (next_row <= cell.y) {
            offsets.push_back(emitted);
            ++next_row;
            if (offsets.size() == offsets.capacity()) flush_offsets();
        }
        float value = static_cast<float>(cell.count);
        require(std::isfinite(value) && value > 0, "SCALE count cannot be represented as float32");
        entries.push_back({cell.x, value});
        ++emitted;
        if (entries.size() == entries.capacity()) flush_entries();
    }
    while (next_row <= bins) {
        offsets.push_back(emitted);
        ++next_row;
        if (offsets.size() == offsets.capacity()) flush_offsets();
    }
    flush_offsets(); flush_entries();
    require(emitted == directed.records, "CSR entry count mismatch");
    index.flush(); data.flush(); index.close(); data.close();
    atomic_rename(result.index + ".tmp", result.index);
    atomic_rename(result.data + ".tmp", result.data);
    std::remove(directed.path.c_str());
    return result;
}

class CsrMatrix {
  public:
    CsrMatrix(const CsrFiles &files, int threads)
        : bins_(files.bins), threads_(std::max(1, threads)),
          index_(files.index, CSR_HEADER_BYTES + uint64_t(files.bins + 1) * 8),
          data_(files.data, CSR_HEADER_BYTES + files.entries * sizeof(CsrEntry)) {
        offsets_ = reinterpret_cast<const uint64_t *>(index_.bytes() + CSR_HEADER_BYTES);
        values_ = reinterpret_cast<const CsrEntry *>(data_.bytes() + CSR_HEADER_BYTES);
        threads_ = std::min<int>(threads_, std::max<uint32_t>(1, bins_));
        if (files.entries < static_cast<uint64_t>(threads_) * 1024) threads_ = 1;
        row_boundaries_.resize(static_cast<size_t>(threads_) + 1);
        row_boundaries_.front() = 0;
        row_boundaries_.back() = bins_;
        for (int thread = 1; thread < threads_; ++thread) {
            const uint64_t target = files.entries * static_cast<uint64_t>(thread) / threads_;
            const uint64_t *at = std::lower_bound(offsets_, offsets_ + bins_ + 1, target);
            row_boundaries_[static_cast<size_t>(thread)] =
                static_cast<uint32_t>(std::min<uint64_t>(bins_, at - offsets_));
        }
        if (threads_ > 1) pool_.reset(new ThreadPool(threads_));
    }
    uint32_t bins() const { return bins_; }
    uint64_t row_entries(uint32_t row) const { return offsets_[row + 1] - offsets_[row]; }
    void multiply(const double *vector, double *result) {
        auto work = [&](uint32_t begin, uint32_t end) {
            for (uint32_t row = begin; row < end; ++row) {
                double sum = 0;
                for (uint64_t at = offsets_[row]; at < offsets_[row + 1]; ++at)
                    sum += double(values_[at].value) * vector[values_[at].column];
                result[row] = sum;
            }
        };
        if (!pool_) { work(0, bins_); return; }
        std::vector<std::future<void>> futures;
        futures.reserve(static_cast<size_t>(threads_));
        for (int thread = 0; thread < threads_; ++thread) {
            uint32_t begin = row_boundaries_[static_cast<size_t>(thread)];
            uint32_t end = row_boundaries_[static_cast<size_t>(thread + 1)];
            futures.push_back(pool_->submit(work, begin, end));
        }
        for (auto &future : futures) future.get();
    }
  private:
    uint32_t bins_;
    int threads_;
    Mapping index_, data_;
    const uint64_t *offsets_ = nullptr;
    const CsrEntry *values_ = nullptr;
    std::vector<uint32_t> row_boundaries_;
    std::unique_ptr<ThreadPool> pool_;
};

using RowMask = std::vector<uint8_t>;

double balanced_error(const std::vector<double> &b, const RowMask &bad,
                      CsrMatrix &matrix, double tolerance,
                      std::vector<double> &product, std::vector<int> &offenders) {
    matrix.multiply(b.data(), product.data());
    offenders.clear();
    double error = 0;
    for (uint32_t i = 0; i < b.size(); ++i) {
        if (bad[i]) continue;
        double e = std::fabs(product[i] * b[i] - 1.0);
        error = std::max(error, e);
        if (e > tolerance) offenders.push_back(static_cast<int>(i));
    }
    return error;
}

uint64_t excluded_at(const std::vector<std::pair<int, uint64_t>> &distribution, int cutoff) {
    uint64_t excluded = 0;
    for (const auto &entry : distribution) {
        if (entry.first >= cutoff) break;
        excluded += entry.second;
    }
    return excluded;
}

int upward_cutoff(const std::vector<std::pair<int, uint64_t>> &distribution,
                  int current, int requested, int bound) {
    int target = std::min(requested, bound);
    if (target <= current) return current;
    if (excluded_at(distribution, target) != excluded_at(distribution, current)) return target;
    for (const auto &entry : distribution) {
        if (entry.first < target) continue;
        if (entry.first + 1 <= bound) return entry.first + 1;
        break;
    }
    return current;
}

int hinted_cutoff(const std::vector<std::pair<int, uint64_t>> &distribution,
                   uint64_t positive, double fraction, int bound) {
    if (!(fraction > 0) || !positive) return 1;
    const uint64_t wanted = static_cast<uint64_t>(std::ceil(fraction * positive));
    uint64_t cumulative = 0;
    for (const auto &entry : distribution) {
        cumulative += entry.second;
        if (cumulative >= wanted) return std::max(1, std::min(bound, entry.first + 1));
    }
    return bound;
}

bool initialize_warm(const ScaleOptions &options, uint32_t resolution,
                     const std::vector<double> &raw_vc, RowMask &warm_valid,
                     std::vector<double> &seed) {
    const std::vector<float> *coarse = options.warm_norm;
    if (!coarse || coarse->empty() || options.warm_resolution <= resolution ||
        !(options.warm_vc_exponent >= 0) || !std::isfinite(options.warm_vc_exponent))
        return false;
    const uint32_t k = static_cast<uint32_t>(seed.size());
    std::vector<double> coverage_mean(coarse->size(), 0);
    size_t active_group = std::numeric_limits<size_t>::max();
    long double group_sum = 0;
    uint64_t group_count = 0;
    auto finish_group = [&]() {
        if (active_group < coverage_mean.size() && group_count)
            coverage_mean[active_group] = static_cast<double>(group_sum / group_count);
    };
    for (uint32_t i = 0; i < k; ++i) {
        const uint64_t center = uint64_t(i) * resolution + resolution / 2;
        const size_t q = static_cast<size_t>(center / options.warm_resolution);
        if (q != active_group) {
            finish_group();
            active_group = q;
            group_sum = 0;
            group_count = 0;
        }
        if (q < coarse->size() && raw_vc[i] > 0 && std::isfinite(raw_vc[i])) {
            group_sum += raw_vc[i];
            ++group_count;
        }
    }
    finish_group();
    long double warm_log = 0, baseline_log = 0;
    uint64_t valid = 0, positive = 0;
    for (uint32_t i = 0; i < k; ++i) {
        if (!(raw_vc[i] > 0) || !std::isfinite(raw_vc[i])) continue;
        ++positive;
        const uint64_t begin = uint64_t(i) * resolution;
        const uint64_t end = begin + resolution;
        const size_t first = static_cast<size_t>(begin / options.warm_resolution);
        const size_t last = static_cast<size_t>((end - 1) / options.warm_resolution);
        long double weighted_log = 0, weight = 0;
        for (size_t q = first; q <= last && q < coarse->size(); ++q) {
            const float value = (*coarse)[q];
            if (!(value > 0) || !std::isfinite(value)) continue;
            const uint64_t coarse_begin = uint64_t(q) * options.warm_resolution;
            const uint64_t coarse_end = coarse_begin + options.warm_resolution;
            const uint64_t overlap_begin = std::max(begin, coarse_begin);
            const uint64_t overlap_end = std::min(end, coarse_end);
            if (overlap_end <= overlap_begin) continue;
            const uint64_t overlap = overlap_end - overlap_begin;
            weighted_log += static_cast<long double>(overlap) * std::log(value);
            weight += overlap;
        }
        if (!(weight > 0)) continue;
        const size_t q = static_cast<size_t>((begin + resolution / 2) /
                                             options.warm_resolution);
        double correction = 1.0;
        if (q < coverage_mean.size() && coverage_mean[q] > 0)
            correction = std::pow(coverage_mean[q] / raw_vc[i],
                                  options.warm_vc_exponent);
        const double lifted_norm = std::exp(static_cast<double>(weighted_log / weight));
        const double value = correction / lifted_norm;
        if (!(value > 0) || !std::isfinite(value)) continue;
        seed[i] = value;
        warm_valid[i] = 1;
        warm_log += std::log(value);
        baseline_log += std::log(std::sqrt(raw_vc[i]));
        ++valid;
    }
    if (!valid || valid * 2 < positive) return false;
    const double multiplier = std::exp(static_cast<double>((baseline_log - warm_log) / valid));
    for (uint32_t i = 0; i < k; ++i)
        if (warm_valid[i]) seed[i] *= multiplier;
    return true;
}

ScaleResult balance(CsrMatrix &matrix, uint32_t resolution, const ScaleOptions &options) {
    const uint32_t k = matrix.bins();
    ScaleResult result;
    std::vector<int> nonzeros(k);
    std::unordered_map<int, uint64_t> histogram;
    uint64_t positive_count = 0;
    long double sum = 0, sum_squared = 0;
    for (uint32_t i = 0; i < k; ++i) {
        nonzeros[i] = static_cast<int>(std::min<uint64_t>(INT_MAX, matrix.row_entries(i)));
        if (nonzeros[i]) {
            ++histogram[nonzeros[i]];
            ++positive_count;
            sum += nonzeros[i];
            sum_squared += static_cast<long double>(nonzeros[i]) * nonzeros[i];
        }
    }
    if (!positive_count) { result.reason = "matrix has no nonempty rows"; return result; }
    std::vector<std::pair<int, uint64_t>> distribution(histogram.begin(), histogram.end());
    std::sort(distribution.begin(), distribution.end());
    const uint64_t percentile_index = static_cast<uint64_t>(positive_count * 0.2);
    uint64_t cumulative = 0;
    int percentile_cutoff = distribution.front().first;
    for (const auto &entry : distribution) {
        cumulative += entry.second;
        if (cumulative > percentile_index) { percentile_cutoff = entry.first; break; }
    }
    const long double mean = sum / positive_count;
    const long double variance = std::max<long double>(0, sum_squared / positive_count - mean * mean);
    const int zscore_cutoff = std::max(1, static_cast<int>(std::ceil(mean - std::sqrt(variance))));
    const int bound = std::max(percentile_cutoff, zscore_cutoff);
    result.percentile_cutoff = percentile_cutoff;
    result.zscore_cutoff = zscore_cutoff;
    int low = 1, low0 = 1;
    double low_converged = 1000, low_diverged = 0;
    bool converged_cutoff = false, diverged_cutoff = false;
    bool accepted = false, erez_retry = true, row_rescue = true;
    double error_converged = options.tolerance;
    RowMask bad(k), bad_converged(k), warm_valid(k), changing(k);
    std::vector<double> current(k), row(k), column(k), dr(k), dc(k), one(k), b(k);
    std::vector<double> b_converged(k);
    std::vector<double> computed_raw_vc;
    for (uint32_t i = 0; i < k; ++i) one[i] = nonzeros[i] ? 1.0 : 0.0;
    if (!options.raw_vc || options.raw_vc->size() != k) {
        computed_raw_vc.resize(k);
        matrix.multiply(one.data(), computed_raw_vc.data());
    }
    const std::vector<double> &raw_vc = computed_raw_vc.empty()
        ? *options.raw_vc : computed_raw_vc;
    auto reset_from_mask = [&](bool allow_row_rescue) {
        bool removed;
        do {
            for (uint32_t i = 0; i < k; ++i) one[i] = bad[i] ? 0.0 : 1.0;
            matrix.multiply(one.data(), row.data());
            removed = false;
            for (uint32_t i = 0; i < k; ++i) {
                if (!bad[i] && !(row[i] > 0)) { bad[i] = 1; removed = true; }
            }
        } while (removed);
        for (uint32_t i = 0; i < k; ++i) {
            one[i] = bad[i] ? 0.0 : 1.0;
            dr[i] = dc[i] = one[i];
            row[i] *= dr[i];
        }
        row_rescue = allow_row_rescue;
    };
    auto reset_cutoff = [&]() {
        for (uint32_t i = 0; i < k; ++i) {
            bad[i] = nonzeros[i] < low;
        }
        reset_from_mask(true);
    };
    auto initialize_baseline = [&]() {
        std::fill(bad.begin(), bad.end(), uint8_t(0));
        for (uint32_t i = 0; i < k; ++i) {
            bad[i] = !nonzeros[i];
            dr[i] = dc[i] = bad[i] ? 0.0 : std::sqrt(raw_vc[i]);
            current[i] = dr[i];
            row[i] = raw_vc[i] * dr[i];
            one[i] = bad[i] ? 0.0 : 1.0;
        }
        row_rescue = true;
    };
    bool warm_active = initialize_warm(options, resolution, raw_vc, warm_valid, dr);
    if (warm_active) {
        for (uint32_t i = 0; i < k; ++i) {
            bad[i] = !nonzeros[i];
            if (!warm_valid[i]) dr[i] = bad[i] ? 0.0 : std::sqrt(raw_vc[i]);
            dc[i] = dr[i];
            current[i] = dr[i];
            one[i] = bad[i] ? 0.0 : 1.0;
        }
        matrix.multiply(dc.data(), row.data());
        for (uint32_t i = 0; i < k; ++i) row[i] *= dr[i];
        result.used_warm_start = true;
    } else {
        initialize_baseline();
    }
    bool warm_retry_available = warm_active;
    const int coarse_hint = hinted_cutoff(distribution, positive_count,
                                           options.warm_excluded_fraction, bound);
    bool cutoff_hint_available = coarse_hint > 1;
    double change = 10 * (1 + options.tolerance);
    int attempt_iterations = 0;
    const int per_attempt_limit = std::min(1000, options.max_iterations);
    const int warm_attempt_limit = std::min(25, per_attempt_limit);
    const double erez_fraction = 1.0e-4;
    std::vector<double> report(static_cast<size_t>(options.max_iterations) + 10);
    while (change > options.tolerance && result.iterations < options.max_iterations) {
        ++result.iterations;
        ++attempt_iterations;
        for (uint32_t i = 0; i < k; ++i) dr[i] = !bad[i] && row[i] ? dr[i] / row[i] : 0;
        matrix.multiply(dr.data(), column.data());
        for (uint32_t i = 0; i < k; ++i) {
            column[i] *= dc[i];
            dc[i] = !bad[i] && column[i] ? dc[i] / column[i] : 0;
        }
        matrix.multiply(dc.data(), row.data());
        for (uint32_t i = 0; i < k; ++i) row[i] *= dr[i];
        change = 0;
        int changing_rows = 0;
        for (uint32_t i = 0; i < k; ++i) {
            b[i] = std::sqrt(dr[i] * dc[i]);
            changing[i] = 0;
            if (!bad[i]) {
                double delta = std::fabs((b[i] - current[i]) /
                                         (b[i] + current[i] + 1e-300));
                change = std::max(change, delta);
                if (delta > options.tolerance) { ++changing_rows; changing[i] = 1; }
            }
            current[i] = b[i];
        }
        report[result.iterations - 1] = change;

        bool row_sum_failed = false;
        std::vector<int> offenders;
        if (change < options.tolerance) {
            double row_error = balanced_error(b, bad, matrix, options.row_sum_tolerance,
                                              column, offenders);
            row_sum_failed = !(row_error <= options.row_sum_tolerance);
        }
        if (change < options.tolerance && !row_sum_failed) {
            error_converged = change;
            low_converged = low;
            erez_retry = true;
            if (low <= low0) { accepted = true; break; }
            converged_cutoff = true;
            b_converged = b;
            bad_converged = bad;
            if (diverged_cutoff) {
                if (low_converged - low_diverged <= 1) { accepted = true; break; }
                low = static_cast<int>((low_converged + low_diverged) / 2);
            } else {
                low = static_cast<int>(low_converged / 2);
            }
            attempt_iterations = 0;
            change = 10;
            reset_cutoff();
            continue;
        }

        if (row_sum_failed && warm_retry_available && low == 1) {
            warm_retry_available = false;
            warm_active = false;
            result.used_warm_start = false;
            attempt_iterations = 0;
            change = 10;
            initialize_baseline();
            continue;
        }

        if (row_sum_failed && row_rescue && !offenders.empty()) {
            for (int row_id : offenders) {
                bad[static_cast<size_t>(row_id)] = 1;
                one[static_cast<size_t>(row_id)] = 0;
            }
            row_rescue = false;
            attempt_iterations = 0;
            change = 10;
            reset_from_mask(false);
            continue;
        }

        // Match the in-memory implementation's convergence-rate test.
        if (!row_sum_failed) {
            const int attempt_limit = warm_retry_available ? warm_attempt_limit
                                                            : per_attempt_limit;
            if (attempt_iterations <= 5 && attempt_iterations < attempt_limit) continue;
            int previous_report = result.iterations - 6;
            if (attempt_iterations < attempt_limit && previous_report >= 0 &&
                report[result.iterations - 1] * (1 + options.delta) < report[previous_report])
                continue;
        }

        if (warm_retry_available && low == 1) {
            warm_retry_available = false;
            warm_active = false;
            result.used_warm_start = false;
            attempt_iterations = 0;
            change = 10;
            initialize_baseline();
            continue;
        }

        diverged_cutoff = true;
        low_diverged = low;
        bool retry_same_mask = false;
        if (converged_cutoff) {
            if (low_converged - low_diverged <= 1) {
                b = b_converged;
                bad = bad_converged;
                change = error_converged;
                accepted = true;
                break;
            } else if (!row_sum_failed &&
                       static_cast<double>(changing_rows) / positive_count < erez_fraction &&
                       erez_retry) {
                for (uint32_t i = 0; i < k; ++i) {
                    if (bad[i]) continue;
                    if (changing[i]) { bad[i] = 1; one[i] = 0; }
                }
                erez_retry = false;
                retry_same_mask = true;
            } else {
                low = static_cast<int>((low_diverged + low_converged) / 2);
                erez_retry = true;
            }
        } else if (!row_sum_failed &&
                   static_cast<double>(changing_rows) / positive_count < erez_fraction &&
                   erez_retry) {
            for (uint32_t i = 0; i < k; ++i) {
                if (bad[i]) continue;
                if (changing[i]) { bad[i] = 1; one[i] = 0; }
            }
            erez_retry = false;
            retry_same_mask = true;
        } else {
            int next = cutoff_hint_available
                ? upward_cutoff(distribution, low, coarse_hint, bound)
                : upward_cutoff(distribution, low, std::min(2 * low, bound), bound);
            cutoff_hint_available = false;
            if (next <= low) break;
            low = next;
            erez_retry = true;
        }
        if (!retry_same_mask) {
            for (uint32_t i = 0; i < k; ++i) {
                bad[i] = nonzeros[i] < low;
                one[i] = bad[i] ? 0 : 1;
            }
        }
        attempt_iterations = 0;
        change = 10;
        reset_from_mask(true);
        if (low > bound) break;
    }
    if (!accepted) { result.reason = "SCALE did not converge"; return result; }
    result.cutoff = low;
    uint64_t excluded = 0;
    for (uint32_t i = 0; i < k; ++i) excluded += bad[i] && nonzeros[i];
    result.excluded_fraction = static_cast<double>(excluded) / positive_count;

    // Post-process in the existing solver buffers, then correct the final
    // float32 vector in the mapped CSR. This replaces a separate cell-file pass.
    for (uint32_t i = 0; i < k; ++i) one[i] = bad[i] ? 0.0 : 1.0;
    matrix.multiply(one.data(), row.data());
    long double raw = 0;
    for (uint32_t i = 0; i < k; ++i) raw += static_cast<long double>(row[i]) * one[i];
    for (uint32_t i = 0; i < k; ++i) dc[i] = bad[i] ? 0.0 : b[i];
    matrix.multiply(dc.data(), column.data());
    long double normalized = 0;
    for (uint32_t i = 0; i < k; ++i)
        normalized += static_cast<long double>(column[i]) * dc[i];
    const double scale = raw > 0 && normalized > 0
        ? std::sqrt(static_cast<double>(normalized / raw)) : 1.0;
    result.values.resize(k);
    bool valid = false;
    for (uint32_t i = 0; i < k; ++i) {
        result.values[i] = bad[i] ? std::numeric_limits<float>::quiet_NaN()
                                  : static_cast<float>(scale / b[i]);
        valid = valid || (result.values[i] > 0 && std::isfinite(result.values[i]));
    }
    if (valid && raw > 0) {
        for (uint32_t i = 0; i < k; ++i) {
            const float value = result.values[i];
            dc[i] = value > 0 && std::isfinite(value) ? 1.0 / value : 0.0;
        }
        matrix.multiply(dc.data(), column.data());
        long double float_normalized = 0;
        for (uint32_t i = 0; i < k; ++i)
            float_normalized += static_cast<long double>(column[i]) * dc[i];
        if (float_normalized > 0) {
            const float correction = static_cast<float>(
                std::sqrt(static_cast<double>(float_normalized / raw)));
            for (float &value : result.values)
                if (value > 0 && std::isfinite(value)) value *= correction;
        }
    }
    result.success = valid;
    result.reason = valid ? "converged" : "SCALE produced no finite rows";
    return result;
}

} // namespace

ScaleResult scale_cis(const RunInfo &cells, uint32_t bins, const ScaleOptions &options) {
    require(cells.chr1 == cells.chr2, "SCALE requires a cis matrix");
    require(options.threads > 0, "SCALE thread count must be positive");
    if (!cells.records) {
        ScaleResult empty;
        empty.reason = "matrix has no occupied cells";
        return empty;
    }
    make_directory(options.temporary_directory);
    const auto build_begin = std::chrono::steady_clock::now();
    CsrFiles files = build_csr(cells, bins, options);
    const auto build_end = std::chrono::steady_clock::now();
    try {
        ScaleResult result;
        {
            CsrMatrix matrix(files, options.threads);
            result = balance(matrix, cells.resolution, options);
        }
        const auto balance_end = std::chrono::steady_clock::now();
        std::cerr << "  SCALE CSR "
                  << std::chrono::duration<double>(build_end - build_begin).count()
                  << "s, balance "
                  << std::chrono::duration<double>(balance_end - build_end).count()
                  << "s; degree cutoffs percentile=" << result.percentile_cutoff
                  << " zscore=" << result.zscore_cutoff
                  << " selected=" << result.cutoff << '\n';
        std::remove(files.index.c_str());
        std::remove(files.data.c_str());
        return result;
    } catch (...) {
        std::remove(files.index.c_str());
        std::remove(files.data.c_str());
        throw;
    }
}

} // namespace hic10large
