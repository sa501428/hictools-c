#include "scale.h"

#include "common/thread_pool.h"
#include "io.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <future>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
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

void emit_directed_chunk(std::vector<CellRecord> &chunk, std::vector<RunInfo> &runs,
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
    RunReader input(cells);
    CellRecord cell;
    while (input.next(cell)) {
        require(cell.x < bins && cell.y < bins && cell.x <= cell.y,
                "SCALE cell outside cis matrix");
        // Run records are sorted by (y,x).  For directed CSR, y is the row and
        // x is the column, hence the deliberate field order below.
        chunk.push_back({cell.y, cell.x, cell.count});
        if (cell.x != cell.y) chunk.push_back({cell.x, cell.y, cell.count});
        if (chunk.size() + 2 > chunk.capacity())
            emit_directed_chunk(chunk, runs, cells, options);
    }
    emit_directed_chunk(chunk, runs, cells, options);
    require(!runs.empty(), "cannot balance an empty matrix");
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
          data_(files.data, CSR_HEADER_BYTES + files.entries * sizeof(CsrEntry)),
          pool_(threads_ > 1 ? new ThreadPool(threads_) : nullptr) {
        offsets_ = reinterpret_cast<const uint64_t *>(index_.bytes() + CSR_HEADER_BYTES);
        values_ = reinterpret_cast<const CsrEntry *>(data_.bytes() + CSR_HEADER_BYTES);
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
        for (int thread = 0; thread < threads_; ++thread) {
            uint32_t begin = static_cast<uint32_t>(uint64_t(bins_) * thread / threads_);
            uint32_t end = static_cast<uint32_t>(uint64_t(bins_) * (thread + 1) / threads_);
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
    std::unique_ptr<ThreadPool> pool_;
};

double balanced_error(const std::vector<double> &b, const std::vector<int> &bad,
                      CsrMatrix &matrix, double tolerance,
                      std::vector<int> *offenders = nullptr) {
    std::vector<double> product(b.size());
    matrix.multiply(b.data(), product.data());
    if (offenders) offenders->clear();
    double error = 0;
    for (uint32_t i = 0; i < b.size(); ++i) {
        if (bad[i]) continue;
        double e = std::fabs(product[i] * b[i] - 1.0);
        error = std::max(error, e);
        if (offenders && e > tolerance) offenders->push_back(static_cast<int>(i));
    }
    return error;
}

void postprocess(CsrMatrix &matrix, std::vector<double> &b) {
    const uint32_t k = matrix.bins();
    std::vector<double> one(k), work(k), vector(k);
    for (uint32_t i = 0; i < k; ++i) one[i] = std::isnan(b[i]) ? 0.0 : 1.0;
    matrix.multiply(one.data(), work.data());
    long double raw = 0;
    for (uint32_t i = 0; i < k; ++i) raw += static_cast<long double>(work[i]) * one[i];
    for (uint32_t i = 0; i < k; ++i) vector[i] = std::isnan(b[i]) ? 0.0 : b[i];
    matrix.multiply(vector.data(), work.data());
    long double normalized = 0;
    for (uint32_t i = 0; i < k; ++i)
        normalized += static_cast<long double>(work[i]) * vector[i];
    if (!(raw > 0) || !(normalized > 0)) return;
    double scale = std::sqrt(static_cast<double>(normalized / raw));
    for (double &value : b)
        if (!std::isnan(value)) value = scale / value;
}

ScaleResult balance(CsrMatrix &matrix, const ScaleOptions &options) {
    const uint32_t k = matrix.bins();
    ScaleResult result;
    std::vector<int> nonzeros(k), bad(k), bad_converged(k);
    std::vector<int> positive;
    positive.reserve(k);
    for (uint32_t i = 0; i < k; ++i) {
        nonzeros[i] = static_cast<int>(std::min<uint64_t>(INT_MAX, matrix.row_entries(i)));
        if (nonzeros[i]) positive.push_back(nonzeros[i]);
    }
    if (positive.empty()) { result.reason = "matrix has no nonempty rows"; return result; }
    std::sort(positive.begin(), positive.end());
    long double mean = std::accumulate(positive.begin(), positive.end(), 0.0L) / positive.size();
    long double variance = 0;
    for (int value : positive) { long double d = value - mean; variance += d * d; }
    variance /= positive.size();
    int bound = std::max(positive[static_cast<size_t>(positive.size() * 0.2)],
                         std::max(1, static_cast<int>(std::ceil(mean - std::sqrt(variance)))));
    int low = 1, low0 = 1;
    double low_converged = 1000, low_diverged = 0;
    bool converged_cutoff = false, diverged_cutoff = false;
    bool accepted = false, erez_retry = true, row_rescue = true;
    double error_converged = options.tolerance;
    std::vector<double> current(k), row(k), column(k), dr(k), dc(k), one(k), b(k), previous(k);
    std::vector<double> b_converged(k);
    std::vector<double> raw_vc(k);
    for (uint32_t i = 0; i < k; ++i) one[i] = nonzeros[i] ? 1.0 : 0.0;
    matrix.multiply(one.data(), raw_vc.data());
    auto reset = [&]() {
        for (uint32_t i = 0; i < k; ++i) {
            bad[i] = nonzeros[i] < low;
            one[i] = bad[i] ? 0.0 : 1.0;
            dr[i] = dc[i] = one[i];
        }
        matrix.multiply(dc.data(), row.data());
        for (uint32_t i = 0; i < k; ++i) row[i] *= dr[i];
        row_rescue = true;
    };
    // Preserve the existing first-attempt sqrt(VC) initialization.
    for (uint32_t i = 0; i < k; ++i) {
        bad[i] = !nonzeros[i];
        dr[i] = dc[i] = bad[i] ? 0.0 : std::sqrt(raw_vc[i]);
        current[i] = dr[i];
        row[i] = raw_vc[i] * dr[i];
    }
    double change = 10 * (1 + options.tolerance);
    int attempt_iterations = 0;
    const int per_attempt_limit = std::min(1000, options.max_iterations);
    const double erez_fraction = 1.0e-4;
    std::vector<double> report(static_cast<size_t>(options.max_iterations) + 10);
    while (change > options.tolerance && attempt_iterations < per_attempt_limit &&
           result.iterations < options.max_iterations) {
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
            if (!bad[i]) {
                double delta = std::fabs((b[i] - current[i]) /
                                         (b[i] + current[i] + 1e-300));
                change = std::max(change, delta);
                if (delta > options.tolerance) ++changing_rows;
            }
            previous[i] = current[i];
            current[i] = b[i];
        }
        report[result.iterations - 1] = change;

        bool row_sum_failed = false;
        std::vector<int> offenders;
        if (change < options.tolerance) {
            double row_error = balanced_error(b, bad, matrix, options.row_sum_tolerance,
                                              &offenders);
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
            reset();
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
            for (uint32_t i = 0; i < k; ++i) dr[i] = dc[i] = one[i];
            matrix.multiply(dc.data(), row.data());
            for (uint32_t i = 0; i < k; ++i) row[i] *= dr[i];
            continue;
        }

        // Match the in-memory implementation's convergence-rate test.
        if (!row_sum_failed) {
            if (attempt_iterations <= 5) continue;
            int previous_report = result.iterations - 6;
            if (previous_report >= 0 &&
                report[result.iterations - 1] * (1 + options.delta) < report[previous_report])
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
                       static_cast<double>(changing_rows) / positive.size() < erez_fraction &&
                       erez_retry) {
                for (uint32_t i = 0; i < k; ++i) {
                    if (bad[i]) continue;
                    double delta = std::fabs((b[i] - previous[i]) /
                                             (b[i] + previous[i] + 1e-300));
                    if (delta > options.tolerance) { bad[i] = 1; one[i] = 0; }
                }
                erez_retry = false;
                retry_same_mask = true;
            } else {
                low = static_cast<int>((low_diverged + low_converged) / 2);
                erez_retry = true;
            }
        } else if (!row_sum_failed &&
                   static_cast<double>(changing_rows) / positive.size() < erez_fraction &&
                   erez_retry) {
            for (uint32_t i = 0; i < k; ++i) {
                if (bad[i]) continue;
                double delta = std::fabs((b[i] - previous[i]) /
                                         (b[i] + previous[i] + 1e-300));
                if (delta > options.tolerance) { bad[i] = 1; one[i] = 0; }
            }
            erez_retry = false;
            retry_same_mask = true;
        } else {
            int next = std::min(2 * low, bound);
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
        row_rescue = true;
        attempt_iterations = 0;
        change = 10;
        for (uint32_t i = 0; i < k; ++i) dr[i] = dc[i] = one[i];
        matrix.multiply(dc.data(), row.data());
        for (uint32_t i = 0; i < k; ++i) row[i] *= dr[i];
        if (low > bound) break;
    }
    if (!accepted) { result.reason = "SCALE did not converge"; return result; }
    for (uint32_t i = 0; i < k; ++i)
        if (bad[i]) b[i] = std::numeric_limits<double>::quiet_NaN();
    postprocess(matrix, b);
    result.values.resize(k);
    bool valid = false;
    for (uint32_t i = 0; i < k; ++i) {
        result.values[i] = static_cast<float>(b[i]);
        valid = valid || (result.values[i] > 0 && std::isfinite(result.values[i]));
    }
    result.success = valid;
    result.reason = valid ? "converged" : "SCALE produced no finite rows";
    return result;
}

} // namespace

ScaleResult scale_cis(const RunInfo &cells, uint32_t bins, const ScaleOptions &options) {
    require(cells.chr1 == cells.chr2, "SCALE requires a cis matrix");
    require(options.threads > 0, "SCALE thread count must be positive");
    if (!cells.records) return {false, 0, "matrix has no occupied cells", {}};
    make_directory(options.temporary_directory);
    CsrFiles files = build_csr(cells, bins, options);
    try {
        ScaleResult result;
        {
            CsrMatrix matrix(files, options.threads);
            result = balance(matrix, options);
        }
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
