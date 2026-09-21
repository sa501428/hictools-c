#include "sort.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <unistd.h>

namespace hic10large {
namespace {

uint64_t key(const CellRecord &record) {
    return (uint64_t(record.y) << 32) | record.x;
}

std::string numbered(const std::string &directory, const std::string &prefix,
                     const char *kind, uint64_t number) {
    std::ostringstream out;
    out << prefix << '-' << kind << '-' << std::setfill('0') << std::setw(8) << number << ".h10r";
    return join_path(directory, out.str());
}

void emit_chunk(std::vector<CellRecord> &chunk, std::vector<RunInfo> &runs,
                const std::string &directory, const std::string &prefix,
                uint32_t a, uint32_t b, uint32_t resolution) {
    if (chunk.empty()) return;
    radix_sort_cells(chunk);
    aggregate_sorted_cells(chunk);
    std::string path = numbered(directory, prefix, "run", runs.size());
    RunWriter writer(path, a, b, resolution);
    for (const auto &record : chunk) writer.add(record);
    runs.push_back(writer.finish());
    chunk.clear();
}

// A coarse output row is usually sparse even when the chromosome has hundreds
// of millions of 1 bp bins.  Keep only touched columns instead of allocating a
// chromosome-width uint64_t array for every rollup task.
class SparseRow {
  public:
    void add(uint32_t key, uint64_t value) {
        if (slots_.empty()) rehash(16);
        if ((used_ + 1) * 10 >= slots_.size() * 7) rehash(slots_.size() * 2);
        size_t at = locate(key);
        if (slots_[at].key == EMPTY) {
            slots_[at] = {key, value};
            touched_.push_back(key);
            ++used_;
        } else {
            require(value <= UINT64_MAX - slots_[at].value,
                    "cell count exceeds uint64 during streaming rollup");
            slots_[at].value += value;
        }
    }

    template <class Emit> void flush(Emit emit) {
        std::sort(touched_.begin(), touched_.end());
        std::vector<size_t> occupied;
        occupied.reserve(touched_.size());
        for (uint32_t key : touched_) {
            size_t at = locate(key);
            require(slots_[at].key == key && slots_[at].value,
                    "invalid sparse rollup row");
            emit(key, slots_[at].value);
            occupied.push_back(at);
        }
        // Open-addressed probe chains must remain intact until all lookups for
        // this row are complete.
        for (size_t at : occupied) slots_[at] = {};
        touched_.clear();
        used_ = 0;
    }

  private:
    static constexpr uint32_t EMPTY = UINT32_MAX;
    struct Slot {
        uint32_t key = EMPTY;
        uint64_t value = 0;
    };
    size_t locate(uint32_t key) const {
        size_t at = (uint64_t(key) * 11400714819323198485ULL) & (slots_.size() - 1);
        while (slots_[at].key != EMPTY && slots_[at].key != key)
            at = (at + 1) & (slots_.size() - 1);
        return at;
    }
    void rehash(size_t size) {
        std::vector<Slot> old = std::move(slots_);
        slots_.assign(size, {});
        for (const Slot &slot : old) {
            if (slot.key == EMPTY) continue;
            size_t at = locate(slot.key);
            slots_[at] = slot;
        }
    }
    std::vector<Slot> slots_;
    std::vector<uint32_t> touched_;
    size_t used_ = 0;
};

} // namespace

std::array<unsigned char, RUN_HEADER_BYTES> encode_run_header(const RunInfo &r) {
    std::array<unsigned char, RUN_HEADER_BYTES> b{};
    std::memcpy(b.data(), "H10R", 4);
    put_u32(b.data() + 4, RUN_VERSION);
    put_u32(b.data() + 8, r.chr1);
    put_u32(b.data() + 12, r.chr2);
    put_u32(b.data() + 16, r.resolution);
    put_u32(b.data() + 20, sizeof(CellRecord));
    put_u64(b.data() + 24, r.records);
    put_u64(b.data() + 32, r.checksum);
    return b;
}

RunInfo decode_run_header(const unsigned char *b, size_t size, const std::string &path) {
    require(size >= RUN_HEADER_BYTES && std::memcmp(b, "H10R", 4) == 0,
            "invalid run magic in " + path);
    require(get_u32(b + 4) == RUN_VERSION && get_u32(b + 20) == sizeof(CellRecord),
            "unsupported run format in " + path);
    RunInfo r;
    r.path = path;
    r.chr1 = get_u32(b + 8);
    r.chr2 = get_u32(b + 12);
    r.resolution = get_u32(b + 16);
    r.records = get_u64(b + 24);
    r.checksum = get_u64(b + 32);
    require(r.resolution && r.chr1 <= r.chr2, "invalid run metadata in " + path);
    return r;
}

RunInfo inspect_run(const std::string &path, bool verify_checksum) {
    File input(path, "rb");
    std::array<unsigned char, RUN_HEADER_BYTES> bytes{};
    input.read(bytes.data(), bytes.size());
    RunInfo r = decode_run_header(bytes.data(), bytes.size(), path);
    require(file_size(path) == RUN_HEADER_BYTES + r.records * sizeof(CellRecord),
            "run length mismatch in " + path);
    if (verify_checksum) {
        std::vector<CellRecord> buffer(65536);
        uint64_t checksum = 14695981039346656037ULL, remaining = r.records;
        uint64_t previous = 0;
        bool first = true;
        while (remaining) {
            size_t n = static_cast<size_t>(std::min<uint64_t>(buffer.size(), remaining));
            input.read(buffer.data(), n * sizeof(CellRecord));
            checksum = fnv1a(buffer.data(), n * sizeof(CellRecord), checksum);
            for (size_t i = 0; i < n; ++i) {
                require(buffer[i].count, "zero count in " + path);
                uint64_t k = key(buffer[i]);
                require(first || previous < k, "duplicate or unordered record in " + path);
                previous = k;
                first = false;
            }
            remaining -= n;
        }
        require(checksum == r.checksum, "run checksum mismatch in " + path);
    }
    return r;
}

void radix_sort_cells(std::vector<CellRecord> &records) {
    if (records.size() < 2) return;
    std::vector<CellRecord> temporary(records.size());
    std::vector<CellRecord> *source = &records, *target = &temporary;
    for (unsigned pass = 0; pass < 8; ++pass) {
        std::array<size_t, 256> offsets{};
        const unsigned shift = pass * 8;
        for (const auto &record : *source)
            ++offsets[(key(record) >> shift) & 255];
        size_t total = 0;
        for (size_t &value : offsets) {
            size_t count = value;
            value = total;
            total += count;
        }
        for (const auto &record : *source)
            (*target)[offsets[(key(record) >> shift) & 255]++] = record;
        std::swap(source, target);
    }
    if (source != &records) records.swap(*source);
}

void aggregate_sorted_cells(std::vector<CellRecord> &records) {
    size_t output = 0;
    for (const auto &record : records) {
        if (!record.count) continue;
        if (output && records[output - 1].x == record.x && records[output - 1].y == record.y) {
            require(record.count <= UINT64_MAX - records[output - 1].count,
                    "cell count exceeds uint64");
            records[output - 1].count += record.count;
        } else {
            records[output++] = record;
        }
    }
    records.resize(output);
}

RunWriter::RunWriter(const std::string &path, uint32_t a, uint32_t b, uint32_t resolution)
    : final_(path), temporary_(path + ".tmp-" + std::to_string(getpid())),
      info_{path, a, b, resolution, 0,
      14695981039346656037ULL}, output_(temporary_, "w+b") {
    auto header = encode_run_header(info_);
    output_.write(header.data(), header.size());
    buffer_.reserve(65536);
}

RunWriter::~RunWriter() {
    if (!finished_) std::remove(temporary_.c_str());
}

void RunWriter::add(const CellRecord &record) {
    require(!finished_ && record.count, "invalid record written to run");
    buffer_.push_back(record);
    info_.checksum = fnv1a(&record, sizeof(record), info_.checksum);
    ++info_.records;
    if (buffer_.size() == buffer_.capacity()) flush();
}

void RunWriter::flush() {
    if (buffer_.empty()) return;
    output_.write(buffer_.data(), buffer_.size() * sizeof(CellRecord));
    buffer_.clear();
}

RunInfo RunWriter::finish() {
    require(!finished_, "run writer finished twice");
    finished_ = true;
    flush();
    auto header = encode_run_header(info_);
    output_.seek(0);
    output_.write(header.data(), header.size());
    output_.flush();
    output_.close();
    atomic_rename(temporary_, final_);
    return info_;
}

RunReader::RunReader(const RunInfo &run, size_t buffer_records)
    : info_(run), input_(run.path, "rb"), buffer_(std::max<size_t>(1, buffer_records)),
      remaining_(run.records) {
    std::array<unsigned char, RUN_HEADER_BYTES> bytes{};
    input_.read(bytes.data(), bytes.size());
    RunInfo actual = decode_run_header(bytes.data(), bytes.size(), run.path);
    require(actual.chr1 == run.chr1 && actual.chr2 == run.chr2 &&
                actual.resolution == run.resolution && actual.records == run.records &&
                actual.checksum == run.checksum,
            "run changed after inspection: " + run.path);
}

bool RunReader::next(CellRecord &record) {
    if (at_ == available_) {
        if (!remaining_) {
            if (!verified_) {
                require(checksum_ == info_.checksum, "run checksum mismatch: " + info_.path);
                verified_ = true;
            }
            return false;
        }
        available_ = static_cast<size_t>(std::min<uint64_t>(buffer_.size(), remaining_));
        input_.read(buffer_.data(), available_ * sizeof(CellRecord));
        checksum_ = fnv1a(buffer_.data(), available_ * sizeof(CellRecord), checksum_);
        remaining_ -= available_;
        at_ = 0;
    }
    record = buffer_[at_++];
    return true;
}

RunInfo merge_runs(const std::vector<RunInfo> &runs, const std::string &output,
                   const std::function<void(const CellRecord &)> &tap) {
    require(!runs.empty(), "cannot merge zero runs");
    uint32_t a = runs.front().chr1, b = runs.front().chr2, resolution = runs.front().resolution;
    struct Head { CellRecord record; size_t run; };
    auto later = [](const Head &x, const Head &y) {
        uint64_t kx = key(x.record), ky = key(y.record);
        return kx != ky ? kx > ky : x.run > y.run;
    };
    std::vector<std::unique_ptr<RunReader>> readers;
    readers.reserve(runs.size());
    std::priority_queue<Head, std::vector<Head>, decltype(later)> heap(later);
    for (size_t i = 0; i < runs.size(); ++i) {
        require(runs[i].chr1 == a && runs[i].chr2 == b && runs[i].resolution == resolution,
                "incompatible runs in merge");
        readers.emplace_back(new RunReader(runs[i]));
        CellRecord record;
        if (readers.back()->next(record)) heap.push({record, i});
    }
    RunWriter writer(output, a, b, resolution);
    while (!heap.empty()) {
        Head head = heap.top();
        heap.pop();
        CellRecord total = head.record;
        CellRecord next;
        if (readers[head.run]->next(next)) heap.push({next, head.run});
        while (!heap.empty() && key(heap.top().record) == key(total)) {
            Head duplicate = heap.top();
            heap.pop();
            require(duplicate.record.count <= UINT64_MAX - total.count,
                    "cell count exceeds uint64 during merge");
            total.count += duplicate.record.count;
            if (readers[duplicate.run]->next(next)) heap.push({next, duplicate.run});
        }
        writer.add(total);
        if (tap) tap(total);
    }
    return writer.finish();
}

RunInfo bounded_merge(std::vector<RunInfo> runs, const std::string &directory,
                      const std::string &prefix, size_t fan_in,
                      const std::function<void(const CellRecord &)> &tap,
                      bool remove_inputs) {
    require(!runs.empty() && fan_in >= 2, "invalid bounded merge");
    uint64_t pass = 0;
    while (runs.size() > fan_in) {
        std::vector<RunInfo> next;
        for (size_t begin = 0; begin < runs.size(); begin += fan_in) {
            size_t end = std::min(runs.size(), begin + fan_in);
            std::vector<RunInfo> group(runs.begin() + begin, runs.begin() + end);
            std::ostringstream name;
            name << prefix << "-merge-" << pass << '-' << (begin / fan_in) << ".h10r";
            next.push_back(merge_runs(group, join_path(directory, name.str())));
            if (remove_inputs || pass > 0)
                for (const auto &run : group) std::remove(run.path.c_str());
        }
        runs.swap(next);
        ++pass;
    }
    RunInfo result = merge_runs(runs, join_path(directory, prefix + "-cells.h10r"), tap);
    if (remove_inputs || pass > 0)
        for (const auto &run : runs) std::remove(run.path.c_str());
    return result;
}

std::vector<RunInfo> sort_staged_shard(const std::string &shard_path,
                                       const std::string &run_directory,
                                       uint32_t output_resolution,
                                       uint64_t memory_bytes,
                                       const std::string &prefix) {
    File input(shard_path, "rb");
    std::array<unsigned char, STAGED_HEADER_BYTES> bytes{};
    input.read(bytes.data(), bytes.size());
    StagedHeader h = decode_staged_header(bytes.data(), bytes.size());
    require(output_resolution >= h.source_resolution &&
                output_resolution % h.source_resolution == 0,
            "output resolution is not a multiple of staged HBS resolution");
    uint32_t factor = output_resolution / h.source_resolution;
    uint64_t capacity = memory_bytes / (2 * sizeof(CellRecord));
    require(capacity >= 1024, "sort memory is too small");
    capacity = std::min<uint64_t>(capacity, std::numeric_limits<size_t>::max());
    std::vector<CellRecord> chunk;
    chunk.reserve(static_cast<size_t>(capacity));
    std::vector<RunInfo> runs;
    std::vector<StagedRecord> buffer(65536);
    uint64_t remaining = h.records;
    uint64_t checksum = 14695981039346656037ULL;
    while (remaining) {
        size_t n = static_cast<size_t>(std::min<uint64_t>(buffer.size(), remaining));
        input.read(buffer.data(), n * sizeof(StagedRecord));
        checksum = fnv1a(buffer.data(), n * sizeof(StagedRecord), checksum);
        for (size_t i = 0; i < n; ++i) {
            chunk.push_back({buffer[i].x / factor, buffer[i].y / factor, buffer[i].count});
            if (chunk.size() == chunk.capacity())
                emit_chunk(chunk, runs, run_directory, prefix, h.chr1, h.chr2,
                           output_resolution);
        }
        remaining -= n;
    }
    require(checksum == h.checksum, "staged shard checksum mismatch: " + shard_path);
    emit_chunk(chunk, runs, run_directory, prefix, h.chr1, h.chr2, output_resolution);
    require(!runs.empty(), "staged shard contains no nonzero records");
    return runs;
}

RunInfo rollup_cell_file(const RunInfo &source,
                         const std::string &output_directory,
                         uint32_t output_resolution,
                         const std::string &prefix) {
    require(output_resolution > source.resolution && output_resolution % source.resolution == 0,
            "invalid rollup resolution");
    uint32_t factor = output_resolution / source.resolution;
    RunWriter output(join_path(output_directory, prefix + "-cells.h10r"),
                     source.chr1, source.chr2, output_resolution);
    RunReader reader(source);
    CellRecord record;
    SparseRow row;
    uint32_t active_row = 0;
    bool have_row = false;
    auto flush_row = [&]() {
        if (!have_row) return;
        row.flush([&](uint32_t x, uint64_t count) {
            output.add({x, active_row, count});
        });
    };
    while (reader.next(record)) {
        const uint32_t y = record.y / factor;
        const uint32_t x = record.x / factor;
        if (!have_row || y != active_row) {
            require(!have_row || y > active_row, "source rows are not ordered during rollup");
            flush_row();
            active_row = y;
            have_row = true;
        }
        row.add(x, record.count);
    }
    flush_row();
    return output.finish();
}

} // namespace hic10large
