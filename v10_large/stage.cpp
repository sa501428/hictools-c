#include "stage.h"

#include "io.h"
#include "manifest.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <unistd.h>
#include <zlib.h>

namespace hic10large {
namespace {

unsigned __int128 add_weight(unsigned __int128 total, uint64_t value) {
    const unsigned __int128 maximum = ~static_cast<unsigned __int128>(0);
    require(total <= maximum - value, "uint128 staged weight overflow");
    return total + value;
}

class GzipInput {
  public:
    explicit GzipInput(const std::string &path) : path_(path) {
        file_ = gzopen(path.c_str(), "rb");
        require(file_ != nullptr, "cannot open gzip input " + path);
        gzbuffer(file_, 8 * 1024 * 1024);
        buffer_.resize(8 * 1024 * 1024);
    }
    ~GzipInput() {
        if (file_) gzclose(file_);
    }
    bool read(void *destination, size_t size, bool eof_allowed = false) {
        auto out = static_cast<unsigned char *>(destination);
        size_t copied = 0;
        while (copied < size) {
            if (at_ == available_) {
                int got = gzread(file_, buffer_.data(), static_cast<unsigned>(buffer_.size()));
                int code = Z_OK;
                const char *message = gzerror(file_, &code);
                require(got >= 0 && (code == Z_OK || code == Z_STREAM_END),
                        "gzip read failed for " + path_ + ": " + (message ? message : "unknown"));
                at_ = 0;
                available_ = got > 0 ? static_cast<size_t>(got) : 0;
                if (!available_) {
                    if (!copied && eof_allowed && gzeof(file_)) return false;
                    fail("truncated HBS record in " + path_);
                }
            }
            size_t n = std::min(size - copied, available_ - at_);
            std::memcpy(out + copied, buffer_.data() + at_, n);
            at_ += n;
            copied += n;
        }
        return true;
    }
    void close() {
        if (!file_) return;
        gzFile current = file_;
        file_ = nullptr;
        require(gzclose(current) == Z_OK, "gzip checksum/close failed for " + path_);
    }

  private:
    std::string path_;
    gzFile file_ = nullptr;
    std::vector<unsigned char> buffer_;
    size_t at_ = 0, available_ = 0;
};

std::string shard_name(uint32_t a, uint32_t b, uint32_t part) {
    std::ostringstream out;
    out << "pair-" << std::setfill('0') << std::setw(5) << a << '-' << std::setw(5) << b
        << '-' << std::setw(6) << part << ".h10s";
    return out.str();
}

class ShardWriter {
  public:
    ShardWriter(const std::string &directory, uint32_t resolution, uint32_t a, uint32_t b,
                uint32_t part, uint64_t first_record)
        : directory_(directory), header_{resolution, a, b, part, 0, 14695981039346656037ULL},
          first_record_(first_record), relative_(shard_name(a, b, part)),
          temporary_(join_path(directory, relative_ + ".tmp-" + std::to_string(getpid()))),
          output_(temporary_, "w+b") {
        auto bytes = encode_staged_header(header_);
        output_.write(bytes.data(), bytes.size());
        pending_.reserve(1 << 18);
    }
    ~ShardWriter() {
        if (!finished_) std::remove(temporary_.c_str());
    }
    void add(const StagedRecord &record) {
        require(header_.records < UINT64_MAX, "staged shard record count overflow");
        pending_.push_back(record);
        header_.checksum = fnv1a(&record, sizeof(record), header_.checksum);
        ++header_.records;
        weight_ = add_weight(weight_, record.count);
        if (pending_.size() == pending_.capacity()) flush_records();
    }
    ShardInfo finish() {
        require(!finished_, "staged shard writer finished twice");
        flush_records();
        auto bytes = encode_staged_header(header_);
        output_.seek(0);
        output_.write(bytes.data(), bytes.size());
        output_.flush();
        output_.close();
        std::string final = join_path(directory_, relative_);
        atomic_rename(temporary_, final);
        finished_ = true;
        ShardInfo result;
        result.chr1 = header_.chr1;
        result.chr2 = header_.chr2;
        result.part = header_.part;
        result.first_record = first_record_;
        result.records = header_.records;
        result.bytes = STAGED_HEADER_BYTES + header_.records * sizeof(StagedRecord);
        result.checksum = header_.checksum;
        result.weight = weight_;
        result.path = relative_;
        return result;
    }

  private:
    void flush_records() {
        if (pending_.empty()) return;
        output_.write(pending_.data(), pending_.size() * sizeof(StagedRecord));
        pending_.clear();
    }
    std::string directory_;
    StagedHeader header_;
    uint64_t first_record_ = 0;
    unsigned __int128 weight_ = 0;
    std::string relative_, temporary_;
    File output_;
    std::vector<StagedRecord> pending_;
    bool finished_ = false;
};

uint64_t source_fingerprint(const std::string &path) {
    File input(path, "rb");
    std::vector<unsigned char> buffer(1024 * 1024);
    uint64_t size = file_size(path);
    uint64_t hash = fnv1a(&size, sizeof(size));
    size_t first = static_cast<size_t>(std::min<uint64_t>(buffer.size(), size));
    if (first) {
        input.read(buffer.data(), first);
        hash = fnv1a(buffer.data(), first, hash);
    }
    if (size > first) {
        size_t last = static_cast<size_t>(std::min<uint64_t>(buffer.size(), size - first));
        input.seek(size - last);
        input.read(buffer.data(), last);
        hash = fnv1a(buffer.data(), last, hash);
    }
    return hash;
}

uint64_t complete_fingerprint(const StageManifest &manifest) {
    uint64_t hash = manifest.source_fingerprint;
    hash = fnv1a(&manifest.source_resolution, sizeof(manifest.source_resolution), hash);
    for (const auto &chromosome : manifest.chromosomes) {
        hash = fnv1a(chromosome.name.data(), chromosome.name.size(), hash);
        hash = fnv1a(&chromosome.length, sizeof(chromosome.length), hash);
    }
    for (const auto &shard : manifest.shards) {
        hash = fnv1a(&shard.chr1, sizeof(shard.chr1), hash);
        hash = fnv1a(&shard.chr2, sizeof(shard.chr2), hash);
        hash = fnv1a(&shard.part, sizeof(shard.part), hash);
        hash = fnv1a(&shard.records, sizeof(shard.records), hash);
        hash = fnv1a(&shard.checksum, sizeof(shard.checksum), hash);
    }
    return hash;
}

} // namespace

void stage_hbs(const std::string &input, const std::string &directory, const StageOptions &options) {
    require(options.chunk_records > 0, "--chunk-records must be positive");
    require(!path_exists(join_path(directory, "stage.manifest")),
            "stage.manifest already exists in " + directory);
    make_directory(directory);

    StageManifest manifest;
    manifest.source = input;
    manifest.source_bytes = file_size(input);
    manifest.source_fingerprint = source_fingerprint(input);

    GzipInput gzip(input);
    std::array<unsigned char, 20> fixed{};
    gzip.read(fixed.data(), fixed.size());
    require(std::memcmp(fixed.data(), "HICBS\0\r\n", 8) == 0, "invalid HBS magic");
    require(get_u16(fixed.data() + 8) == 1, "unsupported HBS version");
    require(get_u16(fixed.data() + 10) == 0, "unsupported HBS flags");
    manifest.source_resolution = get_u32(fixed.data() + 12);
    uint32_t chromosome_count = get_u32(fixed.data() + 16);
    require(manifest.source_resolution && manifest.source_resolution <= INT32_MAX,
            "invalid HBS resolution");
    require(chromosome_count <= 65536, "too many HBS chromosomes");
    manifest.chromosomes.reserve(chromosome_count);
    std::set<std::string> names;
    uint64_t header_bytes = fixed.size();
    for (uint32_t i = 0; i < chromosome_count; ++i) {
        unsigned char width[2];
        gzip.read(width, sizeof(width));
        uint16_t n = get_u16(width);
        header_bytes += uint64_t(n) + 10;
        require(n && n <= 4096 && header_bytes <= 16 * 1024 * 1024,
                "invalid HBS chromosome name/header length");
        std::string name(n, '\0');
        gzip.read(&name[0], n);
        unsigned char length[8];
        gzip.read(length, sizeof(length));
        uint64_t bases = get_u64(length);
        require(name.find('\0') == std::string::npos && names.insert(name).second && bases,
                "invalid HBS chromosome table");
        uint64_t bins = bases / manifest.source_resolution +
                        (bases % manifest.source_resolution != 0);
        require(bins <= UINT32_MAX,
                "HBS source resolution requires more than uint32 bins for " + name);
        manifest.chromosomes.push_back({std::move(name), bases});
    }

    std::pair<uint32_t, uint32_t> active{UINT32_MAX, UINT32_MAX};
    std::set<std::pair<uint32_t, uint32_t>> completed;
    std::unique_ptr<ShardWriter> shard;
    PairInfo pair;
    uint32_t part = 0;
    uint64_t shard_records = 0;
    auto finish_shard = [&]() {
        if (!shard) return;
        ShardInfo info = shard->finish();
        pair.shards.push_back(manifest.shards.size());
        manifest.shards.push_back(std::move(info));
        shard.reset();
        shard_records = 0;
    };
    auto finish_pair = [&]() {
        if (active.first == UINT32_MAX) return;
        finish_shard();
        manifest.pairs.push_back(std::move(pair));
        pair = PairInfo{};
        completed.insert(active);
    };

    auto started = std::chrono::steady_clock::now();
    std::array<unsigned char, 14> raw{};
    while (gzip.read(raw.data(), raw.size(), true)) {
        uint32_t a = get_u16(raw.data());
        uint32_t x = get_u32(raw.data() + 2);
        uint32_t b = get_u16(raw.data() + 6);
        uint32_t y = get_u32(raw.data() + 8);
        uint64_t count = get_u16(raw.data() + 12);
        if (count == 65535) {
            unsigned char wide[8];
            gzip.read(wide, sizeof(wide));
            count = get_u64(wide);
            require(count >= 65535, "noncanonical HBS escaped count");
        }
        require(a < manifest.chromosomes.size() && b < manifest.chromosomes.size(),
                "HBS chromosome ID outside table");
        uint64_t xbins = manifest.chromosomes[a].length / manifest.source_resolution +
                         (manifest.chromosomes[a].length % manifest.source_resolution != 0);
        uint64_t ybins = manifest.chromosomes[b].length / manifest.source_resolution +
                         (manifest.chromosomes[b].length % manifest.source_resolution != 0);
        // HBS permits the numeric endpoint only when it is exactly a bin start;
        // fold that redundant endpoint into the final real V10 bin.
        if (x == xbins && manifest.chromosomes[a].length % manifest.source_resolution == 0)
            --x;
        if (y == ybins && manifest.chromosomes[b].length % manifest.source_resolution == 0)
            --y;
        require(x < xbins && y < ybins, "HBS bin index outside chromosome");
        if (a > b || (a == b && x > y)) {
            std::swap(a, b);
            std::swap(x, y);
        }
        std::pair<uint32_t, uint32_t> key{a, b};
        if (key != active) {
            finish_pair();
            require(!completed.count(key), "chromosome pair is not contiguous in HBS input");
            active = key;
            part = 0;
            pair.chr1 = a;
            pair.chr2 = b;
            pair.first_record = manifest.total_records;
        }
        if (!shard || shard_records >= options.chunk_records) {
            finish_shard();
            shard.reset(new ShardWriter(directory, manifest.source_resolution, a, b, part++,
                                        manifest.total_records));
        }
        StagedRecord record{x, y, count};
        require(shard_records < UINT64_MAX && pair.records < UINT64_MAX &&
                    manifest.total_records < UINT64_MAX,
                "HBS record count exceeds uint64");
        shard->add(record);
        ++shard_records;
        ++pair.records;
        pair.weight = add_weight(pair.weight, count);
        ++manifest.total_records;
        manifest.total_weight = add_weight(manifest.total_weight, count);
        if (options.progress_records && manifest.total_records % options.progress_records == 0) {
            double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
            std::cerr << "Staged " << manifest.total_records << " records ("
                      << std::fixed << std::setprecision(1)
                      << manifest.total_records / std::max(1.0, seconds) / 1e6
                      << " M records/s)\n";
        }
    }
    finish_pair();
    gzip.close();
    // Fold every staged shard checksum into the identity. This retains the
    // cheap compressed-file fingerprint while covering middle-of-file content
    // without a second pass over the input.
    manifest.source_fingerprint = complete_fingerprint(manifest);
    write_stage_manifest(manifest, join_path(directory, "stage.manifest"));
    std::cerr << "HBS staging complete: " << manifest.total_records << " records, "
              << manifest.pairs.size() << " chromosome pairs, " << manifest.shards.size()
              << " independently readable shards\n";
}

void inspect_stage(const std::string &manifest_path, bool verify_shards) {
    StageManifest m = read_stage_manifest(manifest_path);
    std::string directory = manifest_path.substr(0, manifest_path.find_last_of('/'));
    if (directory == manifest_path) directory = ".";
    std::cout << "source\t" << m.source << "\nresolution\t" << m.source_resolution
              << "\nrecords\t" << m.total_records << "\nweight\t"
              << decimal_u128(m.total_weight) << "\nchromosomes\t" << m.chromosomes.size()
              << "\npairs\t" << m.pairs.size() << "\nshards\t" << m.shards.size() << "\n";
    for (const auto &p : m.pairs)
        std::cout << "pair\t" << m.chromosomes[p.chr1].name << '\t'
                  << m.chromosomes[p.chr2].name << '\t' << p.first_record << '\t'
                  << p.first_record + p.records << '\t' << p.records << '\t'
                  << decimal_u128(p.weight) << '\t' << p.shards.size() << "\n";
    for (const auto &s : m.shards)
        std::cout << "shard\t" << s.chr1 << '\t' << s.chr2 << '\t' << s.part << '\t'
                  << s.first_record << '\t' << s.first_record + s.records << '\t'
                  << s.records << '\t' << s.path << "\n";
    if (!verify_shards) return;
    for (const auto &s : m.shards) {
        std::string path = join_path(directory, s.path);
        require(file_size(path) == s.bytes, "staged shard size mismatch: " + path);
        File input(path, "rb");
        std::array<unsigned char, STAGED_HEADER_BYTES> bytes{};
        input.read(bytes.data(), bytes.size());
        StagedHeader h = decode_staged_header(bytes.data(), bytes.size());
        require(h.source_resolution == m.source_resolution && h.chr1 == s.chr1 &&
                    h.chr2 == s.chr2 && h.part == s.part && h.records == s.records &&
                    h.checksum == s.checksum,
                "staged shard header mismatch: " + path);
        std::vector<StagedRecord> records(1 << 18);
        uint64_t checksum = 14695981039346656037ULL, remaining = h.records;
        unsigned __int128 weight = 0;
        uint64_t xbins = m.chromosomes[s.chr1].length / m.source_resolution +
                         (m.chromosomes[s.chr1].length % m.source_resolution != 0);
        uint64_t ybins = m.chromosomes[s.chr2].length / m.source_resolution +
                         (m.chromosomes[s.chr2].length % m.source_resolution != 0);
        while (remaining) {
            size_t n = static_cast<size_t>(std::min<uint64_t>(records.size(), remaining));
            input.read(records.data(), n * sizeof(StagedRecord));
            checksum = fnv1a(records.data(), n * sizeof(StagedRecord), checksum);
            for (size_t i = 0; i < n; ++i) {
                require(records[i].x < xbins && records[i].y < ybins &&
                            (s.chr1 != s.chr2 || records[i].x <= records[i].y),
                        "invalid staged record geometry: " + path);
                weight = add_weight(weight, records[i].count);
            }
            remaining -= n;
        }
        require(checksum == h.checksum, "staged shard checksum mismatch: " + path);
        require(weight == s.weight, "staged shard weight mismatch: " + path);
    }
    std::cout << "verification\tok\n";
}

} // namespace hic10large
