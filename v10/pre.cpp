#include "pre.h"
#include "common/thread_pool.h"
#include "common/expected_value.h"
#include "common/little_endian.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <future>
#include <memory>
#include <set>
#include <sys/stat.h>
#include <unistd.h>

namespace hic10 {
namespace {
struct Spool {
    FILE *file = nullptr;
    std::string path;
    explicit Spool(const std::string &directory, const std::string &prefix = "hic-v10-pairs-") {
        std::string pattern = directory + "/" + prefix + "XXXXXX";
        std::vector<char> name(pattern.begin(), pattern.end());
        name.push_back(0);
        int fd = mkstemp(name.data());
        check(fd >= 0, "cannot create pair spool in " + directory);
        path = name.data();
        file = fdopen(fd, "w+b");
        if (!file) {
            ::close(fd);
            std::remove(path.c_str());
            throw std::runtime_error("V10: cannot open pair spool");
        }
    }
    ~Spool() {
        if (file)
            std::fclose(file);
        if (!path.empty())
            std::remove(path.c_str());
    }
    void close() {
        if (!file)
            return;
        check(std::fflush(file) == 0, "cannot flush spool");
        auto current = file;
        file = nullptr;
        check(std::fclose(current) == 0, "cannot close spool");
    }
    void add(uint32_t x, uint32_t y, const AlignmentPair& pair) {
        unsigned char bytes[16];
        uint64_t fields[3] = {x, y, pair.has_exact_count ? pair.exact_count : bits(pair.score)};
        const unsigned widths[3] = {4, 4, pair.has_exact_count ? 8u : 4u};
        size_t offset = 0;
        for (unsigned field = 0; field < 3; ++field)
            for (unsigned lane = 0; lane < widths[field]; ++lane)
                bytes[offset++] = static_cast<unsigned char>(fields[field] >> (8 * lane));
        check(std::fwrite(bytes, 1, offset, file) == offset, "cannot write pair spool");
    }
};

struct PairJob {
    uint32_t chr1 = 0, chr2 = 0;
    uint64_t records = 0;
    bool scores = false, exact = false;
    std::shared_ptr<Spool> input;
};

struct ResolutionRange {
    uint64_t offset = 0, cells = 0;
};

struct PreparedPair {
    uint32_t chr1 = 0, chr2 = 0;
    uint64_t records = 0;
    bool scores = false;
    std::shared_ptr<Spool> matrices;
    std::vector<ResolutionRange> ranges;
};

uint64_t read_u64(FILE *file, const char *message) {
    unsigned char bytes[8];
    check(std::fread(bytes, 1, sizeof(bytes), file) == sizeof(bytes), message);
    uint64_t value = 0;
    for (unsigned i = 0; i < 8; ++i)
        value |= uint64_t(bytes[i]) << (8 * i);
    return value;
}

PreparedPair prepare_pair(PairJob job, const std::vector<Resolution> &resolutions,
                          const std::string &tmpDir) {
    std::unique_ptr<FILE, decltype(&fclose)> input(fopen(job.input->path.c_str(), "rb"), &fclose);
    check(bool(input), "cannot reopen pair spool");
    auto matrices = std::make_shared<Spool>(tmpDir, "hic-v10-matrices-");
    PreparedPair prepared;
    prepared.chr1 = job.chr1;
    prepared.chr2 = job.chr2;
    prepared.records = job.records;
    prepared.scores = job.scores;
    prepared.matrices = matrices;
    prepared.ranges.resize(resolutions.size());

    struct Acc {
        uint64_t count = 0;
        double score = 0;
        uint32_t firstBits = 0;
        uint64_t n = 0;
    };
    for (uint32_t ri = 0; ri < resolutions.size(); ++ri) {
        if (resolutions[ri].mode)
            continue;
        uint32_t bin = resolutions[ri].bin;
        check(fseeko(input.get(), 0, SEEK_SET) == 0, "cannot rewind pair spool");
        std::map<std::pair<uint32_t, uint32_t>, Acc> cells;
        for (uint64_t i = 0; i < job.records; ++i) {
            uint32_t x = static_cast<uint32_t>(fread_int32(input.get()));
            uint32_t y = static_cast<uint32_t>(fread_int32(input.get()));
            uint64_t stored = 0;
            if (job.exact) {
                stored = read_u64(input.get(), "truncated count spool");
            } else {
                stored = bits(fread_float(input.get()));
            }
            uint32_t word = static_cast<uint32_t>(stored);
            float value;
            std::memcpy(&value, &word, sizeof(value));
            if (job.exact)
                value = static_cast<float>(stored);
            auto &a = cells[{y / bin, x / bin}];
            if (job.scores) {
                if (!a.n)
                    a.firstBits = bits(value);
                a.score += job.exact ? static_cast<double>(stored) : double(value);
                ++a.n;
            } else {
                check(job.exact ||
                          (std::isfinite(value) && value >= 0 && std::floor(value) == value &&
                           double(value) < std::ldexp(1.0, 64)),
                      "invalid count weight");
                a.count = plus(a.count, job.exact ? stored : static_cast<uint64_t>(value));
            }
        }

        auto offset = ftello(matrices->file);
        check(offset >= 0, "cannot determine prepared matrix offset");
        uint64_t count = 0;
        Bytes encoded;
        encoded.reserve(16);
        for (auto it = cells.begin(); it != cells.end();) {
            const auto &entry = *it;
            uint64_t value = job.scores
                                 ? (entry.second.n == 1
                                        ? entry.second.firstBits
                                        : bits(static_cast<float>(entry.second.score)))
                                 : entry.second.count;
            if (job.scores || value) {
                encoded.clear();
                put(encoded, entry.first.second, 4);
                put(encoded, entry.first.first, 4);
                put(encoded, value, 8);
                check(std::fwrite(encoded.data(), 1, encoded.size(), matrices->file) ==
                          encoded.size(),
                      "cannot write prepared matrix spool");
                ++count;
            }
            it = cells.erase(it);
        }
        prepared.ranges[ri] = {static_cast<uint64_t>(offset), count};
    }
    matrices->close();
    return prepared;
}
} // namespace
void pre(const std::string &input, const std::string &output, const std::string &genomeSpec,
         const Options &options, const PreOptions &p) {
    struct stat in{}, out{};
    check(stat(input.c_str(), &in) == 0, "cannot stat pair input");
    check(stat(output.c_str(), &out) != 0 || in.st_dev != out.st_dev || in.st_ino != out.st_ino,
          "input and output must differ");
    Genome genome = Genome::from_spec(genomeSpec);
    Header h;
    h.genome = p.genome.empty() && Genome::is_builtin(genomeSpec) ? genomeSpec : p.genome;
    h.attributes = {{"software", "hictools-c hic_v10"}, {"hicFileScaling", "1.0"}};
    // V10 has no special whole-genome coordinate system. The direct producer
    // advertises real chromosomes only; conversion preserves legacy All data.
    for (auto c : genome.chromosomes_without_all())
        h.chromosomes.push_back({c->name, static_cast<uint64_t>(c->length), {}});
    auto resolutions = p.resolutions;
    if (resolutions.empty())
        for (auto r : DEFAULT_BP_RESOLUTIONS)
            resolutions.push_back(r);
    std::sort(resolutions.begin(), resolutions.end());
    check(std::adjacent_find(resolutions.begin(), resolutions.end()) == resolutions.end(),
          "duplicate resolution");
    std::vector<std::unique_ptr<ExpectedValueCalculation>> ev;
    for (auto r : resolutions) {
        check(r > 0 && r <= INT32_MAX, "invalid resolution");
        h.resolutions[0].push_back({r});
        ev.emplace_back(new ExpectedValueCalculation(genome, static_cast<int>(r)));
    }
    check(options.threads > 0 && options.threads <= 256, "invalid writer thread count");
    auto pool = std::make_shared<ThreadPool>(options.threads);
    Options writerOptions = options;
    // Direct preprocessing defines derived targets by exact source aggregation;
    // unlike V9 conversion, there is no pre-existing target matrix to verify.
    writerOptions.verifyDerived = false;
    Writer writer(output, h, writerOptions, pool);
    const auto matrixResolutions = writer.header().resolutions[0];
    auto iterator = open_pair_iterator(input, genome, p.format);
    if (auto source = iterator->source_resolution())
        for (auto r : resolutions)
            check(r % source == 0, "HBS output resolutions must be multiples of the source resolution");
    AlignmentPair pair;
    std::pair<int32_t, int32_t> active{-1, -1};
    std::set<std::pair<int32_t, int32_t>> seen;
    uint64_t records = 0, total = 0, kept = 0;
    bool scores = options.scores, expectedValid = true;
    bool exactSpool = false;
    std::shared_ptr<Spool> spool;
    std::deque<std::future<PreparedPair>> pending;
    const size_t readAhead = p.readAhead ? p.readAhead : options.threads;
    check(readAhead > 0 && readAhead <= 256, "invalid chromosome-pair read-ahead count");

    auto consume = [&]() {
        PreparedPair prepared = pending.front().get();
        pending.pop_front();
        std::unique_ptr<FILE, decltype(&fclose)> matrices(
            fopen(prepared.matrices->path.c_str(), "rb"), &fclose);
        check(bool(matrices), "cannot reopen prepared matrix spool");
        std::fprintf(stderr, "Writing %s x %s (%llu pairs)\n",
                     genome.at(prepared.chr1 + 1).name.c_str(),
                     genome.at(prepared.chr2 + 1).name.c_str(),
                     static_cast<unsigned long long>(prepared.records));
        writer.matrix(prepared.chr1, prepared.chr2, [&](uint8_t unit, uint32_t ri) {
            check(!unit && ri < prepared.ranges.size(), "direct pre supports BP only");
            const auto range = prepared.ranges[ri];
            check(range.cells <= std::numeric_limits<size_t>::max(),
                  "prepared matrix is too large for this platform");
            check(fseeko(matrices.get(), static_cast<off_t>(range.offset), SEEK_SET) == 0,
                  "cannot seek prepared matrix spool");
            Matrix matrix;
            matrix.scores = prepared.scores;
            matrix.cells.reserve(static_cast<size_t>(range.cells));
            for (uint64_t i = 0; i < range.cells; ++i) {
                uint32_t x = static_cast<uint32_t>(fread_int32(matrices.get()));
                uint32_t y = static_cast<uint32_t>(fread_int32(matrices.get()));
                matrix.cells.push_back(
                    {x, y, read_u64(matrices.get(), "truncated prepared matrix spool")});
            }
            return matrix;
        });
    };

    auto submit = [&]() {
        if (active.first < 0)
            return;
        spool->close();
        PairJob job;
        job.chr1 = static_cast<uint32_t>(active.first - 1);
        job.chr2 = static_cast<uint32_t>(active.second - 1);
        job.records = records;
        job.scores = scores;
        job.exact = exactSpool;
        job.input = std::move(spool);
        pending.push_back(pool->submit(
            [job = std::move(job), matrixResolutions, directory = p.tmpDir]() mutable {
                return prepare_pair(std::move(job), matrixResolutions, directory);
            }));
        records = 0;
        scores = options.scores;
        if (pending.size() >= readAhead)
            consume();
    };
    while (iterator->next(pair)) {
        ++total;
        if (!pair.valid() || pair.chr1 <= 0 || pair.chr2 <= 0)
            continue;
        if (std::min(pair.mapq1, pair.mapq2) < p.mapq)
            continue;
        if (p.intra && pair.chr1 != pair.chr2)
            continue;
        if (p.nearDiagonal && pair.chr1 == pair.chr2 &&
            std::abs(int64_t(pair.pos1) - pair.pos2) > 10000000)
            continue;
        check(pair.pos1 >= 0 && pair.pos2 >= 0 && pair.pos1 < genome.at(pair.chr1).length &&
                  pair.pos2 < genome.at(pair.chr2).length,
              "pair position outside chromosome (V10 coordinates are half-open)");
        if (pair.chr1 > pair.chr2 || (pair.chr1 == pair.chr2 && pair.pos1 > pair.pos2)) {
            std::swap(pair.chr1, pair.chr2);
            std::swap(pair.pos1, pair.pos2);
        }
        std::pair<int32_t, int32_t> key{pair.chr1, pair.chr2};
        if (key != active) {
            submit();
            check(seen.insert(key).second, "chromosome pair is not contiguous in input");
            active = key;
            spool = std::make_shared<Spool>(p.tmpDir);
        }
        float value = pair.score;
        if (!records) exactSpool = pair.has_exact_count;
        check(exactSpool == pair.has_exact_count, "mixed count/score encodings within input pair");
        if (!pair.has_exact_count && (!std::isfinite(value) || value < 0 || std::signbit(value) ||
            std::floor(value) != value || double(value) >= std::ldexp(1.0, 64) || value == 0))
            scores = true;
        spool->add(pair.pos1, pair.pos2, pair);
        ++records;
        ++kept;
        if (!std::isfinite(value) || value < 0)
            expectedValid = false;
        if (pair.chr1 == pair.chr2 && std::isfinite(value) && value > 0)
            for (size_t i = 0; i < resolutions.size(); ++i)
                ev[i]->add_distance(pair.chr1, pair.pos1 / resolutions[i],
                                    pair.pos2 / resolutions[i], pair.has_exact_count ? static_cast<double>(pair.exact_count) : value);
    }
    iterator->close();
    submit();
    while (!pending.empty())
        consume();
    std::vector<Vector> vectors;
    if (expectedValid)
        for (uint32_t ri = 0; ri < ev.size(); ++ri)
            if (ev[ri]->has_data()) {
                ev[ri]->compute_density();
                Vector v;
                v.kind = 1;
                v.ri = ri;
                for (auto value : ev[ri]->density())
                    v.values.push_back(bits(value));
                uint64_t n = 0;
                for (uint32_t c = 0; c < h.chromosomes.size(); ++c)
                    n = std::max(n, h.bins(c, 0, ri));
                v.values.resize(n, 0x7fc00000);
                for (auto scale : ev[ri]->chr_scale_factors())
                    v.scales.emplace(scale.first - 1, bits(scale.second));
                vectors.push_back(std::move(v));
            }
    writer.finish(vectors);
    std::fprintf(stderr, "V10 complete: %llu / %llu pairs retained\n",
                 static_cast<unsigned long long>(kept), static_cast<unsigned long long>(total));
}
} // namespace hic10
