#include "pre.h"
#include "common/expected_value.h"
#include "common/little_endian.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
#include <sys/stat.h>
#include <unistd.h>

namespace hic10 {
namespace {
struct Spool {
    FILE *file = nullptr;
    std::string path;
    explicit Spool(const std::string &directory) {
        std::string pattern = directory + "/hic-v10-pairs-XXXXXX";
        std::vector<char> name(pattern.begin(), pattern.end());
        name.push_back(0);
        int fd = mkstemp(name.data());
        check(fd >= 0, "cannot create pair spool in " + directory);
        path = name.data();
        file = fdopen(fd, "w+b");
        if (!file) {
            close(fd);
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
    void reset() {
        check(std::fflush(file) == 0 && ftruncate(fileno(file), 0) == 0 &&
                  fseeko(file, 0, SEEK_SET) == 0,
              "cannot reset pair spool");
    }
    void add(uint32_t x, uint32_t y, float score) {
        Bytes b;
        put(b, x, 4);
        put(b, y, 4);
        put(b, bits(score), 4);
        check(std::fwrite(b.data(), 1, b.size(), file) == b.size(), "cannot write pair spool");
    }
};
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
    Writer writer(output, h, options);
    Spool spool(p.tmpDir);
    auto iterator = open_pair_iterator(input, genome, p.format);
    AlignmentPair pair;
    std::pair<int32_t, int32_t> active{-1, -1};
    std::set<std::pair<int32_t, int32_t>> seen;
    uint64_t records = 0, total = 0, kept = 0;
    bool scores = options.scores, expectedValid = true;
    auto flush = [&]() {
        if (active.first < 0)
            return;
        check(std::fflush(spool.file) == 0, "cannot flush pair spool");
        std::fprintf(stderr, "Writing %s x %s (%llu pairs)\n", genome.at(active.first).name.c_str(),
                     genome.at(active.second).name.c_str(),
                     static_cast<unsigned long long>(records));
        writer.matrix(active.first - 1, active.second - 1, [&](uint8_t unit, uint32_t ri) {
            check(!unit, "direct pre supports BP only");
            check(fseeko(spool.file, 0, SEEK_SET) == 0, "cannot rewind spool");
            uint32_t bin = resolutions[ri];
            struct Acc {
                uint64_t count = 0;
                double score = 0;
                uint32_t firstBits = 0;
                uint64_t n = 0;
            };
            std::map<std::pair<uint32_t, uint32_t>, Acc> cells;
            for (uint64_t i = 0; i < records; ++i) {
                uint32_t x = static_cast<uint32_t>(fread_int32(spool.file)),
                         y = static_cast<uint32_t>(fread_int32(spool.file));
                float value = fread_float(spool.file);
                auto &a = cells[{y / bin, x / bin}];
                if (scores) {
                    if (!a.n)
                        a.firstBits = bits(value);
                    a.score += double(value);
                    ++a.n;
                } else {
                    check(std::isfinite(value) && value >= 0 && std::floor(value) == value &&
                              double(value) < std::ldexp(1.0, 64),
                          "invalid count weight");
                    a.count = plus(a.count, static_cast<uint64_t>(value));
                }
            }
            Matrix m;
            m.scores = scores;
            m.cells.reserve(cells.size());
            for (auto e : cells) {
                uint64_t value = scores
                                     ? (e.second.n == 1 ? e.second.firstBits
                                                        : bits(static_cast<float>(e.second.score)))
                                     : e.second.count;
                if (scores || value)
                    m.cells.push_back({e.first.second, e.first.first, value});
            }
            return m;
        });
        spool.reset();
        records = 0;
        scores = options.scores;
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
            flush();
            check(seen.insert(key).second, "chromosome pair is not contiguous in input");
            active = key;
        }
        float value = pair.score;
        if (!std::isfinite(value) || value < 0 || std::signbit(value) ||
            std::floor(value) != value || double(value) >= std::ldexp(1.0, 64) || value == 0)
            scores = true;
        spool.add(pair.pos1, pair.pos2, value);
        ++records;
        ++kept;
        if (!std::isfinite(value) || value < 0)
            expectedValid = false;
        if (pair.chr1 == pair.chr2 && std::isfinite(value) && value > 0)
            for (size_t i = 0; i < resolutions.size(); ++i)
                ev[i]->add_distance(pair.chr1, pair.pos1 / resolutions[i],
                                    pair.pos2 / resolutions[i], value);
    }
    iterator->close();
    flush();
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
