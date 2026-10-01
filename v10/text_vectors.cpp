#include "text_vectors.h"
#include "expected_vector.h"
#include "vector_updater.h"
#include "v10_large/vectors.h"
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <tuple>
#include <unistd.h>

namespace hic10 {
namespace {
using Bundle = std::tuple<uint32_t, uint8_t, uint32_t>;
using Key = std::tuple<uint32_t, uint32_t, uint8_t, uint32_t>;
struct Supplied {
    uint8_t unit = 0;
    hic10large::VectorInfo info;
};
struct Workspace {
    std::string path;
    explicit Workspace(const std::string &parent) {
        hic10large::make_directory(parent);
        auto pattern = hic10large::join_path(parent, "hic-v10-norm-input-XXXXXX");
        std::vector<char> name(pattern.begin(), pattern.end()); name.push_back(0);
        check(mkdtemp(name.data()) != nullptr, "cannot create text normalization workspace");
        path = name.data();
    }
    ~Workspace() { std::error_code error; std::filesystem::remove_all(path, error); }
};
uint32_t positive_integer(const std::string &text) {
    check(!text.empty() && text.find_first_not_of("0123456789") == std::string::npos,
          "invalid normalization bin size " + text);
    uint32_t bin = narrow(std::stoull(text));
    check(bin > 0, "normalization bin size must be positive");
    return bin;
}
uint32_t word(const std::string &token) {
    if (token.rfind("bits:", 0) == 0) {
        auto hex = token.substr(5);
        check(hex.size() == 8 && hex.find_first_not_of("0123456789abcdefABCDEF") == std::string::npos,
              "exact float word requires bits:XXXXXXXX");
        return static_cast<uint32_t>(std::stoul(hex, nullptr, 16));
    }
    errno = 0;
    char *end = nullptr;
    float value = std::strtof(token.c_str(), &end);
    check(end != token.c_str() && end == token.c_str() + token.size(),
          "invalid normalization value " + token);
    check(!(errno == ERANGE && !std::isfinite(value)), "normalization value overflows float32");
    return bits(value);
}
bool real_chromosome(const Chromosome &chr) {
    return chr.name != "ALL" && chr.name != "All" && chr.name != "all";
}
std::vector<Supplied> parse_vectors(const std::string &path, Header &header,
                                     const std::string &directory) {
    std::ifstream input(path);
    check(bool(input), "cannot open normalization text " + path);
    std::vector<Supplied> vectors;
    std::set<Key> keys;
    std::unique_ptr<hic10large::VectorFileWriter> writer;
    Supplied active;
    uint64_t count = 0, expected = 0, line_number = 0;
    bool version = false;
    std::string line;
    while (std::getline(input, line)) {
        ++line_number;
        auto comment = line.find('#');
        if (comment != std::string::npos) line.resize(comment);
        std::istringstream row(line);
        std::string token;
        if (!(row >> token)) continue;
        try {
            if (!version) {
                std::string revision, extra;
                check(token == "HIC_NORM_VECTORS" && bool(row >> revision) && revision == "1" &&
                          !(row >> extra), "first content line must be HIC_NORM_VECTORS 1");
                version = true;
            } else if (token == "vector") {
                check(!writer, "missing end before next vector");
                std::string name, chromosome, unit, bin, extra;
                check(bool(row >> std::quoted(name) >> std::quoted(chromosome) >> unit >> bin) &&
                          !(row >> extra), "vector requires TYPE CHROMOSOME UNIT BIN_SIZE");
                check(!name.empty() && name != "NONE", "normalization name must be nonempty and cannot be NONE");
                auto chr = std::find_if(header.chromosomes.begin(), header.chromosomes.end(),
                                       [&](const Chromosome &c) { return c.name == chromosome; });
                check(chr != header.chromosomes.end() && real_chromosome(*chr),
                      "unknown or overview chromosome " + chromosome);
                check(unit == "BP" || unit == "FRAG", "normalization unit must be BP or FRAG");
                auto norm = std::find(header.norms.begin(), header.norms.end(), name);
                active = {};
                active.unit = unit == "FRAG";
                active.info.norm = narrow(norm - header.norms.begin());
                if (norm == header.norms.end()) header.norms.push_back(name);
                active.info.chr = narrow(chr - header.chromosomes.begin());
                active.info.resolution = positive_integer(bin);
                active.info.ri = header.resolution(active.unit, active.info.resolution);
                check(keys.insert({active.info.norm, active.info.chr, active.unit, active.info.ri}).second,
                      "duplicate normalization vector");
                expected = header.bins(active.info.chr, active.unit, active.info.ri);
                count = 0;
                writer.reset(new hic10large::VectorFileWriter(
                    hic10large::join_path(directory, "vector-" + std::to_string(vectors.size()) + ".h10w"), active.info));
            } else if (token == "end") {
                std::string extra;
                check(writer && !(row >> extra), "unexpected end or extra fields after end");
                check(count == expected, "normalization vector length " + std::to_string(count) +
                      " differs from chromosome bin count " + std::to_string(expected));
                active.info = writer->finish(); writer.reset();
                vectors.push_back(std::move(active));
            } else {
                check(bool(writer), "value outside a vector block");
                do {
                    check(count < expected, "too many values in normalization vector");
                    writer->add(word(token)); ++count;
                } while (row >> token);
            }
        } catch (const std::exception &e) {
            throw std::runtime_error(path + ":" + std::to_string(line_number) + ": " + e.what());
        }
    }
    check(input.eof(), "cannot read normalization text " + path);
    check(version && !writer && !vectors.empty(), "empty normalization text or missing final end");
    return vectors;
}
struct Accumulator {
    std::vector<long double> actual;
    std::map<uint32_t, long double> observed;
    std::map<uint32_t, uint32_t> bins;
    explicit Accumulator(uint32_t maximum) : actual(maximum, 0) {}
    void add(uint32_t chr, uint32_t n, uint32_t x, uint32_t y, long double value) {
        check(x <= y && y < n, "expected-value contact outside chromosome bins");
        if (!std::isfinite(value) || value <= 0) return;
        bins[chr] = n;
        actual[y - x] += value;
        observed[chr] += value;
    }
    Vector finish(uint8_t kind, uint32_t norm, uint8_t unit, uint32_t ri) {
        return finish_expected(actual, observed, bins, kind, norm, unit, ri, true);
    }
};
} // namespace
void add_text_vectors_v10(const std::string &path, const std::string &text_path,
                         int level, const std::string &temporary_directory,
                         const ExpectedCellStream &stream) {
    Reader reader(path);
    Header header = reader.header();
    Workspace work(temporary_directory);
    auto supplied = parse_vectors(text_path, header, work.path);
    std::set<Bundle> retained;
    for (const auto &e : reader.vector_entries())
        if (e.kind != 1) retained.insert({e.norm, e.unit, e.ri});
    std::map<Bundle, std::vector<Supplied>> bundles;
    for (auto &v : supplied) {
        auto key = Bundle{v.info.norm, v.unit, v.info.ri};
        if (retained.count(key)) {
            std::cerr << "Preserving existing " << header.norms[v.info.norm] << ' '
                      << (v.unit ? "FRAG " : "BP ") << v.info.resolution << '\n';
        } else bundles[key].push_back(std::move(v));
    }
    if (bundles.empty()) return;
    std::set<uint32_t> cis;
    for (auto key : reader.matrices())
        if (key.chr1 == key.chr2 && real_chromosome(header.chromosomes[key.chr1])) cis.insert(key.chr1);
    auto cells = [&](uint32_t chr, uint8_t unit, uint32_t ri, const ExpectedCellConsumer &consume) {
        if (!cis.count(chr)) return;
        if (stream) stream(reader, chr, unit, ri, consume);
        else {
            Matrix matrix = reader.matrix(chr, chr, unit, ri);
            for (auto cell : matrix.cells) {
                long double value = matrix.scores ? floating(static_cast<uint32_t>(cell.value))
                                                 : static_cast<long double>(cell.value);
                consume(cell.x, cell.y, value);
            }
        }
    };
    VectorOutput output(reader, header, path, level);
    std::set<std::pair<uint8_t, uint32_t>> raw_done;
    for (const auto &bundle : bundles) {
        uint32_t norm = std::get<0>(bundle.first), ri = std::get<2>(bundle.first);
        uint8_t unit = std::get<1>(bundle.first);
        uint32_t maximum = 0;
        for (uint32_t chr = 0; chr < header.chromosomes.size(); ++chr)
            maximum = std::max(maximum, narrow(header.bins(chr, unit, ri)));
        Accumulator expected(maximum);
        for (const auto &v : bundle.second) {
            hic10large::VectorFileReader input(v.info);
            auto values = input.read(0, narrow(v.info.words));
            const uint32_t chr = v.info.chr, n = narrow(v.info.words);
            cells(chr, unit, ri, [&](uint32_t x, uint32_t y, long double value) {
                check(x < n && y < n, "contact outside supplied normalization vector");
                float a = floating(values[x]), b = floating(values[y]);
                if (a > 0 && b > 0 && std::isfinite(a) && std::isfinite(b))
                    expected.add(chr, n, x, y, value / (static_cast<long double>(a) * b));
            });
            Vector vector;
            vector.norm = norm; vector.chr = chr; vector.unit = unit; vector.ri = ri;
            vector.values = std::move(values);
            output.add(vector);
        }
        output.add(expected.finish(2, norm, unit, ri));
        // Raw expected requires all available cis chromosomes, including those
        // omitted from this user-supplied normalization.
        if (!output.contains(1, 0, 0, unit, ri) && raw_done.insert({unit, ri}).second) {
            Accumulator raw(maximum);
            for (uint32_t chr : cis) {
                uint32_t n = narrow(header.bins(chr, unit, ri));
                cells(chr, unit, ri, [&](uint32_t x, uint32_t y, long double value) {
                    raw.add(chr, n, x, y, value);
                });
            }
            output.add(raw.finish(1, 0, unit, ri));
        }
        std::cerr << "Added supplied " << header.norms[norm] << ' ' << (unit ? "FRAG " : "BP ")
                  << header.resolutions[unit][ri].bin << " with normalized expected values\n";
    }
    output.finish(header);
}
} // namespace hic10
