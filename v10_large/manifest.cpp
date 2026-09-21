#include "manifest.h"

#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <unistd.h>

namespace hic10large {

namespace {
unsigned __int128 add128(unsigned __int128 a, unsigned __int128 b) {
    const unsigned __int128 maximum = ~static_cast<unsigned __int128>(0);
    require(a <= maximum - b, "uint128 manifest total overflow");
    return a + b;
}
}

void write_stage_manifest(const StageManifest &m, const std::string &path) {
    std::string temporary = path + ".tmp-" + std::to_string(static_cast<unsigned long long>(getpid()));
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    require(bool(out), "cannot create manifest " + temporary);
    out << "HIC_V10_LARGE_STAGE 1\n";
    out << "source " << std::quoted(m.source) << ' ' << m.source_bytes << ' '
        << hex64(m.source_fingerprint) << "\n";
    out << "resolution " << m.source_resolution << "\n";
    out << "totals " << m.total_records << ' ' << decimal_u128(m.total_weight) << "\n";
    out << "chromosomes " << m.chromosomes.size() << "\n";
    for (size_t i = 0; i < m.chromosomes.size(); ++i)
        out << "chrom " << i << ' ' << m.chromosomes[i].length << ' '
            << std::quoted(m.chromosomes[i].name) << "\n";
    out << "pairs " << m.pairs.size() << "\n";
    for (const auto &p : m.pairs)
        out << "pair " << p.chr1 << ' ' << p.chr2 << ' ' << p.first_record << ' '
            << p.records << ' ' << decimal_u128(p.weight) << ' ' << p.shards.size() << "\n";
    out << "shards " << m.shards.size() << "\n";
    for (const auto &s : m.shards)
        out << "shard " << s.chr1 << ' ' << s.chr2 << ' ' << s.part << ' '
            << s.first_record << ' ' << s.records << ' ' << s.bytes << ' '
            << hex64(s.checksum) << ' ' << decimal_u128(s.weight) << ' '
            << std::quoted(s.path) << "\n";
    out << "end\n";
    out.flush();
    require(bool(out), "cannot write manifest " + temporary);
    out.close();
    atomic_rename(temporary, path);
}

StageManifest read_stage_manifest(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    require(bool(in), "cannot open stage manifest " + path);
    std::string magic;
    unsigned version = 0;
    in >> magic >> version;
    require(magic == "HIC_V10_LARGE_STAGE" && version == 1, "invalid stage manifest");
    StageManifest m;
    std::string tag;
    size_t expected_chromosomes = 0, expected_pairs = 0, expected_shards = 0;
    std::vector<size_t> expected_pair_shards;
    while (in >> tag) {
        if (tag == "source") {
            std::string hash;
            in >> std::quoted(m.source) >> m.source_bytes >> hash;
            m.source_fingerprint = unhex64(hash);
        } else if (tag == "resolution") {
            in >> m.source_resolution;
        } else if (tag == "totals") {
            std::string weight;
            in >> m.total_records >> weight;
            m.total_weight = parse_u128(weight);
        } else if (tag == "chromosomes") {
            in >> expected_chromosomes;
            m.chromosomes.resize(expected_chromosomes);
        } else if (tag == "chrom") {
            size_t id = 0;
            in >> id;
            require(id < m.chromosomes.size(), "chromosome ID outside manifest table");
            in >> m.chromosomes[id].length >> std::quoted(m.chromosomes[id].name);
        } else if (tag == "pairs") {
            in >> expected_pairs;
            m.pairs.reserve(expected_pairs);
        } else if (tag == "pair") {
            PairInfo p;
            std::string weight;
            size_t shard_count = 0;
            in >> p.chr1 >> p.chr2 >> p.first_record >> p.records >> weight >> shard_count;
            p.weight = parse_u128(weight);
            p.shards.reserve(shard_count);
            m.pairs.push_back(std::move(p));
            expected_pair_shards.push_back(shard_count);
        } else if (tag == "shards") {
            in >> expected_shards;
            m.shards.reserve(expected_shards);
        } else if (tag == "shard") {
            ShardInfo s;
            std::string checksum, weight;
            in >> s.chr1 >> s.chr2 >> s.part >> s.first_record >> s.records >> s.bytes
               >> checksum >> weight >> std::quoted(s.path);
            s.checksum = unhex64(checksum);
            s.weight = parse_u128(weight);
            m.shards.push_back(std::move(s));
        } else if (tag == "end") {
            break;
        } else {
            fail("unknown stage manifest field " + tag);
        }
        require(bool(in), "truncated stage manifest " + path);
    }
    require(m.source_resolution && m.chromosomes.size() == expected_chromosomes &&
                m.pairs.size() == expected_pairs && m.shards.size() == expected_shards,
            "incomplete stage manifest " + path);
    std::set<std::string> chromosome_names;
    for (const auto &chromosome : m.chromosomes)
        require(chromosome.length && !chromosome.name.empty() &&
                    chromosome_names.insert(chromosome.name).second,
                "invalid chromosome table in stage manifest");
    std::map<std::pair<uint32_t, uint32_t>, size_t> pair_index;
    for (size_t i = 0; i < m.pairs.size(); ++i) {
        require(m.pairs[i].chr1 < m.chromosomes.size() &&
                    m.pairs[i].chr2 < m.chromosomes.size() &&
                    m.pairs[i].chr1 <= m.pairs[i].chr2,
                "invalid chromosome pair in stage manifest");
        auto key = std::make_pair(m.pairs[i].chr1, m.pairs[i].chr2);
        require(pair_index.emplace(key, i).second, "duplicate pair in stage manifest");
    }
    for (size_t i = 0; i < m.shards.size(); ++i) {
        require(m.shards[i].records <=
                    (UINT64_MAX - STAGED_HEADER_BYTES) / sizeof(StagedRecord) &&
                    m.shards[i].chr1 < m.chromosomes.size() &&
                    m.shards[i].chr2 < m.chromosomes.size() &&
                    m.shards[i].chr1 <= m.shards[i].chr2 && m.shards[i].records &&
                    m.shards[i].bytes == STAGED_HEADER_BYTES +
                                         m.shards[i].records * sizeof(StagedRecord) &&
                    base_name(m.shards[i].path) == m.shards[i].path,
                "invalid shard metadata in stage manifest");
        auto it = pair_index.find({m.shards[i].chr1, m.shards[i].chr2});
        require(it != pair_index.end(), "shard references unknown pair");
        m.pairs[it->second].shards.push_back(i);
    }
    uint64_t next_record = 0;
    unsigned __int128 total_weight = 0;
    for (size_t pair_id = 0; pair_id < m.pairs.size(); ++pair_id) {
        const auto &p = m.pairs[pair_id];
        require(!p.shards.empty() && p.shards.size() == expected_pair_shards[pair_id] &&
                    p.first_record == next_record,
                "invalid pair shard range in stage manifest");
        uint64_t pair_records = 0;
        unsigned __int128 pair_weight = 0;
        uint32_t expected_part = 0;
        for (size_t shard_id : p.shards) {
            const auto &s = m.shards[shard_id];
            require(s.part == expected_part++ && s.first_record == next_record,
                    "noncontiguous shard range in stage manifest");
            require(s.records <= UINT64_MAX - pair_records &&
                        s.records <= UINT64_MAX - next_record,
                    "stage record total overflow");
            pair_records += s.records;
            next_record += s.records;
            pair_weight = add128(pair_weight, s.weight);
        }
        require(pair_records == p.records && pair_weight == p.weight,
                "pair totals disagree with shards");
        total_weight = add128(total_weight, p.weight);
    }
    require(next_record == m.total_records && total_weight == m.total_weight,
            "stage totals disagree with pairs");
    return m;
}

} // namespace hic10large
