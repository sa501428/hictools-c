#include "genome.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <cstdio>

// --------------------------------------------------------------------------
// Static factory methods
// --------------------------------------------------------------------------

namespace {
using ChromList = std::vector<std::pair<std::string, int64_t>>;

const ChromList* builtin_chromosomes(const std::string& raw_id) {
    std::string id = raw_id;
    std::transform(id.begin(), id.end(), id.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    static const ChromList hg19 = {
        {"1",249250621},{"2",243199373},{"3",198022430},{"4",191154276},
        {"5",180915260},{"6",171115067},{"7",159138663},{"8",146364022},
        {"9",141213431},{"10",135534747},{"11",135006516},{"12",133851895},
        {"13",115169878},{"14",107349540},{"15",102531392},{"16",90354753},
        {"17",81195210},{"18",78077248},{"19",59128983},{"20",63025520},
        {"21",48129895},{"22",51304566},{"X",155270560},{"Y",59373566},
        {"MT",16569}
    };
    static const ChromList hg38 = {
        {"chr1",248956422},{"chr2",242193529},{"chr3",198295559},{"chr4",190214555},
        {"chr5",181538259},{"chr6",170805979},{"chr7",159345973},{"chr8",145138636},
        {"chr9",138394717},{"chr10",133797422},{"chr11",135086622},{"chr12",133275309},
        {"chr13",114364328},{"chr14",107043718},{"chr15",101991189},{"chr16",90338345},
        {"chr17",83257441},{"chr18",80373285},{"chr19",58617616},{"chr20",64444167},
        {"chr21",46709983},{"chr22",50818468},{"chrX",156040895},{"chrY",57227415},
        {"chrM",16569}
    };
    static const ChromList mm9 = {
        {"chr1",197195432},{"chr2",181748087},{"chr3",159599783},{"chr4",155630120},
        {"chr5",152537259},{"chr6",149517037},{"chr7",152524553},{"chr8",131738871},
        {"chr9",124076172},{"chr10",129993255},{"chr11",121843856},{"chr12",121257530},
        {"chr13",120284312},{"chr14",125194864},{"chr15",103494974},{"chr16",98319150},
        {"chr17",95272651},{"chr18",90772031},{"chr19",61342430},{"chrX",166650296},
        {"chrY",15902555},{"chrM",16300}
    };
    static const ChromList mm10 = {
        {"chr1",195471971},{"chr2",182113224},{"chr3",160039680},{"chr4",156508116},
        {"chr5",151834684},{"chr6",149736546},{"chr7",145441459},{"chr8",129401213},
        {"chr9",124595110},{"chr10",130694993},{"chr11",122082543},{"chr12",120129022},
        {"chr13",120421639},{"chr14",124902244},{"chr15",104043685},{"chr16",98207768},
        {"chr17",94987271},{"chr18",90702639},{"chr19",61431566},{"chrM",16299},
        {"chrX",171031299},{"chrY",91744698}
    };

    if (id == "hg19") return &hg19;
    if (id == "hg38") return &hg38;
    if (id == "mm9") return &mm9;
    if (id == "mm10") return &mm10;
    return nullptr;
}
}

Genome Genome::from_chrom_sizes(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open chrom.sizes file: " + path);

    std::vector<std::pair<std::string,int64_t>> chroms;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string name;
        int64_t length;
        if (!(ss >> name >> length)) continue;
        chroms.emplace_back(name, length);
    }
    if (chroms.empty()) throw std::runtime_error("No chromosomes found in: " + path);
    return from_list(chroms);
}

bool Genome::is_builtin(const std::string& id) {
    return builtin_chromosomes(id) != nullptr;
}

Genome Genome::from_spec(const std::string& path_or_id) {
    if (const auto* chroms = builtin_chromosomes(path_or_id)) return from_list(*chroms);
    return from_chrom_sizes(path_or_id);
}

Genome Genome::from_list(const std::vector<std::pair<std::string,int64_t>>& chroms) {
    Genome g;
    int64_t total_len = 0;
    for (auto& [name, len] : chroms) total_len += len;

    // Index 0: synthetic "All" chromosome (whole-genome view)
    Chromosome all;
    all.name   = "ALL";
    all.index  = 0;
    // The synthetic All chromosome is expressed in kilobases in Juicer.
    // Whole-genome contact positions are likewise converted to kb below.
    all.length = total_len / 1000;
    g.chromosomes_.push_back(all);

    int idx = 1;
    for (auto& [name, len] : chroms) {
        Chromosome c;
        c.name   = name;
        c.index  = idx++;
        c.length = len;
        g.chromosomes_.push_back(c);
    }
    g.build_index();
    return g;
}

// --------------------------------------------------------------------------
// Name cleaning and lookup
// --------------------------------------------------------------------------

std::string Genome::clean_name(const std::string& name) {
    std::string s = name;
    // lowercase
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

void Genome::build_index() {
    name_to_idx_.clear();
    chr_offsets_.resize(chromosomes_.size(), 0);

    int64_t offset = 0;
    for (auto& c : chromosomes_) {
        // Insert clean name (lower-case, with 'chr' if present)
        name_to_idx_[clean_name(c.name)] = c.index;

        // Also insert without 'chr' prefix and with 'chr' prefix for flexibility
        std::string cn = clean_name(c.name);
        if (cn.size() > 3 && cn.substr(0, 3) == "chr") {
            name_to_idx_[cn.substr(3)] = c.index;
        } else if (cn != "all") {
            name_to_idx_["chr" + cn] = c.index;
        }

        chr_offsets_[c.index] = offset;
        if (c.name != "ALL") offset += c.length;
    }
}

int Genome::index_of(const std::string& name) const {
    auto it = name_to_idx_.find(clean_name(name));
    if (it != name_to_idx_.end()) return it->second;
    return -1;
}

std::vector<const Chromosome*> Genome::chromosomes_without_all() const {
    std::vector<const Chromosome*> result;
    for (auto& c : chromosomes_) {
        if (c.name != "ALL") result.push_back(&c);
    }
    return result;
}

int64_t Genome::genome_position(int chr_idx, int64_t pos) const {
    if (chr_idx < 0 || chr_idx >= static_cast<int>(chr_offsets_.size()))
        return -1;
    return (chr_offsets_[chr_idx] + pos) / 1000;
}
