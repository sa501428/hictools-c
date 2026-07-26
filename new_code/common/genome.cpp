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
