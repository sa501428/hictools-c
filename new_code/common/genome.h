#pragma once
// Genome / chromosome utilities.
// Loads chromosome sizes from a chrom.sizes file or a genome ID, and provides
// helpers used throughout the pipeline.

#include "hic_file_def.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>

class Genome {
public:
    // Load from a two-column chrom.sizes file (name\tlength per line).
    // Chromosomes are stored in file order; index 0 is the first chromosome.
    // A synthetic "All" (whole-genome) chromosome is prepended at index 0,
    // with length = sum of all chromosome lengths / 1000 (Juicer stores All in kb).
    static Genome from_chrom_sizes(const std::string& path);

    // Convenience: build from an explicit list (name, length) — used in tests.
    static Genome from_list(const std::vector<std::pair<std::string,int64_t>>& chroms);

    // Number of chromosomes INCLUDING the synthetic "All" entry at index 0.
    int size() const { return static_cast<int>(chromosomes_.size()); }

    // Chromosome by 0-based index.
    const Chromosome& at(int idx) const { return chromosomes_.at(idx); }

    // Lookup by name (case-sensitive, also tries adding/removing "chr" prefix).
    // Returns -1 if not found.
    int index_of(const std::string& name) const;

    // Clean a chromosome name: lowercase, remove leading "chr" prefix for matching.
    static std::string clean_name(const std::string& name);

    // All chromosomes (including "All" at [0]).
    const std::vector<Chromosome>& chromosomes() const { return chromosomes_; }

    // Chromosomes without the synthetic "All" entry.
    std::vector<const Chromosome*> chromosomes_without_all() const;

    // Whole-genome position: offset of the start of chr (index) in the virtual
    // concatenated genome used for the GW matrix, expressed in kb.
    int64_t genome_position(int chr_idx, int64_t pos) const;

private:
    std::vector<Chromosome>                  chromosomes_;   // [0] = "ALL"
    std::unordered_map<std::string, int>     name_to_idx_;   // clean name → index
    std::vector<int64_t>                     chr_offsets_;   // for genome_position()

    void build_index();
};
