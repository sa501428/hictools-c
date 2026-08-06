#pragma once
// Input contact-pair parsers:
//   1. .mnd / merged_nodups (Juicer standard)
//   2. Short format: chr1 pos1 chr2 pos2
//   3. 4DN .pairs format (with ## header section)
//   4. Juicer .bin and .bn little-endian binary formats
//
// All parsers emit AlignmentPair records.  The caller decides what to do with
// contacts that don't pass MAPQ or chr-filter criteria.

#include "genome.h"
#include <string>
#include <memory>
#include <functional>

// Raw alignment pair as read from input
struct AlignmentPair {
    int32_t chr1   = -1;
    int32_t chr2   = -1;
    int32_t pos1   = 0;
    int32_t pos2   = 0;
    int8_t  strand1 = 0;  // 0 = forward, 1 = reverse
    int8_t  strand2 = 0;
    int32_t frag1  = 0;
    int32_t frag2  = 0;
    // Java uses 1000 as the sentinel for formats that do not carry MAPQ.
    int32_t mapq1  = 1000;
    int32_t mapq2  = 1000;
    float   score  = 1.0f;

    bool valid() const { return chr1 >= 0 && chr2 >= 0; }
    // Is this a contig-to-contig pair (filtered out by the java code)?
    bool is_contig_pair() const { return false; } // we skip unmapped chrs at parse time
};

// Callback type: called for each successfully parsed pair
using PairCallback = std::function<void(const AlignmentPair&)>;

// Detect format from file extension / first line
enum class InputFormat {
    MND,    // merged_nodups (Juicer): strand chr pos frag strand chr pos frag [mapq mapq score]
    SHORT,  // chr1 pos1 chr2 pos2
    PAIRS,  // 4DN .pairs: readID chr1 pos1 chr2 pos2 strand1 strand2 [...]
    BIN,    // Juicer .bin: strands, chromosome indices, positions, fragments
    BN,     // Juicer .bn: chromosome indices, positions, score
    AUTO,   // detect from content
};

// Base class for pair iterators
class PairIterator {
public:
    virtual ~PairIterator() = default;
    virtual bool next(AlignmentPair& out) = 0;
    virtual void close() = 0;
};

// Factory: open the appropriate parser based on format (AUTO = detect from header/extension).
// genome is used to map chromosome names to indices.
// mapq_threshold: pairs with min(mapq1, mapq2) < threshold are skipped by caller (or parser).
std::unique_ptr<PairIterator> open_pair_iterator(
    const std::string& path,
    const Genome& genome,
    InputFormat fmt = InputFormat::AUTO
);

// Detect format from first few bytes / file extension
InputFormat detect_format(const std::string& path);
