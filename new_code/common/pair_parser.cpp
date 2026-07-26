#include "pair_parser.h"
#include "little_endian.h"
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <vector>

// ============================================================
//  Helpers
// ============================================================

static bool ends_with(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) return false;
    return s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// Open a file, supporting .gz via pipe to zcat/gzip
static FILE* open_maybe_gz(const std::string& path) {
    if (ends_with(path, ".gz")) {
        std::string cmd = "gzip -dc " + path;
        FILE* f = popen(cmd.c_str(), "r");
        if (!f) throw std::runtime_error("Cannot decompress: " + path);
        return f;
    }
    FILE* f = fopen(path.c_str(), "r");
    if (!f) throw std::runtime_error("Cannot open file: " + path);
    return f;
}

// Read one tab/space-delimited field from a line buffer
// Returns pointer past the consumed field, or nullptr on end of string.
static const char* next_field(const char* p, char* field, size_t max_len) {
    while (*p == ' ' || *p == '\t') p++;
    if (*p == '\0' || *p == '\n' || *p == '\r') return nullptr;
    size_t i = 0;
    while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') {
        if (i + 1 < max_len) field[i++] = *p;
        p++;
    }
    field[i] = '\0';
    return p;
}

// ============================================================
//  Format detection
// ============================================================

InputFormat detect_format(const std::string& path) {
    // Extension-based detection first
    std::string p = path;
    if (ends_with(p, ".bin")) return InputFormat::BIN;
    if (ends_with(p, ".bn")) return InputFormat::BN;
    if (ends_with(p, ".gz")) p = p.substr(0, p.size() - 3);

    // Text formats may be gzip-compressed.
    if (ends_with(p, ".pairs")) return InputFormat::PAIRS;
    if (ends_with(p, ".mnd")) return InputFormat::MND;
    // .txt is ambiguous — fall through to content detection

    // Peek at first non-comment line
    FILE* f = open_maybe_gz(path);
    char line[4096];
    InputFormat fmt = InputFormat::MND; // default
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') {
            // 4DN .pairs files start with "## pairs format"
            if (strncmp(line, "## pairs", 8) == 0) { fmt = InputFormat::PAIRS; break; }
            continue;
        }
        // Count whitespace-separated fields
        int fields = 0;
        const char* q = line;
        char tmp[256];
        while ((q = next_field(q, tmp, sizeof(tmp))) != nullptr) fields++;
        // Extra-short format has 4-5 fields: chr1 pos1 chr2 pos2 [score]
        // MND/medium/long have 8-16 fields
        if (fields == 4 || fields == 5) fmt = InputFormat::SHORT;
        break;
    }
    if (ends_with(path, ".gz")) pclose(f);
    else fclose(f);
    return fmt;
}

// ============================================================
//  MND (merged_nodups) parser
//  Format: strand chr pos frag strand chr pos frag [score mapq1 mapq2]
//          or (old): strand chr pos frag strand chr pos frag mapq1 mapq2
// ============================================================

class MNDIterator : public PairIterator {
public:
    MNDIterator(const std::string& path, const Genome& genome)
        : genome_(genome), is_pipe_(ends_with(path, ".gz")) {
        file_ = open_maybe_gz(path);
    }

    ~MNDIterator() override { close(); }

    void close() override {
        if (file_) {
            if (is_pipe_) pclose(file_);
            else fclose(file_);
            file_ = nullptr;
        }
    }

    bool next(AlignmentPair& out) override {
        char line[4096];
        while (fgets(line, sizeof(line), file_)) {
            if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
            if (parse_mnd_line(line, out)) return true;
        }
        return false;
    }

private:
    const Genome& genome_;
    FILE* file_ = nullptr;
    bool  is_pipe_;

    bool parse_mnd_line(const char* line, AlignmentPair& out) {
        // Handles four related formats (all tab/space separated):
        //
        // Short (8):           str1 chr1 pos1 frag1 str2 chr2 pos2 frag2
        // Short w/ score (9):  str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 score
        // Old MND (10):        str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 mapq2
        // New MND (11):        str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 mapq2 score
        // Long/Juicer (16):    str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 cigar1 seq1 mapq2 ...
        // Medium (11):         readname str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 mapq2
        //
        // Medium format is detected by col0 NOT being a single strand character.

        char fields[17][256];
        int nf = 0;
        const char* p = line;
        while (nf < 17 && (p = next_field(p, fields[nf], 256)) != nullptr) nf++;
        if (nf < 8) return false;

        // Detect medium format: col0 is a readname, not a strand char (single 0/1/+/-)
        int o = 0; // field offset (1 for medium format)
        {
            char c = fields[0][0];
            bool is_strand = (c == '0' || c == '1' || c == '+' || c == '-') && fields[0][1] == '\0';
            if (!is_strand) o = 1;
        }
        if (nf - o < 8) return false;

        // strand
        out.strand1 = (fields[o+0][0] == '1' || fields[o+0][0] == '-') ? 1 : 0;
        out.strand2 = (fields[o+4][0] == '1' || fields[o+4][0] == '-') ? 1 : 0;

        // chromosomes
        out.chr1 = genome_.index_of(fields[o+1]);
        out.chr2 = genome_.index_of(fields[o+5]);
        if (out.chr1 < 0 || out.chr2 < 0) return false;

        // Juicer MND coordinates are already zero-based.  In particular,
        // merged_nodups.txt is written with the same genomic coordinates that
        // the Java preprocessor bins directly.
        out.pos1  = std::max(0, (int)std::atol(fields[o+2]));
        out.frag1 = std::atoi(fields[o+3]);
        out.pos2  = std::max(0, (int)std::atol(fields[o+6]));
        out.frag2 = std::atoi(fields[o+7]);

        out.mapq1 = 1000;
        out.mapq2 = 1000;
        out.score = 1.0f;

        int eff = nf - o; // effective field count after readname offset
        if (eff == 9 && o == 0) {
            // Short with score: 9 fields, score is col8
            out.score = (float)std::atof(fields[8]);
        } else if (eff == 16 && o == 0) {
            // Long/Juicer:
            // str chr pos frag str chr pos frag mapq1 cigar1 seq1 mapq2 ...
            out.mapq1 = std::atoi(fields[o+8]);
            out.mapq2 = std::atoi(fields[o+11]);
        } else if (eff >= 10) {
            // Old/new short MND and medium formats have adjacent MAPQs.
            out.mapq1 = std::atoi(fields[o+8]);
            out.mapq2 = std::atoi(fields[o+9]);
            // Score at col10 only when exactly 11 fields and no readname offset
            // (nf==16 long format: col10 is sequence, not score)
            if (eff == 11 && o == 0) {
                out.score = (float)std::atof(fields[10]);
            }
        }

        return true;
    }
};

// ============================================================
//  Short format parser
//  Format: chr1 pos1 chr2 pos2
// ============================================================

class ShortIterator : public PairIterator {
public:
    ShortIterator(const std::string& path, const Genome& genome)
        : genome_(genome), is_pipe_(ends_with(path, ".gz")) {
        file_ = open_maybe_gz(path);
    }

    ~ShortIterator() override { close(); }

    void close() override {
        if (file_) {
            if (is_pipe_) pclose(file_);
            else fclose(file_);
            file_ = nullptr;
        }
    }

    bool next(AlignmentPair& out) override {
        char line[4096];
        while (fgets(line, sizeof(line), file_)) {
            if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
            if (parse_short_line(line, out)) return true;
        }
        return false;
    }

private:
    const Genome& genome_;
    FILE* file_ = nullptr;
    bool  is_pipe_;

    bool parse_short_line(const char* line, AlignmentPair& out) {
        // Extra-short format: chr1 pos1 chr2 pos2 [score]
        char f[5][256];
        int nf = 0;
        const char* p = line;
        while (nf < 5 && (p = next_field(p, f[nf], 256)) != nullptr) nf++;
        if (nf < 4) return false;

        out.chr1 = genome_.index_of(f[0]);
        out.chr2 = genome_.index_of(f[2]);
        if (out.chr1 < 0 || out.chr2 < 0) return false;

        // Match Java's SUPER_SHORT parser: coordinates are consumed as-is.
        out.pos1   = std::max(0, (int)std::atol(f[1]));
        out.pos2   = std::max(0, (int)std::atol(f[3]));
        out.strand1 = out.strand2 = 0;
        out.frag1 = out.frag2 = 0;
        out.mapq1 = out.mapq2 = 1000;
        out.score  = (nf >= 5) ? (float)std::atof(f[4]) : 1.0f;
        return true;
    }
};

// ============================================================
//  4DN .pairs format parser
//  Header: lines starting with #
//  Data:   readID chr1 pos1 chr2 pos2 strand1 strand2 [rest]
// ============================================================

class PairsIterator : public PairIterator {
public:
    PairsIterator(const std::string& path, const Genome& genome)
        : genome_(genome), is_pipe_(ends_with(path, ".gz")) {
        file_ = open_maybe_gz(path);
        // Parse header: scan #columns: to locate optional fragment and MAPQ columns.
        char line[4096];
        while (fgets(line, sizeof(line), file_)) {
            if (line[0] != '#') {
                first_line_ = line;
                has_first_ = true;
                break;
            }
            // Check for columns header: "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 ..."
            if (strncmp(line, "#columns:", 9) == 0) {
                const char* q = line + 9;
                char col[64];
                int ci = 0;
                frag1_col_ = -1;
                frag2_col_ = -1;
                mapq1_col_ = -1;
                mapq2_col_ = -1;
                while ((q = next_field(q, col, sizeof(col))) != nullptr) {
                    if (strcmp(col, "frag1") == 0) frag1_col_ = ci;
                    if (strcmp(col, "frag2") == 0) frag2_col_ = ci;
                    if (strcmp(col, "mapq1") == 0) mapq1_col_ = ci;
                    if (strcmp(col, "mapq2") == 0) mapq2_col_ = ci;
                    ci++;
                }
            }
        }
    }

    ~PairsIterator() override { close(); }

    void close() override {
        if (file_) {
            if (is_pipe_) pclose(file_);
            else fclose(file_);
            file_ = nullptr;
        }
    }

    bool next(AlignmentPair& out) override {
        char line[4096];
        if (has_first_) {
            has_first_ = false;
            return parse_pairs_line(first_line_.c_str(), out);
        }
        while (fgets(line, sizeof(line), file_)) {
            if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
            if (parse_pairs_line(line, out)) return true;
        }
        return false;
    }

private:
    const Genome& genome_;
    FILE* file_ = nullptr;
    bool  is_pipe_;
    bool  has_first_ = false;
    std::string first_line_;
    int   frag1_col_ = -1; // column index (0-based) of frag1 in data lines, -1 if absent
    int   frag2_col_ = -1;
    int   mapq1_col_ = -1;
    int   mapq2_col_ = -1;

    bool parse_pairs_line(const char* line, AlignmentPair& out) {
        // readID chr1 pos1 chr2 pos2 strand1 strand2 [extra columns...]
        // Maximum columns we ever need: frag1/frag2 can appear at any position > 6.
        // Read enough to cover whatever frag_col_ indices we found.
        int max_optional_col = std::max(
            std::max(frag1_col_, frag2_col_),
            std::max(mapq1_col_, mapq2_col_));
        int max_cols = std::max(8, max_optional_col + 1);
        // Cap for safety
        if (max_cols > 32) max_cols = 32;

        std::vector<std::string> f;
        f.reserve(max_cols);
        char tmp[256];
        const char* p = line;
        while ((int)f.size() < max_cols && (p = next_field(p, tmp, sizeof(tmp))) != nullptr)
            f.emplace_back(tmp);

        if ((int)f.size() < 5) return false;

        out.chr1 = genome_.index_of(f[1].c_str());
        out.chr2 = genome_.index_of(f[3].c_str());
        if (out.chr1 < 0 || out.chr2 < 0) return false;

        // Match Java's DCIC parser: coordinates are consumed as-is.
        out.pos1 = std::max(0, (int)std::atol(f[2].c_str()));
        out.pos2 = std::max(0, (int)std::atol(f[4].c_str()));

        out.strand1 = ((int)f.size() >= 6 && f[5][0] == '-') ? 1 : 0;
        out.strand2 = ((int)f.size() >= 7 && f[6][0] == '-') ? 1 : 0;

        // frag1/frag2: use #columns-specified positions if available, else 0/1
        // (Java: "set frag1=0 and frag2=1 so that no reads are discarded")
        bool has_fragments = frag1_col_ >= 0 && frag2_col_ >= 0
                          && frag1_col_ < (int)f.size() && frag2_col_ < (int)f.size();
        out.frag1 = has_fragments ? std::atoi(f[frag1_col_].c_str()) : 0;
        out.frag2 = has_fragments ? std::atoi(f[frag2_col_].c_str()) : 1;

        // Match Java: MAPQ filtering is enabled only when both columns exist.
        bool has_mapqs = mapq1_col_ >= 0 && mapq2_col_ >= 0
                      && mapq1_col_ < (int)f.size() && mapq2_col_ < (int)f.size();
        out.mapq1 = has_mapqs ? std::atoi(f[mapq1_col_].c_str()) : 1000;
        out.mapq2 = has_mapqs ? std::atoi(f[mapq2_col_].c_str()) : 1000;
        out.score = 1.0f;
        return true;
    }
};

// ============================================================
//  Juicer binary pair formats
// ============================================================

class BinaryPairIterator : public PairIterator {
public:
    BinaryPairIterator(const std::string& path, const Genome& genome, bool short_format)
        : genome_(genome), short_format_(short_format), path_(path) {
        file_ = fopen(path.c_str(), "rb");
        if (!file_) throw std::runtime_error("Cannot open binary pair file: " + path);
    }

    ~BinaryPairIterator() override { close(); }

    void close() override {
        if (file_) {
            fclose(file_);
            file_ = nullptr;
        }
    }

    bool next(AlignmentPair& out) override {
        const size_t record_size = short_format_ ? 20 : 26;
        uint8_t record[26];
        size_t n = fread(record, 1, record_size, file_);
        if (n == 0) {
            if (ferror(file_)) throw std::runtime_error("Error reading binary pair file: " + path_);
            return false;
        }
        if (n != record_size) {
            throw std::runtime_error("Truncated binary pair record in: " + path_);
        }

        auto get_i32 = [](const uint8_t* p) {
            uint32_t bits = static_cast<uint32_t>(p[0])
                          | (static_cast<uint32_t>(p[1]) << 8)
                          | (static_cast<uint32_t>(p[2]) << 16)
                          | (static_cast<uint32_t>(p[3]) << 24);
            return static_cast<int32_t>(bits);
        };

        if (short_format_) {
            out.chr1 = get_i32(record);
            out.pos1 = get_i32(record + 4);
            out.chr2 = get_i32(record + 8);
            out.pos2 = get_i32(record + 12);
            uint32_t score_bits = static_cast<uint32_t>(record[16])
                                | (static_cast<uint32_t>(record[17]) << 8)
                                | (static_cast<uint32_t>(record[18]) << 16)
                                | (static_cast<uint32_t>(record[19]) << 24);
            std::memcpy(&out.score, &score_bits, sizeof(out.score));
            out.strand1 = 0;
            out.strand2 = 1;
            out.frag1 = 0;
            out.frag2 = 1;
        } else {
            // Java stores true for forward strand; C++ stores 0 for forward.
            out.strand1 = record[0] ? 0 : 1;
            out.chr1 = get_i32(record + 1);
            out.pos1 = get_i32(record + 5);
            out.frag1 = get_i32(record + 9);
            out.strand2 = record[13] ? 0 : 1;
            out.chr2 = get_i32(record + 14);
            out.pos2 = get_i32(record + 18);
            out.frag2 = get_i32(record + 22);
            out.score = 1.0f;
        }
        out.mapq1 = out.mapq2 = 1000;

        // Binary chromosome IDs are 1-based indices into the supplied genome,
        // exactly as in Java's BinPairIterator.
        if (out.chr1 <= 0 || out.chr1 >= genome_.size()
                || out.chr2 <= 0 || out.chr2 >= genome_.size()) {
            throw std::runtime_error(
                "Binary pair chromosome index is incompatible with the supplied genome: " + path_);
        }
        return true;
    }

private:
    const Genome& genome_;
    bool short_format_;
    std::string path_;
    FILE* file_ = nullptr;
};

// ============================================================
//  Factory
// ============================================================

std::unique_ptr<PairIterator> open_pair_iterator(
    const std::string& path,
    const Genome& genome,
    InputFormat fmt)
{
    if (fmt == InputFormat::AUTO) fmt = detect_format(path);

    switch (fmt) {
        case InputFormat::MND:
            return std::make_unique<MNDIterator>(path, genome);
        case InputFormat::SHORT:
            return std::make_unique<ShortIterator>(path, genome);
        case InputFormat::PAIRS:
            return std::make_unique<PairsIterator>(path, genome);
        case InputFormat::BIN:
            return std::make_unique<BinaryPairIterator>(path, genome, false);
        case InputFormat::BN:
            return std::make_unique<BinaryPairIterator>(path, genome, true);
        default:
            throw std::runtime_error("Unknown input format");
    }
}
