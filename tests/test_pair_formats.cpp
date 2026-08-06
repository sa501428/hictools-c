#include "common/genome.h"
#include "common/pair_parser.h"
#include "common/little_endian.h"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static void check_pairs_mapq(const fs::path& path, const Genome& genome) {
    {
        std::ofstream out(path);
        out << "## pairs format v1.0\n"
            << "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 pair_type mapq2 mapq1 frag2 frag1\n"
            << "r1 chr1 100 chr2 200 + - UU 17 42 8 7\n";
    }

    require(detect_format(path.string()) == InputFormat::PAIRS, "failed to detect .pairs");
    auto iter = open_pair_iterator(path.string(), genome);
    AlignmentPair pair;
    require(iter->next(pair), "missing .pairs record");
    require(pair.chr1 == 1 && pair.chr2 == 2, "wrong .pairs chromosomes");
    require(pair.pos1 == 100 && pair.pos2 == 200, "wrong .pairs positions");
    require(pair.strand1 == 0 && pair.strand2 == 1, "wrong .pairs strands");
    require(pair.mapq1 == 42 && pair.mapq2 == 17, "wrong .pairs MAPQs");
    require(pair.frag1 == 7 && pair.frag2 == 8, "wrong .pairs fragments");
    require(!iter->next(pair), "extra .pairs record");
}

static void check_bin(const fs::path& path, const Genome& genome) {
    FILE* out = fopen(path.c_str(), "wb");
    require(out != nullptr, "failed to create .bin");
    // Java boolean true means forward.
    fwrite_byte(out, 1);
    fwrite_int32(out, 1);
    fwrite_int32(out, 123);
    fwrite_int32(out, 4);
    fwrite_byte(out, 0);
    fwrite_int32(out, 2);
    fwrite_int32(out, 456);
    fwrite_int32(out, 5);
    fclose(out);

    require(detect_format(path.string()) == InputFormat::BIN, "failed to detect .bin");
    auto iter = open_pair_iterator(path.string(), genome);
    AlignmentPair pair;
    require(iter->next(pair), "missing .bin record");
    require(pair.chr1 == 1 && pair.chr2 == 2, "wrong .bin chromosomes");
    require(pair.pos1 == 123 && pair.pos2 == 456, "wrong .bin positions");
    require(pair.frag1 == 4 && pair.frag2 == 5, "wrong .bin fragments");
    require(pair.strand1 == 0 && pair.strand2 == 1, "wrong .bin strands");
    require(pair.mapq1 == 1000 && pair.mapq2 == 1000, "wrong .bin MAPQ sentinel");
    require(!iter->next(pair), "extra .bin record");
}

static void check_bn(const fs::path& path, const Genome& genome) {
    FILE* out = fopen(path.c_str(), "wb");
    require(out != nullptr, "failed to create .bn");
    fwrite_int32(out, 2);
    fwrite_int32(out, 789);
    fwrite_int32(out, 1);
    fwrite_int32(out, 987);
    fwrite_float(out, 2.5f);
    fclose(out);

    require(detect_format(path.string()) == InputFormat::BN, "failed to detect .bn");
    auto iter = open_pair_iterator(path.string(), genome);
    AlignmentPair pair;
    require(iter->next(pair), "missing .bn record");
    require(pair.chr1 == 2 && pair.chr2 == 1, "wrong .bn chromosomes");
    require(pair.pos1 == 789 && pair.pos2 == 987, "wrong .bn positions");
    require(pair.frag1 == 0 && pair.frag2 == 1, "wrong .bn fragment defaults");
    require(std::fabs(pair.score - 2.5f) < 1e-6f, "wrong .bn score");
    require(pair.mapq1 == 1000 && pair.mapq2 == 1000, "wrong .bn MAPQ sentinel");
    require(!iter->next(pair), "extra .bn record");
}

int main() {
    const Genome genome = Genome::from_list({{"chr1", 1000}, {"chr2", 2000}});
    const fs::path dir = fs::temp_directory_path() / "hictools_pair_format_test";
    fs::create_directories(dir);

    check_pairs_mapq(dir / "input.pairs", genome);
    check_bin(dir / "input.bin", genome);
    check_bn(dir / "input.bn", genome);

    require(Genome::is_builtin("hg19"), "hg19 not built in");
    require(Genome::is_builtin("HG38"), "hg38 lookup not case insensitive");
    require(Genome::is_builtin("mm9"), "mm9 not built in");
    require(Genome::is_builtin("mm10"), "mm10 not built in");
    require(!Genome::is_builtin("not-a-genome"), "accepted unknown built-in genome");

    const Genome hg38 = Genome::from_spec("hg38");
    require(hg38.size() == 26, "wrong hg38 chromosome count");
    require(hg38.index_of("chr1") == 1, "wrong hg38 chr1 index");
    require(hg38.index_of("1") == 1, "hg38 chr-prefix alias failed");
    require(hg38.at(1).length == 248956422, "wrong hg38 chr1 length");
    require(hg38.at(25).name == "chrM", "wrong hg38 mitochondrial chromosome");

    const Genome hg19 = Genome::from_spec("hg19");
    require(hg19.size() == 26, "wrong hg19 chromosome count");
    require(hg19.index_of("chr1") == 1, "wrong hg19 chr1 index");
    require(hg19.at(25).name == "MT", "wrong hg19 mitochondrial chromosome");

    const Genome mm9 = Genome::from_spec("mm9");
    require(mm9.size() == 23, "wrong mm9 chromosome count");
    require(mm9.at(1).length == 197195432, "wrong mm9 chr1 length");
    require(mm9.at(22).name == "chrM", "wrong mm9 mitochondrial chromosome");

    const Genome mm10 = Genome::from_spec("mm10");
    require(mm10.size() == 23, "wrong mm10 chromosome count");
    require(mm10.at(1).length == 195471971, "wrong mm10 chr1 length");
    require(mm10.at(20).name == "chrM", "wrong mm10 chromosome order");

    fs::remove_all(dir);
    return 0;
}
