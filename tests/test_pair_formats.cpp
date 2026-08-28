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
#include <zlib.h>

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

static void check_hbs(const fs::path& path, const Genome& genome) {
    // Independent byte fixture: IDs deliberately differ from the genome order.
    std::vector<uint8_t> data{'H','I','C','B','S',0,13,10};
    auto put = [&](uint64_t n, unsigned bytes) {
        for (unsigned i = 0; i < bytes; ++i) data.push_back(static_cast<uint8_t>(n >> (8*i)));
    };
    put(1, 2); put(0, 2); put(10, 4); put(2, 4);
    for (auto entry : std::vector<std::pair<std::string, uint64_t>>{{"chr2",2000}, {"1",1000}}) {
        put(entry.first.size(), 2);
        data.insert(data.end(), entry.first.begin(), entry.first.end());
        put(entry.second, 8);
    }
    const size_t headerSize = data.size();
    const std::vector<uint64_t> counts{0, 1, 65534, 65535, 65536, (uint64_t{1}<<24)+1,
                                       (uint64_t{1}<<53)+1, UINT64_MAX};
    for (auto count : counts) {
        put(0, 2); put(12, 4); put(1, 2); put(34, 4);
        put(count < 65535 ? count : 65535, 2);
        if (count >= 65535) put(count, 8);
    }
    auto save = [&](const std::vector<uint8_t>& bytes) {
        gzFile f = gzopen(path.c_str(), "wb");
        require(f != nullptr, "cannot create HBS fixture");
        require(gzwrite(f, bytes.data(), bytes.size()) == static_cast<int>(bytes.size()), "cannot write HBS fixture");
        require(gzclose(f) == Z_OK, "cannot close HBS fixture");
    };
    save(data);
    require(detect_format(path.string()) == InputFormat::HBS, "HBS detection");
    auto iterator = open_pair_iterator(path.string(), genome);
    require(iterator->source_resolution() == 10, "HBS resolution");
    AlignmentPair pair;
    for (auto count : counts) {
        require(iterator->next(pair), "missing HBS record");
        require(pair.chr1 == 2 && pair.chr2 == 1, "HBS name mapping");
        require(pair.pos1 == 120 && pair.pos2 == 340, "HBS bin positions");
        require(pair.has_exact_count && pair.exact_count == count, "HBS lost integer precision");
        require(pair.mapq1 == 1000 && pair.mapq2 == 1000 && pair.frag2 == 1, "HBS defaults");
    }
    require(!iterator->next(pair), "extra HBS record");
    iterator->close();
    auto reject = [&]() {
        bool failed = false;
        try {
            auto bad = open_pair_iterator(path.string(), genome, InputFormat::HBS);
            while (bad->next(pair)) {}
            bad->close();
        } catch (const std::exception&) { failed = true; }
        require(failed, "accepted malformed HBS");
    };
    for (auto mutation : std::vector<std::pair<size_t, uint8_t>>{
             {0, 'X'}, {8, 2}, {10, 1}, {12, 0}, {18, 2}, {20, 0}, {22, 'X'},
             {26, 0}, {headerSize, 2}, {headerSize+2, 255}}) {
        auto bad = data;
        bad[mutation.first] = mutation.second;
        save(bad); reject();
    }
    for (size_t length : {size_t(19), headerSize-1, headerSize+1, data.size()-1}) {
        save(std::vector<uint8_t>(data.begin(), data.begin()+length)); reject();
    }
    auto badEscape = data;
    badEscape[headerSize + 3*14 + 14] = 1;
    badEscape[headerSize + 3*14 + 15] = 0;
    save(badEscape); reject();
    save(data);
    fs::resize_file(path, fs::file_size(path)-4); reject(); // truncated gzip trailer
    save(data);
    {
        std::fstream f(path, std::ios::in | std::ios::out | std::ios::binary);
        f.seekg(-8, std::ios::end);
        char crc = 0; f.get(crc); crc ^= 1;
        f.seekp(-8, std::ios::end); f.put(crc);
    }
    reject();
    { std::ofstream f(path, std::ios::binary); f.write(reinterpret_cast<const char*>(data.data()), data.size()); }
    reject(); // .hbs.gz requires gzip, not zlib's transparent uncompressed mode.
}

int main() {
    const Genome genome = Genome::from_list({{"chr1", 1000}, {"chr2", 2000}});
    const fs::path dir = fs::temp_directory_path() / "hictools_pair_format_test";
    fs::create_directories(dir);

    check_pairs_mapq(dir / "input.pairs", genome);
    check_bin(dir / "input.bin", genome);
    check_bn(dir / "input.bn", genome);
    check_hbs(dir / "input.hbs.gz", genome);

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
