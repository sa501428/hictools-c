#pragma once

#include "io.h"
#include <map>
#include <string>
#include <vector>

namespace hic10large {

struct Chromosome {
    std::string name;
    uint64_t length = 0;
};

struct ShardInfo {
    uint32_t chr1 = 0, chr2 = 0, part = 0;
    uint64_t first_record = 0, records = 0, bytes = 0, checksum = 0;
    unsigned __int128 weight = 0;
    std::string path;
};

struct PairInfo {
    uint32_t chr1 = 0, chr2 = 0;
    uint64_t first_record = 0, records = 0;
    unsigned __int128 weight = 0;
    std::vector<size_t> shards;
};

struct StageManifest {
    std::string source;
    uint64_t source_bytes = 0;
    uint64_t source_fingerprint = 0;
    uint32_t source_resolution = 0;
    uint64_t total_records = 0;
    unsigned __int128 total_weight = 0;
    std::vector<Chromosome> chromosomes;
    std::vector<ShardInfo> shards;
    std::vector<PairInfo> pairs;
};

void write_stage_manifest(const StageManifest &manifest, const std::string &path);
StageManifest read_stage_manifest(const std::string &path);

} // namespace hic10large
