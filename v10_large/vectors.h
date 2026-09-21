#pragma once

#include "io.h"
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace hic10large {

constexpr uint32_t NO_ID = UINT32_MAX;

struct VectorInfo {
    uint8_t kind = 0;
    uint32_t norm = NO_ID, chr = NO_ID, ri = 0, resolution = 0;
    uint64_t words = 0, checksum = 0;
    std::string path;
    std::map<uint32_t, uint32_t> scales;
};

struct VectorManifest {
    uint64_t source_fingerprint = 0;
    std::vector<std::string> norms;
    std::vector<VectorInfo> vectors;
};

class VectorFileWriter {
  public:
    VectorFileWriter(const std::string &path, VectorInfo info);
    ~VectorFileWriter();
    void add(uint32_t word);
    void add(const std::vector<uint32_t> &words);
    VectorInfo finish();
  private:
    void flush();
    std::string final_, temporary_;
    VectorInfo info_;
    File output_;
    std::vector<uint32_t> buffer_;
    bool finished_ = false;
};

class VectorFileReader {
  public:
    explicit VectorFileReader(const VectorInfo &info);
    std::vector<uint32_t> read(uint64_t begin, uint32_t count);
  private:
    VectorInfo info_;
    File input_;
    uint64_t next_word_ = 0;
    uint64_t checksum_ = 14695981039346656037ULL;
};

void write_vector_manifest(const VectorManifest &manifest, const std::string &path);
VectorManifest read_vector_manifest(const std::string &path, bool inspect_files = true);

uint32_t float_bits(float value);
float bits_float(uint32_t value);

} // namespace hic10large
