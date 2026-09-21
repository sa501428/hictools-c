#include "vectors.h"

#include <array>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <unistd.h>

namespace hic10large {
namespace {
constexpr size_t VECTOR_HEADER_BYTES = 64;

std::array<unsigned char, VECTOR_HEADER_BYTES> header(const VectorInfo &v) {
    std::array<unsigned char, VECTOR_HEADER_BYTES> b{};
    std::memcpy(b.data(), "H10W", 4);
    put_u32(b.data() + 4, 1);
    put_u32(b.data() + 8, v.kind);
    put_u32(b.data() + 12, v.norm);
    put_u32(b.data() + 16, v.chr);
    put_u32(b.data() + 20, v.ri);
    put_u32(b.data() + 24, v.resolution);
    put_u32(b.data() + 28, 4);
    put_u64(b.data() + 32, v.words);
    put_u64(b.data() + 40, v.checksum);
    return b;
}
}

uint32_t float_bits(float value) {
    uint32_t word;
    std::memcpy(&word, &value, sizeof(word));
    return word;
}
float bits_float(uint32_t value) {
    float word;
    std::memcpy(&word, &value, sizeof(word));
    return word;
}

VectorFileWriter::VectorFileWriter(const std::string &path, VectorInfo info)
    : final_(path), temporary_(path + ".tmp-" + std::to_string(getpid())),
      info_(std::move(info)),
      output_(temporary_, "w+b") {
    info_.path = path;
    info_.words = 0;
    info_.checksum = 14695981039346656037ULL;
    auto b = header(info_);
    output_.write(b.data(), b.size());
    buffer_.reserve(65536);
}
VectorFileWriter::~VectorFileWriter() {
    if (!finished_) std::remove(temporary_.c_str());
}
void VectorFileWriter::add(uint32_t word) {
    buffer_.push_back(word);
    info_.checksum = fnv1a(&word, sizeof(word), info_.checksum);
    ++info_.words;
    if (buffer_.size() == buffer_.capacity()) flush();
}
void VectorFileWriter::add(const std::vector<uint32_t> &words) {
    for (uint32_t word : words) add(word);
}
void VectorFileWriter::flush() {
    if (buffer_.empty()) return;
    output_.write(buffer_.data(), buffer_.size() * sizeof(uint32_t));
    buffer_.clear();
}
VectorInfo VectorFileWriter::finish() {
    require(!finished_, "vector writer finished twice");
    flush();
    auto b = header(info_);
    output_.seek(0);
    output_.write(b.data(), b.size());
    output_.flush();
    output_.close();
    atomic_rename(temporary_, final_);
    finished_ = true;
    return info_;
}

VectorFileReader::VectorFileReader(const VectorInfo &info) : info_(info), input_(info.path, "rb") {
    std::array<unsigned char, VECTOR_HEADER_BYTES> b{};
    input_.read(b.data(), b.size());
    require(std::memcmp(b.data(), "H10W", 4) == 0 && get_u32(b.data() + 4) == 1 &&
                get_u32(b.data() + 8) == info.kind && get_u32(b.data() + 12) == info.norm &&
                get_u32(b.data() + 16) == info.chr && get_u32(b.data() + 20) == info.ri &&
                get_u32(b.data() + 24) == info.resolution && get_u32(b.data() + 28) == 4 &&
                get_u64(b.data() + 32) == info.words && get_u64(b.data() + 40) == info.checksum,
            "invalid vector sidecar " + info.path);
}
std::vector<uint32_t> VectorFileReader::read(uint64_t begin, uint32_t count) {
    require(begin <= info_.words && count <= info_.words - begin, "vector sidecar read outside range");
    require(begin == next_word_, "vector sidecar reads must be sequential");
    input_.seek(VECTOR_HEADER_BYTES + begin * 4);
    std::vector<uint32_t> result(count);
    if (count) input_.read(result.data(), uint64_t(count) * 4);
    checksum_ = fnv1a(result.data(), uint64_t(count) * 4, checksum_);
    next_word_ += count;
    if (next_word_ == info_.words)
        require(checksum_ == info_.checksum, "vector sidecar checksum mismatch: " + info_.path);
    return result;
}

void write_vector_manifest(const VectorManifest &m, const std::string &path) {
    std::string temporary = path + ".tmp-" + std::to_string(getpid());
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    require(bool(out), "cannot create vector manifest");
    out << "HIC_V10_LARGE_VECTORS 1\n";
    out << "source " << hex64(m.source_fingerprint) << "\n";
    out << "norms " << m.norms.size();
    for (const auto &norm : m.norms) out << ' ' << std::quoted(norm);
    out << "\n";
    out << "vectors " << m.vectors.size() << "\n";
    for (const auto &v : m.vectors) {
        out << "vector " << unsigned(v.kind) << ' ' << v.norm << ' ' << v.chr << ' '
            << v.ri << ' ' << v.resolution << ' ' << v.words << ' ' << hex64(v.checksum)
            << ' ' << std::quoted(v.path) << ' ' << v.scales.size();
        for (const auto &scale : v.scales)
            out << ' ' << scale.first << ' ' << hex64(scale.second);
        out << "\n";
    }
    out << "end\n";
    out.flush();
    require(bool(out), "cannot write vector manifest");
    out.close();
    atomic_rename(temporary, path);
}

VectorManifest read_vector_manifest(const std::string &path, bool inspect_files) {
    std::ifstream in(path, std::ios::binary);
    require(bool(in), "cannot open vector manifest " + path);
    std::string magic;
    unsigned version = 0;
    in >> magic >> version;
    require(magic == "HIC_V10_LARGE_VECTORS" && version == 1, "invalid vector manifest");
    VectorManifest result;
    std::string tag;
    size_t expected = 0;
    while (in >> tag) {
        if (tag == "source") {
            std::string hash; in >> hash; result.source_fingerprint = unhex64(hash);
        } else if (tag == "norms") {
            size_t count = 0; in >> count; result.norms.resize(count);
            for (auto &norm : result.norms) in >> std::quoted(norm);
        } else if (tag == "vectors") {
            in >> expected; result.vectors.reserve(expected);
        } else if (tag == "vector") {
            VectorInfo v;
            unsigned kind = 0;
            std::string checksum;
            size_t scales = 0;
            in >> kind >> v.norm >> v.chr >> v.ri >> v.resolution >> v.words >> checksum
               >> std::quoted(v.path) >> scales;
            require(kind <= 2, "invalid vector kind");
            v.kind = static_cast<uint8_t>(kind);
            v.checksum = unhex64(checksum);
            for (size_t i = 0; i < scales; ++i) {
                uint32_t chr = 0;
                std::string word;
                in >> chr >> word;
                require(v.scales.emplace(chr, static_cast<uint32_t>(unhex64(word))).second,
                        "duplicate vector scale");
            }
            if (inspect_files) {
                VectorFileReader reader(v);
                require(file_size(v.path) == VECTOR_HEADER_BYTES + v.words * 4,
                        "vector sidecar length mismatch");
            }
            result.vectors.push_back(std::move(v));
        } else if (tag == "end") {
            break;
        } else fail("unknown vector manifest field " + tag);
        require(bool(in), "truncated vector manifest");
    }
    require(result.vectors.size() == expected, "incomplete vector manifest");
    return result;
}

} // namespace hic10large
