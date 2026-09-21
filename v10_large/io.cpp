#include "io.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

namespace hic10large {

[[noreturn]] void fail(const std::string &message) {
    throw std::runtime_error("hic_v10_large: " + message);
}

void require(bool condition, const std::string &message) {
    if (!condition)
        fail(message);
}

uint64_t parse_size(const std::string &text) {
    require(!text.empty(), "empty size");
    size_t used = 0;
    long double value = std::stold(text, &used);
    require(std::isfinite(value) && value >= 0, "invalid size: " + text);
    long double multiplier = 1;
    std::string suffix = text.substr(used);
    if (suffix == "K" || suffix == "KB") multiplier = 1000.0L;
    else if (suffix == "M" || suffix == "MB") multiplier = 1000.0L * 1000;
    else if (suffix == "G" || suffix == "GB") multiplier = 1000.0L * 1000 * 1000;
    else if (suffix == "T" || suffix == "TB") multiplier = 1000.0L * 1000 * 1000 * 1000;
    else if (suffix == "KiB") multiplier = 1024.0L;
    else if (suffix == "MiB") multiplier = 1024.0L * 1024;
    else if (suffix == "GiB") multiplier = 1024.0L * 1024 * 1024;
    else if (suffix == "TiB") multiplier = 1024.0L * 1024 * 1024 * 1024;
    else require(suffix.empty(), "unknown size suffix: " + suffix);
    value *= multiplier;
    require(value <= std::numeric_limits<uint64_t>::max(), "size is too large: " + text);
    return static_cast<uint64_t>(value);
}

std::string decimal_u128(unsigned __int128 value) {
    if (!value)
        return "0";
    std::string result;
    while (value) {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    std::reverse(result.begin(), result.end());
    return result;
}

unsigned __int128 parse_u128(const std::string &text) {
    require(!text.empty(), "empty uint128");
    unsigned __int128 result = 0;
    const unsigned __int128 maximum = ~static_cast<unsigned __int128>(0);
    for (char c : text) {
        require(c >= '0' && c <= '9', "invalid uint128: " + text);
        unsigned digit = static_cast<unsigned>(c - '0');
        require(result <= (maximum - digit) / 10, "uint128 overflow: " + text);
        result = result * 10 + digit;
    }
    return result;
}

uint64_t fnv1a(const void *data, size_t size, uint64_t seed) {
    auto p = static_cast<const unsigned char *>(data);
    uint64_t value = seed;
    for (size_t i = 0; i < size; ++i) {
        value ^= p[i];
        value *= 1099511628211ULL;
    }
    return value;
}

std::string hex64(uint64_t value) {
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16) << value;
    return out.str();
}

uint64_t unhex64(const std::string &value) {
    require(value.size() == 16, "invalid 64-bit hexadecimal value");
    uint64_t result = 0;
    for (char c : value) {
        unsigned digit = c >= '0' && c <= '9' ? c - '0' :
                         c >= 'a' && c <= 'f' ? c - 'a' + 10 :
                         c >= 'A' && c <= 'F' ? c - 'A' + 10 : 16;
        require(digit < 16, "invalid 64-bit hexadecimal value");
        result = (result << 4) | digit;
    }
    return result;
}

void put_u16(unsigned char *p, uint16_t v) {
    p[0] = static_cast<unsigned char>(v);
    p[1] = static_cast<unsigned char>(v >> 8);
}
void put_u32(unsigned char *p, uint32_t v) {
    for (unsigned i = 0; i < 4; ++i) p[i] = static_cast<unsigned char>(v >> (8 * i));
}
void put_u64(unsigned char *p, uint64_t v) {
    for (unsigned i = 0; i < 8; ++i) p[i] = static_cast<unsigned char>(v >> (8 * i));
}
uint16_t get_u16(const unsigned char *p) {
    return static_cast<uint16_t>(p[0] | (uint16_t(p[1]) << 8));
}
uint32_t get_u32(const unsigned char *p) {
    uint32_t v = 0;
    for (unsigned i = 0; i < 4; ++i) v |= uint32_t(p[i]) << (8 * i);
    return v;
}
uint64_t get_u64(const unsigned char *p) {
    uint64_t v = 0;
    for (unsigned i = 0; i < 8; ++i) v |= uint64_t(p[i]) << (8 * i);
    return v;
}

File::File(const std::string &path, const char *mode) : path_(path) {
    file_ = std::fopen(path.c_str(), mode);
    require(file_ != nullptr, "cannot open " + path + ": " + std::strerror(errno));
    std::setvbuf(file_, nullptr, _IOFBF, 8 * 1024 * 1024);
}
File::~File() {
    if (file_)
        std::fclose(file_);
}
File::File(File &&other) noexcept : file_(other.file_), path_(std::move(other.path_)) {
    other.file_ = nullptr;
}
File &File::operator=(File &&other) noexcept {
    if (this != &other) {
        if (file_) std::fclose(file_);
        file_ = other.file_;
        path_ = std::move(other.path_);
        other.file_ = nullptr;
    }
    return *this;
}
void File::close() {
    if (!file_) return;
    FILE *f = file_;
    file_ = nullptr;
    require(std::fclose(f) == 0, "cannot close " + path_ + ": " + std::strerror(errno));
}
void File::flush() {
    require(file_ && std::fflush(file_) == 0, "cannot flush " + path_);
}
uint64_t File::tell() const {
    require(file_, "tell on closed file");
    off_t value = ftello(file_);
    require(value >= 0, "cannot determine position in " + path_);
    return static_cast<uint64_t>(value);
}
void File::seek(uint64_t offset) {
    require(file_ && offset <= static_cast<uint64_t>(std::numeric_limits<off_t>::max()) &&
                fseeko(file_, static_cast<off_t>(offset), SEEK_SET) == 0,
            "cannot seek " + path_);
}
void File::write(const void *data, size_t size) {
    require(file_ && std::fwrite(data, 1, size, file_) == size, "cannot write " + path_);
}
bool File::read(void *data, size_t size, bool eof_allowed) {
    require(file_, "read on closed file");
    size_t got = std::fread(data, 1, size, file_);
    if (!got && eof_allowed && std::feof(file_)) return false;
    require(got == size, "truncated or unreadable file " + path_);
    return true;
}

void make_directory(const std::string &path) {
    if (::mkdir(path.c_str(), 0775) == 0) return;
    require(errno == EEXIST, "cannot create directory " + path + ": " + std::strerror(errno));
    struct stat s{};
    require(stat(path.c_str(), &s) == 0 && S_ISDIR(s.st_mode), path + " is not a directory");
}
bool path_exists(const std::string &path) {
    struct stat s{};
    return stat(path.c_str(), &s) == 0;
}
uint64_t file_size(const std::string &path) {
    struct stat s{};
    require(stat(path.c_str(), &s) == 0 && s.st_size >= 0, "cannot stat " + path);
    return static_cast<uint64_t>(s.st_size);
}
std::string join_path(const std::string &directory, const std::string &name) {
    return directory.empty() || directory.back() == '/' ? directory + name : directory + "/" + name;
}
std::string base_name(const std::string &path) {
    size_t slash = path.find_last_of('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}
void atomic_rename(const std::string &from, const std::string &to) {
    require(std::rename(from.c_str(), to.c_str()) == 0,
            "cannot publish " + to + ": " + std::strerror(errno));
}

std::array<unsigned char, STAGED_HEADER_BYTES> encode_staged_header(const StagedHeader &h) {
    std::array<unsigned char, STAGED_HEADER_BYTES> b{};
    std::memcpy(b.data(), "H10S", 4);
    put_u32(b.data() + 4, STAGED_VERSION);
    put_u32(b.data() + 8, h.source_resolution);
    put_u32(b.data() + 12, h.chr1);
    put_u32(b.data() + 16, h.chr2);
    put_u32(b.data() + 20, h.part);
    put_u64(b.data() + 24, h.records);
    put_u64(b.data() + 32, h.checksum);
    put_u32(b.data() + 40, sizeof(StagedRecord));
    return b;
}

StagedHeader decode_staged_header(const unsigned char *b, size_t size) {
    require(size >= STAGED_HEADER_BYTES && std::memcmp(b, "H10S", 4) == 0,
            "invalid staged shard magic");
    require(get_u32(b + 4) == STAGED_VERSION, "unsupported staged shard version");
    require(get_u32(b + 40) == sizeof(StagedRecord), "incompatible staged record width");
    StagedHeader h;
    h.source_resolution = get_u32(b + 8);
    h.chr1 = get_u32(b + 12);
    h.chr2 = get_u32(b + 16);
    h.part = get_u32(b + 20);
    h.records = get_u64(b + 24);
    h.checksum = get_u64(b + 32);
    return h;
}

} // namespace hic10large
