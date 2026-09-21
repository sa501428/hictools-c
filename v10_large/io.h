#pragma once

#include <array>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace hic10large {

[[noreturn]] void fail(const std::string &message);
void require(bool condition, const std::string &message);

uint64_t parse_size(const std::string &text);
std::string decimal_u128(unsigned __int128 value);
unsigned __int128 parse_u128(const std::string &text);

uint64_t fnv1a(const void *data, size_t size,
               uint64_t seed = 14695981039346656037ULL);
std::string hex64(uint64_t value);
uint64_t unhex64(const std::string &value);

void put_u16(unsigned char *p, uint16_t value);
void put_u32(unsigned char *p, uint32_t value);
void put_u64(unsigned char *p, uint64_t value);
uint16_t get_u16(const unsigned char *p);
uint32_t get_u32(const unsigned char *p);
uint64_t get_u64(const unsigned char *p);

class File {
  public:
    File() = default;
    File(const std::string &path, const char *mode);
    ~File();
    File(const File &) = delete;
    File &operator=(const File &) = delete;
    File(File &&other) noexcept;
    File &operator=(File &&other) noexcept;

    FILE *get() const { return file_; }
    const std::string &path() const { return path_; }
    explicit operator bool() const { return file_ != nullptr; }
    void close();
    void flush();
    uint64_t tell() const;
    void seek(uint64_t offset);
    void write(const void *data, size_t size);
    bool read(void *data, size_t size, bool eof_allowed = false);

  private:
    FILE *file_ = nullptr;
    std::string path_;
};

void make_directory(const std::string &path);
bool path_exists(const std::string &path);
uint64_t file_size(const std::string &path);
std::string join_path(const std::string &directory, const std::string &name);
std::string base_name(const std::string &path);
void atomic_rename(const std::string &from, const std::string &to);

struct StagedRecord {
    uint32_t x = 0;
    uint32_t y = 0;
    uint64_t count = 0;
};
static_assert(sizeof(StagedRecord) == 16, "staged records must remain compact");

constexpr uint32_t STAGED_VERSION = 1;
constexpr size_t STAGED_HEADER_BYTES = 64;

struct StagedHeader {
    uint32_t source_resolution = 0;
    uint32_t chr1 = 0;
    uint32_t chr2 = 0;
    uint32_t part = 0;
    uint64_t records = 0;
    uint64_t checksum = 0;
};

std::array<unsigned char, STAGED_HEADER_BYTES> encode_staged_header(const StagedHeader &header);
StagedHeader decode_staged_header(const unsigned char *bytes, size_t size);

} // namespace hic10large
