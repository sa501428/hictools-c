#pragma once
// Portable little-endian binary I/O helpers.
// HiC format uses little-endian for all multi-byte values (same as Java DataOutputStream on LE platforms).

#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>

// ---- Write helpers (to a byte buffer) ----

inline void write_byte(std::vector<uint8_t>& buf, uint8_t v) {
    buf.push_back(v);
}

inline void write_int16(std::vector<uint8_t>& buf, int16_t v) {
    buf.push_back(static_cast<uint8_t>(v & 0xFF));
    buf.push_back(static_cast<uint8_t>((v >> 8) & 0xFF));
}

inline void write_int32(std::vector<uint8_t>& buf, int32_t v) {
    buf.push_back(static_cast<uint8_t>(v & 0xFF));
    buf.push_back(static_cast<uint8_t>((v >> 8) & 0xFF));
    buf.push_back(static_cast<uint8_t>((v >> 16) & 0xFF));
    buf.push_back(static_cast<uint8_t>((v >> 24) & 0xFF));
}

inline void write_int64(std::vector<uint8_t>& buf, int64_t v) {
    for (int i = 0; i < 8; i++) {
        buf.push_back(static_cast<uint8_t>((v >> (8 * i)) & 0xFF));
    }
}

inline void write_float(std::vector<uint8_t>& buf, float v) {
    uint32_t bits;
    std::memcpy(&bits, &v, 4);
    buf.push_back(static_cast<uint8_t>(bits & 0xFF));
    buf.push_back(static_cast<uint8_t>((bits >> 8) & 0xFF));
    buf.push_back(static_cast<uint8_t>((bits >> 16) & 0xFF));
    buf.push_back(static_cast<uint8_t>((bits >> 24) & 0xFF));
}

// Null-terminated string
inline void write_string(std::vector<uint8_t>& buf, const std::string& s) {
    for (char c : s) buf.push_back(static_cast<uint8_t>(c));
    buf.push_back(0);
}

// ---- Write helpers (direct to FILE* for streaming output) ----

inline void fwrite_byte(FILE* f, uint8_t v) {
    fwrite(&v, 1, 1, f);
}

inline void fwrite_int16(FILE* f, int16_t v) {
    uint8_t b[2] = { static_cast<uint8_t>(v & 0xFF), static_cast<uint8_t>((v >> 8) & 0xFF) };
    fwrite(b, 1, 2, f);
}

inline void fwrite_int32(FILE* f, int32_t v) {
    uint8_t b[4];
    for (int i = 0; i < 4; i++) b[i] = static_cast<uint8_t>((v >> (8*i)) & 0xFF);
    fwrite(b, 1, 4, f);
}

inline void fwrite_int64(FILE* f, int64_t v) {
    uint8_t b[8];
    for (int i = 0; i < 8; i++) b[i] = static_cast<uint8_t>((v >> (8*i)) & 0xFF);
    fwrite(b, 1, 8, f);
}

inline void fwrite_float(FILE* f, float v) {
    uint32_t bits;
    std::memcpy(&bits, &v, 4);
    uint8_t b[4];
    for (int i = 0; i < 4; i++) b[i] = static_cast<uint8_t>((bits >> (8*i)) & 0xFF);
    fwrite(b, 1, 4, f);
}

inline void fwrite_string(FILE* f, const std::string& s) {
    fwrite(s.c_str(), 1, s.size(), f);
    uint8_t nul = 0;
    fwrite(&nul, 1, 1, f);
}

// ---- Read helpers (direct from FILE*) ----

inline uint8_t fread_byte(FILE* f) {
    uint8_t v;
    if (fread(&v, 1, 1, f) != 1) throw std::runtime_error("Unexpected EOF reading byte");
    return v;
}

inline int16_t fread_int16(FILE* f) {
    uint8_t b[2];
    if (fread(b, 1, 2, f) != 2) throw std::runtime_error("Unexpected EOF reading int16");
    return static_cast<int16_t>(b[0] | (b[1] << 8));
}

inline int32_t fread_int32(FILE* f) {
    uint8_t b[4];
    if (fread(b, 1, 4, f) != 4) throw std::runtime_error("Unexpected EOF reading int32");
    return static_cast<int32_t>(
        (uint32_t)b[0] | ((uint32_t)b[1] << 8) | ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24));
}

inline int64_t fread_int64(FILE* f) {
    uint8_t b[8];
    if (fread(b, 1, 8, f) != 8) throw std::runtime_error("Unexpected EOF reading int64");
    int64_t v = 0;
    for (int i = 0; i < 8; i++) v |= (int64_t)b[i] << (8*i);
    return v;
}

inline float fread_float(FILE* f) {
    uint8_t b[4];
    if (fread(b, 1, 4, f) != 4) throw std::runtime_error("Unexpected EOF reading float");
    uint32_t bits = (uint32_t)b[0] | ((uint32_t)b[1] << 8) | ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
    float v;
    std::memcpy(&v, &bits, 4);
    return v;
}

inline std::string fread_string(FILE* f) {
    std::string s;
    int c;
    while ((c = fgetc(f)) != 0 && c != EOF) s += static_cast<char>(c);
    return s;
}

// ---- Back-patch helpers (seek + overwrite) ----

inline void patch_int64(FILE* f, int64_t file_pos, int64_t value) {
    long saved = ftell(f);
    fseek(f, static_cast<long>(file_pos), SEEK_SET);
    fwrite_int64(f, value);
    fseek(f, saved, SEEK_SET);
}

inline void patch_int32(FILE* f, int64_t file_pos, int32_t value) {
    long saved = ftell(f);
    fseek(f, static_cast<long>(file_pos), SEEK_SET);
    fwrite_int32(f, value);
    fseek(f, saved, SEEK_SET);
}
