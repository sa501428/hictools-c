#include "repack.h"
#include <algorithm>

namespace hic10 {
namespace {
void write_bytes(FILE *output, const Bytes &bytes) {
    check(std::fwrite(bytes.data(), 1, bytes.size(), output) == bytes.size(),
          "cannot write repacked V10 file");
}
uint64_t position(FILE *output) {
    auto result = ftello(output);
    check(result >= 0, "cannot determine repacked V10 position");
    return result;
}
void patch_wide(FILE *output, uint64_t at, uint64_t value) {
    uint64_t saved = position(output);
    check(fseeko(output, at, SEEK_SET) == 0, "cannot seek repacked V10 file");
    Bytes field;
    put(field, value, 8);
    write_bytes(output, field);
    check(fseeko(output, saved, SEEK_SET) == 0, "cannot restore repacked V10 position");
}
uint32_t get_word(const Bytes &bytes, size_t at) {
    check(at <= bytes.size() && 4 <= bytes.size() - at, "truncated V10 u32 field");
    uint32_t value = 0;
    for (unsigned i = 0; i < 4; ++i)
        value |= uint32_t(bytes[at + i]) << (8 * i);
    return value;
}
uint64_t get_wide(const Bytes &bytes, size_t at) {
    check(at <= bytes.size() && 8 <= bytes.size() - at, "truncated V10 u64 field");
    uint64_t value = 0;
    for (unsigned i = 0; i < 8; ++i)
        value |= uint64_t(bytes[at + i]) << (8 * i);
    return value;
}
void set_wide(Bytes &bytes, size_t at, uint64_t value) {
    check(at <= bytes.size() && 8 <= bytes.size() - at, "truncated V10 u64 field");
    for (unsigned i = 0; i < 8; ++i)
        bytes[at + i] = static_cast<uint8_t>(value >> (8 * i));
}
uint64_t shifted(uint64_t value, int64_t delta) {
    if (delta >= 0)
        return plus(value, static_cast<uint64_t>(delta));
    uint64_t magnitude = static_cast<uint64_t>(-(delta + 1)) + 1;
    check(value >= magnitude, "V10 relocation underflow");
    return value - magnitude;
}
Bytes serialize_header(const Header &header) {
    Bytes bytes;
    magic(bytes, "HIC\0");
    put(bytes, 10, 4);
    bytes.resize(88, 0); // Includes reserved footer/NVI/EVI/NEVI locator space.
    str(bytes, header.genome);
    put(bytes, header.attributes.size(), 4);
    for (const auto &attribute : header.attributes) {
        str(bytes, attribute.first);
        str(bytes, attribute.second);
    }
    put(bytes, header.chromosomes.size(), 4);
    for (const auto &chromosome : header.chromosomes) {
        str(bytes, chromosome.name);
        put(bytes, chromosome.length, 8);
    }
    for (const auto &resolutions : header.resolutions) {
        put(bytes, resolutions.size(), 4);
        for (const auto &resolution : resolutions) {
            put(bytes, resolution.bin, 4);
            put(bytes, resolution.mode, 1);
            put(bytes, resolution.aggregation, 1);
            put(bytes, 0, 2);
            put(bytes, resolution.source, 4);
        }
    }
    if (!header.resolutions[1].empty())
        for (const auto &chromosome : header.chromosomes) {
            put(bytes, chromosome.sites.size(), 4);
            for (uint64_t site : chromosome.sites)
                put(bytes, site, 8);
        }
    put(bytes, header.norms.size(), 4);
    for (const auto &name : header.norms)
        str(bytes, name);
    set_wide(bytes, 8, bytes.size());
    return bytes;
}
} // namespace

Bytes repack_matrix_prefix(FILE *output, Reader &reader, const Header &header) {
    Bytes new_header = serialize_header(header);
    check(new_header.size() <= INT64_MAX && reader.header_length() <= INT64_MAX,
          "V10 header too large to relocate");
    int64_t delta = static_cast<int64_t>(new_header.size()) -
                    static_cast<int64_t>(reader.header_length());
    write_bytes(output, new_header);

    uint64_t vector_start = reader.vector_data_start();
    check(reader.header_length() <= vector_start, "V10 vectors overlap the header");
    constexpr uint64_t chunk = 8 * 1024 * 1024;
    for (uint64_t at = reader.header_length(); at < vector_start; at += chunk)
        write_bytes(output, reader.read_bytes(at, std::min(chunk, vector_start - at)));

    for (uint64_t field : reader.matrix_relocation_fields()) {
        check(field >= reader.header_length() && field + 8 <= vector_start,
              "V10 relocation field outside matrix section");
        Bytes old = reader.read_bytes(field, 8);
        patch_wide(output, shifted(field, delta), shifted(get_wide(old, 0), delta));
    }

    Bytes footer = reader.read_bytes(reader.footer().position, reader.footer().length);
    check(footer.size() >= 24 && std::equal(footer.begin(), footer.begin() + 4, "H10F") &&
              get_wide(footer, 8) == footer.size(),
          "invalid V10 footer during repack");
    uint32_t matrices = get_word(footer, 16);
    check(matrices == reader.matrices().size() && footer.size() == 24 + uint64_t(matrices) * 24,
          "invalid V10 footer matrix count");
    for (uint32_t i = 0; i < matrices; ++i) {
        size_t field = 24 + uint64_t(i) * 24 + 8;
        set_wide(footer, field, shifted(get_wide(footer, field), delta));
    }
    return footer;
}
} // namespace hic10
