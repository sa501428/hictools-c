#pragma once

#include "pair_parser.h"
#include <climits>
#include <cstring>
#include <fstream>
#include <set>
#include <zlib.h>

// HBS v1: independently implemented reader; see HBS_FORMAT.md.
class HbsIterator : public PairIterator {
public:
    HbsIterator(const std::string& path, const Genome& genome) : genome_(genome) {
        std::ifstream raw(path, std::ios::binary);
        unsigned char magic[2]{};
        raw.read(reinterpret_cast<char*>(magic), 2);
        require(raw && magic[0] == 0x1f && magic[1] == 0x8b, "input is not gzip");
        file_.reset(gzopen(path.c_str(), "rb"));
        require(bool(file_), "cannot open input");
        gzbuffer(file_.get(), 128 * 1024);
        unsigned char header[20];
        read(header, sizeof(header));
        require(std::memcmp(header, "HICBS\0\r\n", 8) == 0, "invalid magic");
        require(number(header + 8, 2) == 1, "unsupported version");
        require(number(header + 10, 2) == 0, "unsupported flags");
        resolution_ = static_cast<uint32_t>(number(header + 12, 4));
        require(resolution_ && resolution_ <= INT32_MAX, "invalid resolution");
        const auto n = number(header + 16, 4);
        require(n <= 65536, "too many chromosomes");
        std::set<int> used;
        size_t headerBytes = 20;
        for (uint64_t i = 0; i < n; ++i) {
            unsigned char word[8];
            read(word, 2);
            const auto length = number(word, 2);
            headerBytes += 10 + length;
            require(length && length <= 4096 && headerBytes <= 16 * 1024 * 1024, "invalid name/header length");
            std::string name(static_cast<size_t>(length), '\0');
            read(&name[0], static_cast<unsigned>(length));
            require(name.find('\0') == std::string::npos, "NUL in chromosome name");
            read(word, 8);
            const auto bp = number(word, 8);
            const int id = genome.index_of(name);
            require(id > 0, "unknown or synthetic chromosome: " + name);
            require(used.insert(id).second, "duplicate chromosome: " + name);
            require(bp && bp == static_cast<uint64_t>(genome.at(id).length), "chromosome length mismatch: " + name);
            ids_.push_back(id);
        }
    }

    uint32_t source_resolution() const override { return resolution_; }

    bool next(AlignmentPair& out) override {
        unsigned char record[14];
        if (!read(record, sizeof(record), true)) return false;
        auto a = number(record, 2), b = number(record + 6, 2);
        require(a < ids_.size() && b < ids_.size(), "chromosome ID outside table");
        auto x = number(record + 2, 4) * resolution_;
        auto y = number(record + 8, 4) * resolution_;
        require(x < static_cast<uint64_t>(genome_.at(ids_[a]).length) &&
                y < static_cast<uint64_t>(genome_.at(ids_[b]).length), "bin start outside chromosome");
        require(x <= INT32_MAX && y <= INT32_MAX, "position exceeds builder's int32 coordinate range");
        auto count = number(record + 12, 2);
        if (count == 65535) {
            unsigned char wide[8];
            read(wide, sizeof(wide));
            count = number(wide, 8);
            require(count >= 65535, "noncanonical escaped count");
        }
        out = AlignmentPair{};
        out.chr1 = ids_[a]; out.chr2 = ids_[b];
        out.pos1 = static_cast<int32_t>(x); out.pos2 = static_cast<int32_t>(y);
        out.frag2 = 1;
        out.has_exact_count = true;
        out.exact_count = count;
        out.score = static_cast<float>(count); // V9's existing float storage path.
        return true;
    }

    void close() override {
        if (file_) require(gzclose(file_.release()) == Z_OK, "gzip close failed");
    }

private:
    static void require(bool condition, const std::string& message) {
        if (!condition) throw std::runtime_error("HBS: " + message);
    }
    static uint64_t number(const unsigned char* b, unsigned size) {
        uint64_t n = 0;
        for (unsigned i = 0; i < size; ++i) n |= uint64_t(b[i]) << (8 * i);
        return n;
    }
    bool read(void* destination, unsigned size, bool eofAllowed = false) {
        require(bool(file_), "read after close");
        int got = gzread(file_.get(), destination, size);
        int error = Z_OK;
        const char* message = gzerror(file_.get(), &error);
        require(error == Z_OK || error == Z_STREAM_END, std::string("gzip read failed: ") + message);
        if (!got && eofAllowed && gzeof(file_.get())) return false;
        require(got == static_cast<int>(size), "truncated header or record");
        return true;
    }
    const Genome& genome_;
    std::unique_ptr<gzFile_s, decltype(&gzclose)> file_{nullptr, gzclose};
    std::vector<int> ids_;
    uint32_t resolution_ = 0;
};
