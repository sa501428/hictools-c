#include "distributed_addnorm.h"
#include "addnorm.h"
#include "vectors.h"
#include "v10/text_vectors.h"
#include "v10/vector_updater.h"
#include "v10/expected_vector.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <tuple>
#include <fcntl.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zstd.h>

namespace hic10large {
namespace {
namespace fs = std::filesystem;
using hic10::Header;
using hic10::Reader;
using hic10::StoredVectorChunk;
using hic10::SuppliedTextVector;
constexpr uint64_t SEED = 14695981039346656037ULL;
uint64_t hash_file(const std::string &path) {
    File input(path, "rb");
    std::vector<char> buffer(1024 * 1024);
    uint64_t h = SEED;
    for (;;) {
        size_t n = std::fread(buffer.data(), 1, buffer.size(), input.get());
        h = fnv1a(buffer.data(), n, h);
        if (n < buffer.size()) {
            require(!std::ferror(input.get()), "cannot checksum " + path);
            return h;
        }
    }
}
void save(const std::string &path, const std::string &text) {
    std::string temporary = path + ".tmp-" + std::to_string(getpid());
    struct Cleanup {
        std::string path;
        ~Cleanup() { std::remove(path.c_str()); }
    } cleanup{temporary};
    File out(temporary, "wb");
    auto full = text + "checksum " + hex64(fnv1a(text.data(), text.size())) + "\n";
    out.write(full.data(), full.size());
    out.flush();
    require(fsync(fileno(out.get())) == 0, "cannot sync task manifest");
    out.close();
    atomic_rename(temporary, path);
}
std::string load(const std::string &path) {
    require(file_size(path) <= 512ULL * 1024 * 1024, "task manifest exceeds limit");
    std::ifstream in(path, std::ios::binary);
    require(bool(in), "missing task result: " + path);
    std::string text((std::istreambuf_iterator<char>(in)), {});
    auto pos = text.rfind("checksum ");
    require(pos != std::string::npos && (pos == 0 || text[pos - 1] == '\n') &&
                text.substr(pos) == "checksum " + hex64(fnv1a(text.data(), pos)) + "\n",
            "task manifest checksum mismatch: " + path);
    text.resize(pos);
    return text;
}
void end(std::istream &in) {
    std::string tag;
    require(bool(in >> tag) && tag == "end" && !(in >> tag), "truncated/trailing task manifest");
}
struct Lock {
    int fd;
    explicit Lock(const std::string &path) {
        fd = open((path + ".lock").c_str(), O_CREAT | O_RDWR, 0600);
        require(fd >= 0, "cannot create task lock");
        if (flock(fd, LOCK_EX | LOCK_NB) != 0) {
            close(fd);
            fail("task is already running: " + path);
        }
    }
    ~Lock() { close(fd); }
};
struct Scratch {
    std::string path;
    explicit Scratch(const std::string &parent) {
        make_directory(parent);
        auto pattern = join_path(parent, "hic-v10-import-XXXXXX");
        std::vector<char> bytes(pattern.begin(), pattern.end());
        bytes.push_back(0);
        require(mkdtemp(bytes.data()) != nullptr, "cannot create distributed addnorm scratch");
        path = bytes.data();
    }
    ~Scratch() {
        std::error_code e;
        fs::remove_all(path, e);
    }
};
struct Identity {
    uint64_t inode = 0, bytes = 0;
    int64_t mt = 0, mn = 0, ct = 0, cn = 0;
    bool operator==(const Identity &o) const {
        return std::tie(inode, bytes, mt, mn, ct, cn) ==
               std::tie(o.inode, o.bytes, o.mt, o.mn, o.ct, o.cn);
    }
};
Identity identity(const std::string &path) {
    struct stat s{};
    require(stat(path.c_str(), &s) == 0 && S_ISREG(s.st_mode), "cannot stat input " + path);
    Identity i;
    i.inode = s.st_ino;
    i.bytes = s.st_size;
#ifdef __APPLE__
    i.mt = s.st_mtimespec.tv_sec;
    i.mn = s.st_mtimespec.tv_nsec;
    i.ct = s.st_ctimespec.tv_sec;
    i.cn = s.st_ctimespec.tv_nsec;
#else
    i.mt = s.st_mtim.tv_sec;
    i.mn = s.st_mtim.tv_nsec;
    i.ct = s.st_ctim.tv_sec;
    i.cn = s.st_ctim.tv_nsec;
#endif
    return i;
}
void write_identity(std::ostream &out, const Identity &i) {
    out << i.inode << ' ' << i.bytes << ' ' << i.mt << ' ' << i.mn << ' ' << i.ct << ' ' << i.cn;
}
Identity read_identity(std::istream &in) {
    Identity i;
    in >> i.inode >> i.bytes >> i.mt >> i.mn >> i.ct >> i.cn;
    require(bool(in), "truncated input identity");
    return i;
}
uint64_t header_hash(Reader &r) {
    auto b = r.read_bytes(0, r.header_length());
    return fnv1a(b.data(), b.size());
}
bool real(const std::string &name) { return name != "ALL" && name != "All" && name != "all"; }
struct Bundle {
    uint8_t kind = 2, unit = 0;
    uint32_t norm = 0, ri = 0;
};
struct CellTask {
    uint32_t chr = 0, ri = 0, rows = 0;
    uint8_t unit = 0;
    std::vector<size_t> bundles;
};
struct ExpectedTask {
    size_t cell = 0;
    uint32_t row = 0;
    std::vector<size_t> bundles;
};
struct Slot {
    uint8_t unit = 0;
    VectorInfo info;
    std::vector<uint64_t> hashes;
};
struct CompressionTask {
    size_t slot = 0;
    uint64_t first = 0, last = 0;
};
struct Plan {
    std::string source, directory;
    Identity original;
    uint64_t original_header = 0, signature = 0;
    AddnormPlanOptions options;
    std::vector<std::pair<std::string, Identity>> inputs;
    std::vector<std::string> norms;
    std::vector<std::pair<std::string, std::string>> attributes;
    std::vector<SuppliedTextVector> supplied;
    std::vector<Bundle> bundles;
    std::vector<CellTask> cells;
    std::vector<ExpectedTask> expected;
    std::vector<Slot> slots;
    std::vector<CompressionTask> compression;
    std::unique_ptr<Reader> reader;
    Header header;
    std::string path(const std::string &name) const { return join_path(directory, name); }
    void verify_source() const {
        require(identity(source) == original,
                "v10 input changed since plan-addnorm; use a new plan");
    }
};
void write_info(std::ostream &out, const VectorInfo &i) {
    out << unsigned(i.kind) << ' ' << i.norm << ' ' << i.chr << ' ' << i.ri << ' ' << i.resolution
        << ' ' << i.words << ' ' << hex64(i.checksum) << ' ' << std::quoted(i.path);
}
VectorInfo read_info(std::istream &in) {
    VectorInfo i;
    unsigned kind;
    std::string h;
    in >> kind >> i.norm >> i.chr >> i.ri >> i.resolution >> i.words >> h >> std::quoted(i.path);
    require(bool(in) && kind <= 2, "invalid staged vector info");
    i.kind = kind;
    i.checksum = unhex64(h);
    return i;
}
std::vector<uint64_t> chunk_hashes(const VectorInfo &i) {
    VectorFileReader input(i);
    std::vector<uint64_t> hashes;
    for (uint64_t first = 0; first < i.words; first += 65536) {
        auto words =
            input.read(first, static_cast<uint32_t>(std::min<uint64_t>(65536, i.words - first)));
        hashes.push_back(fnv1a(words.data(), words.size() * 4));
    }
    if (!i.words)
        input.read(0, 0);
    return hashes;
}
void write_hashes(std::ostream &out, const std::vector<uint64_t> &hashes) {
    out << ' ' << hashes.size();
    for (auto h : hashes)
        out << ' ' << hex64(h);
}
std::vector<uint64_t> read_hashes(std::istream &in, uint64_t words) {
    uint64_t n;
    in >> n;
    require(bool(in) && n == hic10::ceil_div(words, 65536) && n <= UINT32_MAX,
            "invalid vector chunk hash count");
    std::vector<uint64_t> hashes;
    for (uint64_t j = 0; j < n; ++j) {
        std::string h;
        in >> h;
        require(bool(in), "truncated chunk hashes");
        hashes.push_back(unhex64(h));
    }
    return hashes;
}
void expand(Plan &p) {
    p.verify_source();
    p.reader.reset(new Reader(p.source));
    require(header_hash(*p.reader) == p.original_header, "v10 input header changed");
    p.header = p.reader->header();
    require(p.norms.size() >= p.header.norms.size() &&
                std::equal(p.header.norms.begin(), p.header.norms.end(), p.norms.begin()),
            "normalization IDs changed in plan");
    p.header.norms = p.norms;
    p.header.attributes.insert(p.header.attributes.end(), p.attributes.begin(), p.attributes.end());
    using Key = std::tuple<uint32_t, uint8_t, uint32_t>;
    std::map<Key, std::vector<size_t>> cells;
    for (const auto &v : p.supplied) {
        require(v.info.norm < p.norms.size() && v.info.chr < p.header.chromosomes.size() &&
                    v.unit <= 1 && p.header.resolution(v.unit, v.info.resolution) == v.info.ri &&
                    p.header.bins(v.info.chr, v.unit, v.info.ri) == v.info.words,
                "staged norm geometry mismatch");
    }
    for (size_t b = 0; b < p.bundles.size(); ++b) {
        const auto &bundle = p.bundles[b];
        require((bundle.kind == 1 || bundle.kind == 2) && bundle.unit <= 1 &&
                    bundle.ri < p.header.resolutions[bundle.unit].size() &&
                    bundle.norm < p.norms.size(),
                "invalid expected bundle");
        if (bundle.kind == 2) {
            for (const auto &v : p.supplied)
                if (v.info.norm == bundle.norm && v.unit == bundle.unit && v.info.ri == bundle.ri)
                    cells[{v.info.chr, bundle.unit, bundle.ri}].push_back(b);
        } else
            for (auto key : p.reader->matrices())
                if (key.chr1 == key.chr2 && real(p.header.chromosomes[key.chr1].name))
                    cells[{key.chr1, bundle.unit, bundle.ri}].push_back(b);
    }
    for (auto &item : cells) {
        CellTask c;
        std::tie(c.chr, c.unit, c.ri) = item.first;
        c.bundles = std::move(item.second);
        c.rows = hic10::narrow(std::max<uint64_t>(
            1, hic10::ceil_div(p.header.bins(c.chr, c.unit, c.ri), p.options.rows_per_task)));
        size_t id = p.cells.size();
        p.cells.push_back(c);
        for (uint32_t row = 0; row < c.rows; ++row)
            for (size_t begin = 0; begin < c.bundles.size(); begin += p.options.norms_per_task) {
                ExpectedTask task;
                task.cell = id;
                task.row = row;
                auto last = std::min(c.bundles.size(), begin + p.options.norms_per_task);
                task.bundles.assign(c.bundles.begin() + begin, c.bundles.begin() + last);
                p.expected.push_back(std::move(task));
            }
    }
    for (const auto &v : p.supplied)
        p.slots.push_back({v.unit, v.info, {}});
    for (size_t b = 0; b < p.bundles.size(); ++b) {
        const auto &bundle = p.bundles[b];
        Slot slot;
        slot.unit = bundle.unit;
        slot.info.kind = bundle.kind;
        slot.info.norm = bundle.norm;
        slot.info.chr = NO_ID;
        slot.info.ri = bundle.ri;
        slot.info.resolution = p.header.resolutions[bundle.unit][bundle.ri].bin;
        for (uint32_t chr = 0; chr < p.header.chromosomes.size(); ++chr)
            slot.info.words = std::max(slot.info.words, p.header.bins(chr, bundle.unit, bundle.ri));
        slot.info.path = p.path("expected-" + std::to_string(b) + ".h10w");
        p.slots.push_back(slot);
    }
    for (size_t slot = 0; slot < p.slots.size(); ++slot) {
        uint64_t chunks = hic10::ceil_div(p.slots[slot].info.words, 65536);
        for (uint64_t first = 0; first < chunks; first += p.options.chunks_per_task)
            p.compression.push_back(
                {slot, first, std::min(chunks, first + p.options.chunks_per_task)});
    }
}
Plan read_plan(const std::string &path, bool expand_tasks = true) {
    auto text = load(path);
    Plan p;
    p.directory = fs::absolute(fs::path(path).parent_path()).string();
    p.signature = fnv1a(text.data(), text.size());
    std::istringstream in(text);
    std::string magic, hash;
    unsigned version;
    in >> magic >> version >> std::quoted(p.source);
    require(magic == "HIC_V10_ADDNORM_PLAN" && version == 1, "invalid addnorm plan");
    p.original = read_identity(in);
    in >> hash;
    p.original_header = unhex64(hash);
    in >> p.options.rows_per_task >> p.options.chunks_per_task >> p.options.norms_per_task >>
        p.options.level;
    require(bool(in) && p.options.rows_per_task && p.options.chunks_per_task &&
                p.options.norms_per_task,
            "invalid plan options");
    size_t n;
    in >> n;
    require(bool(in) && n <= text.size(), "invalid input count");
    for (size_t j = 0; j < n; ++j) {
        std::string path;
        in >> std::quoted(path);
        p.inputs.push_back({path, read_identity(in)});
    }
    in >> n;
    require(bool(in) && n <= text.size(), "invalid norm count");
    for (size_t j = 0; j < n; ++j) {
        std::string norm;
        in >> std::quoted(norm);
        p.norms.push_back(norm);
    }
    in >> n;
    require(bool(in) && n <= text.size(), "invalid attribute count");
    for (size_t j = 0; j < n; ++j) {
        std::string key, value;
        in >> std::quoted(key) >> std::quoted(value);
        p.attributes.emplace_back(key, value);
    }
    in >> n;
    require(bool(in) && n <= text.size(), "invalid supplied vector count");
    for (size_t j = 0; j < n; ++j) {
        SuppliedTextVector v;
        unsigned unit;
        in >> unit;
        require(unit <= 1, "invalid vector unit");
        v.unit = unit;
        v.info = read_info(in);
        v.info.path = p.path(v.info.path);
        p.supplied.push_back(v);
    }
    // Hashes for normalization chunks are prepared once, avoiding a full scan
    // by every parallel compression job.
    std::vector<std::vector<uint64_t>> hashes;
    for (const auto &v : p.supplied)
        hashes.push_back(read_hashes(in, v.info.words));
    in >> n;
    require(bool(in) && n <= text.size(), "invalid bundle count");
    for (size_t j = 0; j < n; ++j) {
        Bundle b;
        unsigned kind, unit;
        in >> kind >> unit >> b.norm >> b.ri;
        require(kind <= 2 && unit <= 1, "invalid bundle");
        b.kind = kind;
        b.unit = unit;
        p.bundles.push_back(b);
    }
    end(in);
    if (expand_tasks) {
        expand(p);
        for (size_t j = 0; j < hashes.size(); ++j)
            p.slots[j].hashes = std::move(hashes[j]);
    }
    return p;
}
std::string serialize(const Plan &p, const std::vector<std::vector<uint64_t>> &hashes) {
    std::ostringstream out;
    out << "HIC_V10_ADDNORM_PLAN 1\n" << std::quoted(p.source) << ' ';
    write_identity(out, p.original);
    out << ' ' << hex64(p.original_header) << '\n'
        << p.options.rows_per_task << ' ' << p.options.chunks_per_task << ' '
        << p.options.norms_per_task << ' ' << p.options.level << '\n';
    out << p.inputs.size() << '\n';
    for (const auto &i : p.inputs) {
        out << std::quoted(i.first) << ' ';
        write_identity(out, i.second);
        out << '\n';
    }
    out << p.norms.size();
    for (const auto &norm : p.norms)
        out << ' ' << std::quoted(norm);
    out << '\n';
    out << p.attributes.size() << '\n';
    for (const auto &a : p.attributes)
        out << std::quoted(a.first) << ' ' << std::quoted(a.second) << '\n';
    out << p.supplied.size() << '\n';
    for (const auto &v : p.supplied) {
        out << unsigned(v.unit) << ' ';
        write_info(out, v.info);
        out << '\n';
    }
    for (const auto &h : hashes) {
        write_hashes(out, h);
        out << '\n';
    }
    out << p.bundles.size() << '\n';
    for (auto b : p.bundles)
        out << unsigned(b.kind) << ' ' << unsigned(b.unit) << ' ' << b.norm << ' ' << b.ri << '\n';
    out << "end\n";
    return out.str();
}
std::string result_header(const Plan &p, const std::string &phase, size_t task) {
    return "HIC_V10_ADDNORM_RESULT 1 " + hex64(p.signature) + " " + phase + " " +
           std::to_string(task) + "\n";
}
std::istringstream result(const Plan &p, const std::string &path, const std::string &phase,
                          size_t task) {
    std::istringstream in(load(path));
    std::string magic, hash, name;
    unsigned version;
    size_t id;
    in >> magic >> version >> hash >> name >> id;
    require(magic == "HIC_V10_ADDNORM_RESULT" && version == 1 && unhex64(hash) == p.signature &&
                name == phase && id == task,
            "task result belongs to another plan/task");
    return in;
}
std::string cell_path(const Plan &p, size_t id, uint32_t row) {
    return p.path("cells-" + std::to_string(id) + "-row-" + std::to_string(row) + ".h10r");
}
std::vector<RunInfo> cell_results(const Plan &p, size_t id) {
    auto in = result(p, p.path("cells-" + std::to_string(id) + ".done"), "cells", id);
    size_t n;
    in >> n;
    require(n == p.cells.at(id).rows, "missing cell row shards");
    std::vector<RunInfo> runs;
    const auto &c = p.cells[id];
    for (size_t j = 0; j < n; ++j) {
        RunInfo r;
        std::string h;
        in >> r.records >> h;
        require(bool(in), "truncated cell result");
        r.checksum = unhex64(h);
        r.path = cell_path(p, id, j);
        r.chr1 = r.chr2 = c.chr;
        r.resolution = p.header.resolutions[c.unit][c.ri].bin;
        auto actual = inspect_run(r.path, false);
        require(actual.checksum == r.checksum && actual.records == r.records &&
                    actual.chr1 == r.chr1 && actual.chr2 == r.chr2 &&
                    actual.resolution == r.resolution,
                "cell sidecar identity mismatch");
        runs.push_back(r);
    }
    end(in);
    return runs;
}
void cells(Plan &p, size_t id, const NormalizeOptions &options) {
    const auto &task = p.cells.at(id);
    auto done = p.path("cells-" + std::to_string(id) + ".done");
    Lock lock(done);
    if (path_exists(done)) {
        for (auto r : cell_results(p, id))
            inspect_run(r.path, true);
        return;
    }
    require(options.memory_bytes >= 32768 && options.scale_options.merge_fan_in >= 2,
            "invalid cell sort memory/fan-in");
    Scratch work(options.temporary_directory.empty() ? p.directory : options.temporary_directory);
    auto r = p.header.resolutions[task.unit][task.ri];
    RunInfo target;
    bool cis = false;
    for (auto key : p.reader->matrices())
        cis |= key.chr1 == task.chr && key.chr2 == task.chr;
    if (cis) {
        uint32_t source = r.mode ? r.source : task.ri;
        auto run =
            extract_addnorm_cis_run(*p.reader, task.chr, task.unit, source, work.path, options);
        if (!run.path.empty())
            target = r.mode ? rollup_cell_file(run, work.path, r.bin, "target") : run;
    }
    std::unique_ptr<RunReader> input;
    if (!target.path.empty())
        input.reset(new RunReader(target));
    CellRecord cell;
    bool have = input && input->next(cell);
    std::ostringstream out;
    out << result_header(p, "cells", id) << task.rows << '\n';
    for (uint32_t row = 0; row < task.rows; ++row) {
        RunWriter output(cell_path(p, id, row), task.chr, task.chr, r.bin);
        uint64_t lo = uint64_t(row) * p.options.rows_per_task, hi = lo + p.options.rows_per_task;
        while (have && cell.y < hi) {
            require(cell.y >= lo, "unordered cell run");
            output.add(cell);
            have = input->next(cell);
        }
        auto info = output.finish();
        {
            File durable(info.path, "rb");
            require(fsync(fileno(durable.get())) == 0, "cannot sync cell shard");
        }
        out << info.records << ' ' << hex64(info.checksum) << '\n';
    }
    require(!have, "unconsumed cell rows");
    p.verify_source();
    out << "end\n";
    save(done, out.str());
}
// Portable lossless positive long-double representation: 120-bit significand
// and signed binary exponent. Padding/ABI bytes never enter durable files.
void put_real(unsigned char *b, long double x) {
    require(std::isfinite(x) && x >= 0 && std::numeric_limits<long double>::digits <= 120,
            "invalid expected accumulator");
    int exponent = 0;
    auto mantissa = static_cast<unsigned __int128>(std::ldexp(std::frexp(x, &exponent), 120));
    put_u32(b, static_cast<uint32_t>(exponent));
    put_u64(b + 4, static_cast<uint64_t>(mantissa));
    put_u64(b + 12, static_cast<uint64_t>(mantissa >> 64));
}
long double get_real(const unsigned char *b) {
    auto mantissa = (static_cast<unsigned __int128>(get_u64(b + 12)) << 64) | get_u64(b + 4);
    int64_t exponent = int32_t(get_u32(b));
    require(exponent >= std::numeric_limits<long double>::min_exponent - 120 &&
                exponent <= std::numeric_limits<long double>::max_exponent + 120 &&
                mantissa < (static_cast<unsigned __int128>(1) << 120),
            "invalid partial exponent/significand");
    auto value = std::ldexp(static_cast<long double>(mantissa), static_cast<int>(exponent - 120));
    require(std::isfinite(value) && value >= 0, "invalid expected partial");
    return value;
}
struct Hist {
    std::vector<long double> actual;
    long double observed = 0;
    explicit Hist(uint32_t n) : actual(n, 0) {}
};
std::string hist_path(const Plan &p, size_t task, size_t b) {
    return p.path("partial-" + std::to_string(task) + "-b" + std::to_string(b) + ".h10e");
}
void write_hist(const std::string &path, const Hist &h) {
    auto temp = path + ".tmp-" + std::to_string(getpid());
    struct Cleanup {
        std::string path;
        ~Cleanup() { std::remove(path.c_str()); }
    } cleanup{temp};
    File out(temp, "wb");
    std::array<unsigned char, 32> head{};
    std::memcpy(head.data(), "H10E", 4);
    put_u32(head.data() + 4, 1);
    put_u32(head.data() + 8, h.actual.size());
    put_real(head.data() + 12, h.observed);
    out.write(head.data(), head.size());
    for (uint32_t d = 0; d < h.actual.size(); ++d)
        if (h.actual[d] > 0) {
            std::array<unsigned char, 24> v{};
            put_u32(v.data(), d);
            put_real(v.data() + 4, h.actual[d]);
            out.write(v.data(), v.size());
        }
    out.flush();
    require(fsync(fileno(out.get())) == 0, "cannot sync expected partial");
    out.close();
    atomic_rename(temp, path);
}
std::map<size_t, std::pair<uint64_t, uint64_t>> expected_results(const Plan &p, size_t id) {
    auto in = result(p, p.path("expected-" + std::to_string(id) + ".done"), "expected", id);
    size_t n;
    in >> n;
    require(n == p.expected.at(id).bundles.size(), "expected partial count mismatch");
    std::map<size_t, std::pair<uint64_t, uint64_t>> outputs;
    for (auto b : p.expected[id].bundles) {
        uint64_t bytes;
        std::string h;
        in >> bytes >> h;
        require(bool(in), "truncated expected result");
        outputs[b] = {bytes, unhex64(h)};
    }
    end(in);
    return outputs;
}
void expected(Plan &p, size_t id) {
    const auto &task = p.expected.at(id);
    const auto &c = p.cells[task.cell];
    auto done = p.path("expected-" + std::to_string(id) + ".done");
    Lock lock(done);
    if (path_exists(done)) {
        for (auto item : expected_results(p, id)) {
            auto path = hist_path(p, id, item.first);
            require(file_size(path) == item.second.first && hash_file(path) == item.second.second,
                    "expected partial checksum mismatch");
        }
        return;
    }
    auto runs = cell_results(p, task.cell);
    RunReader input(runs[task.row]);
    uint32_t n = hic10::narrow(p.header.bins(c.chr, c.unit, c.ri));
    std::vector<Hist> hist;
    std::vector<std::vector<uint32_t>> norms;
    for (auto b : task.bundles) {
        hist.emplace_back(n);
        std::vector<uint32_t> values;
        if (p.bundles[b].kind == 2) {
            auto v = std::find_if(
                p.supplied.begin(), p.supplied.end(), [&](const SuppliedTextVector &v) {
                    return v.unit == c.unit && v.info.ri == c.ri && v.info.chr == c.chr &&
                           v.info.norm == p.bundles[b].norm;
                });
            require(v != p.supplied.end(), "missing supplied norm for expected task");
            VectorFileReader reader(v->info);
            values = reader.read(0, n);
        }
        norms.push_back(std::move(values));
    }
    CellRecord cell;
    uint64_t lo = uint64_t(task.row) * p.options.rows_per_task, hi = lo + p.options.rows_per_task;
    while (input.next(cell)) {
        require(cell.x <= cell.y && cell.y < n && cell.y >= lo && cell.y < hi,
                "contact outside expected shard");
        for (size_t j = 0; j < task.bundles.size(); ++j) {
            long double value = cell.count;
            if (!norms[j].empty()) {
                float a = bits_float(norms[j][cell.x]), b = bits_float(norms[j][cell.y]);
                if (!(a > 0 && b > 0 && std::isfinite(a) && std::isfinite(b)))
                    continue;
                value /= static_cast<long double>(a) * b;
            }
            if (std::isfinite(value) && value > 0) {
                hist[j].actual[cell.y - cell.x] += value;
                hist[j].observed += value;
            }
        }
    }
    std::ostringstream out;
    out << result_header(p, "expected", id) << task.bundles.size() << '\n';
    for (size_t j = 0; j < task.bundles.size(); ++j) {
        auto path = hist_path(p, id, task.bundles[j]);
        write_hist(path, hist[j]);
        out << file_size(path) << ' ' << hex64(hash_file(path)) << '\n';
    }
    p.verify_source();
    out << "end\n";
    save(done, out.str());
}
std::string vector_done(const Plan &p, size_t slot) {
    return p.path("vector-" + std::to_string(slot) + ".done");
}
Slot slot_result(const Plan &p, size_t id) {
    if (id < p.supplied.size())
        return p.slots.at(id);
    auto in = result(p, vector_done(p, id), "reduce", id - p.supplied.size());
    Slot s;
    s.unit = p.slots.at(id).unit;
    s.info = read_info(in);
    s.hashes = read_hashes(in, s.info.words);
    const auto &wanted = p.slots[id].info;
    require(s.info.kind == wanted.kind && s.info.norm == wanted.norm && s.info.chr == wanted.chr &&
                s.info.ri == wanted.ri && s.info.resolution == wanted.resolution &&
                s.info.words == wanted.words && s.info.path == wanted.path,
            "expected vector identity mismatch");
    size_t n;
    in >> n;
    require(n <= p.header.chromosomes.size(), "invalid expected scales");
    for (size_t j = 0; j < n; ++j) {
        uint32_t chr;
        std::string bits;
        in >> chr >> bits;
        auto word = unhex64(bits);
        require(chr < p.header.chromosomes.size() && word <= UINT32_MAX &&
                    s.info.scales.emplace(chr, word).second,
                "invalid/duplicate expected scale");
    }
    end(in);
    VectorFileReader check(s.info);
    require(file_size(s.info.path) == 64 + s.info.words * 4, "expected sidecar size mismatch");
    return s;
}
void reduce(Plan &p, size_t b) {
    auto id = p.supplied.size() + b;
    auto done = vector_done(p, id);
    Lock lock(done);
    if (path_exists(done)) {
        auto s = slot_result(p, id);
        require(chunk_hashes(s.info) == s.hashes, "expected sidecar checksum mismatch");
        return;
    }
    const auto &bundle = p.bundles.at(b);
    uint32_t n = hic10::narrow(p.slots[id].info.words);
    std::vector<long double> actual(n, 0);
    std::map<uint32_t, long double> observed;
    std::map<uint32_t, uint32_t> bins;
    for (size_t task = 0; task < p.expected.size(); ++task) {
        const auto &e = p.expected[task];
        if (std::find(e.bundles.begin(), e.bundles.end(), b) == e.bundles.end())
            continue;
        auto outputs = expected_results(p, task);
        auto info = outputs.at(b);
        auto path = hist_path(p, task, b);
        require(file_size(path) == info.first && hash_file(path) == info.second,
                "expected partial checksum mismatch");
        const auto &c = p.cells[e.cell];
        File input(path, "rb");
        std::array<unsigned char, 32> head{};
        input.read(head.data(), head.size());
        uint32_t length = hic10::narrow(p.header.bins(c.chr, c.unit, c.ri));
        require(std::memcmp(head.data(), "H10E", 4) == 0 && get_u32(head.data() + 4) == 1 &&
                    get_u32(head.data() + 8) == length && info.first >= 32 &&
                    (info.first - 32) % 24 == 0,
                "invalid expected partial header");
        auto total = get_real(head.data() + 12);
        if (total > 0) {
            observed[c.chr] += total;
            bins[c.chr] = length;
        }
        uint32_t previous = 0;
        bool have = false;
        for (uint64_t j = 0; j < (info.first - 32) / 24; ++j) {
            std::array<unsigned char, 24> v{};
            input.read(v.data(), v.size());
            auto d = get_u32(v.data());
            require(d < length && (!have || d > previous), "invalid expected distance");
            previous = d;
            have = true;
            actual[d] += get_real(v.data() + 4);
        }
    }
    auto v = hic10::finish_expected(actual, observed, bins, bundle.kind, bundle.norm, bundle.unit,
                                    bundle.ri, true);
    auto info = p.slots[id].info;
    info.scales = v.scales;
    VectorFileWriter writer(info.path, info);
    writer.add(v.values);
    info = writer.finish();
    {
        File durable(info.path, "rb");
        require(fsync(fileno(durable.get())) == 0, "cannot sync expected vector");
    }
    auto hashes = chunk_hashes(info);
    std::ostringstream out;
    out << result_header(p, "reduce", b);
    write_info(out, info);
    write_hashes(out, hashes);
    out << ' ' << info.scales.size();
    for (auto s : info.scales)
        out << ' ' << s.first << ' ' << hex64(s.second);
    out << "\nend\n";
    p.verify_source();
    save(done, out.str());
}
std::string pack_path(const Plan &p, size_t id) {
    return p.path("compressed-" + std::to_string(id) + ".pack");
}
std::vector<StoredVectorChunk> compression_result(const Plan &p, size_t id) {
    auto in = result(p, p.path("compressed-" + std::to_string(id) + ".done"), "compress", id);
    const auto &task = p.compression.at(id);
    size_t n;
    in >> n;
    require(n == task.last - task.first, "compressed fragment count mismatch");
    std::vector<StoredVectorChunk> result;
    uint64_t next = 0;
    auto size = file_size(pack_path(p, id));
    for (uint64_t j = task.first; j < task.last; ++j) {
        StoredVectorChunk c;
        unsigned t;
        std::string h;
        in >> c.begin >> c.count >> c.offset >> c.stored >> t >> h;
        require(bool(in) && c.begin == j * 65536 &&
                    c.count == std::min<uint64_t>(65536, p.slots[task.slot].info.words - c.begin) &&
                    c.offset == next && c.stored > 16 && t <= 2 && c.offset <= size &&
                    c.stored <= size - c.offset,
                "invalid compressed result coverage");
        c.transform = t;
        c.checksum = unhex64(h);
        c.path = pack_path(p, id);
        next += c.stored;
        result.push_back(c);
    }
    require(next == size, "compressed fragment size mismatch");
    end(in);
    return result;
}
void compress(Plan &p, size_t id) {
    const auto &task = p.compression.at(id);
    auto done = p.path("compressed-" + std::to_string(id) + ".done");
    Lock lock(done);
    auto slot = slot_result(p, task.slot);
    if (path_exists(done)) {
        File input(pack_path(p, id), "rb");
        for (auto c : compression_result(p, id)) {
            std::vector<unsigned char> b(c.stored);
            input.read(b.data(), b.size());
            require(fnv1a(b.data(), b.size()) == c.checksum,
                    "compressed fragment checksum mismatch");
        }
        return;
    }
    VectorFileReader validate(slot.info);
    require(file_size(slot.info.path) == 64 + slot.info.words * 4, "vector sidecar size mismatch");
    File input(slot.info.path, "rb");
    auto path = pack_path(p, id), temp = path + ".tmp-" + std::to_string(getpid());
    struct Cleanup {
        std::string path;
        ~Cleanup() { std::remove(path.c_str()); }
    } cleanup{temp};
    File output(temp, "wb");
    std::ostringstream out;
    out << result_header(p, "compress", id) << task.last - task.first << '\n';
    for (uint64_t chunk = task.first; chunk < task.last; ++chunk) {
        uint64_t begin = chunk * 65536;
        uint32_t count = static_cast<uint32_t>(std::min<uint64_t>(65536, slot.info.words - begin));
        std::vector<uint32_t> words(count);
        input.seek(64 + begin * 4);
        input.read(words.data(), words.size() * 4);
        require(fnv1a(words.data(), words.size() * 4) == slot.hashes.at(chunk),
                "normalization/expected chunk checksum mismatch");
        hic10::Bytes best;
        uint8_t transform = 0;
        for (uint8_t t = 0; t < 3; ++t) {
            hic10::Bytes raw;
            raw.reserve(count * 4);
            if (t == 1) {
                for (unsigned lane = 0; lane < 4; ++lane)
                    for (auto w : words)
                        hic10::put(raw, w >> (8 * lane), 1);
            } else
                for (uint32_t j = 0; j < count; ++j)
                    hic10::put(raw, words[j] ^ (t == 2 && j ? words[j - 1] : 0), 4);
            hic10::Bytes frame(ZSTD_compressBound(raw.size()));
            auto size =
                ZSTD_compress(frame.data(), frame.size(), raw.data(), raw.size(), p.options.level);
            require(!ZSTD_isError(size), "vector compression failed");
            frame.resize(size);
            if (best.empty() || frame.size() < best.size()) {
                best = std::move(frame);
                transform = t;
            }
        }
        hic10::Bytes bytes;
        hic10::magic(bytes, "H10V");
        hic10::put(bytes, 1, 1);
        hic10::put(bytes, transform, 1);
        hic10::put(bytes, 0, 2);
        hic10::put(bytes, uint64_t(count) * 4, 4);
        hic10::put(bytes, count, 4);
        hic10::append(bytes, best);
        out << begin << ' ' << count << ' ' << output.tell() << ' ' << bytes.size() << ' '
            << unsigned(transform) << ' ' << hex64(fnv1a(bytes.data(), bytes.size())) << '\n';
        output.write(bytes.data(), bytes.size());
    }
    output.flush();
    require(fsync(fileno(output.get())) == 0, "cannot sync compressed fragments");
    output.close();
    p.verify_source();
    atomic_rename(temp, path);
    out << "end\n";
    save(done, out.str());
}
std::vector<size_t> phase_ids(const Plan &p, const std::string &phase) {
    size_t count = phase == "cells"      ? p.cells.size()
                   : phase == "expected" ? p.expected.size()
                   : phase == "reduce"   ? p.bundles.size()
                                         : p.compression.size();
    require(phase == "cells" || phase == "expected" || phase == "reduce" ||
                phase == "compress-norm" || phase == "compress-expected",
            "unknown addnorm phase " + phase);
    std::vector<size_t> ids;
    for (size_t j = 0; j < count; ++j)
        if (phase == "cells" || phase == "expected" || phase == "reduce" ||
            (p.compression[j].slot < p.supplied.size()) == (phase == "compress-norm"))
            ids.push_back(j);
    return ids;
}
std::string shell_quote(const std::string &s) {
    std::string result = "'";
    for (char c : s)
        result += c == '\'' ? "'\\''" : std::string(1, c);
    return result + "'";
}
bool completed(const Plan &p) {
    auto path = p.path("assembled.done");
    if (!path_exists(path))
        return false;
    auto in = result(p, path, "assemble", 0);
    uint64_t inode, bytes;
    int64_t mt, mn;
    std::string hash;
    in >> inode >> bytes >> mt >> mn >> hash;
    end(in);
    auto actual = identity(p.source);
    if (actual.inode != inode || actual.bytes != bytes || actual.mt != mt || actual.mn != mn)
        return false;
    Reader reader(p.source);
    return header_hash(reader) == unhex64(hash);
}
} // namespace

void plan_addnorm(const std::string &source, const std::vector<std::string> &texts,
                  const std::string &directory, const AddnormPlanOptions &options) {
    require(!texts.empty() && options.rows_per_task && options.chunks_per_task &&
                options.norms_per_task,
            "plan-addnorm needs norm files and positive task sizes");
    require(options.level >= ZSTD_minCLevel() && options.level <= ZSTD_maxCLevel(),
            "unsupported Zstandard compression level");
    auto final = fs::absolute(directory).lexically_normal().string(),
         plan_path = join_path(final, "addnorm.plan");
    std::vector<std::pair<std::string, Identity>> inputs;
    for (const auto &path : texts)
        inputs.push_back({fs::absolute(path).lexically_normal().string(), identity(path)});
    if (path_exists(plan_path)) {
        auto p = read_plan(plan_path, false);
        require(p.source == fs::absolute(source).lexically_normal().string() &&
                    p.inputs == inputs && p.options.rows_per_task == options.rows_per_task &&
                    p.options.chunks_per_task == options.chunks_per_task &&
                    p.options.norms_per_task == options.norms_per_task &&
                    p.options.level == options.level,
                "existing plan has different inputs/options; use another directory");
        p.verify_source();
        std::cout << plan_path << '\n';
        return;
    }
    require(
        !path_exists(final),
        "addnorm work directory already exists without a completed plan; use another directory");
    make_directory(fs::path(final).parent_path().string());
    Scratch staging(fs::path(final).parent_path().string());
    Plan p;
    p.source = fs::absolute(source).lexically_normal().string();
    p.directory = final;
    p.original = identity(p.source);
    p.inputs = inputs;
    p.options = options;
    Reader reader(p.source);
    p.original_header = header_hash(reader);
    Header header = reader.header();
    using BKey = std::tuple<uint32_t, uint8_t, uint32_t>;
    std::set<BKey> retained;
    std::set<std::pair<uint8_t, uint32_t>> raw;
    for (auto v : reader.vector_entries()) {
        if (v.kind == 1)
            raw.insert({v.unit, v.ri});
        else
            retained.insert({v.norm, v.unit, v.ri});
    }
    std::set<std::tuple<uint32_t, uint32_t, uint8_t, uint32_t>> keys;
    std::set<BKey> new_bundles;
    std::vector<std::vector<uint64_t>> hashes;
    for (size_t j = 0; j < inputs.size(); ++j) {
        auto sub = join_path(staging.path, "input-" + std::to_string(j));
        make_directory(sub);
        auto vectors = hic10::stage_text_vectors(inputs[j].first, header, sub);
        for (auto &v : vectors) {
            require(keys.insert({v.info.norm, v.info.chr, v.unit, v.info.ri}).second,
                    "duplicate supplied normalization across norm files");
            BKey key{v.info.norm, v.unit, v.info.ri};
            if (retained.count(key)) {
                std::cerr << "Preserving existing " << header.norms[v.info.norm] << ' '
                          << (v.unit ? "FRAG" : "BP") << ' ' << v.info.resolution << '\n';
                continue;
            }
            if (!v.source_geometry.empty())
                p.attributes.emplace_back(
                    "hictools.import.vector.0." + std::to_string(v.info.norm) + "." +
                        std::to_string(v.info.chr) + "." + std::to_string(v.unit) + "." +
                        std::to_string(v.info.ri),
                    v.source_geometry);
            {
                File durable(v.info.path, "rb");
                require(fsync(fileno(durable.get())) == 0, "cannot sync supplied vector");
            }
            hashes.push_back(chunk_hashes(v.info));
            v.info.path = fs::relative(v.info.path, staging.path).string();
            p.supplied.push_back(v);
            new_bundles.insert(key);
        }
    }
    p.norms = header.norms;
    std::set<std::pair<uint8_t, uint32_t>> missing;
    for (auto key : new_bundles) {
        auto [norm, unit, ri] = key;
        p.bundles.push_back({2, unit, norm, ri});
        if (!raw.count({unit, ri}))
            missing.insert({unit, ri});
    }
    for (auto key : missing)
        p.bundles.push_back({1, key.first, 0, key.second});
    p.verify_source();
    for (auto i : inputs)
        require(identity(i.first) == i.second, "normalization text changed during preparation");
    save(join_path(staging.path, "addnorm.plan"), serialize(p, hashes));
    atomic_rename(staging.path, final);
    std::cout << plan_path << '\n';
}
void print_addnorm_tasks(const std::string &path, const std::string &phase, bool count_only,
                         bool dependencies) {
    require(!(count_only && dependencies), "--count and --dependencies are mutually exclusive");
    auto p = read_plan(path);
    auto ids = phase_ids(p, phase);
    if (count_only) {
        std::cout << ids.size() << '\n';
        return;
    }
    if (dependencies) {
        // Ordinals in this machine-readable DAG use the same phase-local IDs
        // as addnorm-run. Compression's internal ID is deliberately hidden.
        std::vector<std::vector<size_t>> expected_by_bundle(p.bundles.size());
        if (phase == "reduce")
            for (size_t task = 0; task < p.expected.size(); ++task)
                for (auto b : p.expected[task].bundles)
                    expected_by_bundle[b].push_back(task);
        for (size_t ordinal = 0; ordinal < ids.size(); ++ordinal) {
            auto id = ids[ordinal];
            std::cout << ordinal;
            if (phase == "expected")
                std::cout << " cells:" << p.expected[id].cell;
            else if (phase == "reduce")
                for (auto task : expected_by_bundle[id])
                    std::cout << " expected:" << task;
            else if (phase == "compress-expected")
                std::cout << " reduce:" << p.compression[id].slot - p.supplied.size();
            std::cout << '\n';
        }
        return;
    }
    for (size_t ordinal = 0; ordinal < ids.size(); ++ordinal)
        std::cout << "hic_v10_large addnorm-run " << shell_quote(fs::absolute(path).string()) << ' '
                  << phase << ' ' << ordinal << '\n';
}
void run_addnorm_task(const std::string &path, const std::string &phase, uint64_t ordinal,
                      const NormalizeOptions &options) {
    auto p = read_plan(path);
    auto ids = phase_ids(p, phase);
    require(ordinal < ids.size(), "addnorm task index outside phase");
    auto id = ids[ordinal];
    if (phase == "cells")
        cells(p, id, options);
    else if (phase == "expected")
        expected(p, id);
    else if (phase == "reduce")
        reduce(p, id);
    else
        compress(p, id);
}
void assemble_addnorm(const std::string &path, uint32_t threads) {
    require(threads > 0, "assembly threads must be positive");
    auto initial = read_plan(path, false);
    Lock plan_lock(initial.path("assemble"));
    Lock target_lock(initial.source + ".distributed-addnorm");
    if (completed(initial)) {
        std::cerr << "Distributed addnorm already assembled\n";
        return;
    }
    auto p = read_plan(path);
    std::vector<std::vector<StoredVectorChunk>> chunks(p.slots.size());
    // Require every dependency, including otherwise empty expected vectors,
    // before opening a replacement file or copying the large raw-data section.
    for (size_t slot = 0; slot < p.slots.size(); ++slot)
        slot_result(p, slot);
    for (size_t id = 0; id < p.compression.size(); ++id) {
        auto c = compression_result(p, id);
        auto &out = chunks[p.compression[id].slot];
        out.insert(out.end(), c.begin(), c.end());
    }
    auto publish_marker = [&](const std::string &output) {
        p.verify_source();
        Reader ready(output);
        auto info = identity(output);
        std::ostringstream out;
        out << result_header(p, "assemble", 0) << info.inode << ' ' << info.bytes << ' ' << info.mt
            << ' ' << info.mn << ' ' << hex64(header_hash(ready)) << "\nend\n";
        save(p.path("assembled.done"), out.str());
        p.verify_source();
    };
    if (p.slots.empty()) {
        publish_marker(p.source);
        return;
    }
    hic10::VectorOutput output(*p.reader, p.header, p.source, p.options.level, threads);
    for (size_t id = 0; id < p.slots.size(); ++id) {
        auto s = slot_result(p, id);
        hic10::Vector v;
        v.kind = s.info.kind;
        v.unit = s.unit;
        v.norm = s.info.norm;
        v.chr = s.info.kind ? 0 : s.info.chr;
        v.ri = s.info.ri;
        v.scales = s.info.scales;
        output.add_precompressed(v, s.info.words, chunks[id]);
    }
    output.finish(p.header, publish_marker);
    std::cerr << "Distributed normalization vectors assembled: " << p.source << '\n';
}
} // namespace hic10large
