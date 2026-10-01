#include "vector_import.h"
#include "vector_updater.h"
#include "v10_large/vectors.h"
#include <set>
#include <tuple>

namespace hic10 {
void add_vector_manifest_v10(const std::string &path, const std::string &manifest_path, int level) {
    auto manifest = hic10large::read_vector_manifest(manifest_path);
    Reader reader(path);
    Header header = reader.header();
    std::vector<uint32_t> ids;
    std::set<std::string> names;
    for (const auto &name : manifest.norms) {
        check(names.insert(name).second, "duplicate imported normalization name");
        check(!name.empty() && name != "NONE", "invalid imported normalization name");
        auto found = std::find(header.norms.begin(), header.norms.end(), name);
        if (found == header.norms.end()) {
            ids.push_back(narrow(header.norms.size()));
            header.norms.push_back(name);
        } else ids.push_back(narrow(found - header.norms.begin()));
    }
    auto existing = reader.vector_entries();
    // Preserve type/resolution bundles together: mixing old norms with a new
    // expected vector computed from other norms would change O/E semantics.
    std::set<std::pair<uint32_t, uint32_t>> retained;
    for (const auto &e : existing)
        if (e.kind != 1 && e.unit == 0) retained.insert({e.norm, e.ri});
    VectorOutput output(reader, header, path, level);
    std::set<std::tuple<uint8_t, uint32_t, uint32_t, uint32_t>> keys;
    for (const auto &info : manifest.vectors) {
        check(keys.insert({info.kind, info.kind == 1 ? 0 : info.norm,
                           info.kind == 0 ? info.chr : 0, info.resolution}).second,
              "duplicate imported vector key");
        Vector v;
        v.kind = info.kind;
        check(v.kind <= 2, "invalid imported vector kind");
        v.ri = header.resolution(0, info.resolution);
        if (v.kind != 1) {
            check(info.norm < ids.size(), "invalid imported normalization ID");
            v.norm = ids[info.norm];
            if (retained.count({v.norm, v.ri})) continue;
        }
        if (v.kind == 0) {
            check(info.chr < header.chromosomes.size(), "imported chromosome outside header");
            v.chr = info.chr;
            check(info.words == header.bins(v.chr, 0, v.ri), "imported normalization length mismatch");
        }
        for (auto scale : info.scales)
            check(scale.first < header.chromosomes.size(), "imported expected scale outside header");
        v.scales = info.scales;
        v.streamed_values = info.words;
        auto input = std::make_shared<hic10large::VectorFileReader>(info);
        v.loader = [input](uint64_t begin, uint32_t count) { return input->read(begin, count); };
        output.add(v);
    }
    output.finish(header);
}
}
