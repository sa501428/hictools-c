#include "v10/format.h"

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <unistd.h>

int main() {
    // hg38 chr1 and chr2 at 10 bp: the previous fixed 256-bin grid overflowed
    // while assigning a block near the ends of the chromosomes.
    const uint64_t chr1Bins = hic10::ceil_div(248956422, 10);
    const uint64_t chr2Bins = hic10::ceil_div(242193529, 10);
    const uint32_t blockBins = hic10::safe_block_bin_count(chr1Bins, chr2Bins, 256);
    assert(blockBins == 380);

    const uint64_t columns = hic10::ceil_div(chr1Bins, blockBins);
    const uint64_t rows = hic10::ceil_div(chr2Bins, blockBins);
    assert(columns * rows <= uint64_t{UINT32_MAX} + 1);
    assert((rows - 1) * columns + (columns - 1) <= UINT32_MAX);

    // Existing coarser/smaller matrices retain the user-requested block scale.
    assert(hic10::safe_block_bin_count(1000, 750, 256) == 256);

    // Exercise the real writer at the last block in the hg38-scale grid. This
    // was the call site that previously threw "uint32 overflow".
    char output[] = "/tmp/hic-v10-geometry-XXXXXX";
    int fd = mkstemp(output);
    assert(fd >= 0);
    close(fd);
    std::remove(output);

    hic10::Header header;
    header.genome = "hg38";
    header.chromosomes = {{"chr1", 248956422, {}}, {"chr2", 242193529, {}}};
    header.resolutions[0] = {{10}};
    {
        hic10::Writer writer(output, header, {});
        writer.matrix(0, 1, [&](uint8_t unit, uint32_t resolution) {
            assert(unit == 0 && resolution == 0);
            hic10::Matrix matrix;
            matrix.cells.push_back(
                {static_cast<uint32_t>(chr1Bins - 1), static_cast<uint32_t>(chr2Bins - 1), 1});
            return matrix;
        });
        writer.finish({});
    }
    std::remove(output);
    return 0;
}
