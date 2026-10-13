#include "v10/reader.h"
#include <iostream>
using namespace hic10;
int main(int argc, char **argv) {
    try {
        check(argc == 2, "expected fixture path");
        Reader r(argv[1]);
        r.validate_matrix_storage();
        uint64_t batches = 0, largest = 0;
        r.relocate_matrix_pointers(37, [&](uint64_t at, const Bytes &bytes) {
            check(at >= r.header_length() + 37, "bad shifted offset");
            check(bytes.size() <= 1024 * 1024, "unbounded pointer batch");
            largest = std::max<uint64_t>(largest, bytes.size()); ++batches;
        });
        r.assert_unchanged();
        std::cout << batches << ' ' << largest << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
