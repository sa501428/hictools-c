#!/usr/bin/env python3
"""Regenerate the cross-reader kilobase All fixture: EXE OUTPUT.hic."""
import gzip
from pathlib import Path
import struct
import subprocess
import sys
import tempfile


def make(executable, output):
    with tempfile.TemporaryDirectory() as directory:
        p = Path(directory)
        sizes = p/'sizes'; sizes.write_text('chrA\t95000\nchrB\t75000\n')
        data = bytearray(b'HICBS\0\r\n') + struct.pack('<HHII', 1, 0, 1, 2)
        for name, length in [(b'chrA', 95000), (b'chrB', 75000)]:
            data += struct.pack('<H', len(name)) + name + struct.pack('<Q', length)
        for a, x, b, y, count in [(0, 14000, 0, 25000, (1 << 53) + 17),
                                  (0, 94000, 1, 74000, 13), (1, 0, 1, 74000, 19)]:
            data += struct.pack('<HIHIH', a, x, b, y, min(count, 65535))
            if count >= 65535: data += struct.pack('<Q', count)
        source = p/'source.hbs.gz'; source.write_bytes(gzip.compress(data))
        subprocess.run([str(executable), 'pre', '-r', '1,2,5,1000,10000',
                        '--all-bins', '10', str(source), str(output), str(sizes)], check=True)


if __name__ == '__main__': make(sys.argv[1], sys.argv[2])
