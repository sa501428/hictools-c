#!/usr/bin/env python3
"""Supplied vectors: independent expected arithmetic and atomic input validation."""
import gzip
import math
from pathlib import Path
import struct
import sys
import tempfile
sys.dont_write_bytecode = True
from inspect_v10 import Hic
from test_additive_norm import run, preserve
from v9_fixture import make as v9_fixture


def f32(value):
    return struct.unpack('<f', struct.pack('<f', value))[0]


def bits(value):
    return struct.unpack('<I', struct.pack('<f', value))[0]


def floats(words):
    return [struct.unpack('<f', struct.pack('<I', word))[0] for word in words]


def source(path):
    names = [(b'chr1', 50), (b'chr2', 30)]
    data = bytearray(b'HICBS\0\r\n') + struct.pack('<HHII', 1, 0, 1, 2)
    for name, length in names:
        data += struct.pack('<H', len(name)) + name + struct.pack('<Q', length)
    for chr, bins in enumerate((5, 3)):
        for y in range(bins):
            for x in range(y + 1):
                data += struct.pack('<HIHIHQ', chr, x*10, chr, y*10, 65535, 100000)
    data += struct.pack('<HIHIHQ', 0, 0, 1, 10, 65535, 100000)
    path.write_bytes(gzip.compress(data))


def block(name, chr, bin, values, unit='BP'):
    return f'vector {name} {chr} {unit} {bin}\n' + '\n'.join(map(str, values)) + '\nend\n'


def validate_expected(h, name, bin, divisors):
    size = max((length+bin-1)//bin for _, length in h.chroms)
    actual = [0.0]*size
    observed = {}
    for chr, norm in divisors.items():
        total = 0.0
        for x, y, value in h.records(chr, chr, bin):
            if h.matrices[chr, chr, 0, bin]['type']:
                value = floats([value])[0]
            if not (norm[x] > 0 and norm[y] > 0 and math.isfinite(norm[x]) and math.isfinite(norm[y])):
                continue
            value /= norm[x]*norm[y]
            actual[y-x] += value
            total += value
        if total: observed[chr] = total
    words, scales = h.vectors[2, name, None, 0, bin]
    expected = []
    support = max((len(divisors[chr]) for chr in observed), default=0)
    for d in range(size):
        possible = sum(max(0, len(divisors[chr])-d) for chr in observed)
        # Every supported distance exceeds the smoothing count threshold;
        # its expected is independently the normalized sum / possible pairs.
        if possible:
            assert actual[d] > 400
            numerator = actual[d]
            # The shared smoother retains a two-distance window at the final
            # supported distance when it cannot advance its right boundary.
            if d == support-1 and d > 0:
                numerator += actual[d-1]
                possible += sum(max(0, len(divisors[chr])-(d-1)) for chr in observed)
            expected.append(f32(numerator/possible))
        else: expected.append(float('nan'))
    for a, b in zip(floats(words), expected):
        assert math.isnan(a) if math.isnan(b) else math.isclose(a, b, rel_tol=2e-6), (a, b)
    expected_scales = {chr: f32(sum((len(divisors[chr])-d)*expected[d]
                                  for d in range(len(divisors[chr])))/total)
                       for chr, total in observed.items()}
    assert dict(scales) == {chr: bits(value) for chr, value in expected_scales.items()}


def main():
    regular, large = map(Path, sys.argv[1:3])
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        hbs = root/'source.hbs.gz'; source(hbs)
        chrom = root/'chrom.sizes'; chrom.write_text('chr1\t50\nchr2\t30\n')
        raw = root/'raw.hic'
        run([regular, 'pre', '-r', '10,20', hbs, raw, chrom])
        text = root/'vectors.txt'
        text.write_text('# user divisors, not inverse weights\nHIC_NORM_VECTORS 1\n' +
                        block('CUSTOM', 'chr1', 10, [1, 2, 3, 4, 5]) +
                        block('CUSTOM', 'chr2', 10, [2, 2, 2]) +
                        block('CUSTOM', 'chr1', 20, [1, 2, 3]) +
                        block('CUSTOM', 'chr2', 20, [2, 4]) +
                        block('SUBSET', 'chr1', 20, [3, 2, 1]) +
                        block('MASKED', 'chr1', 10, ['0', 'nan', 'bits:7fc01234', 'inf', '-0']) +
                        block('EXACT', 'chr1', 10, ['1.25', 'bits:3f800001', '2e-1', '-inf', 'bits:80000000']))
        outputs = []
        for tool in (regular, large):
            path = root/f'{tool.name}.hic'; path.write_bytes(raw.read_bytes())
            run([regular, 'addnorm', '--norm', 'VC', path])
            before = Hic(path)
            opts = ['--memory', '32KiB', '--fan-in', '2'] if tool == large else []
            run([tool, 'addnorm', '--norm-file', text, '--tmp', root, *opts, path])
            after = Hic(path); preserve(before, after)
            assert after.norms == ['VC', 'CUSTOM', 'SUBSET', 'MASKED', 'EXACT']
            assert after.vectors[0, 'CUSTOM', 0, 0, 10][0] == list(map(bits, [1, 2, 3, 4, 5]))
            assert after.vectors[0, 'CUSTOM', 1, 0, 20][0] == list(map(bits, [2, 4]))
            validate_expected(after, 'CUSTOM', 10, {0: [1, 2, 3, 4, 5], 1: [2, 2, 2]})
            validate_expected(after, 'CUSTOM', 20, {0: [1, 2, 3], 1: [2, 4]})
            validate_expected(after, 'SUBSET', 20, {0: [3, 2, 1]})
            assert after.vectors[0, 'MASKED', 0, 0, 10][0] == [0, 0x7fc00000, 0x7fc01234, 0x7f800000, 0x80000000]
            assert after.vectors[2, 'MASKED', None, 0, 10] == ([0x7fc00000]*5, [])
            assert after.vectors[0, 'EXACT', 0, 0, 10][0] == [bits(1.25), 0x3f800001, bits(0.2), 0xff800000, 0x80000000]
            before_bytes = path.read_bytes()
            run([tool, 'add-norm', '--norm-file', text, '--tmp', root, *opts, path])
            assert path.read_bytes() == before_bytes
            outputs.append(after)
            # Missing raw expected is rebuilt at target resolutions.
            missing = root/f'missing-{tool.name}.hic'
            data = bytearray(raw.read_bytes()); struct.pack_into('<QQ', data, 48, 0, 0); missing.write_bytes(data)
            run([tool, 'addnorm', '--norm-file', text, '--tmp', root, *opts, missing])
            built = Hic(missing)
            for bin in (10, 20):
                assert built.vectors[1, None, None, 0, bin] == Hic(raw).vectors[1, None, None, 0, bin]
            bad = root/'bad.txt'
            valid = 'HIC_NORM_VECTORS 1\n' + block('NEW', 'chr1', 10, [1]*5)
            malformed = [valid + block('NEW', 'chr1', 10, [1]*5),
                         'HIC_NORM_VECTORS 2\n',
                         'HIC_NORM_VECTORS 1\n' + block('NEW', 'unknown', 10, [1]*5),
                         'HIC_NORM_VECTORS 1\n' + block('NONE', 'chr1', 10, [1]*5),
                         'HIC_NORM_VECTORS 1\n' + block('NEW', 'chr1', 15, [1]*4),
                         'HIC_NORM_VECTORS 1\n' + block('NEW', 'chr1', 10, [1]*4),
                         'HIC_NORM_VECTORS 1\n' + block('NEW', 'chr1', 10, [1]*6),
                         valid.replace('1\nend', 'broken\nend'), valid[:-4],
                         valid.replace('1\nend', '1e100\nend'),
                         valid.replace('1\nend', 'bits:123\nend')]
            for invalid in malformed:
                bad.write_text(invalid)
                result = run([tool, 'addnorm', '--norm-file', bad, '--tmp', root, *opts, path], ok=False)
                assert path.read_bytes() == before_bytes
                assert not list(root.glob('*.addnorm.*'))
                assert not list(root.glob('hic-v10-norm-input-*'))
                assert not list(root.glob('hic-v10-addnorm-*'))
            run([tool, 'addnorm', '--norm', 'VC', '--norm-file', text, path], ok=False)
            run([tool, 'addnorm', '--vectors', text, '--norm-file', text, path], ok=False)
        assert outputs[0].vectors == outputs[1].vectors
        # FRAG vectors use the file's fragment-bin count. Both count paths work.
        v9, frag = root/'frag.v9.hic', root/'frag.hic'
        v9_fixture(v9, True, True, False, True, frag=True)
        run([regular, 'convert', v9, frag])
        h = Hic(frag); count = len(h.sites[0]) + 1
        text.write_text('HIC_NORM_VECTORS 1\n' + block('FRAG_CUSTOM', 'chr1', 1, [2]*count, 'FRAG'))
        for tool in (regular, large):
            path = root/f'frag-{tool.name}.hic'; path.write_bytes(frag.read_bytes())
            run([tool, 'addnorm', '--norm-file', text, '--tmp', root, path])
            after = Hic(path); preserve(h, after)
            assert after.vectors[0, 'FRAG_CUSTOM', 0, 1, 1][0] == [bits(2)]*count
            assert after.vectors[2, 'FRAG_CUSTOM', None, 1, 1][0]
        # The regular importer also handles score matrices.
        scores = root/'scores.hic'
        run([regular, 'pre', '--scores', '-r', '10', hbs, scores, chrom])
        text.write_text('HIC_NORM_VECTORS 1\n' + block('SCORE_CUSTOM', 'chr1', 10, [1]*5))
        before = Hic(scores)
        run([regular, 'addnorm', '--norm-file', text, '--tmp', root, scores])
        preserve(before, Hic(scores))
        validate_expected(Hic(scores), 'SCORE_CUSTOM', 10, {0: [1]*5})
        failed_score = root/'large-score.hic'; failed_score.write_bytes(before.data)
        run([large, 'addnorm', '--norm-file', text, '--tmp', root, failed_score], ok=False)
        assert failed_score.read_bytes() == before.data
        assert not list(root.glob('*.addnorm.*'))
        assert not list(root.glob('hic-v10-addnorm-*'))
        # A vector spanning multiple 64K chunks, with external-sort spills for
        # derived expected calculation. Neither importer may rescale it.
        chrom.write_text('chr1\t700010\nchr2\t30\n')
        pairs = root/'many.txt'
        pairs.write_text(''.join(f'chr1 {x*10} chr1 {y*10}\n'
                                 for y in range(60) for x in range(y+1)))
        long_raw = root/'long-raw.hic'
        run([regular, 'pre', '-r', '10,20', pairs, long_raw, chrom])
        text.write_text('HIC_NORM_VECTORS 1\n' + block('LONG', 'chr1', 10, [1.25]*70001)
                        + block('LONG', 'chr1', 20, [2]*35001))
        long_vectors = []
        for tool in (regular, large):
            path = root/f'long-{tool.name}.hic'; path.write_bytes(long_raw.read_bytes())
            opts = ['--memory', '32KiB', '--fan-in', '2'] if tool == large else []
            run([tool, 'addnorm', '--norm-file', text, '--tmp', root, *opts, path])
            after = Hic(path); preserve(Hic(long_raw), after)
            assert after.vectors[0, 'LONG', 0, 0, 10][0] == [bits(1.25)]*70001
            assert after.vectors[0, 'LONG', 0, 0, 20][0] == [bits(2)]*35001
            long_vectors.append(after.vectors)
        assert long_vectors[0] == long_vectors[1]
    print('Text vector import: exact divisors, independently checked NEVI/scales, derived/FRAG/score data, masks, preservation, failure safety passed')

if __name__ == '__main__': main()
