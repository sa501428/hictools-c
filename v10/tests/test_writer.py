#!/usr/bin/env python3
import gzip
import math
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
sys.dont_write_bytecode = True
from inspect_v10 import Hic, Cursor
from v9_fixture import make as v9_fixture

def run(args, ok=True):
    p = subprocess.run(list(map(str, args)), capture_output=True, text=True)
    assert (p.returncode == 0) == ok, (args, p.stdout, p.stderr)
    return p.stdout

def main():
    v10, v9, addnorm = sys.argv[1:4]
    straw = sys.argv[4] if len(sys.argv) > 4 else None
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp)
        chrom = p/'chrom.sizes'; chrom.write_text('chr1\t10000\nchr2\t7500\n')
        rows = [('chr1', 100, 'chr1', 250), ('chr1', 100, 'chr1', 250),
                ('chr1', 300, 'chr1', 450), ('chr1', 9999, 'chr1', 9999),
                ('chr1', 150, 'chr2', 300), ('chr2', 50, 'chr2', 250)]
        layouts = {
            'extra.txt': '\n'.join(f'{a} {x} {b} {y}' for a, x, b, y in rows),
            'short.mnd': '\n'.join(f'0 {a} {x} 0 1 {b} {y} 1' for a, x, b, y in rows),
            'long.mnd': '\n'.join(f'0 {a} {x} 0 1 {b} {y} 1 60 50M A 60 50M A r1 r2' for a, x, b, y in rows),
            'input.pairs': '#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 mapq1 mapq2\n' + '\n'.join(f'r {a} {x} {b} {y} + - 60 60' for a, x, b, y in rows),
        }
        expected100 = [(1, 2, 2), (3, 4, 1), (99, 99, 1)]
        for name, text in layouts.items():
            path = p/name; path.write_text(text+'\n')
            out = p/(name+'.hic')
            run([v10, 'pre', '-r', '100,200', '-q', '30', '--derive', '200:100', path, out, chrom])
            h = Hic(out)
            assert h.norms == []
            assert h.vector_locs[0] == (0, 0) and h.vector_locs[2] == (0, 0)
            assert h.vector_locs[1][0] and h.vector_locs[1][1]
            assert h.records(0, 0, 100) == expected100
            assert h.records(0, 0, 200) == [(0, 1, 2), (1, 2, 1), (49, 49, 1)]
            assert h.records(0, 1, 100) == [(1, 3, 1)]
            assert h.records(1, 1, 100) == [(0, 2, 1)]
            assert h.vectors[1, None, None, 0, 100][0]
            if straw:
                assert '100\t200\t2' in run([straw, 'observed', 'NONE', out, 'chr1', 'chr1', 'BP', 100])
                assert '300\t100\t1' in run([straw, 'observed', 'NONE', out, 'chr2', 'chr1', 'BP', 100])
        # V10 addnorm computes vectors for both physical and derived resolutions.
        # It also rebuilds raw expected and normalized expected indexes, without
        # touching the matrix section. A second run replaces rather than grows it.
        normalized = p/'input.pairs.hic'
        before = Hic(normalized)
        before_records = {res: before.records(0, 0, res) for res in (100, 200)}
        run([v10, 'addnorm', '--no-scale', normalized])
        normalized_size = normalized.stat().st_size
        h = Hic(normalized)
        assert h.norms == ['VC', 'VC_SQRT']
        assert all(position and length for position, length in h.vector_locs)
        def float_words(words):
            return [struct.unpack('<f', struct.pack('<I', word))[0] for word in words]
        def scaled_vc(records, bins):
            vc = [0.0] * bins
            for x, y, value in records:
                vc[x] += value
                if x != y: vc[y] += value
            raw = normed = 0.0
            for x, y, value in records:
                if vc[x] <= 0 or vc[y] <= 0: continue
                multiple = 1 if x == y else 2
                raw += multiple * value
                normed += multiple * value / (vc[x] * vc[y])
            factor = math.sqrt(normed / raw)
            return [value * factor for value in vc]
        for res, bins in ((100, 100), (200, 50)):
            assert h.records(0, 0, res) == before_records[res]
            expected_vc = scaled_vc(before_records[res], bins)
            actual_vc = float_words(h.vectors[0, 'VC', 0, 0, res][0])
            assert len(actual_vc) == bins
            for actual, expected in zip(actual_vc, expected_vc):
                assert math.isclose(actual, expected, rel_tol=1e-6, abs_tol=1e-8)
            assert len(h.vectors[0, 'VC_SQRT', 0, 0, res][0]) == bins
            assert len(h.vectors[1, None, None, 0, res][0]) == bins
            assert len(h.vectors[2, 'VC', None, 0, res][0]) == bins
            assert len(h.vectors[2, 'VC_SQRT', None, 0, res][0]) == bins
            assert len(h.vectors[1, None, None, 0, res][1]) == 2
            assert len(h.vectors[2, 'VC', None, 0, res][1]) == 2
        run([v10, 'addnorm', '--no-scale', normalized])
        assert normalized.stat().st_size == normalized_size
        Hic(normalized)
        run([v10, 'addnorm', '--no-vc-sqrt', '--no-scale', normalized])
        h = Hic(normalized)
        assert h.norms == ['VC']
        assert h.records(0, 0, 200) == before_records[200]
        run([v10, 'addnorm', '--no-scale', normalized])
        assert Hic(normalized).norms == ['VC', 'VC_SQRT']
        raw_only = p/'raw-only.hic'
        raw_only.write_bytes(normalized.read_bytes())
        run([v10, 'addnorm', '--no-vc', '--no-vc-sqrt', '--no-scale', raw_only])
        raw = Hic(raw_only)
        assert raw.norms == []
        assert raw.vector_locs[0] == (0, 0) and raw.vector_locs[2] == (0, 0)
        assert raw.vector_locs[1][0] and raw.vector_locs[1][1]
        if straw:
            assert run([straw, 'observed', 'VC', normalized, 'chr1', 'chr1', 'BP', 200])
            assert run([straw, 'expected', 'VC', normalized, 'chr1', 'chr1', 'BP', 200])
            assert run([straw, 'oe', 'VC_SQRT', normalized, 'chr1', 'chr1', 'BP', 200])
        scale_input = p/'scale.txt'
        scale_input.write_text(''.join(
            f'chr1 {i*100} chr1 {(i+d)*100}\n'
            for i in range(50) for d in range(4)))
        scale_hic = p/'scale.hic'
        run([v10, 'pre', '-r', '100,200', '--derive', '200:100',
             scale_input, scale_hic, chrom])
        run([v10, 'addnorm', '--no-vc', '--no-vc-sqrt', '-t', '2', scale_hic])
        h = Hic(scale_hic)
        assert h.norms == ['SCALE']
        for res, bins in ((100, 100), (200, 50)):
            assert len(h.vectors[0, 'SCALE', 0, 0, res][0]) == bins
            assert len(h.vectors[1, None, None, 0, res][0]) == bins
            assert len(h.vectors[2, 'SCALE', None, 0, res][0]) == bins
        gz = p/'input.pairs.gz'
        with gzip.open(gz, 'wt') as f: f.write(layouts['input.pairs']+'\n')
        run([v10, 'pre', '-r', '100', gz, p/'gzip.hic', chrom])
        assert Hic(p/'gzip.hic').records(0, 0, 100) == expected100
        # New V10 counter never passes accumulated integer counts through float.
        weights = p/'weights.txt'; weights.write_text('chr1 100 chr1 250 16777216\nchr1 100 chr1 250 1\n')
        run([v10, 'pre', '-r', '100', weights, p/'counts.hic', chrom])
        assert Hic(p/'counts.hic').records(0, 0, 100) == [(1, 2, 16777217)]
        # Independent HBS fixture: uint16 IDs, bins, escaped exact uint64 counts.
        hbs = p/'counts.hbs.gz'
        def hbs_bytes(records):
            data = b'HICBS\0\r\n' + struct.pack('<HHII', 1, 0, 100, 2)
            for name, length in [(b'chr2', 7500), (b'chr1', 10000)]:
                data += struct.pack('<H', len(name)) + name + struct.pack('<Q', length)
            for a, x, b, y, count in records:
                data += struct.pack('<HIHIH', a, x, b, y, min(count, 65535))
                if count >= 65535: data += struct.pack('<Q', count)
            return data
        records = [(1, 1, 1, 2, (1 << 53)+1), (1, 1, 1, 2, 1),
                   (1, 3, 1, 4, 65534), (1, 5, 1, 6, 65535),
                   (1, 0, 0, 1, 2)]
        hbs.write_bytes(gzip.compress(hbs_bytes(records)))
        run([v10, 'pre', '-r', '100,200', hbs, p/'hbs.hic', chrom])
        assert Hic(p/'hbs.hic').records(0, 0, 100) == [(1, 2, (1<<53)+2), (3, 4, 65534), (5, 6, 65535)]
        assert Hic(p/'hbs.hic').records(0, 1, 100) == [(0, 1, 2)]
        run([v9, '-r', '100,200', hbs, p/'hbs.v9.hic', chrom])
        run([v10, 'convert', p/'hbs.v9.hic', p/'hbs-converted.hic'])
        assert Hic(p/'hbs-converted.hic').records(1, 1, 100) == [(1, 2, 1<<53), (3, 4, 65534), (5, 6, 65535)]
        # Explicit score storage remains a requested, lossy conversion.
        run([v10, 'pre', '--scores', '-f', 'hbs', '-r', '100', hbs, p/'hbs-score.hic', chrom])
        assert Hic(p/'hbs-score.hic').records(0, 0, 100)[0][2] == struct.unpack('<I', struct.pack('<f', (1<<53)+2))[0]
        for cmd in [[v9], [v10, 'pre']]:
            run(cmd + ['-r', '50', hbs, p/'bad.hic', chrom], ok=False)
            run(cmd + ['-r', '150', hbs, p/'bad.hic', chrom], ok=False)
        hbs.write_bytes(gzip.compress(hbs_bytes([(1, 1, 1, 2, (1<<64)-1)])))
        run([v10, 'pre', '-r', '100', hbs, p/'hbs-max.hic', chrom])
        assert Hic(p/'hbs-max.hic').records(0, 0, 100) == [(1, 2, (1<<64)-1)]
        original_hbs_output = (p/'hbs-max.hic').read_bytes()
        hbs.write_bytes(gzip.compress(hbs_bytes([(1, 1, 1, 2, (1<<64)-1), (1, 1, 1, 2, 1)])))
        run([v10, 'pre', '-r', '100', hbs, p/'hbs-max.hic', chrom], ok=False)
        assert (p/'hbs-max.hic').read_bytes() == original_hbs_output
        hbs.write_bytes(gzip.compress(hbs_bytes(records))[:-4])
        for cmd in [[v9], [v10, 'pre']]:
            run(cmd + ['-r', '100', hbs, p/'hbs-max.hic', chrom], ok=False)
            assert (p/'hbs-max.hic').read_bytes() == original_hbs_output
        weights.write_text('chr1 100 chr1 250 -0.0\nchr1 300 chr1 450 1.25\n')
        run([v10, 'pre', '-r', '100', weights, p/'scores.hic', chrom])
        assert Hic(p/'scores.hic').records(0, 0, 100) == [(1, 2, 0x80000000), (3, 4, 0x3fa00000)]
        # Binary input uses the same parser as V9.
        binary = p/'input.bn'
        binary.write_bytes(b''.join(struct.pack('<iiiif', int(a[-1]), x, int(b[-1]), y, 1) for a, x, b, y in rows))
        run([v10, 'pre', '-r', '100', binary, p/'binary.hic', chrom])
        assert Hic(p/'binary.hic').records(0, 0, 100) == expected100
        original = p/'v9.hic'; converted = p/'converted.hic'
        run([v9, '-r', '100,200', p/'extra.txt', original, chrom])
        assert struct.unpack_from('<I', original.read_bytes(), 4)[0] == 9
        run([addnorm, '--no-scale', original])
        run([v10, 'convert', '--derive', '200:100', original, converted])
        h = Hic(converted)
        assert h.records(1, 1, 100) == expected100
        assert h.records(1, 1, 200) == [(0, 1, 2), (1, 2, 1), (49, 49, 1)]
        assert h.chroms[0][0] == 'ALL' and h.records(0, 0, 1)
        # Compare source norm words, including signed zeros/NaNs, without float conversion.
        data = original.read_bytes(); r = Cursor(data); r.take(16); r.string(); nvi, length = r.unpack('QQ')
        r = Cursor(data[nvi:nvi+length]); entries = []
        for _ in range(r.u32()):
            norm, chr, unit, bin, pos, size = r.string(), r.u32(), r.string(), r.u32(), r.u64(), r.u64()
            entries.append((norm, chr, unit, bin, pos, size))
        for norm, chr, unit, bin, pos, size in entries:
            c = Cursor(data[pos:pos+size]); count = c.u64(); words = list(c.unpack('I'*count))
            new = h.vectors[0, norm, chr, 0, bin][0]
            assert words[:len(new)] == new[:len(words)]
            if len(words) > len(new): assert any(k.startswith('hictools.v9.vector.') for k, _ in h.attributes)
        if straw:
            for norm in ('NONE', 'VC', 'VC_SQRT'):
                for res in (100, 200):
                    left = run([straw, 'observed', norm, original, 'chr1', 'chr1', 'BP', res])
                    right = run([straw, 'observed', norm, converted, 'chr1', 'chr1', 'BP', res])
                    assert sorted(left.splitlines()) == sorted(right.splitlines()), (norm, res, left, right)
                    for matrix_type in ('oe', 'expected'):
                        left = run([straw, matrix_type, norm, original, 'chr1', 'chr1', 'BP', res])
                        right = run([straw, matrix_type, norm, converted, 'chr1', 'chr1', 'BP', res])
                        left = sorted([line.split() for line in left.splitlines()])
                        right = sorted([line.split() for line in right.splitlines()])
                        assert len(left) == len(right)
                        for l, r in zip(left, right):
                            assert l[:2] == r[:2] and math.isclose(float(l[2]), float(r[2]), rel_tol=1e-6)

        for x_int in (False, True):
            for y_int in (False, True):
                for floating in (False, True):
                    for dense in (False, True):
                        source, target = p/'fixture.v9.hic', p/'fixture.v10.hic'
                        expected, words = v9_fixture(source, x_int, y_int, floating, dense, frag=True)
                        run([v10, 'convert', source, target])
                        fixture = Hic(target)
                        assert fixture.attributes[:2] == [('unknown', 'first'), ('unknown', 'second')]
                        assert fixture.records(0, 0, 10) == expected
                        assert fixture.records(0, 0, 1, unit=1) == expected
                        assert fixture.vectors[0, 'VC', 0, 0, 10][0] == words
                        if x_int and y_int and not floating and dense:
                            run([v10, 'addnorm', '--no-scale', target])
                            fixture = Hic(target)
                            assert fixture.vectors[0, 'VC', 0, 1, 1]
                            assert fixture.vectors[0, 'VC_SQRT', 0, 1, 1]
                            assert fixture.vectors[1, None, None, 1, 1]
                            assert fixture.vectors[2, 'VC', None, 1, 1]
                        if straw:
                            assert len(run([straw, 'observed', 'NONE', target, 'chr1', 'chr1', 'FRAG', 1]).splitlines()) == 3
        # Many tiny pages force multiple checkpoint groups.
        bigchrom = p/'big.sizes'; bigchrom.write_text('chr1\t100000\n')
        many = p/'many.txt'; many.write_text(''.join(f'chr1 {i*20} chr1 {i*20+10}\n' for i in range(3000)))
        serial_many = p/'many-serial.hic'
        run([v10, 'pre', '-t', '1', '-r', '10', '--block-bins', '1',
             '--page-bytes', '1024', many, serial_many, bigchrom])
        run([v10, 'pre', '-t', '4', '-r', '10', '--block-bins', '1',
             '--page-bytes', '1024', many, p/'many.hic', bigchrom])
        assert serial_many.read_bytes() == (p/'many.hic').read_bytes()
        many_hic = Hic(p/'many.hic')
        assert len(many_hic.pages) > 64 and len(many_hic.records(0, 0, 10)) == 3000
        run([v10, 'addnorm', '--no-scale', p/'many.hic'])
        many_hic = Hic(p/'many.hic')
        assert many_hic.norms == ['VC', 'VC_SQRT']
        assert len(many_hic.pages) > 64 and len(many_hic.records(0, 0, 10)) == 3000
        if straw:
            assert len(run([straw, 'observed', 'NONE', p/'many.hic', 'chr1:40000:40100', 'chr1:40000:40100', 'BP', 10]).splitlines()) == 5
        # Error handling is transactional and never overwrites an existing file.
        saved = converted.read_bytes()
        for args in [
            ['convert', converted, converted],
            ['convert', converted, p/'no.hic'],
            ['convert', '--derive', '200:75', original, converted],
            ['pre', '-r', '0', p/'extra.txt', converted, chrom],
            ['pre', '--wat', p/'extra.txt', converted, chrom],
        ]: run([v10, *args], ok=False)
        unsorted = p/'unsorted.txt'; unsorted.write_text('chr1 100 chr1 200\nchr1 100 chr2 200\nchr1 200 chr1 400\n')
        run([v10, 'pre', '-r', '100', unsorted, converted, chrom], ok=False)
        assert converted.read_bytes() == saved
        assert not list(p.glob('*.tmp.*'))
        damaged = p/'damaged.hic'; damaged.write_bytes(original.read_bytes()[:-10])
        run([v10, 'convert', damaged, converted], ok=False)
        assert converted.read_bytes() == saved
        v9_saved = original.read_bytes()
        run([v10, 'addnorm', original], ok=False)
        assert original.read_bytes() == v9_saved
    print('V10 writer/addnorm: matrices, derived norms, expected vectors, formats, counts, checkpoints, failure safety passed')
if __name__ == '__main__': main()
