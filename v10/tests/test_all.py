#!/usr/bin/env python3
"""Stored All construction, repair, and independent wire-format validation."""
import gzip
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
sys.dont_write_bytecode = True
from inspect_v10 import Hic
from v9_fixture import make as v9_fixture


def run(*args, ok=True):
    p = subprocess.run(list(map(str, args)), capture_output=True, text=True)
    assert (p.returncode == 0) == ok, (args, p.stdout, p.stderr)
    return p


def hbs(path, records):
    names = [(b'chrA', 95), (b'chrB', 75)]
    data = bytearray(b'HICBS\0\r\n') + struct.pack('<HHII', 1, 0, 1, len(names))
    for name, length in names:
        data += struct.pack('<H', len(name)) + name + struct.pack('<Q', length)
    for a, x, b, y, count in records:
        data += struct.pack('<HIHIH', a, x, b, y, min(count, 65535))
        if count >= 65535: data += struct.pack('<Q', count)
    path.write_bytes(gzip.compress(data))


def main():
    regular, large = sys.argv[1:3]
    with tempfile.TemporaryDirectory() as directory:
        p = Path(directory)
        sizes = p/'sizes'; sizes.write_text('chrA\t95\nchrB\t75\n')
        pairs = p/'pairs'
        pairs.write_text('chrA 0 chrA 9\nchrA 14 chrA 25\nchrA 94 chrA 94\nchrA 14 chrB 6\nchrA 94 chrB 74\nchrB 0 chrB 74\n')
        original = p/'original.hic'
        built = p/'built.hic'
        common = ['-r', '1,2,5,10,100', '--all-source-resolution', '10']
        run(regular, 'pre', *common, '--no-all', pairs, original, sizes)
        run(regular, 'addnorm', '--no-scale', original)
        before = Hic(original)
        run(regular, 'pre', *common, pairs, built, sizes)
        direct = Hic(built)
        assert direct.chroms == before.chroms + [('All', 1)]
        assert direct.res == before.res
        assert dict(direct.attributes)['hicAll.chromosomeBins'] == '10,8'
        expected = [(0, 0, 6)]
        assert sorted(direct.records(2, 2, 1)) == sorted(expected)
        for r in [1, 2, 5, 10, 100]:
            assert sum(c for _, _, c in direct.records(2, 2, r)) == 6
            assert (0, 2, 0, r) not in direct.matrices
            assert (1, 2, 0, r) not in direct.matrices
        repaired = p/'repaired.hic'
        run(regular, 'add-all', '--all-source-resolution', '10', original, repaired)
        after = Hic(repaired)
        assert after.chroms == direct.chroms and after.res == before.res
        assert after.vectors == before.vectors
        for key, matrix in before.matrices.items(): assert after.matrices[key] == matrix
        # Complete original compressed block frames are copied byte-for-byte.
        for _, size, pos in before.blocks:
            assert before.data[pos:pos+size] in after.data
        assert sorted(after.records(2, 2, 1)) == sorted(expected)
        snapshot = repaired.read_bytes()
        run(regular, 'add-all', '--in-place', repaired)
        assert repaired.read_bytes() == snapshot
        run(regular, 'add-all', '--all-bins', '0', original, p/'invalid.hic', ok=False)
        assert original.read_bytes() == before.data
        run(regular, 'add-all', '--in-place', '--all-source-resolution', '10', original)
        assert Hic(original).vectors == before.vectors
        # A derived coarse preference must select an actual prepared anchor.
        # Otherwise the direct writer would silently construct an empty All.
        run(regular, 'pre', '-r', '1,2,5', '--all-source-resolution', '5', pairs, p/'derived-source.hic', sizes)
        derived_source = Hic(p/'derived-source.hic')
        assert dict(derived_source.attributes)['hicAll.sourceResolution'] == '1'
        assert sum(v for _, _, v in derived_source.records(2, 2, 1)) == 6
        # Tight bin budget and tiny genomes preserve vector lengths.
        run(regular, 'pre', *common, '--all-bins', '2', pairs, p/'two.hic', sizes)
        two = Hic(p/'two.hic'); assert two.chroms[-1] == ('All', 1)
        assert sum(v for _, _, v in two.records(2, 2, 1)) == 6
        # Coarse-only large genomes must also fit desktop synthetic boundaries.
        big_sizes = p/'big.sizes'; big_sizes.write_text('chrA\t5000000000\nchrB\t4000000000\n')
        big_pairs = p/'big.pairs'; big_pairs.write_text('chrA 0 chrA 0\n')
        run(regular, 'pre', '-f', 'short', '-r', '5000000,10000000', big_pairs, p/'big.hic', big_sizes)
        big = Hic(p/'big.hic')
        assert big.chroms[-1][1] <= (1 << 31) - 1
        assert sum(v for _, _, v in big.records(2, 2, 5000000)) == 1
        # A standalone score + FRAG file also repairs without changing original data.
        legacy = p/'legacy.hic'; v9_fixture(legacy, floating=True, frag=True)
        frag = p/'frag.hic'
        run(regular, 'convert', '--no-all', legacy, frag)
        old = Hic(frag)
        run(regular, 'add-all', '--in-place', frag)
        new = Hic(frag)
        assert new.vectors == old.vectors and new.sites[-1] == []
        for key, matrix in old.matrices.items(): assert new.matrices[key] == matrix
        assert new.records(len(old.chroms), len(old.chroms), 1, unit=1) == []
        # The out-of-core and distributed paths retain exact counts above 2^53.
        count = (1 << 53) + 17
        source = p/'source.hbs.gz'
        hbs(source, [(0, 14, 0, 25, count), (0, 14, 1, 6, 13), (1, 0, 1, 74, 19)])
        stage, build, parts = p/'stage', p/'build', p/'parts'
        run(large, 'stage', source, stage)
        run(large, 'build-cells', '-r', '1,2,5,10,100', '--memory', '1MiB', stage/'stage.manifest', build)
        run(large, 'write', '--all-source-resolution', '10', '--memory', '1MiB', stage/'stage.manifest', build/'build.manifest', p/'large.hic')
        l = Hic(p/'large.hic')
        assert l.chroms[-1] == ('All', 1)
        assert (0, 0, count + 32) in l.records(2, 2, 1)
        for i in range(3):
            run(large, 'write-pair', '--memory', '1MiB', stage/'stage.manifest', build/'build.manifest', parts, i)
        run(large, 'merge-pairs', '--all-source-resolution', '10', stage/'stage.manifest', build/'build.manifest', parts, p/'merged.hic')
        merged = Hic(p/'merged.hic')
        for key in l.matrices: assert merged.matrices[key] == l.matrices[key]
        run(large, 'merge-pairs', '--no-all', stage/'stage.manifest', build/'build.manifest', parts, p/'without.hic')
        assert len(Hic(p/'without.hic').chroms) == 2
        # V9-style kb span, unpadded boundaries, sparse finest coordinates, and
        # replacement of an existing overview without renumbering chromosomes.
        kb_sizes = p/'kb.sizes'; kb_sizes.write_text('chrA\t95000\nchrB\t75000\n')
        kb_pairs = p/'kb.contacts'
        kb_pairs.write_text('chrA 14000 chrA 25000\nchrA 16999 chrA 17000\nchrA 94000 chrB 74000\nchrB 0 chrB 74000\n')
        kb = p/'kb.hic'
        run(regular, 'pre', '-r', '1,10000', '--all-source-resolution', '1', kb_pairs, kb, kb_sizes)
        initial = Hic(kb)
        assert initial.chroms[-1] == ('All', 170)
        assert dict(initial.attributes)['hicAll.version'] == '2'
        assert initial.records(2, 2, 1) == [(14, 25, 1), (16, 17, 1), (94, 169, 1), (95, 169, 1)]
        rebuilt = p/'kb.rebuilt.hic'
        run(regular, 'add-all', '--rebuild-all', '--all-bins', '10', '--all-source-resolution', '1', kb, rebuilt)
        changed = Hic(rebuilt)
        assert changed.chroms == initial.chroms
        assert dict(changed.attributes)['hicAll.storageStride'] == '17'
        assert dict(changed.attributes)['hicAll.coordinateScale'] == '1000'
        assert dict(changed.attributes)['hicAll.binCount'] == '10'
        assert changed.records(2, 2, 1) == [(0, 17, 2), (85, 153, 2)]
        assert len({key[:2] for key in changed.matrices}) == len({key[:2] for key in initial.matrices})
        assert changed.vectors == initial.vectors
        for key, matrix in initial.matrices.items():
            if key[0] != 2: assert changed.matrices[key] == matrix
        # A shorter metadata header is padded safely during full-body shifting.
        run(regular, 'add-all', '--rebuild-all', '--all-bins', '10', '--in-place', kb)
        assert Hic(kb).chroms == initial.chroms
        assert sum(c for _, _, c in Hic(kb).records(2, 2, 1)) == 4
        # Existing convention-one files migrate from their real matrices.
        old_fixture = Path(__file__).resolve().parents[3]/'straw'/'C++'/'tests'/'pseudo-all.hic'
        if old_fixture.exists():
            migrated = p/'migrated.hic'
            run(regular, 'add-all', '--rebuild-all', old_fixture, migrated)
            old, new = Hic(old_fixture), Hic(migrated)
            assert new.chroms[-1] == ('All', 1)
            assert new.vectors == old.vectors
            assert sum(c for _, _, c in new.records(2, 2, 1)) == (1 << 53) + 49
        fallback_sizes = p/'fallback.sizes'
        fallback_sizes.write_text('chrA\t1500000000\nchrB\t1500000000\n')
        fallback_pairs = p/'fallback.contacts'; fallback_pairs.write_text('chrA 0 chrB 0\n')
        run(regular, 'pre', '-r', '1000,2500000', fallback_pairs, p/'fallback.hic', fallback_sizes)
        fallback = Hic(p/'fallback.hic')
        assert dict(fallback.attributes)['hicAll.sourceResolution'] == '2500000'
        assert fallback.chroms[-1] == ('All', 3000000)
        if len(sys.argv) > 4 and sys.argv[4]:
            v9, straw = sys.argv[3:5]
            legacy_kb = p/'kb.v9.hic'
            run(v9, '-r', '1,10000', kb_pairs, legacy_kb, kb_sizes)
            coarse_legacy = p/'coarse.v9.hic'
            run(v9, '-r', '1000,10000', kb_pairs, coarse_legacy, kb_sizes)
            preserved = p/'preserved-legacy.v10.hic'
            run(regular, 'convert', coarse_legacy, preserved)
            assert Hic(preserved).chroms[0][0].lower() == 'all'
            run(regular, 'add-all', '--rebuild-all', '--in-place', preserved)
            migrated_legacy = Hic(preserved)
            assert migrated_legacy.chroms[0][1] == 170
            assert dict(migrated_legacy.attributes)['hicAll.sourceResolution'] == '1000'
            assert sum(c for _, _, c in migrated_legacy.records(0, 0, 1)) == 4
            raw = run(straw, 'observed', 'NONE', legacy_kb, 'ALL', 'ALL', 'BP', '1').stdout
            legacy_records = sorted((int(x), int(y), int(float(v))) for x,y,v in
                (line.split() for line in raw.splitlines() if line.strip()))
            assert sorted(initial.records(2, 2, 1)) == legacy_records
    print('Scaled All: default writers, exact counts, levels, no trans pairs, atomic repair, vectors, FRAG, and distributed assembly passed')

if __name__ == '__main__': main()
