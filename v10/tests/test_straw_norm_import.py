#!/usr/bin/env python3
"""Straw v9/v10 exact exports imported additively by both v10 executables."""
from pathlib import Path
import struct
import sys
import tempfile
sys.dont_write_bytecode = True
from inspect_v10 import Hic
from test_additive_norm import run, preserve
from v9_fixture import make, p, s

WEIRD = 'CUSTOM #/ "quoted"'


def add_norms(path):
    data = bytearray(path.read_bytes())
    patch = 16+len('fixture')+1
    original, nvi_len = struct.unpack_from('<QQ', data, patch)
    # The base VC is retained; multiple resolutions/units are represented by
    # BP10/BP20 and FRAG1. Every custom norm is unknown to built-in norm algorithms.
    nvi = len(data)
    entries = [('RU', 'BP', 10, [0x3f800001, 0x40000000, 0x40400000, 0x40800000,
                             0x40a00000, 0x40c00000, 0x40e00000, 0x41000000, 0x7fc01234]),
               ('RU', 'FRAG', 1, [0x40000000]*9),
               ('NDSCALE', 'BP', 10, [0x40000000]*4),
               ('NDSCALE', 'FRAG', 1, [0x40400000]*8),
               (WEIRD, 'BP', 10, [0x80000000, 0x7f801234, 0x7fc05678, 0xffc05678,
                              0x7f800000, 0xff800000, 0, 0x3f800001]),
               ('EMPTY', 'BP', 10, []),
               ('RU', 'BP', 20, [0x40000000]*5),
               ('NDSCALE', 'BP', 20, [0x40400000]*3)]
    index = bytearray(p('I', len(entries)+1))
    index += data[original+4:original+nvi_len]
    patches = []
    for norm, unit, bin, words in entries:
        index += s(norm)+p('I', 0)+s(unit)+p('I', bin)
        patches.append(len(index)); index += p('QQ', 0, 8+len(words)*4)
    payloads = bytearray()
    for patch_vector, (_, _, _, words) in zip(patches, entries):
        struct.pack_into('<Q', index, patch_vector, nvi+len(index)+len(payloads))
        payloads += p('Q', len(words))+b''.join(p('I', w) for w in words)
    data += index+payloads
    struct.pack_into('<QQ', data, patch, nvi, len(index))
    path.write_bytes(data)
    return {(norm, unit, bin): words for norm, unit, bin, words in entries}


def main():
    straw, regular, large = map(Path, sys.argv[1:4])
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        source = root/'source.v9.hic'; make(source, frag=True, bp_resolutions=(10, 20))
        target = root/'target.v10.hic'; run([regular, 'convert', '--materialize', '20', source, target])
        words = add_norms(source)
        original_bytes = source.read_bytes()
        dump = root/'norms'; run([straw, 'dump-norms', source, '--output-dir', dump])
        paths = list(dump.glob('*.norm.txt')); assert len(paths) == 5
        assert source.read_bytes() == original_bytes
        outputs = []
        for tool in (regular, large):
            path = root/f'{tool.name}.hic'; path.write_bytes(target.read_bytes())
            # A preexisting computed type and v9 VC must both survive every
            # custom import. Large import operates independently of build data.
            run([regular, 'addnorm', '--norm', 'SCALE', path])
            before = Hic(path)
            for text in sorted(paths):
                run([tool, 'addnorm', '--norm-file', text, '--tmp', root, path])
            after = Hic(path); preserve(before, after, allow_new_attributes=True)
            assert set(after.norms) == {'VC', 'SCALE', 'RU', 'NDSCALE', WEIRD, 'EMPTY'}
            for (name, unit, bin), original in words.items():
                u = int(unit == 'FRAG')
                count = 8 if bin in (1, 10) else 4
                want = (original+[0x7fc00000]*count)[:count]
                assert after.vectors[0, name, 0, u, bin][0] == want
                assert len(after.vectors[2, name, None, u, bin][0]) == count
                if len(original) != count:
                    ri = int(bin == 20); norm = after.norms.index(name)
                    key = f'hictools.import.vector.0.{norm}.0.{u}.{ri}'
                    value = str(len(original))+':'+''.join(f'{w:08x}' for w in original[count:])
                    assert (key, value) in after.attributes
            outputs.append(after)
            saved = path.read_bytes()
            for text in paths: run([tool, 'addnorm', '--norm-file', text, path])
            assert path.read_bytes() == saved
            # A v10 exact dump also feeds addnorm without converting through
            # double or losing signaling/quiet NaN payloads.
            v10dump = root/f'{tool.name}-dump'
            run([straw, 'dump-norms', path, '--norm', WEIRD, '--norm', 'RU', '--output-dir', v10dump])
            clean = root/f'{tool.name}-again.hic'; clean.write_bytes(target.read_bytes())
            for text in v10dump.glob('*.norm.txt'): run([tool, 'addnorm', '--norm-file', text, clean])
            imported = Hic(clean)
            for key, value in after.vectors.items():
                if key[0] == 0 and key[1] in ('RU', WEIRD): assert imported.vectors[key] == value
            # Source geometry is explicitly opt-in and validated. A malformed
            # export never modifies the target, even after valid staged blocks.
            bad = root/'bad.txt'
            text = 'HIC_NORM_VECTORS 1\nvector NEW chr1 BP 10\n'
            failures = ['source-length 9\n'+'1\n'*8+'end\n',
                        'source-length -1\nend\n',
                        'source-length 0\nsource-length 0\nend\n',
                        '1\nsource-length 8\n'+'1\n'*7+'end\n',
                        'source-length 18446744073709551616\nend\n',
                        'source-length 0 extra\nend\n',
                        'source-length 0\n1\nend\n']
            for invalid in failures:
                bad.write_text(text+invalid)
                run([tool, 'addnorm', '--norm-file', bad, path], ok=False)
                assert path.read_bytes() == saved
        assert outputs[0].vectors == outputs[1].vectors
    print('Straw normalization import round trip passed')

if __name__ == '__main__': main()
