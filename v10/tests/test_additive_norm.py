#!/usr/bin/env python3
"""Independent regression for additive writers, exact preservation, and failure safety."""
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
sys.dont_write_bytecode = True
from inspect_v10 import Hic


def run(args, ok=True):
    result = subprocess.run(list(map(str, args)), capture_output=True, text=True)
    assert (result.returncode == 0) == ok, (args, result.stdout, result.stderr)
    return result


def vector_objects(h):
    """Return all directory records/chunks and absolute-pointer fields."""
    objects, fields = {}, {}
    def add(pos, size):
        assert pos not in objects
        objects[pos] = size
        fields[pos] = []
    footer, flen = struct.unpack_from('<QQ', h.data, 16)
    add(footer, flen)
    n = struct.unpack_from('<I', h.data, footer + 16)[0]
    for i in range(n):
        at = footer + 24 + 24*i
        pos, size = struct.unpack_from('<QQ', h.data, at + 8)
        add(pos, size); fields[footer].append(at + 8 - footer)
        nr = struct.unpack_from('<I', h.data, pos + 16)[0]
        for ri in range(nr):
            pointer = pos + 24 + 76*ri + 52
            ip, il = struct.unpack_from('<QQ', h.data, pointer)
            if not il: continue
            add(ip, il); fields[pos].append(pointer - pos)
            nb = struct.unpack_from('<I', h.data, ip + 16)[0]
            for j in range(nb):
                pointer = ip + 24 + 16*j + 8
                bp = struct.unpack_from('<Q', h.data, pointer)[0]
                bl = struct.unpack_from('<I', h.data, pointer - 4)[0]
                add(bp, bl); fields[ip].append(pointer - ip)
    for kind, (pos, size) in enumerate(h.vector_locs):
        if not size: continue
        add(pos, size)
        n = struct.unpack_from('<I', h.data, pos + 8)[0]
        at = pos + 16
        for _ in range(n):
            length = struct.unpack_from('<I', h.data, at)[0]
            base = at + 4 + (4 if kind != 1 else 0) + (4 if kind == 0 else 0)
            nc = struct.unpack_from('<I', h.data, base + 24)[0]
            desc = base + 28
            if kind:
                ns = struct.unpack_from('<I', h.data, desc)[0]
                desc += 8 + 8*ns
            for j in range(nc):
                pointer = desc + 32*j + 16
                vp, vl = struct.unpack_from('<QI', h.data, pointer)
                add(vp, vl); fields[pos].append(pointer - pos)
            at += length
    return objects, fields


def reorder(path):
    """Reverse every physical record: footer first, vectors before matrix data."""
    h = Hic(path)
    header = struct.unpack_from('<Q', h.data, 8)[0]
    objects, fields = vector_objects(h)
    output = bytearray(h.data[:header]); relocated = {}
    for pos in sorted(objects, reverse=True):
        relocated[pos] = len(output)
        output += h.data[pos:pos + objects[pos]]
    for pos, offsets in fields.items():
        for offset in offsets:
            old = struct.unpack_from('<Q', h.data, pos + offset)[0]
            struct.pack_into('<Q', output, relocated[pos] + offset, relocated[old])
    for offset in (16, 32, 48, 64):
        old = struct.unpack_from('<Q', h.data, offset)[0]
        if old: struct.pack_into('<Q', output, offset, relocated[old])
    path.write_bytes(output)
    assert Hic(path).vectors == h.vectors


def snapshot(h):
    objects, _ = vector_objects(h)
    chunks = sorted(h.data[pos:pos + size] for pos, size in objects.items()
                    if h.data[pos:pos+4] == b'H10V')
    blocks = sorted(h.data[pos:pos + size] for pos, size in objects.items()
                    if h.data[pos:pos+4] == b'H10B')
    records = {key: h.records(key[0], key[1], key[3], key[2]) for key in h.matrices}
    return chunks, blocks, records


def preserve(before, after):
    assert after.norms[:len(before.norms)] == before.norms
    for key, value in before.vectors.items(): assert after.vectors[key] == value, key
    old_chunks, old_blocks, old_records = snapshot(before)
    chunks, blocks, records = snapshot(after)
    assert all(chunk in chunks for chunk in old_chunks)
    assert blocks == old_blocks and records == old_records
    assert before.attributes == after.attributes
    assert before.chroms == after.chroms and before.res == after.res


def custom_manifest(root, h):
    vectors = []
    for ri, resolution in enumerate(h.res[0]):
        bin = resolution[0]
        count = (h.chroms[0][1] + bin - 1)//bin
        words = ([0x3f800000, 0x80000000, 0x7fc01234, 0x7f800000, 0xff800000] * ((count+4)//5))[:count]
        for kind in (0, 2):
            chr = 0 if kind == 0 else 0xffffffff
            payload = struct.pack('<' + 'I'*count, *words)
            checksum = 14695981039346656037
            for byte in payload: checksum = ((checksum ^ byte)*1099511628211) & ((1<<64)-1)
            path = root/f'custom-{kind}-{ri}.h10w'
            header = bytearray(64)
            struct.pack_into('<4s7I2Q', header, 0, b'H10W', 1, kind, 0, chr, ri, bin, 4, count, checksum)
            path.write_bytes(header + payload)
            scales = '' if kind == 0 else ' 0 0000000080000000'
            vectors.append(f'vector {kind} 0 {chr} {ri} {bin} {count} {checksum:016x} "{path}" {int(kind==2)}{scales}')
    manifest = root/'custom.manifest'
    manifest.write_text('HIC_V10_LARGE_VECTORS 1\nsource 0000000000000000\nnorms 1 "CUSTOM"\n'
                        f'vectors {len(vectors)}\n' + '\n'.join(vectors) + '\nend\n')
    return manifest


def main():
    regular, large = map(Path, sys.argv[1:3])
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        chrom = root/'chrom.sizes'; chrom.write_text('chr1\t700010\nchr2\t30\n')
        pairs = root/'pairs.txt'
        pairs.write_text(''.join(f'chr1 {i*10} chr1 {(i+d)*10}\n' for i in range(600) for d in range(4))
                         + 'chr1 0 chr2 10\nchr2 0 chr2 10\n')
        raw = root/'raw.hic'
        run([regular, 'pre', '-r', '10,20', pairs, raw, chrom])
        for first, second in ((regular, large), (large, regular)):
            path = root/f'{first.name}.hic'; path.write_bytes(raw.read_bytes())
            path.chmod(0o640)
            before = Hic(path)
            run([first, 'add-norm', '--norm', 'VC', path])
            vc = Hic(path); preserve(before, vc); assert vc.norms == ['VC']
            reorder(path)
            before = Hic(path)
            opts = ['--memory', '32KiB', '--fan-in', '2', '--tmp', root] if second == large else []
            run([second, 'addnorm', '--norm', 'VC_SQRT', *opts, path])
            after = Hic(path); preserve(before, after)
            assert after.norms == ['VC', 'VC_SQRT']
            before = after
            run([second, 'addnorm', '--norm', 'SCALE', '-t', '2', *opts, path])
            after = Hic(path); preserve(before, after)
            assert after.norms == ['VC', 'VC_SQRT', 'SCALE']
            assert any(key[1] == 'SCALE' for key in after.vectors)
            size = path.stat().st_size
            before_bytes = path.read_bytes()
            run([regular, 'addnorm', '--norm', 'VC', path])
            preserve(after, Hic(path)); assert path.stat().st_size == size
            assert path.read_bytes() == before_bytes
            assert path.stat().st_mode & 0o777 == 0o640
            manifest = custom_manifest(root, after)
            run([second, 'addnorm', '--vectors', manifest, path])
            custom = Hic(path); preserve(after, custom)
            assert custom.norms == ['VC', 'VC_SQRT', 'SCALE', 'CUSTOM']
            for bin in (10, 20):
                assert custom.vectors[0, 'CUSTOM', 0, 0, bin][0][:5] == [0x3f800000, 0x80000000, 0x7fc01234, 0x7f800000, 0xff800000]
                assert custom.vectors[2, 'CUSTOM', None, 0, bin][1] == [(0, 0x80000000)]
            # Recompute a known type after importing an unknown type.
            run([regular, 'addnorm', '--norm', 'VC,VC_SQRT', path])
            preserve(custom, Hic(path))
            # Checksum failures occur after staging but must leave input untouched.
            bad = root/'bad.hic'; bad.write_bytes(raw.read_bytes())
            sidecar = root/'custom-0-0.h10w'
            damaged = bytearray(sidecar.read_bytes()); damaged[-1] ^= 1; sidecar.write_bytes(damaged)
            for tool in (regular, large):
                before_bytes = bad.read_bytes()
                run([tool, 'addnorm', '--vectors', manifest, bad], ok=False)
                assert bad.read_bytes() == before_bytes
                assert not list(root.glob('*.addnorm.*'))
                run([tool, 'addnorm', '--norm', 'UNKNOWN', bad], ok=False)
                assert bad.read_bytes() == before_bytes
            assert not list(root.glob('hic-v10-addnorm-*'))
        # Unknown imported types remain intact while another algorithm is added.
        path = root/'custom-first.hic'; path.write_bytes(raw.read_bytes())
        manifest = custom_manifest(root, Hic(path))
        run([regular, 'addnorm', '--vectors', manifest, path])
        before = Hic(path)
        run([regular, 'addnorm', '--norm', 'VC', path])
        after = Hic(path); preserve(before, after)
        run([large, 'addnorm', '--norm', 'VC_SQRT', '--memory', '32KiB',
             '--fan-in', '2', '--tmp', root, path])
        preserve(after, Hic(path))
        # Failing during disk-backed extraction also cleans the workspace.
        scored = root/'scores.hic'
        run([regular, 'pre', '--scores', '-r', '10', pairs, scored, chrom])
        before_bytes = scored.read_bytes()
        run([large, 'addnorm', '--norm', 'VC', '--tmp', root, scored], ok=False)
        assert scored.read_bytes() == before_bytes
        assert not list(root.glob('hic-v10-addnorm-*'))
    print('Additive V10 normalization: both tools, exact vectors/chunks/blocks, custom types, physical ordering, atomic failures passed')

if __name__ == '__main__': main()
