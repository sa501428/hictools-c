"""Small independent decoder used by writer tests (not a production reader)."""
import ctypes
import ctypes.util
import math
import struct
from pathlib import Path
Z = ctypes.CDLL(ctypes.util.find_library('zstd') or '/opt/homebrew/lib/libzstd.dylib')
Z.ZSTD_decompress.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_size_t]
Z.ZSTD_decompress.restype = ctypes.c_size_t

def decompress(data, size):
    out = ctypes.create_string_buffer(size)
    n = Z.ZSTD_decompress(out, size, data, len(data))
    assert n == size
    return out.raw

class Cursor:
    def __init__(self, b): self.b, self.at = b, 0
    def take(self, n):
        assert n >= 0 and self.at+n <= len(self.b)
        out = self.b[self.at:self.at+n]; self.at += n
        return out
    def unpack(self, fmt): return struct.unpack('<'+fmt, self.take(struct.calcsize('<'+fmt)))
    def u32(self): return self.unpack('I')[0]
    def u64(self): return self.unpack('Q')[0]
    def string(self):
        end = self.b.index(0, self.at)
        return self.take(end-self.at+1)[:-1].decode()
    def var(self):
        v = 0
        for i in range(10):
            b = self.take(1)[0]
            v |= (b & 127) << (i*7)
            if not b & 128:
                assert i == 0 or b
                assert v < (1 << 64)
                return v
        raise AssertionError('bad varint')
    def done(self): assert self.at == len(self.b), (self.at, len(self.b))

def block(data):
    c = Cursor(data)
    ver, rep, mode, score, flags = c.unpack('5B3x')
    x, y, w, h, n, np, nv = c.unpack('4IQ2I')
    assert ver == 1 and n > 0
    p, v = Cursor(c.take(np)), Cursor(c.take(nv)); c.done()
    if rep == 0:
        positions, last = [], 0
        for i in range(n):
            delta = p.var(); last = last+delta if i else delta
            positions.append(last)
        p.done()
    elif rep == 1 or score:
        positions = [i for i in range(w*h) if p.b[i//8] & (1 << (i%8))]
        assert len(positions) == n
    else: positions = list(range(w*h))
    slots = w*h if rep == 2 else n
    scalar = v.u32 if score else v.var
    if mode == 0: values = [scalar()] * slots
    elif mode == 1:
        default, ne = scalar(), v.var(); ordinals, last = [], 0
        for i in range(ne):
            delta = v.var(); last = last+delta if i else delta; ordinals.append(last)
        values = [default] * slots
        for i in ordinals: values[i] = scalar()
    else: values = [scalar() for _ in range(slots)]
    v.done()
    if rep == 2:
        present = set(positions)
        result = [(x+i%w, y+i//w, a) for i, a in enumerate(values) if (i in present if score else a != 0)]
    else: result = [(x+p%w, y+p//w, a) for p, a in zip(positions, values)]
    assert len(result) == n
    return score, result

class Hic:
    def __init__(self, path):
        self.data = Path(path).read_bytes(); c = Cursor(self.data)
        assert c.take(4) == b'HIC\0' and c.u32() == 10
        headerlen, footer, footerlen, nvi, nvil, evi, evil, nevi, nevil, flags = c.unpack('10Q')
        assert flags == 0
        self.vector_locs = [(nvi, nvil), (evi, evil), (nevi, nevil)]
        self.genome = c.string()
        self.attributes = [(c.string(), c.string()) for _ in range(c.u32())]
        self.chroms = [(c.string(), c.u64()) for _ in range(c.u32())]
        self.res = [[c.unpack('IBBHI') for _ in range(c.u32())] for _ in range(2)]
        if self.res[1]:
            self.sites = [[c.u64() for _ in range(c.u32())] for _ in self.chroms]
        self.norms = [c.string() for _ in range(c.u32())]
        assert c.at == headerlen
        c = Cursor(self.data[footer:footer+footerlen]); assert c.take(4) == b'H10F' and c.u32() == 1 and c.u64() == footerlen
        n, zero = c.unpack('II'); assert zero == 0
        entries = [c.unpack('IIQQ') for _ in range(n)]; c.done()
        self.matrices, self.blocks = {}, []
        for a, b, pos, size in entries:
            c = Cursor(self.data[pos:pos+size]); assert c.take(4) == b'H10M'
            ver, ca, cb, nr, zero = c.unpack('5I'); assert (ver, ca, cb, zero) == (1, a, b, 0)
            assert nr == sum(map(len, self.res))
            for _ in range(nr):
                u, mode, agg, typ, ri, bin, source, grid, total, occupied, sd, pc, B, cols, ip, il, nb, reserved = c.unpack('4B3IB3xQQII2IQQ2I')
                assert (bin, mode, agg, 0, source) == self.res[u][ri]
                assert sd == pc == 0x7fc00000
                assert grid == int(a == b)
                matrix = {'type': typ, 'source': source, 'mode': mode, 'grid': grid,
                          'block_bins': B, 'columns': cols, 'occupied': occupied,
                          'sum': total, 'records': []}
                self.matrices[a, b, u, bin] = matrix
                assert reserved == 0
                if not mode and nb:
                    matrix['records'] = self.read_blocks(ip, il, nb)
                if not mode:
                    assert len(matrix['records']) == occupied
                    if typ == 0: assert sum(v for _, _, v in matrix['records']) == total
            c.done()
        self.vectors = {}
        for kind, pos, length in [(0, nvi, nvil), (1, evi, evil), (2, nevi, nevil)]:
            if not length: continue
            c = Cursor(self.data[pos:pos+length]); assert c.take(4) == [b'NVI0', b'EVI0', b'NEVI'][kind]
            ver, n, zero = c.unpack('III'); assert (ver, zero) == (1, 0)
            for _ in range(n):
                length = c.u32(); e = Cursor(c.take(length-4))
                norm = e.u32() if kind != 1 else None
                chr = e.u32() if kind == 0 else None
                unit, ri, bin, count, nominal, nc = e.unpack('B3xIIQII')
                scales = []
                if kind:
                    ns, zero = e.unpack('II'); assert zero == 0
                    scales = [e.unpack('II') for _ in range(ns)]
                words, next = [], 0
                for _ in range(nc):
                    first, n, transform, codec, flags, pos, stored, raw = e.unpack('QIBBHQII')
                    assert first == next and raw == n*4 and codec == 1 and flags == 0
                    next += n
                    chunk = Cursor(self.data[pos:pos+stored]); assert chunk.take(4) == b'H10V'
                    assert chunk.unpack('BBHII') == (1, transform, 0, raw, n)
                    values = decompress(chunk.take(stored-16), raw)
                    if transform == 1: values = bytes(values[lane*n+i] for i in range(n) for lane in range(4))
                    vals = list(struct.unpack('<'+'I'*n, values))
                    if transform == 2:
                        for i in range(1, n): vals[i] ^= vals[i-1]
                    words += vals
                assert next == count
                e.done()
                self.vectors[kind, self.norms[norm] if norm is not None else None, chr, unit, bin] = (words, scales)
            c.done()
    def read_blocks(self, pos, length, expected_blocks):
        c = Cursor(self.data[pos:pos+length]); assert c.take(4) == b'H10I'
        ver, index_length, n, zero = c.unpack('IQII')
        assert ver == 2 and index_length == length == 24 + n*16
        assert n == expected_blocks and zero == 0
        entries = [c.unpack('IIQ') for _ in range(n)]; c.done()
        assert all(entries[i-1][0] < entries[i][0] for i in range(1, n))
        records = []
        for number, size, pos in entries:
            stored = Cursor(self.data[pos:pos+size]); assert stored.take(4) == b'H10B'
            codec, record_version, flags, raw, embedded_number = stored.unpack('BBHII')
            assert (codec, record_version, flags, embedded_number) == (1, 1, 0, number)
            _, cells = block(decompress(stored.take(size-16), raw)); stored.done()
            records += cells
            self.blocks.append((number, size, pos))
        assert len({(x, y) for x, y, _ in records}) == len(records)
        return sorted(records)
    def records(self, a, b, bin, unit=0):
        m = self.matrices.get((a, b, unit, bin))
        if m is None: return []
        if not m['mode']: return sorted(m['records'])
        src = self.res[unit][m['source']][0]
        counts = {}
        for x, y, value in self.records(a, b, src, unit):
            key = x//(bin//src), y//(bin//src)
            counts[key] = counts.get(key, 0)+value
        return sorted((x, y, v) for (x, y), v in counts.items())
