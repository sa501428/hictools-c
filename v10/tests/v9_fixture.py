"""Independent V9 fixture writer for conversion coverage."""
import struct
import zlib

def p(fmt, *args): return struct.pack('<'+fmt, *args)
def s(value): return value.encode()+b'\0'

def make(path, x_int=False, y_int=False, floating=False, dense=False, frag=False):
    header = bytearray(b'HIC\0'+p('IQ', 9, 0)+s('fixture'))
    nvi_patch = len(header); header += bytes(16)
    header += p('I', 2)+s('unknown')+s('first')+s('unknown')+s('second')
    header += p('I', 1)+s('chr1')+p('Q', 80)+p('II', 1, 10)
    header += p('I', int(frag))
    if frag: header += p('I', 1)+p('I', 7)+b''.join(p('I', i) for i in range(5, 66, 10))
    data = header
    values = [1.25, -0.0, 2.5] if floating else [1, 2, 300]
    if dense:
        # 2 rows x 3 columns: cells 0, 2 and 4 are present.
        raw = p('iii4Bih', 3, 0, 0, int(floating), int(x_int), int(y_int), 2, 6, 3)
        raw += b''.join(p('f' if floating else 'h', v) for v in [values[0], float('nan') if floating else -32768, values[1], float('nan') if floating else -32768, values[2], float('nan') if floating else -32768])
        expected = [(0, 0, values[0]), (2, 0, values[1]), (1, 1, values[2])]
        # Use trans rather than cis to allow x>y in dense fixture.
    else:
        raw = p('iii4B', 3, 0, 0, int(floating), int(x_int), int(y_int), 1)
        raw += p('i' if y_int else 'h', 3)
        for i, v in enumerate(values):
            raw += p('i' if y_int else 'h', i)
            raw += p('i' if x_int else 'h', 1)
            raw += p('i' if x_int else 'h', i)
            raw += p('f' if floating else 'h', v)
        expected = [(i, i, v) for i, v in enumerate(values)]
    # Dense fixture uses a cis-valid 1-column block instead.
    if dense:
        raw = p('iii4Bih', 3, 0, 0, int(floating), int(x_int), int(y_int), 2, 6, 1)
        raw += b''.join(p('f' if floating else 'h', v) for v in [values[0], float('nan') if floating else -32768, values[1], float('nan') if floating else -32768, values[2], float('nan') if floating else -32768])
        expected = [(0, 0, values[0]), (0, 2, values[1]), (0, 4, values[2])]
    compressed = zlib.compress(raw)
    units = [('BP', 10)]+([('FRAG', 1)] if frag else [])
    matrix_pos = len(data)
    meta = bytearray(p('iii', 0, 0, len(units)))
    patches = []
    for unit, bin in units:
        meta += s(unit)+p('i4f4i', 0, 0, 0, 0, 0, bin, 8, 1, 1)
        patches.append(len(meta)+4)
        meta += p('iQi', 0, 0, len(compressed))
    for patch in patches: struct.pack_into('<Q', meta, patch, matrix_pos+len(meta))
    data += meta+compressed
    footer_pos = len(data)
    # Empty expected and normalized-expected lists are legal in V9.
    footer = p('I', 1)+s('0_0')+p('Qi', matrix_pos, len(meta))+p('I', 0)
    data += p('Q', len(footer))+footer+p('I', 0)
    nvi_pos = len(data)
    words = [0x80000000, 0x7fc01234, 0x3f800000, 0x40000000, 0x3f000000, 0x7f800000, 0, 0x40400000]
    nvi = p('I', 1)+s('VC')+p('I', 0)+s('BP')+p('IQQ', 10, 0, 8+4*len(words))
    nvi = bytearray(nvi); struct.pack_into('<Q', nvi, len(nvi)-16, nvi_pos+len(nvi))
    data += nvi+p('Q', len(words))+b''.join(p('I', v) for v in words)
    struct.pack_into('<Q', data, 8, footer_pos)
    struct.pack_into('<QQ', data, nvi_patch, nvi_pos, len(nvi))
    path.write_bytes(data)
    return [(x, y, struct.unpack('<I', p('f', v))[0] if floating else v) for x, y, v in expected], words
