#!/usr/bin/env python3
"""Many-block repair, sparse offset/index checks and Slurm launcher integration.

The largest logical fixture is 5 GiB with only a few KiB allocated; it is read
only by the streaming-reader probe, never copied to a materialized output.
"""
import ctypes
import importlib.util
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
sys.dont_write_bytecode = True
from inspect_v10 import Hic, Z

Z.ZSTD_compressBound.argtypes = [ctypes.c_size_t]
Z.ZSTD_compressBound.restype = ctypes.c_size_t
Z.ZSTD_compress.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int]
Z.ZSTD_compress.restype = ctypes.c_size_t


def pack(fmt, *values):
    return struct.pack('<' + fmt, *values)


def var(value):
    result = bytearray()
    while value >= 128:
        result.append((value & 127) | 128); value >>= 7
    result.append(value)
    return bytes(result)


def block(number, x, y, value):
    values = var(value)
    raw = pack('5B3x4IQ2I', 1, 0, 0, 0, 0, x, y, 1, 1, 1, 1, len(values)) + b'\0' + values
    out = ctypes.create_string_buffer(Z.ZSTD_compressBound(len(raw)))
    n = Z.ZSTD_compress(out, len(out), raw, len(raw), 1)
    return b'H10B' + pack('BBHII', 1, 1, 0, len(raw), number) + out.raw[:n]


def fixture(path, count=131075):
    header = bytearray(b'HIC\0' + pack('I', 10) + b'\0' * 80)
    header += b'test\0' + pack('I', 0) + pack('I', 2)
    header += b'chrA\0' + pack('Q', 512) + b'chrB\0' + pack('Q', 512)
    header += pack('I', 2) + pack('IBBHI', 1, 0, 1, 0, 0xffffffff) + pack('IBBHI', 1000, 0, 1, 0, 0xffffffff)
    header += pack('II', 0, 0)
    struct.pack_into('<Q', header, 8, len(header))
    value = 100000000001
    blocks, entries = bytearray(), bytearray()
    for number in range(count):
        data = block(number, number % 512, number // 512, value)
        entries += pack('IIQ', number, len(data), len(header) + len(blocks))
        blocks += data
    fine_at = len(header) + len(blocks)
    fine = b'H10I' + pack('IQII', 2, 24 + len(entries), count, 0) + entries
    coarse_at = fine_at + len(fine)
    coarse_block = block(0, 0, 0, count * value)
    coarse_index_at = coarse_at + len(coarse_block)
    coarse = b'H10I' + pack('IQII', 2, 40, 1, 0) + pack('IIQ', 0, len(coarse_block), coarse_at)
    def zoom(ri, bin_size, occupied, ip, il, nb, columns):
        return pack('4B3IB3xQQII2IQQ2I', 0, 0, 1, 0, ri, bin_size, 0xffffffff, 0,
                    count * value, occupied, 0x7fc00000, 0x7fc00000, 1, columns, ip, il, nb, 0)
    matrix = b'H10M' + pack('5I', 1, 0, 1, 2, 0)
    matrix += zoom(0, 1, count, fine_at, len(fine), count, 512)
    matrix += zoom(1, 1000, 1, coarse_index_at, len(coarse), 1, 1)
    matrix_at = coarse_index_at + len(coarse)
    footer_at = matrix_at + len(matrix)
    footer = b'H10F' + pack('IQII', 1, 48, 1, 0) + pack('IIQQ', 0, 1, matrix_at, len(matrix))
    struct.pack_into('<QQ', header, 16, footer_at, len(footer))
    path.write_bytes(header + blocks + fine + coarse_block + coarse + matrix + footer)
    return count * value


def run(*args, ok=True, limited=False):
    def limit():
        if sys.platform.startswith('linux'):
            import resource
            resource.setrlimit(resource.RLIMIT_AS, (256 << 20, 256 << 20))
    result = subprocess.run(list(map(str, args)), text=True, capture_output=True,
                            preexec_fn=limit if limited else None)
    assert (result.returncode == 0) == ok, (args, result.stdout, result.stderr)
    return result


def main():
    exe, probe, launcher = map(lambda x: Path(x).resolve(), sys.argv[1:4])
    spec = importlib.util.spec_from_file_location('slurm_all', launcher)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp); original = p / 'many.hic'; expected = fixture(original)
        info = module.inspect(original, 500, 1000, 4)
        assert info['stored_blocks'] == 131076 and info['source_cells'] == 1
        assert info['source_resolution_bp'] == 1000
        assert info['conservative_ram_estimate_bytes'] < 10 * module.GIB
        repaired = p / 'repaired.hic'
        run(exe, 'add-all', '-t', 4, '--all-source-resolution', 1000, original, repaired, limited=True)
        before, after = Hic(original), Hic(repaired)
        assert after.records(0, 1, 1) == before.records(0, 1, 1)
        assert sum(v for _, _, v in after.records(2, 2, 1)) == expected
        assert 'no per-block relocation map' in run(exe, 'add-all', '-t', 1, repaired, p / 'copy.hic').stderr
        run(exe, 'add-all', '-t', 0, original, p / 'bad-threads.hic', ok=False)
        # A bad fine index is rejected without replacing an existing destination.
        bad = p / 'bad.hic'; data = bytearray(original.read_bytes())
        matrix_at = struct.unpack_from('<Q', data, struct.unpack_from('<Q', data, 16)[0] + 32)[0]
        ip = struct.unpack_from('<Q', data, matrix_at + 24 + 52)[0]
        struct.pack_into('<I', data, ip + 24 + 16, 0) # duplicate block number
        bad.write_bytes(data); protected = p / 'protected.hic'; protected.write_bytes(b'keep me')
        run(exe, 'add-all', '--all-source-resolution', 1000, bad, protected, ok=False)
        assert protected.read_bytes() == b'keep me'
        assert not list(p.glob('protected.hic.add-all.*'))
        # Move metadata/indexes beyond 4 GiB using sparse holes. No full repair/copy.
        tiny = p / 'tiny.hic'; fixture(tiny, 1); data = bytearray(tiny.read_bytes())
        footer = struct.unpack_from('<Q', data, 16)[0]
        matrix = struct.unpack_from('<Q', data, footer + 32)[0]
        shift = 5 * (1 << 30)
        moved_matrix = bytearray(data[matrix:footer])
        index_data = []
        for i in range(2):
            field = 24 + i * 76 + 52
            index_at, index_len = struct.unpack_from('<QQ', moved_matrix, field)
            struct.pack_into('<Q', moved_matrix, field, index_at + shift)
            index_data.append((index_at + shift, data[index_at:index_at + index_len]))
        moved_footer = bytearray(data[footer:])
        struct.pack_into('<Q', moved_footer, 32, matrix + shift)
        header_length = struct.unpack_from('<Q', data, 8)[0]
        header = bytearray(data[:header_length]); struct.pack_into('<Q', header, 16, footer + shift)
        sparse = p / 'sparse-offset.hic'
        with sparse.open('wb') as out:
            out.write(header)
            # Source block positions remain valid below the sparse metadata.
            out.write(data[header_length:matrix])
            for at, content in index_data:
                out.seek(at); out.write(content)
            out.seek(matrix + shift); out.write(moved_matrix); out.write(moved_footer)
        run(probe, sparse, limited=True)
        assert sparse.stat().st_blocks * 512 < 1 << 20
        # >512 MiB logical index: traverse the first batch, then diagnose invalid
        # ordering in the second batch, rather than rejecting its total size.
        huge_index = p / 'sparse-index.hic'
        n = (512 * (1 << 20)) // 16 + 1; index_length = 24 + n * 16
        matrix_bytes = bytearray(data[matrix:footer]); large_at = 1 << 20
        struct.pack_into('<QQI', matrix_bytes, 24 + 52, large_at, index_length, n)
        end_matrix = large_at + index_length
        footer_bytes = bytearray(data[footer:]); struct.pack_into('<Q', footer_bytes, 32, end_matrix)
        header = bytearray(data[:header_length]); struct.pack_into('<Q', header, 16, end_matrix + len(matrix_bytes))
        first_block = before.blocks[0]
        with huge_index.open('wb') as out:
            out.write(header); out.write(data[header_length:matrix])
            out.seek(large_at); out.write(b'H10I' + pack('IQII', 2, index_length, n, 0))
            for i in range(65536): out.write(pack('IIQ', i, first_block[1], first_block[2]))
            out.seek(end_matrix); out.write(matrix_bytes); out.write(footer_bytes)
        error = run(probe, huge_index, ok=False, limited=True)
        assert 'unordered V10 block index' in error.stderr, error.stderr
        # Script generation, shell quoting (paths with spaces), worker and mocked sbatch.
        spaced = p / 'space dir'; spaced.mkdir(); source = spaced / 'tiny input.hic'; fixture(source, 1)
        output = spaced / 'tiny output.hic'; work = spaced / 'jobs'
        args = [sys.executable, launcher, source, output, '--hic-v10', exe, '--work-dir', work,
                '--threads', 2, '--all-source-resolution', 1000]
        run(*args, '--dry-run')
        directory = next(work.glob('add-all-*'))
        job_text = (directory / 'repair.sbatch').read_text()
        assert '--run' in job_text
        run('bash', '-n', directory / 'repair.sbatch')
        run('bash', directory / 'repair.sbatch')
        assert Hic(output).chroms[-1][0] == 'All'
        assert json.loads((directory / 'result.json').read_text())['returncode'] == 0
        assert not Path(str(output) + '.add-all.lock').exists()
        assert module.inspect(source, 500, 1000, 4)['conservative_ram_estimate_bytes'] == info['conservative_ram_estimate_bytes']
        # A queued job must not silently run against a different input.
        stale_work = spaced / 'stale jobs'; stale_output = spaced / 'stale output.hic'
        stale_args = [sys.executable, launcher, source, stale_output, '--hic-v10', exe,
                      '--work-dir', stale_work, '--all-source-resolution', 1000, '--dry-run']
        run(*stale_args)
        old_time = source.stat().st_mtime_ns
        os.utime(source, ns=(source.stat().st_atime_ns, old_time + 1000000000))
        stale_dir = next(stale_work.glob('add-all-*'))
        rejected = run('bash', stale_dir / 'repair.sbatch', ok=False)
        assert 'input changed since submission' in rejected.stderr
        assert not stale_output.exists() and not Path(str(stale_output) + '.add-all.lock').exists()
        # Free-space refusal uses metadata only, without creating a huge file.
        from types import SimpleNamespace
        from unittest.mock import patch
        settings = SimpleNamespace(input=source, output=stale_output, all_bins=500,
            all_source_resolution=1000, threads=4, memory='100G', overwrite=False, copy_gbps=0.5)
        with patch.object(module.shutil, 'disk_usage', return_value=SimpleNamespace(free=0)):
            try:
                module.preflight(settings)
                raise AssertionError('expected disk preflight rejection')
            except ValueError as error:
                assert 'free' in str(error)
        run(*args, '--dry-run', ok=False) # no silent overwrite
        run(*args, '--memory', '1G', '--overwrite', '--dry-run', ok=False)
        mock = p / 'sbatch'; mock.write_text('#!/bin/sh\nprintf "12345;testcluster\\n"\n'); mock.chmod(0o755)
        run(*args, '--overwrite', '--sbatch', mock)
        records = [json.loads(x.read_text()) for x in work.glob('*/submission.json')]
        assert any(x.get('job_id') == '12345' for x in records)
    print('Streamed repair: 131k blocks, u64 counts, sparse >4GiB offsets, >512MiB index streaming, atomic failure and Slurm launcher passed')


if __name__ == '__main__': main()
