#!/usr/bin/env python3
"""Preflight and submit a streamed, parallel V10 All repair to one Slurm node.

All paths must be visible on compute nodes. --dry-run writes a reviewable job
and metadata report without submitting it. No contact blocks are read during
preflight. The worker repeats preflight on its allocated node.
"""
import argparse
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
import time

GIB = 1 << 30


def memory_bytes(value):
    match = re.fullmatch(r'(\d+)([KMGT]?)', value.upper())
    if not match or int(match[1]) <= 0:
        raise argparse.ArgumentTypeError('use a positive Slurm memory size such as 100G or 150G')
    # Slurm --mem with no suffix means MiB.
    return int(match[1]) * (1 << {'': 20, 'K': 10, 'M': 20, 'G': 30, 'T': 40}[match[2]])


class Cursor:
    def __init__(self, data):
        self.data, self.at = data, 0

    def take(self, count):
        if count < 0 or count > len(self.data) - self.at:
            raise ValueError('truncated V10 metadata')
        value = self.data[self.at:self.at + count]
        self.at += count
        return value

    def unpack(self, fmt):
        return struct.unpack('<' + fmt, self.take(struct.calcsize('<' + fmt)))

    def word(self):
        return self.unpack('I')[0]

    def string(self):
        end = self.data.index(0, self.at)
        return self.take(end - self.at + 1)[:-1].decode('utf-8')


def inspect(path, bins, preferred, threads, rebuild=False):
    """Read only header, matrix descriptors, footer entries and vector entry sizes."""
    st = path.stat()
    with path.open('rb') as stream:
        def read(at, length):
            if at < 0 or length < 0 or at + length > st.st_size or length > 512 * (1 << 20):
                raise ValueError('invalid or oversized V10 metadata interval')
            stream.seek(at)
            data = stream.read(length)
            if len(data) != length:
                raise ValueError('short V10 metadata read')
            return data
        prefix = read(0, 88)
        if prefix[:8] != b'HIC\0\x0a\0\0\0':
            raise ValueError('input must be a V10 .hic file')
        header_length, footer_at, footer_length, *locators = struct.unpack_from('<9Q', prefix, 8)
        c = Cursor(read(0, header_length)); c.take(88)
        genome = c.string()
        attributes = [(c.string(), c.string()) for _ in range(c.word())]
        chroms = [(c.string(), c.unpack('Q')[0]) for _ in range(c.word())]
        resolutions = [[c.unpack('IBBHI') for _ in range(c.word())] for _ in range(2)]
        if resolutions[1]:
            for _ in chroms:
                c.take(c.word() * 8)
        norms = [c.string() for _ in range(c.word())]
        if c.at != header_length:
            raise ValueError('V10 header length mismatch')
        real_lengths = [length for name, length in chroms if name.lower() != 'all']
        genome_length = sum(real_lengths)
        coordinate_scale = max(1000, (genome_length + max(real_lengths) - 1) // max(real_lengths),
                               (genome_length + 2147483646) // 2147483647)
        pseudo_length = max(1, genome_length // coordinate_scale)
        finest = resolutions[0][0][0] if resolutions[0] else 1
        display_resolution = max(1, (max(1, pseudo_length // bins) + finest - 1) // finest) * finest
        if (pseudo_length + display_resolution - 1) // display_resolution > bins + 1:
            minimum = (pseudo_length + bins) // (bins + 1)
            display_resolution = ((minimum + finest - 1) // finest) * finest
        source_preference = min(preferred, display_resolution * coordinate_scale)
        if not resolutions[0]:
            raise ValueError('All requires BP resolutions')
        def choose(allowed=None):
            eligible = [r[0] for r in resolutions[0] if not r[1] and (allowed is None or r[0] in allowed)]
            if not eligible: raise ValueError('All has no populated BP source')
            below = [bin for bin in eligible if bin <= source_preference]
            above = next((bin for bin in eligible if source_preference <= bin <= display_resolution * coordinate_scale), None)
            if not below: return above or eligible[0]
            base = below[-1]
            return above if above and above * base < source_preference * source_preference else base
        source = choose()
        f = Cursor(read(footer_at, 24))
        if f.take(4) != b'H10F' or f.word() != 1 or f.unpack('Q')[0] != footer_length:
            raise ValueError('invalid V10 footer')
        pairs = f.word(); f.take(4)
        if footer_length != 24 + pairs * 24:
            raise ValueError('invalid V10 matrix directory length')
        nres = sum(map(len, resolutions))
        blocks = index_bytes = source_cells = source_blocks = 0
        all_ids = [i for i, (name, _) in enumerate(chroms) if name.lower() == 'all']
        all_pair = False
        populated = set()
        for first in range(0, pairs, 4096):
            entries = read(footer_at + 24 + first * 24, min(4096, pairs - first) * 24)
            for a, b, at, size in struct.iter_unpack('<IIQQ', entries):
                if a > b or b >= len(chroms):
                    raise ValueError('invalid V10 matrix chromosome ID')
                all_pair |= a == b and a in all_ids
                if size != 24 + nres * 76:
                    raise ValueError('invalid V10 matrix descriptor length')
                m = Cursor(read(at, size))
                if m.take(4) != b'H10M' or m.unpack('5I') != (1, a, b, nres, 0):
                    raise ValueError('invalid V10 matrix descriptor')
                for _ in range(nres):
                    z = m.unpack('4B3IB3xQQII2IQQ2I')
                    unit, mode, _, _, ri = z[:5]
                    occupied, il, nb = z[9], z[15], z[16]
                    blocks += nb; index_bytes += il
                    if unit == 0 and not mode and occupied and a not in all_ids and b not in all_ids:
                        populated.add(z[5])
        if all_ids and not all_pair and not rebuild:
            raise ValueError('existing All header lacks an All/All matrix')
        source = choose(populated or None)
        source_resolution_pairs = {}
        # Repeat only the small directory/descriptors so per-pair fallbacks are
        # counted accurately without retaining O(pairs * resolutions) Python data.
        for first in range(0, pairs, 4096):
            entries = read(footer_at + 24 + first * 24, min(4096, pairs - first) * 24)
            for a, b, at, size in struct.iter_unpack('<IIQQ', entries):
                if a in all_ids or b in all_ids: continue
                m = Cursor(read(at, size)); m.take(24)
                samples = []
                for _ in range(nres):
                    z = m.unpack('4B3IB3xQQII2IQQ2I')
                    if z[0] == 0 and not z[1] and z[9]: samples.append((z[5], z[9], z[16]))
                if not samples: continue
                allowed = {sample[0] for sample in samples}
                selected = source if source in allowed else choose(allowed)
                _, cells, count = next(sample for sample in samples if sample[0] == selected)
                source_cells += cells; source_blocks += count
                source_resolution_pairs[selected] = source_resolution_pairs.get(selected, 0) + 1
        max_entry = vector_entries = 0
        for kind, (at, length) in enumerate(zip(locators[::2], locators[1::2])):
            if not length:
                continue
            v = Cursor(read(at, 16)); v.take(8); count = v.word(); v.take(4)
            offset = 16
            for _ in range(count):
                entry_length = struct.unpack('<I', read(at + offset, 4))[0]
                if entry_length < 4 or entry_length > min(length - offset, 512 * (1 << 20)):
                    raise ValueError('invalid or oversized V10 vector entry')
                if kind == 0 and rebuild:
                    if entry_length < 12: raise ValueError('invalid normalization entry length')
                    norm, chromosome = struct.unpack('<II', read(at + offset + 4, 8))
                    if chromosome in all_ids:
                        raise ValueError('cannot rebuild All while its normalization vectors exist')
                max_entry = max(max_entry, entry_length); offset += entry_length
            if offset != length:
                raise ValueError('invalid V10 vector index length')
            vector_entries += count
    existing_action = 'rebuild' if all_ids and rebuild else 'preserve' if all_ids else 'create'
    reported_length, reported_display, reported_bins = pseudo_length, display_resolution, (pseudo_length + display_resolution - 1) // display_resolution
    reported_scale = coordinate_scale
    if existing_action == 'preserve':
        attrs = dict(attributes)
        reported_length = chroms[all_ids[0]][1]
        reported_bins = int(attrs['hicAll.binCount']) if 'hicAll.binCount' in attrs else None
        reported_display = (int(attrs['hicAll.finestResolution']) * int(attrs.get('hicAll.storageStride', '1'))
                            if 'hicAll.finestResolution' in attrs else None)
        source = int(attrs['hicAll.sourceResolution']) if 'hicAll.sourceResolution' in attrs else None
        reported_scale = int(attrs['hicAll.coordinateScale']) if 'hicAll.coordinateScale' in attrs else None
        source_cells = source_blocks = 0
        source_resolution_pairs = {}
    real = [length for name, length in chroms if name.lower() != 'all']
    source_bins = sum((length + source - 1) // source for length in real) if source else None
    actual_bins = (pseudo_length + display_resolution - 1) // display_resolution
    pixels = actual_bins * (actual_bins + 1) // 2
    # Conservative planning allowance, not a formal bound: worst accepted block
    # decoding, matrix metadata caches, header/footer copies, All and vector entry.
    estimated_ram = (8 * GIB + threads * 8 * (1 << 20) + pairs * (256 + 192 * nres)
                     + header_length * 3 + footer_length * 3 + max_entry * 3
                     + pixels * (128 + 32 * len(resolutions[0])))
    return dict(input=str(path), input_bytes=st.st_size, genome=genome,
                fingerprint=dict(size=st.st_size, mtime_ns=st.st_mtime_ns, inode=st.st_ino),
                chromosomes=len(chroms), matrix_pairs=pairs, resolutions=nres,
                stored_blocks=blocks, block_index_bytes=index_bytes,
                source_resolution_bp=source, source_cells=source_cells,
                source_resolution_pairs=source_resolution_pairs,
                source_blocks=source_blocks, genome_source_bins=source_bins,
                existing_all=bool(all_ids), existing_all_action=existing_action, vector_entries=vector_entries,
                max_vector_entry_bytes=max_entry, conservative_ram_estimate_bytes=estimated_ram,
                all_bin_budget=bins, copy_threads=threads,
                output_pseudo_length=reported_length, coordinate_scale=reported_scale,
                display_resolution=reported_display, output_all_bins=reported_bins,
                warnings=([f'Selected {source} bp source is finer than requested {preferred} bp; '
                           'aggregation may take longer.'] if existing_action != 'preserve' and source < preferred else
                          [f'Selected {source} bp source is coarser than requested {preferred} bp; '
                           'this avoids scanning an excessively fine anchor.'] if existing_action != 'preserve' and source > preferred else []))


def preflight(args):
    report = inspect(args.input, args.all_bins, args.all_source_resolution, args.threads, args.rebuild_all)
    if report['existing_all_action'] == 'preserve':
        report['warnings'].append('Existing All will be preserved; use --rebuild-all to migrate it.')
    if report['conservative_ram_estimate_bytes'] > memory_bytes(args.memory) * 0.8:
        raise ValueError('metadata/block planning estimate exceeds 80% of requested RAM; request more RAM')
    if args.output.exists() and args.output != args.input and not args.overwrite:
        raise ValueError('output already exists; use a different path or explicit --overwrite')
    # Output is temporary beside the final target, also for in-place repair.
    needed = (0 if report['existing_all'] and not args.rebuild_all and args.output == args.input else
              report['input_bytes'] + max(2 * GIB, report['input_bytes'] // 20))
    free = shutil.disk_usage(args.output.parent).free
    if free < needed:
        raise ValueError(f'output filesystem needs about {needed / GIB:.1f} GiB free; has {free / GIB:.1f} GiB')
    report.update(output=str(args.output), requested_memory=args.memory,
                  required_free_space_bytes=needed, available_free_space_bytes=free,
                  ideal_copy_seconds=(0 if report['existing_all_action'] == 'preserve' and args.output == args.input
                                      else report['input_bytes'] / (args.copy_gbps * 1e9)))
    report['warnings'].append('Filesystem free space does not check user/project quotas; verify quota separately.')
    return report


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def worker(args):
    lock = Path(str(args.output) + '.add-all.lock')
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    try:
        with os.fdopen(fd, 'w') as out:
            json.dump(dict(job=os.environ.get('SLURM_JOB_ID'), pid=os.getpid()), out)
        report = preflight(args)
        planned = json.loads((args.work_dir / 'preflight.json').read_text())
        if report['fingerprint'] != planned['fingerprint']:
            raise ValueError('input changed since submission; generate a fresh preflight/job')
        write_json(args.work_dir / 'preflight-worker.json', report)
        print(json.dumps(report, indent=2), flush=True)
        command = [args.hic_v10, 'add-all', '-t', str(args.threads), '--all-bins', str(args.all_bins),
                   '--all-source-resolution', str(args.all_source_resolution)]
        if args.rebuild_all:
            command.append('--rebuild-all')
        if args.output == args.input:
            command += ['--in-place', str(args.input)]
        else:
            command += [str(args.input), str(args.output)]
        if sys.platform.startswith('linux') and Path('/usr/bin/time').exists():
            command = ['/usr/bin/time', '-v', '-o', str(args.work_dir / 'resources.txt'), *command]
        print(shlex.join(command), flush=True)
        begun = time.monotonic()
        result = subprocess.run(command)
        write_json(args.work_dir / 'result.json', dict(returncode=result.returncode,
                   elapsed_seconds=time.monotonic() - begun,
                   output_bytes=args.output.stat().st_size if result.returncode == 0 else None))
        return result.returncode
    finally:
        lock.unlink()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path)
    p.add_argument('output', nargs='?', type=Path)
    p.add_argument('--in-place', action='store_true')
    p.add_argument('--hic-v10', default='hic_v10')
    p.add_argument('--work-dir', type=Path, required=True, help='Shared parent for job scripts, reports and logs')
    p.add_argument('--threads', type=int, default=16, help='Copy workers and Slurm CPUs (1..256; start with 8..16)')
    p.add_argument('--memory', default='100G', help='Slurm total node memory, e.g. 100G or 150G')
    p.add_argument('--time', default='12:00:00', help='Slurm wall time')
    p.add_argument('--partition')
    p.add_argument('--account')
    p.add_argument('--all-bins', type=int, default=500)
    p.add_argument('--all-source-resolution', type=int, default=500000)
    p.add_argument('--copy-gbps', type=float, default=0.5, help='Assumed end-to-end copy GB/s for estimate only')
    p.add_argument('--rebuild-all', action='store_true', help='Reconstruct and replace an existing All')
    p.add_argument('--overwrite', action='store_true', help='Explicitly allow replacing an existing output')
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--sbatch', default='sbatch')
    p.add_argument('--run', action='store_true', help=argparse.SUPPRESS)
    args = p.parse_args(argv)
    try:
        if args.in_place == bool(args.output):
            p.error('choose an output path or --in-place')
        if not 1 <= args.threads <= 256 or not 1 <= args.all_bins <= 4096 or args.all_source_resolution <= 0 or not math.isfinite(args.copy_gbps) or args.copy_gbps <= 0:
            p.error('invalid threads, All settings or copy speed')
        memory_bytes(args.memory)
        args.input = args.input.resolve(strict=True)
        args.output = args.input if args.in_place else args.output.resolve()
        if not args.in_place and args.output == args.input:
            p.error('use --in-place explicitly to replace the input')
        if args.output.exists() and not args.in_place and os.path.samefile(args.input, args.output):
            p.error('output must not alias input')
        executable = shutil.which(args.hic_v10)
        if not executable:
            p.error('hic_v10 not found; build the updated executable or pass --hic-v10')
        args.hic_v10 = str(Path(executable).resolve())
        args.work_dir = args.work_dir.resolve()
        if args.run:
            return worker(args)
        report = preflight(args)
        args.work_dir.mkdir(parents=True, exist_ok=True)
        run_dir = Path(tempfile.mkdtemp(prefix='add-all-', dir=args.work_dir))
        write_json(run_dir / 'preflight.json', report)
        command = [sys.executable, str(Path(__file__).resolve()), str(args.input)]
        command += ['--in-place'] if args.in_place else [str(args.output)]
        command += ['--run', '--work-dir', str(run_dir), '--hic-v10', args.hic_v10,
                    '--threads', str(args.threads), '--memory', args.memory,
                    '--all-bins', str(args.all_bins), '--all-source-resolution', str(args.all_source_resolution),
                    '--copy-gbps', str(args.copy_gbps)]
        if args.rebuild_all:
            command.append('--rebuild-all')
        if args.overwrite:
            command.append('--overwrite')
        job = run_dir / 'repair.sbatch'
        job.write_text('#!/bin/bash\nset -euo pipefail\n' + shlex.join(command) + '\n')
        submit = [args.sbatch, '--parsable', '--job-name=hic-add-all', '--nodes=1', '--ntasks=1',
                  f'--cpus-per-task={args.threads}', f'--mem={args.memory}', f'--time={args.time}',
                  '--output=' + str(run_dir / 'slurm-%j.out'),
                  '--error=' + str(run_dir / 'slurm-%j.err')]
        if args.partition:
            submit.append('--partition=' + args.partition)
        if args.account:
            submit.append('--account=' + args.account)
        submit.append(str(job))
        print(json.dumps(report, indent=2))
        print('Job directory:', run_dir)
        print(shlex.join(submit), flush=True)
        write_json(run_dir / 'submission.json', dict(command=submit, status='dry-run' if args.dry_run else 'submitting'))
        if not args.dry_run:
            accepted = subprocess.run(submit, text=True, capture_output=True, check=True)
            job_id = accepted.stdout.strip().split(';')[0]
            if not job_id.isdigit():
                raise ValueError('unexpected sbatch response; check scheduler before resubmitting: ' + accepted.stdout)
            write_json(run_dir / 'submission.json', dict(command=submit, status='submitted', job_id=job_id))
            print('Submitted Slurm job', job_id)
        return 0
    except (OSError, ValueError, argparse.ArgumentTypeError, subprocess.CalledProcessError) as error:
        p.exit(1, f'Error: {error}\n')


if __name__ == '__main__':
    sys.exit(main())
