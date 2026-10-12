#!/usr/bin/env python3
"""Add stored scaled All using the native V10 reader and atomic repacker."""
import argparse
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Local v10 .hic file')
    parser.add_argument('output', nargs='?', help='Output .hic file; required unless --in-place')
    parser.add_argument('--in-place', action='store_true', help='Replace input atomically after validation')
    parser.add_argument('--hic-v10', default='hic_v10', help='Path to native hic_v10 executable')
    parser.add_argument('--all-bins', type=int, default=500, help='Maximum overview width (1..4096)')
    parser.add_argument('--all-source-resolution', type=int, default=2500000)
    parser.add_argument('--level', type=int, default=6)
    parser.add_argument('-t', type=int, default=4)
    args = parser.parse_args()
    if args.in_place == bool(args.output):
        parser.error('choose an output path or --in-place')
    executable = shutil.which(args.hic_v10)
    if not executable: parser.error('hic_v10 was not found; build with BUILD_V10=ON or pass --hic-v10')
    command = [executable, 'add-all', '--all-bins', str(args.all_bins),
               '--all-source-resolution', str(args.all_source_resolution), '--level', str(args.level), '-t', str(args.t)]
    if args.in_place: command += ['--in-place', args.input]
    else: command += [args.input, args.output]
    subprocess.run(command, check=True)


if __name__ == '__main__': main()
