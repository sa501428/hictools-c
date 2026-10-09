#!/usr/bin/env python3
"""Submit a prepared v10 text-normalization import as Slurm job arrays.

Does not submit jobs with --dry-run. Input/planning remain explicit CLI steps.
The target, staged plan directory, and executable must be visible on nodes.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile

PHASES = ('cells', 'expected', 'reduce', 'compress-norm', 'compress-expected')


def positive(value):
    n = int(value)
    if n <= 0:
        raise argparse.ArgumentTypeError('must be positive')
    return n


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('--hic-v10-large', default='hic_v10_large')
    parser.add_argument('--sbatch', default='sbatch')
    parser.add_argument('--partition')
    parser.add_argument('--account')
    parser.add_argument('--max-jobs', type=positive, default=100,
                        help='concurrent tasks per submitted array (not a global ceiling)')
    parser.add_argument('--array-batch-size', type=positive, default=1000,
                        help='split arrays to fit the site MaxArraySize; indexes restart at zero')
    parser.add_argument('--cell-array-batch-size', type=positive, default=1,
                        help='cell tasks per array; 1 releases each chromosome immediately')
    parser.add_argument('--dependency-mode', choices=('ready', 'phase'), default='ready',
                        help='ready releases each dependency group immediately; phase uses whole-phase barriers')
    parser.add_argument('--job-memory', default='32G', help='Slurm --mem per worker')
    parser.add_argument('--assemble-memory', default='32G')
    parser.add_argument('--assemble-threads', type=positive, default=8)
    parser.add_argument('--sort-memory', default='8GiB', help='cell sort budget; separate from --mem')
    parser.add_argument('--fan-in', type=positive, default=64)
    parser.add_argument('--tmp', help='worker scratch; omit to use the shared plan directory')
    parser.add_argument('--time', default='2-00:00:00')
    parser.add_argument('--assemble-time', default='2-00:00:00')
    parser.add_argument('--log-dir', type=Path)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--resubmit', action='store_true',
                        help='explicitly allow another submission in a work directory with a submission ledger')
    args = parser.parse_args(argv)
    if args.fan_in < 2:
        parser.error('--fan-in must be at least 2')
    plan = args.plan.resolve(strict=True)
    executable = shutil.which(args.hic_v10_large)
    if executable is None:
        parser.error('hic_v10_large executable not found')
    executable = str(Path(executable).resolve())
    task_dependencies = {}
    for phase in PHASES:
        result = subprocess.run([executable, 'addnorm-tasks', '--phase', phase,
                                 '--dependencies', str(plan)],
                                text=True, capture_output=True, check=True)
        tasks = []
        for line in result.stdout.splitlines():
            words = line.split()
            if not words or int(words[0]) != len(tasks):
                raise RuntimeError('invalid task DAG from hic_v10_large')
            dependencies = []
            for word in words[1:]:
                dep_phase, dep_id = word.split(':')
                dependencies.append((dep_phase, int(dep_id)))
            tasks.append(dependencies)
        task_dependencies[phase] = tasks
    logs = (args.log_dir or plan.parent/'slurm-logs').resolve()
    submit_dir = plan.parent/'slurm-submit'
    logs.mkdir(parents=True, exist_ok=True)
    submit_dir.mkdir(exist_ok=True)
    # Hold a local/shared-filesystem lock for the complete submission so two
    # launchers cannot both pass the ledger check and create duplicate DAGs.
    submit_lock = (submit_dir/'submit.lock').open('a')
    try:
        fcntl.flock(submit_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        parser.error('another launcher is submitting this plan')
    ledger = submit_dir/'jobs.json'
    if ledger.exists() and not args.resubmit:
        parser.error(f'{ledger} already exists; inspect it and use --resubmit only when retries are intended')
    run_dir = Path(tempfile.mkdtemp(prefix='submission-', dir=submit_dir))
    if ledger.exists() and not args.dry_run:
        shutil.copy2(ledger, run_dir/'previous-jobs.json')
    job_records = []
    serial = 0
    journal_path = run_dir/'jobs.jsonl'
    journal = journal_path.open('w') if not args.dry_run else None

    def save_ledger(status):
        if args.dry_run:
            return
        temporary = ledger.with_suffix('.tmp')
        with temporary.open('w') as output:
            json.dump({'plan': str(plan), 'jobs': job_records, 'journal': str(journal_path),
                       'status': status}, output, indent=2)
            output.write('\n')
            output.flush()
            os.fsync(output.fileno())
        temporary.replace(ledger)
        shutil.copy2(ledger, run_dir/'jobs.json')

    # Reserve the ledger before any scheduler calls. Each accepted job is
    # appended/fsynced in O(1) journal writes; a final snapshot avoids O(N²)
    # rewriting of a potentially large dependency graph during submission.
    save_ledger('submitting')

    def submit(phase, command, dependencies, array=None, tasks=None):
        nonlocal serial
        serial += 1
        # Keep scripts/logs distinct on resubmission. JSON ledger remains a
        # durable record of successful sbatch calls if a later call fails.
        script = run_dir/f'{phase}-{serial}.sh'
        script.write_text('#!/bin/bash\nset -euo pipefail\n'+command+'\n')
        flags = [args.sbatch, '--parsable', '--job-name', 'hic-addnorm-'+phase,
                 '--cpus-per-task', str(args.assemble_threads if phase == 'assemble' else 1),
                 '--mem', args.assemble_memory if phase == 'assemble' else args.job_memory,
                 '--time', args.assemble_time if phase == 'assemble' else args.time,
                 '--output', str(logs/(phase+'-%A_%a.out')),
                 '--error', str(logs/(phase+'-%A_%a.err'))]
        if args.partition:
            flags += ['--partition', args.partition]
        if args.account:
            flags += ['--account', args.account]
        if dependencies:
            flags += ['--dependency', 'afterok:'+':'.join(dependencies)]
        if array:
            flags += ['--array', array]
        flags.append(str(script))
        print(shlex.join(flags), flush=True)
        if args.dry_run:
            job_id = str(900000+serial)
        else:
            result = subprocess.run(flags, text=True, capture_output=True, check=True)
            job_id = result.stdout.strip().split(';')[0]
            if not job_id.isdigit():
                raise RuntimeError(f'unrecognized sbatch --parsable response: {result.stdout!r}')
        job_records.append({'phase': phase, 'id': job_id, 'array': array, 'dependencies': dependencies,
                            'script': str(script), 'command': flags, 'tasks': tasks})
        if journal:
            journal.write(json.dumps(job_records[-1])+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        return job_id

    task_jobs = {}

    def arrays(phase, phase_dependencies):
        ids = []
        task_jobs[phase] = {}
        # Tasks with the same submitted prerequisites share arrays; they need
        # not have contiguous phase ordinals. The literal Bash lookup maps the
        # compact Slurm index to the actual immutable plan task.
        groups = {}
        for ordinal, prerequisite_tasks in enumerate(task_dependencies[phase]):
            if args.dependency_mode == 'phase':
                prerequisites = tuple(phase_dependencies)
            else:
                prerequisites = tuple(dict.fromkeys(task_jobs[p][i] for p, i in prerequisite_tasks))
            groups.setdefault(prerequisites, []).append(ordinal)
        batch_size = args.cell_array_batch_size if phase == 'cells' else args.array_batch_size
        for prerequisites, ordinals in groups.items():
            for offset in range(0, len(ordinals), batch_size):
                tasks = ordinals[offset:offset+batch_size]
                words = [executable, 'addnorm-run']
                if phase == 'cells':
                    words += ['--memory', args.sort_memory, '--fan-in', str(args.fan_in)]
                    if args.tmp:
                        words += ['--tmp', args.tmp]
                words += [str(plan), phase]
                lookup = 'task_ids=('+' '.join(map(str, tasks))+')\n'
                command = lookup+'exec '+shlex.join(words)+' "${task_ids[$SLURM_ARRAY_TASK_ID]}"'
                job_id = submit(phase, command, list(prerequisites),
                                f'0-{len(tasks)-1}%{args.max_jobs}', tasks)
                ids.append(job_id)
                for ordinal in tasks:
                    task_jobs[phase][ordinal] = job_id
        return ids or phase_dependencies

    submission_complete = False
    try:
        cells = arrays('cells', [])
        norm_compression = arrays('compress-norm', [])  # independent of cis extraction/expectations
        expected = arrays('expected', cells)
        reductions = arrays('reduce', expected)
        expected_compression = arrays('compress-expected', reductions)
        dependencies = list(dict.fromkeys(norm_compression+expected_compression+reductions))
        command = shlex.join([executable, 'addnorm-assemble', '-t', str(args.assemble_threads), str(plan)])
        final = submit('assemble', 'exec '+command, dependencies)
        print(('Dry-run assembly placeholder: ' if args.dry_run else 'Assembly job: ')+final)
        submission_complete = True
    finally:
        if journal:
            journal.close()
        save_ledger('submitted' if submission_complete else 'incomplete')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, subprocess.CalledProcessError, RuntimeError) as error:
        print(f'Error: {error}', file=sys.stderr)
        raise SystemExit(1)
