#!/usr/bin/env python3
"""Concurrent imported-norm workers, exact preservation, and failed/retried tasks."""
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
sys.dont_write_bytecode = True
from inspect_v10 import Hic
from test_text_norm import source, block, floats
from test_additive_norm import run, preserve, reorder
from v9_fixture import make

PHASES = ('cells', 'expected', 'reduce', 'compress-norm', 'compress-expected')


def count(tool, plan, phase):
    return int(run([tool, 'addnorm-tasks', '--count', '--phase', phase, plan]).stdout)


def phase(tool, plan, name, workers=4):
    n = count(tool, plan, name)
    # Reverse order and concurrent subprocesses exercise scheduler-independent
    # task IDs and immutable preparation; thread pools only launch processes.
    with ThreadPoolExecutor(max_workers=workers) as executor:
        list(executor.map(lambda i: run([tool, 'addnorm-run', '--memory', '32KiB',
                                        '--fan-in', '2', plan, name, i]), reversed(range(n))))
    return n


def equivalent(left, right):
    assert left.norms == right.norms
    assert left.attributes == right.attributes
    assert set(left.vectors) == set(right.vectors)
    for key in left.vectors:
        a, sa = left.vectors[key]; b, sb = right.vectors[key]
        assert len(a) == len(b), key
        if key[0] == 0:
            assert a == b, key
        else:
            for x, y in zip(floats(a), floats(b)):
                assert math.isnan(x) if math.isnan(y) else math.isclose(x, y, rel_tol=2e-6, abs_tol=1e-7), (key, x, y)
            sa, sb = dict(sa), dict(sb)
            assert set(sa) == set(sb)
            for chr in sa:
                x, y = floats([sa[chr]])[0], floats([sb[chr]])[0]
                assert math.isclose(x, y, rel_tol=2e-6), (key, chr, x, y)


def prepare(tool, path, files, directory, groups=1, rows=2, chunks=1):
    args = [tool, 'plan-addnorm', '--rows-per-task', rows, '--chunks-per-task', chunks,
            '--norms-per-task', groups]
    for file in files:
        args += ['--norm-file', file]
    run([*args, path, directory])
    plan = directory/'addnorm.plan'
    snapshot = plan.read_bytes()
    run([*args, path, directory])  # same-input preparation is idempotent
    assert plan.read_bytes() == snapshot
    return plan


def execute(tool, plan):
    # Compression of supplied divisors is independent of expected/cell work.
    with ThreadPoolExecutor(max_workers=2) as executor:
        tasks = [executor.submit(phase, tool, plan, name) for name in ('cells', 'compress-norm')]
        for task in tasks: task.result()
    phase(tool, plan, 'expected')
    phase(tool, plan, 'reduce')
    phase(tool, plan, 'compress-expected')


def mock_slurm(root, tool, plan):
    script = Path(__file__).resolve().parents[2]/'scripts/slurm_addnorm_v10_large.py'
    sbatch = root/'sbatch'
    # Record job IDs, then execute the generated scripts in dependency order
    # locally, with Slurm's zero-based array environment. This validates the
    # submission DAG without claiming access to a real cluster.
    sbatch.write_text('#!/usr/bin/env python3\nimport pathlib,sys\np=pathlib.Path(__file__).with_suffix(".counter")\nn=int(p.read_text())+1 if p.exists() else 100\np.write_text(str(n))\nprint(str(n)+";cluster")\n')
    sbatch.chmod(0o755)
    args = [sys.executable, script, plan, '--hic-v10-large', tool, '--sbatch', sbatch,
            '--array-batch-size', '3', '--max-jobs', '2', '--sort-memory', '32KiB', '--fan-in', '2',
            '--assemble-threads', '3']
    run([*args, '--dry-run'])
    assert not (plan.parent/'slurm-submit/jobs.json').exists()
    run(args)
    ledger = plan.parent/'slurm-submit/jobs.json'
    saved = ledger.read_bytes()
    run(args, ok=False)
    assert ledger.read_bytes() == saved
    assert json.loads(saved)['status'] == 'submitted'
    jobs = json.loads(saved)['jobs']; seen = set()
    journal = Path(json.loads(saved)['journal'])
    assert [json.loads(line) for line in journal.read_text().splitlines()] == jobs
    # Each ready group references precisely the scheduler jobs supplying its
    # declared task dependencies; no whole-phase barrier is introduced.
    dependencies = {}
    for name in PHASES:
        rows = run([tool, 'addnorm-tasks', '--phase', name, '--dependencies', plan]).stdout.splitlines()
        dependencies[name] = {int(row.split()[0]): row.split()[1:] for row in rows}
    owners = {(j['phase'], task): j['id'] for j in jobs for task in (j['tasks'] or [])}
    for j in jobs:
        for task in j['tasks'] or []:
            need = dependencies[j['phase']][task]
            wanted = {owners[(p, int(i))] for p, i in (word.split(':') for word in need)}
            assert set(j['dependencies']) == wanted
    cells = [j['id'] for j in jobs if j['phase'] == 'cells']
    assert len(cells) > 1
    assert all(len(j['dependencies']) == 1 for j in jobs if j['phase'] == 'expected')
    scripts = {j['script']: Path(j['script']).read_bytes() for j in jobs}
    run([*args, '--resubmit', '--dependency-mode', 'phase'])
    retry = json.loads(ledger.read_bytes())['jobs']
    assert not (set(scripts) & {j['script'] for j in retry})
    assert all(set(j['dependencies']) == {c['id'] for c in retry if c['phase']=='cells'}
               for j in retry if j['phase']=='expected')
    assert all(Path(p).read_bytes() == data for p, data in scripts.items())
    import os
    for job in jobs:
        assert set(job['dependencies']) <= seen
        command = job['command']
        assert '--parsable' in command
        if job['array']:
            lo, hi = map(int, job['array'].split('%')[0].split('-'))
            for i in range(lo, hi+1):
                env = dict(os.environ, SLURM_ARRAY_TASK_ID=str(i))
                result = subprocess.run(['bash', job['script']], env=env, text=True, capture_output=True)
                assert result.returncode == 0, result.stderr
        else:
            run(['bash', job['script']])
        seen.add(job['id'])
    assert jobs[-1]['phase'] == 'assemble'
    assert all('--dependency' not in j['command'] for j in jobs if j['phase'] == 'compress-norm')
    # A failure midway through submission leaves a usable partial ledger and
    # durable journal, and it cannot silently create a second DAG on retry.
    failing = root/'failing-sbatch'
    failing.write_text('#!/usr/bin/env python3\nimport pathlib,sys\np=pathlib.Path(__file__).with_suffix(".counter")\nn=int(p.read_text())+1 if p.exists() else 1\np.write_text(str(n))\nprint(1000+n)\nsys.exit(1 if n==2 else 0)\n')
    failing.chmod(0o755)
    # Use a fresh plan to retain the immutable input identity after assembly.
    target = root/'failed-submit.hic'; target.write_bytes((root/'base.hic').read_bytes())
    failed_plan = prepare(tool, target, [root/'RU.norm.txt'], root/'failed-submit-work')
    run([sys.executable, script, failed_plan, '--hic-v10-large', tool,
         '--sbatch', failing], ok=False)
    partial = json.loads((failed_plan.parent/'slurm-submit/jobs.json').read_text())
    assert partial['status']=='incomplete' and len(partial['jobs'])==1
    assert [json.loads(line) for line in Path(partial['journal']).read_text().splitlines()] == partial['jobs']
    run([sys.executable, script, failed_plan, '--hic-v10-large', tool,
         '--sbatch', failing], ok=False)
    assert failing.with_suffix('.counter').read_text()=='2'



def main():
    regular, large = map(Path, sys.argv[1:3])
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        hbs = root/'source.hbs.gz'; source(hbs)
        chrom = root/'chrom.sizes'; chrom.write_text('chr1\t50\nchr2\t30\n')
        base = root/'base.hic'; run([regular, 'pre', '-r', '10,20', hbs, base, chrom])
        run([regular, 'addnorm', '--norm', 'VC', base]); reorder(base)
        first = root/'RU.norm.txt'
        first.write_text('HIC_NORM_VECTORS 1\n'+block('RU', 'chr1', 10, [1.25, 2, 3, 4, 5])+ 
                         block('RU', 'chr2', 10, [2]*3)+block('RU', 'chr1', 20, [1,2,3])+block('RU', 'chr2', 20, [2,4]))
        # Explicit legacy geometry and custom quoted names; an existing VC
        # bundle must be skipped, including its potential metadata attributes.
        second = root/'custom.norm.txt'
        second.write_text('HIC_NORM_VECTORS 1\n'+block('"NDSCALE #/quoted"', 'chr1', 10,
                         ['bits:3f800001','bits:80000000','bits:7fc01234','inf','bits:40000000','bits:7fa01234'])
                         .replace('BP 10\n','BP 10\nsource-length 6\n')+
                         block('VC', 'chr1', 10, [1]*6).replace('BP 10\n','BP 10\nsource-length 6\n'))
        combined = root/'combined.txt'; combined.write_text(first.read_text()+second.read_text().split('\n',1)[1])
        baseline = root/'serial.hic'; baseline.write_bytes(base.read_bytes()); run([regular, 'addnorm', '--norm-file', combined, baseline]); reference = Hic(baseline)
        outputs = []
        for grouping in (1, 3):
            target = root/f'parallel-{grouping}.hic'; target.write_bytes(base.read_bytes()); before = Hic(target)
            directory = root/f'work {grouping}'; plan = prepare(large, target, [first,second], directory, grouping)
            untouched = target.read_bytes()
            assert count(large, plan, 'expected') > count(large, plan, 'cells')
            run([large, 'addnorm-assemble', plan], ok=False)
            run([large, 'addnorm-run', plan, 'reduce', 0], ok=False)
            assert target.read_bytes() == untouched
            phase(large, plan, 'compress-norm')
            phase(large, plan, 'cells')
            # An interrupted task may leave sidecars without its done marker.
            phase(large, plan, 'expected')
            result = directory/'expected-0.done'; result.unlink()
            run([large, 'addnorm-run', plan, 'expected', 0])
            partial = next(directory.glob('partial-*-b0.h10e'))
            original = partial.read_bytes(); partial.write_bytes(original[:-1])
            run([large, 'addnorm-run', plan, 'reduce', 0], ok=False)
            assert target.read_bytes() == untouched
            partial.write_bytes(original)
            phase(large, plan, 'reduce')
            phase(large, plan, 'compress-expected')
            snapshot = {p: p.read_bytes() for p in directory.glob('*.done')}
            for name in PHASES: phase(large, plan, name)
            assert snapshot == {p: p.read_bytes() for p in directory.glob('*.done')}
            pack = next(directory.glob('compressed-*.pack')); original = pack.read_bytes()
            pack.write_bytes(original[:-1]+bytes([original[-1]^1]))
            run([large, 'addnorm-assemble', '-t', 4, plan], ok=False)
            assert target.read_bytes() == untouched
            assert not list(root.glob('*.addnorm.*'))
            pack.write_bytes(original)
            run([large, 'addnorm-assemble', '-t', 4 if grouping==1 else 1, plan])
            after = Hic(target); preserve(before, after, allow_new_attributes=True); equivalent(reference, after)
            saved = target.read_bytes(); run([large, 'addnorm-assemble', plan]); assert target.read_bytes() == saved
            if outputs:
                # Single-thread and threaded matrix/vector relocation produce
                # exactly the same file, independent of norm worker grouping.
                assert saved == outputs[0].data
            outputs.append(after)
            # The completion marker must not mask a subsequent target change.
            modified = bytearray(saved); modified[-1] ^= 1; target.write_bytes(modified)
            run([large, 'addnorm-assemble', plan], ok=False)
            assert target.read_bytes() == modified
        equivalent(outputs[0], outputs[1])
        # A full submission script with batched arrays and dependencies is
        # exercised against a mocked scheduler, then scripts run locally.
        target = root/'slurm.hic'; target.write_bytes(base.read_bytes())
        plan = prepare(large, target, [first,second], root/'slurm-work')
        mock_slurm(root, large, plan); equivalent(reference, Hic(target))
        # File mutation invalidates the plan before worker outputs are used.
        target = root/'changed.hic'; target.write_bytes(base.read_bytes())
        plan = prepare(large, target, [first], root/'changed-work')
        data = bytearray(target.read_bytes()); data[-1] ^= 1; target.write_bytes(data)
        run([large, 'addnorm-run', plan, 'cells', 0], ok=False)
        run([large, 'addnorm-assemble', plan], ok=False)
        assert target.read_bytes() == data
        # Supplied text is staged once; compute nodes need only durable sidecars.
        target = root/'staged.hic'; target.write_bytes(base.read_bytes())
        gone = root/'gone.txt'; gone.write_bytes(first.read_bytes())
        plan = prepare(large, target, [gone], root/'staged-work'); gone.unlink()
        execute(large, plan); run([large, 'addnorm-assemble', plan])
        after = Hic(target); preserve(Hic(base), after)
        for key, value in reference.vectors.items():
            if key[0] == 0 and key[1] == 'RU': assert after.vectors[key] == value
        # Bad preparation does not publish a plan or modify the target.
        target = root/'bad-plan.hic'; target.write_bytes(base.read_bytes()); untouched=target.read_bytes()
        run([large, 'plan-addnorm', '--norm-file', first, '--norm-file', first, target, root/'bad-work'], ok=False)
        assert not (root/'bad-work').exists() and target.read_bytes() == untouched
        run([large, 'plan-addnorm', '--rows-per-task', 0, '--norm-file', first, target, root/'bad-work'], ok=False)
        assert not (root/'bad-work').exists()
        scores = root/'scores.hic'; run([regular, 'pre', '--scores', '-r', 10, hbs, scores, chrom])
        text = root/'score.norm.txt'; text.write_text('HIC_NORM_VECTORS 1\n'+block('SCORE', 'chr1', 10, [1]*5))
        plan = prepare(large, scores, [text], root/'score-work'); saved=scores.read_bytes()
        run([large, 'addnorm-run', plan, 'cells', 0], ok=False)
        run([large, 'addnorm-assemble', plan], ok=False); assert scores.read_bytes()==saved
        # Missing raw EVI is rebuilt from all cis chromosomes while normalized
        # expected contributions include only chromosomes with supplied norms.
        target = root/'no-evi.hic'; data = bytearray(base.read_bytes()); struct.pack_into('<QQ',data,48,0,0); target.write_bytes(data)
        plan = prepare(large,target,[first],root/'no-evi-work'); execute(large,plan); run([large,'addnorm-assemble',plan]); after=Hic(target)
        for bin in (10,20): assert after.vectors[1,None,None,0,bin] == Hic(base).vectors[1,None,None,0,bin]
        # FRAG metadata, score rejection, and no-op preservation of bundles.
        v9, frag = root/'frag.v9.hic', root/'frag.hic'; make(v9,frag=True); run([regular,'convert',v9,frag])
        text=root/'frag.txt'; text.write_text('HIC_NORM_VECTORS 1\n'+block('FRAG_CUSTOM','chr1',1,[2]*8,'FRAG'))
        serial=root/'frag-serial.hic'; serial.write_bytes(frag.read_bytes()); run([regular,'addnorm','--norm-file',text,serial])
        plan=prepare(large,frag,[text],root/'frag-work'); execute(large,plan); run([large,'addnorm-assemble',plan]); equivalent(Hic(serial),Hic(frag))
        plan=prepare(large,frag,[text],root/'noop-work'); assert all(count(large,plan,name)==0 for name in PHASES)
        saved=frag.read_bytes();run([large,'addnorm-assemble',plan]);assert frag.read_bytes()==saved
        # Large vectors span multiple independently compressed 64K chunks.
        long_chrom=root/'long.sizes';long_chrom.write_text('chr1\t700010\n')
        pairs=root/'long.pairs';pairs.write_text('chr1 0 chr1 0\nchr1 0 chr1 700000\n')
        long_hic=root/'long.hic';run([regular,'pre','-r',10,pairs,long_hic,long_chrom])
        text=root/'long.txt';text.write_text('HIC_NORM_VECTORS 1\nvector LONG chr1 BP 10\n'+'bits:3f800001\n'*70001+'end\n')
        serial=root/'long-serial.hic';serial.write_bytes(long_hic.read_bytes());run([regular,'addnorm','--norm-file',text,serial])
        plan=prepare(large,long_hic,[text],root/'long-work',rows=40000)
        assert count(large,plan,'compress-norm')==2
        execute(large,plan);run([large,'addnorm-assemble','-t',3,plan]);equivalent(Hic(serial),Hic(long_hic))
    print('Distributed supplied normalization tests passed')

if __name__ == '__main__': main()
