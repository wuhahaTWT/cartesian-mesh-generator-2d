#!/usr/bin/env python3
"""Native CLI relaxation experiments; no independent PDE implementation.

Immutable manifests and per-run records make interrupted batches resumable.
Selection uses training cases only. Limits and different solutions never win.
See docs/DEVELOPMENT_CN.md for controls, normalization and interpretation.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import queue
import random
import shutil
import signal
import statistics
import subprocess
import sys
import tarfile
import threading
import time

# Isolated resource accounting avoids mixing concurrent children's CPU time.
# Handle this before importing NumPy or benchmark helpers to keep launch cheap.
if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == '--measure-native':
    child = subprocess.Popen(sys.argv[3:])
    _, wait_status, usage = os.wait4(child.pid, 0)
    child.returncode = os.waitstatus_to_exitcode(wait_status)
    Path(sys.argv[2]).write_text(json.dumps(dict(userSeconds=usage.ru_utime,
        systemSeconds=usage.ru_stime, maxRssKiB=usage.ru_maxrss)) + '\n')
    sys.exit(child.returncode if child.returncode >= 0 else 128-child.returncode)

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools' / 'flow'))
import laminar
import extract_steady_iterate as iterate
import native_mesh

SCRIPT_SOURCE = Path(__file__).read_bytes()
SCRIPT_SHA = hashlib.sha256(SCRIPT_SOURCE).hexdigest()
THREAD_ENV = {k: '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                              'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']}
METHODS = ['simple', 'simple-consistent']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def freeze(path, value):
    path = Path(path)
    if path.exists():
        if json.loads(path.read_text()) != value:
            raise ValueError(f'Immutable experiment changed: {path}; use a new batch/root')
    else:
        write(path, value)


def compact(summary):
    return {k: v for k, v in summary.items() if k != 'boundaryConditions'}


def initialise(root, flow_cli, mesh_cli):
    root.mkdir(parents=True, exist_ok=True)
    if (root / 'manifest.json').exists():
        raise ValueError('Study already exists; resume with run/select/report')
    inputs = root / 'inputs'
    inputs.mkdir()
    names = [f'{case}.{suffix}' for case in ['annulus-l5', 'annulus-l6',
                                            'annulus-l7', 'dfg20-l7']
             for suffix in ['solver.cm2d', 'boundaries']]
    archive = ROOT / 'artifacts/current/laminar-foundation/failure-inputs.tar.gz'
    with tarfile.open(archive) as bundle:
        for name in names:
            members = [m for m in bundle.getmembers() if Path(m.name).name == name]
            if len(members) != 1 or not members[0].isfile():
                raise ValueError(f'Missing/ambiguous archived input: {name}')
            (inputs / name).write_bytes(bundle.extractfile(members[0]).read())
    cases = {}
    for name, size in [('annulus-l5', 484), ('annulus-l6', 1672),
                       ('annulus-l7', 6208), ('dfg20-l7', 3063)]:
        physical = name.split('-')[0]
        cases[name] = dict(mesh=str((inputs / f'{name}.solver.cm2d').relative_to(root)),
                           boundary=str((inputs / f'{name}.boundaries').relative_to(root)),
                           cells=size, physical=physical, **laminar.controls(physical),
                           convection='limited-linear' if physical == 'dfg20' else 'upwind',
                           closed=physical == 'annulus', length=2 if physical == 'annulus' else .41)
    for case in ['poiseuille', 'cavity100']:
        mesh = laminar.prepare(case, 5, root, mesh_cli, 90)
        if mesh['generation']['returnCode'] != 0:
            raise RuntimeError(f'Mesh generation failed: {mesh}')
        cases[case + '-l5'] = dict(mesh=str(Path(mesh['path']).relative_to(root)),
                                 boundary=str(Path(mesh['boundary']).relative_to(root)),
                                 cells=mesh['cells'], physical=case, **laminar.controls(case),
                                 convection='upwind', closed=case.startswith('cavity'), length=1.)
    for case in cases.values():
        case['hashes'] = {key: sha(root / case[key]) for key in ['mesh', 'boundary']}
    binary = root / 'bin' / 'flow-cli'
    binary.parent.mkdir()
    shutil.copy2(flow_cli, binary)
    manifest = dict(format='cartmesh2d-relaxation-study-v1', cases=cases,
                    sourceCommit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                    sourceDirty=subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True),
                    binary='bin/flow-cli', binarySha256=sha(binary), meshBinarySha256=sha(mesh_cli),
                    archiveSha256=sha(archive), platform=platform.platform(), python=platform.python_version(),
                    numpy=np.__version__, cpus=sorted(os.sched_getaffinity(0)), threadEnv=THREAD_ENV,
                    training=['annulus-l5', 'poiseuille-l5'],
                    holdout=['annulus-l6', 'annulus-l7', 'cavity100-l5', 'dfg20-l7'],
                    alphaU=[.2, .4, .6, .8, 1.], alphaP=[.05, .1, .25, .5, 1.], seed=20261009,
                    controls=dict(tolerance=1e-8, referenceTolerance=1e-10, linearPolicy='adaptive',
                                  convergence='strict', viscousStress='symmetric', pressureCorrections=4,
                                  pressurePreconditioner='ic0', inertia=1, acceleration='none'),
                    agreementThreshold=1e-4,
                    agreementDefinition=dict(velocity='area RMS vector difference / Uref',
                        pressure='area RMS difference / (Uref^2 + nu*Uref/H); subtract area mean only in closed domains',
                        flux='RMS face-flux difference / (Uref*sqrt(mean cell area))',
                        scope='same discrete solution; not physical accuracy certification'),
                    objective='cold-start full native process wall time, including loading and all native exports',
                    selection='require eligible on every training case; minimize geometric mean time relative to coupled; held-out data excluded')
    write(root / 'manifest.json', manifest)
    return manifest


def load(root):
    manifest = json.loads((root / 'manifest.json').read_text())
    if sha(root / manifest['binary']) != manifest['binarySha256']:
        raise ValueError('Frozen solver binary changed')
    for case in manifest['cases'].values():
        for key, digest in case['hashes'].items():
            if sha(root / case[key]) != digest:
                raise ValueError(f'Frozen input changed: {case[key]}')
    return manifest


def fields(prefix):
    return tuple(np.genfromtxt(str(prefix) + suffix, delimiter=',', names=True, ndmin=1)
                 for suffix in ['.cells.csv', '.faces.csv'])


def agreement(prefix, reference, case):
    a, af = fields(prefix)
    b, bf = fields(reference)
    if a.shape != b.shape or af.shape != bf.shape:
        raise ValueError('Field coverage differs from reference')
    for key in ['cell', 'x', 'y', 'area']:
        if not np.array_equal(a[key], b[key]):
            raise ValueError(f'Cell identity/geometry mismatch: {key}')
    for key in ['face', 'owner', 'neighbour']:
        if not np.array_equal(af[key], bf[key]):
            raise ValueError(f'Face identity mismatch: {key}')
    if not all(np.isfinite(a[key]).all() for key in ['u', 'v', 'p']) or not np.isfinite(af['flux']).all():
        return dict(finite=False)
    weights = a['area'] / np.sum(a['area'])
    dp = a['p'] - b['p']
    if case['closed']:
        dp -= np.sum(weights * dp)
    u, nu, height = case['speed'], case['nu'], case['length']
    errors = dict(velocity=float(np.sqrt(np.sum(weights * ((a['u']-b['u'])**2 + (a['v']-b['v'])**2))) / u),
                  pressure=float(np.sqrt(np.sum(weights * dp**2)) / (u*u + nu*u/height)),
                  flux=float(np.sqrt(np.mean((af['flux']-bf['flux'])**2)) / (u*np.sqrt(np.mean(a['area'])))))
    return dict(finite=all(math.isfinite(v) for v in errors.values()), **errors)


def classify(code, timed_out, summary, errors, threshold, log):
    if timed_out:
        return 'time-limit'
    if summary.get('converged'):
        if code != 0 or not summary.get('strictLinearFinal'):
            return 'uncertified'
        if errors is None:
            return 'reference'
        if not errors['finite']:
            return 'nonfinite'
        return 'eligible' if max(errors[k] for k in ['velocity', 'pressure', 'flux']) <= threshold else 'different-solution'
    if code == 1:
        return 'linear-failure' if 'linear' in log.lower() else 'solver-error'
    if summary.get('failureReason'):
        return 'solver-failure'
    return 'iteration-limit' if code == 2 else 'process-error'


def token(value):
    return format(value, '.6g').replace('.', 'p')


def run_one(root, manifest, batch, case_name, method, pair, repeat, budget, core,
            reference=False, initial=None, label=None, linear_policy=None, compare=True, reference_tolerance=None):
    runner = root / 'runners' / (SCRIPT_SHA + '.py')
    runner.parent.mkdir(exist_ok=True)
    if not runner.exists():
        runner.write_bytes(SCRIPT_SOURCE)
    case = manifest['cases'][case_name]
    name = label or f'{case_name}__{method}__u{token(pair[0])}_p{token(pair[1])}__r{repeat}'
    directory = root / 'runs' / batch / name
    directory.mkdir(parents=True, exist_ok=True)
    prefix = directory / 'flow'
    tolerance = manifest['controls']['referenceTolerance' if reference else 'tolerance']
    if reference_tolerance is not None:
        if not reference or not 0 < reference_tolerance < manifest['controls']['tolerance']:
            raise ValueError('Reference override must be positive and tighter than candidate tolerance')
        tolerance = reference_tolerance
    command = [str(root / manifest['binary']), '--mesh', str(root / case['mesh']),
               '--boundary', str(root / case['boundary']), '--case', 'custom',
               '--coupling', method, '--nu', str(case['nu']), '--speed', str(case['speed']),
               '--viscous-stress', 'symmetric', '--convection', case['convection'],
               '--momentum-inertia', '1', '--max-iterations', str(budget['iterations']),
               '--tolerance', str(tolerance), '--profile', '--output', str(prefix),
               '--linear-policy', linear_policy or manifest['controls']['linearPolicy'], '--convergence', 'strict',
               '--pressure-preconditioner', 'ic0', '--pressure-corrections', '4',
               '--velocity-relaxation', str(pair[0]), '--pressure-relaxation', str(pair[1])]
    if initial:
        command += ['--initial-guess', str(initial) + '.cells.csv',
                    '--initial-flux', str(initial) + '.flux.csv']
    spec = dict(case=case_name, method=method, pair=list(pair), repeat=repeat, reference=reference,
                budget=budget, command=command, binarySha256=manifest['binarySha256'],
                inputHashes=case['hashes'], initialHashes=None if initial is None else
                {ext: sha(str(initial) + ext) for ext in ['.cells.csv', '.flux.csv']})
    if not compare:
        spec['compare'] = False
    freeze(directory / 'spec.json', spec)
    record_path = directory / 'record.json'
    if record_path.exists():
        record = json.loads(record_path.read_text())
        if reference and record['status'] == 'reference':
            preserve_reference(root, case_name, prefix)
        return record
    # Never overwrite a partially executed run after interruption.
    if (directory / 'native.log').exists():
        raise RuntimeError(f'Interrupted run retained at {directory}; use a new batch')
    timeout_event = threading.Event()
    started = time.perf_counter()
    process_command = [sys.executable, str(Path(__file__).resolve()), '--measure-native',
                       str(directory / 'resources.txt'), 'taskset', '-c', str(core), *command]
    with (directory / 'native.log').open('w') as stream:
        proc = subprocess.Popen(process_command, stdout=stream, stderr=subprocess.STDOUT,
                                env={**os.environ, **THREAD_ENV}, start_new_session=True)
        def expire():
            if proc.poll() is None:
                timeout_event.set()
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        timer = threading.Timer(budget['seconds'], expire)
        timer.start()
        try:
            proc.wait()
        except BaseException:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            raise
        finally:
            timer.cancel()
            timer.join()
    timed_out = timeout_event.is_set()
    elapsed = time.perf_counter() - started
    log = (directory / 'native.log').read_text()
    summary_path = Path(str(prefix) + '.json')
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    perf_path = Path(str(prefix) + '.performance.json')
    perf = json.loads(perf_path.read_text()) if perf_path.exists() else {}
    resources = {}
    if (directory / 'resources.txt').exists():
        for line in (directory / 'resources.txt').read_text().splitlines():
            if line.startswith('{'):
                resources = json.loads(line)
    ref_prefix = root / 'references' / case_name
    errors = agreement(prefix, ref_prefix, case) if summary.get('converged') and not reference and compare else None
    status = classify(proc.returncode, timed_out, summary, errors, manifest['agreementThreshold'], log)
    if status == 'reference' and not reference:
        status = 'uncompared'
    record = dict(**spec, batch=batch, name=name, prefix=str(prefix.relative_to(root)), status=status,
                  returnCode=proc.returncode, timedOut=timed_out, wallSeconds=elapsed, cpu=core,
                  runnerSha256=SCRIPT_SHA, summary=compact(summary), performance=perf, agreement=errors,
                  resources=resources, logTail=log[-1800:] if not summary.get('converged') else None)
    record['outputHashes'] = {p.name: sha(p) for p in directory.glob('flow.*')}
    write(record_path, record)
    if reference and status == 'reference':
        preserve_reference(root, case_name, prefix)
    print(json.dumps(dict(batch=batch, case=case_name, method=method, pair=pair, status=status,
                          iterations=summary.get('iterations'), seconds=round(elapsed, 3))), flush=True)
    return record


def preserve_reference(root, case_name, prefix):
    target = root / 'references' / case_name
    target.parent.mkdir(exist_ok=True)
    hashes = {}
    for ext in ['.cells.csv', '.faces.csv', '.json']:
        hashes[ext] = sha(str(prefix) + ext)
        dest = Path(str(target) + ext)
        if dest.exists() and sha(dest) != hashes[ext]:
            raise ValueError(f'Frozen reference changed: {dest}')
        if not dest.exists():
            shutil.copy2(str(prefix) + ext, dest)
    freeze(Path(str(target) + '.hashes.json'), hashes)


def run_batch(root, manifest, batch, cases, methods, pairs, budget, repeats, workers, reference=False, linear_policy=None, reference_tolerance=None):
    batch_spec = dict(cases=cases, methods=methods, pairs=pairs, budget=budget,
                      repeats=repeats, workers=workers, reference=reference)
    if linear_policy is not None:
        batch_spec['linearPolicy'] = linear_policy
    if reference_tolerance is not None:
        batch_spec['referenceTolerance'] = reference_tolerance
    freeze(root / 'batches' / (batch + '.json'), batch_spec)
    if not reference:
        for case in cases:
            if not (root / 'references' / (case + '.json')).exists():
                raise ValueError(f'Run a converged reference first: {case}')
            hashes = json.loads((root / 'references' / (case + '.hashes.json')).read_text())
            for ext, digest in hashes.items():
                if sha(root / 'references' / (case + ext)) != digest:
                    raise ValueError(f'Frozen reference changed: {case}{ext}')
    jobs = [(case, method, pair, repeat) for case in cases for method in methods
            for pair in ([(1., 1.)] if method == 'coupled' else pairs) for repeat in range(repeats)]
    random.Random(manifest['seed']).shuffle(jobs)
    cores = queue.Queue()
    for cpu in manifest['cpus'][:workers]:
        cores.put(cpu)
    def job(args):
        core = cores.get()
        try:
            return run_one(root, manifest, batch, *args, budget, core, reference,
                           linear_policy=linear_policy, reference_tolerance=reference_tolerance)
        finally:
            cores.put(core)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(job, jobs))


def records(root, batches=None):
    result = [json.loads(path.read_text()) for path in sorted((root / 'runs').glob('*/*/record.json'))]
    return [r for r in result if batches is None or r['batch'] in batches]


def run_role(record):
    if record['reference']:
        return 'reference'
    if record['batch'].startswith('schedule-'):
        return 'segment'
    if record.get('compare') is False:
        return 'reference-bootstrap'
    return 'warm-diagnostic' if record.get('initialHashes') else 'cold-trial'


def bootstrap_reference(root, manifest, case_name):
    """Tight certification from a converged ordinary solve; retain BOTH costs."""
    base = run_one(root, manifest, 'reference-bootstrap', case_name, 'coupled', [1.,1.], 0,
                   dict(iterations=2000, seconds=600), manifest['cpus'][0], compare=False)
    if base['status'] != 'uncompared':
        raise ValueError('Reference bootstrap failed native strict convergence')
    case = manifest['cases'][case_name]
    prefix = root / base['prefix']
    mesh = native_mesh.read_cm2d(root / case['mesh'])
    cells, faces = iterate.extract(mesh, iterate.rows(Path(str(prefix)+'.cells.csv')),
                                  iterate.rows(Path(str(prefix)+'.faces.csv')))
    initial = root / 'references' / (case_name+'-bootstrap')
    initial.parent.mkdir(exist_ok=True)
    for ext, header, data in [('.cells.csv','cell,x,y,u,v,p',cells),
                              ('.flux.csv','face,owner,neighbour,x,y,flux',faces)]:
        with Path(str(initial)+ext).open('w') as stream:
            writer = csv.writer(stream); writer.writerow(header.split(',')); writer.writerows(data)
    return run_one(root, manifest, 'reference-warm-tight', case_name, 'coupled', [1.,1.], 0,
                   dict(iterations=2000, seconds=600), manifest['cpus'][0], reference=True, initial=initial)


def select(root, manifest, batches, output):
    data = [r for r in records(root, batches) if r['case'] in manifest['training'] and run_role(r)=='cold-trial']
    selected = {}
    for method in METHODS:
        candidates = []
        pairs = sorted({tuple(r['pair']) for r in data if r['method'] == method})
        for pair in pairs:
            costs = {}; valid = True
            for case in manifest['training']:
                rows = [r for r in data if r['case'] == case and r['method'] == method and tuple(r['pair']) == pair]
                if not rows or any(r['status'] != 'eligible' for r in rows):
                    valid = False
                    break
                costs[case] = statistics.median(r['wallSeconds'] for r in rows)
            if valid:
                candidates.append(dict(pair=list(pair), geometricMeanSeconds=statistics.geometric_mean(costs.values()), costs=costs))
        candidates.sort(key=lambda r: r['geometricMeanSeconds'])
        selected[method] = candidates
    result = dict(training=manifest['training'], batches=batches, rankings=selected,
                  note='Only training observations used; common coupled normalization cancels in candidate ranking.')
    freeze(output, result)
    return result


def summarize(root, manifest, destination, data):
    from collections import Counter
    selections = {p.stem: json.loads(p.read_text()) for p in root.glob('*selection.json')}
    locked = selections.get('locked-selection')
    summary = dict(completedNativeRuns=len(data), statusCounts=dict(Counter(r['status'] for r in data)),
                   runRoleCounts=dict(Counter(run_role(r) for r in data)),
                   coldTrialStatusCounts=dict(Counter(r['status'] for r in data if run_role(r)=='cold-trial')),
                   sumFullProcessSeconds=sum(r['wallSeconds'] for r in data),
                   summedCostMeaning='sum of per-run process times; concurrent runs overlap, not project elapsed time',
                   selections=selections,
                   batches={p.stem: json.loads(p.read_text()) for p in (root/'batches').glob('*.json')},
                   schedules=[json.loads(p.read_text()) for p in sorted((root/'schedules').glob('*/*/record.json'))],
                   protocolAmendments=json.loads((root/'protocol-amendments.json').read_text())
                        if (root/'protocol-amendments.json').exists() else [],
                   trainingComparisons=[], validation=[], timingReplicates=[], postSelectionDiagnostics=[], physicalReferenceDiagnostics={})
    by_prefix={r['prefix']:r for r in data}
    for schedule in summary['schedules']:
        runs=[by_prefix[s['run']] for s in schedule['steps']]
        schedule['completeNativeProfiles']=all(bool(r['performance']) for r in runs)
        schedule['knownNativeWork']={k:sum(r['performance'].get(k,0) for r in runs) for k in
            ['coupledKrylovIterations','pressureIterations','momentumIterations','lineSearchTrials']}
    for method in METHODS+['coupled']:
        rows=[r for r in data if r['case']=='annulus-l7' and r['method']==method and r['pair']==[1.,1.]
              and (r['batch']=='timing-fine-repeat' or r['batch'].startswith('validation-'))]
        if rows:
            times=[r['wallSeconds'] for r in rows]
            summary['timingReplicates'].append(dict(case='annulus-l7', method=method, times=times,
                medianSeconds=statistics.median(times), minSeconds=min(times), maxSeconds=max(times),
                nativeWork=[{k:r['performance'].get(k) for k in ['coupledKrylovIterations','pressureIterations']} for r in rows]))
    for case in manifest['training']:
        for method in METHODS+['coupled']:
            rows = [r for r in data if r['case']==case and r['method']==method and
                    (r['batch'].startswith('confirm-') or r['batch'].startswith('repeat-'))]
            for pair in sorted({tuple(r['pair']) for r in rows}):
                group = [r for r in rows if tuple(r['pair']) == pair]
                times = [r['wallSeconds'] for r in group]
                summary['trainingComparisons'].append(dict(case=case, method=method, pair=pair, repeats=len(group),
                    statuses=dict(Counter(r['status'] for r in group)), medianSeconds=statistics.median(times),
                    minSeconds=min(times), maxSeconds=max(times),
                    iterations=sorted({r['summary'].get('iterations') for r in group}, key=lambda x: math.inf if x is None else x),
                    chosen=bool(locked and method in METHODS and list(pair)==locked['rankings'][method][0]['pair'])))
    for r in data:
        if r['batch']=='fine-diagnostic':
            summary['postSelectionDiagnostics'].append(dict(case=r['case'], method=r['method'], pair=r['pair'],
                status=r['status'], iterations=r['summary'].get('iterations'), wallSeconds=r['wallSeconds'],
                momentumResidual=r['summary'].get('momentumResidual'), agreement=r['agreement'],
                scope='post-selection diagnostic, excluded from locked ranking and transfer success rate'))
        if r['batch'].startswith('validation-'):
            summary['validation'].append(dict(case=r['case'], method=r['method'], pair=r['pair'],
                status=r['status'], iterations=r['summary'].get('iterations'), wallSeconds=r['wallSeconds'],
                agreement=r['agreement'], failureReason=r['summary'].get('failureReason'),
                prefix=r['prefix']))
    for name, case in manifest['cases'].items():
        prefix = root/'references'/name
        if Path(str(prefix)+'.json').exists():
            reference = json.loads(Path(str(prefix)+'.json').read_text())
            summary['physicalReferenceDiagnostics'][name] = dict(tolerance=reference['tolerance'],
                metrics=laminar.evaluate(case['physical'], prefix, reference),
                scope='same fixed spatial discretization; diagnostic only, not accuracy qualification')
    write(destination/'summary.json', summary)
    return summary


def report(root, manifest, destination):
    destination.mkdir(parents=True, exist_ok=True)
    data = records(root)
    summary = summarize(root, manifest, destination, data)
    # The archive already contains the complete records; keep this convenient
    # aggregate local rather than commit a second multi-megabyte copy.
    write(root / 'results.json', dict(manifest=manifest, runs=data,
        protocolAmendments=summary['protocolAmendments'],
        referenceTolerances={k:v['tolerance'] for k,v in summary['physicalReferenceDiagnostics'].items()}))
    columns = ['batch', 'runRole', 'case', 'partition', 'cells', 'nu', 'speed', 'convection',
               'initialization', 'iterationBudget', 'secondsBudget',
               'method', 'relaxationActive', 'alphaU', 'alphaP', 'repeat', 'status', 'eligible', 'objectiveSeconds', 'wallSeconds',
               'cpuSeconds', 'maxRssKiB',
               'iterations', 'momentumResidual', 'continuity', 'velocityError', 'pressureError', 'fluxError']
    with (destination / 'results.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator='\n')
        writer.writeheader()
        for r in data:
            row = {k: r[k] for k in columns if k in r}
            row.update(runRole=run_role(r), initialization='same-mesh warm' if r.get('initialHashes') else 'cold',
                       iterationBudget=r['budget']['iterations'], secondsBudget=r['budget']['seconds'])
            active=r['method'] in METHODS
            row.update(relaxationActive=active, alphaU=r['pair'][0] if active else None,
                       alphaP=r['pair'][1] if active else None)
            row.update(eligible=r['status']=='eligible',
                       objectiveSeconds=r['wallSeconds'] if r['status']=='eligible' and run_role(r)=='cold-trial' else None)
            row.update({k: manifest['cases'][r['case']][k] for k in ['cells','nu','speed','convection']})
            row['partition']='training' if r['case'] in manifest['training'] else 'holdout'
            row['cpuSeconds']=r['resources'].get('userSeconds',0)+r['resources'].get('systemSeconds',0) if r['resources'] else None
            row['maxRssKiB']=r['resources'].get('maxRssKiB')
            row.update({k: r['summary'].get(k) for k in ['iterations', 'momentumResidual', 'continuity']})
            row.update({k+'Error': (r['agreement'] or {}).get(k) for k in ['velocity', 'pressure', 'flux']})
            writer.writerow(row)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import LogLocator, NullFormatter
    matplotlib.rcParams['svg.hashsalt']='cartmesh2d-relaxation-study-v1'
    def save_figure(fig, stem):
        for ext in ['png','svg']:
            path=destination/(stem+'.'+ext)
            fig.savefig(path,dpi=160,metadata={'Date':None} if ext=='svg' else None)
            if ext=='svg':
                path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
        plt.close(fig)
    grid = [r for r in data if r['batch'] == 'grid']
    if not grid:
        return
    costs = [r['wallSeconds'] for r in grid if r['status'] == 'eligible']
    fig, axes = plt.subplots(2, len(manifest['training']), figsize=(12, 10), squeeze=False, layout='constrained')
    abbreviations = {'time-limit': 'TIME', 'iteration-limit': 'ITER', 'linear-failure': 'LIN',
                     'solver-failure': 'FAIL', 'solver-error': 'ERR', 'different-solution': 'DIFF'}
    for i, method in enumerate(METHODS):
        for j, case in enumerate(manifest['training']):
            ax = axes[i, j]
            z = np.full((len(manifest['alphaP']), len(manifest['alphaU'])), np.nan)
            rows = [r for r in grid if r['method'] == method and r['case'] == case]
            for r in rows:
                x = manifest['alphaU'].index(r['pair'][0]); y = manifest['alphaP'].index(r['pair'][1])
                if r['status'] == 'eligible':
                    z[y, x] = r['wallSeconds']
                label = f"{r['wallSeconds']:.2f}s\n{r['summary']['iterations']} it" if r['status'] == 'eligible' else abbreviations.get(r['status'], 'FAIL')
                colour = ('white' if LogNorm(min(costs), max(costs))(r['wallSeconds']) < .7 else '#15232b') if r['status'] == 'eligible' else '#963333'
                ax.text(x, y, label, ha='center', va='center', fontsize=9, color=colour)
            cmap = plt.colormaps['viridis'].copy(); cmap.set_bad('#eeeeee')
            im = ax.imshow(z, origin='lower', cmap=cmap, norm=LogNorm(min(costs), max(costs)))
            ax.set_xticks(range(len(manifest['alphaU'])), manifest['alphaU'])
            ax.set_yticks(range(len(manifest['alphaP'])), manifest['alphaP'])
            ax.set(xlabel='Velocity relaxation', ylabel='Pressure relaxation', title=f'{case} ({manifest["cases"][case]["cells"]} cells) | {method}')
            ax.scatter([2], [2], s=3300, facecolors='none', edgecolors='#ff8033', linewidths=2)
    fig.colorbar(im, ax=axes, label='Full process wall time (s), discovery runs')
    grid_budget=summary['batches']['grid']['budget']
    fig.suptitle('Measured 5 x 5 grids: native strict convergence + same-solution gate\n'
                 f'Orange: default (0.6, 0.25); grey: did not qualify within {grid_budget["iterations"]} iterations / {grid_budget["seconds"]:g} seconds\n'
                 'LIN: linear failure | ERR: native error | FAIL: line search failure | ITER: iteration limit', fontsize=12)
    save_figure(fig,'parameter-maps')
    validation = summary['validation']
    if validation:
        fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
        for ax, case in zip(axes.flat, manifest['holdout']):
            rows = [r for r in validation if r['case']==case]
            for i, r in enumerate(rows):
                if r['status']=='eligible':
                    ax.barh(i, r['wallSeconds'], color='#316e9a')
                    ax.text(r['wallSeconds']*1.05, i, f"{r['wallSeconds']:.2f}s / {r['iterations']} it", va='center', fontsize=9)
                else:
                    ax.text(.03, i, r['status'], transform=ax.get_yaxis_transform(), va='center', color='#a03636')
            names={'simple':'SIMPLE','simple-consistent':'consistent'}
            ax.set_yticks(range(len(rows)), ['coupled' if r['method']=='coupled' else
                f"{names[r['method']]} {r['pair'][0]:g}/{r['pair'][1]:g}" for r in rows])
            ax.invert_yaxis(); ax.set_xscale('log')
            ax.xaxis.set_major_locator(LogLocator(base=10))
            ax.xaxis.set_minor_formatter(NullFormatter())
            valid_times=[r['wallSeconds'] for r in rows if r['status']=='eligible']
            if valid_times:
                ax.set_xlim(min(valid_times)*.6, max(valid_times)*5)
            ax.set(title=f'{case} | {manifest["cases"][case]["cells"]} cells', xlabel='Full process time (s), log scale')
            ax.grid(axis='x', alpha=.2)
        fig.suptitle('Held-out transfer: locked training choices, same native strict + field gates\nSerial cold starts; one run per condition; failures have no speedup score\nconsistent = simple-consistent; coupled does not use relaxation factors', fontsize=12)
        save_figure(fig,'heldout-costs')
    if any(r['policy']=='adaptive' for r in summary['schedules']):
        adaptive = [r for r in summary['schedules'] if r['policy']=='adaptive']
        fig, axes = plt.subplots(len(adaptive), 1, figsize=(10, max(3, 2.7*len(adaptive))), squeeze=False, layout='constrained')
        for ax, result in zip(axes.flat, adaptive):
            steps=result['steps']; x=np.cumsum([s['chargedIterations'] for s in steps])
            edges=np.r_[0,x]
            ax.step(edges, [s['pair'][0] for s in steps]+[steps[-1]['pair'][0]], where='post', label='velocity factor')
            ax.step(edges, [s['pair'][1] for s in steps]+[steps[-1]['pair'][1]], where='post', label='pressure factor')
            rejected=[i for i,s in enumerate(steps) if not s['accepted']]
            ax.scatter(x[rejected], [steps[i]['pair'][0] for i in rejected], marker='x', color='red', label='rejected segment')
            ax.set(title=f"{result['case']} | {result['method']} | {result['status']} | {result['wallSeconds']:.2f}s total",
                   xlabel='Charged iterations (includes rejected work)', ylabel='Factor', ylim=(0,1.07))
            ax.legend(loc='upper left', ncol=3, fontsize=8); ax.grid(alpha=.2)
        fig.suptitle('Predeclared residual-feedback heuristic through segmented native restarts', fontsize=13)
        save_figure(fig,'adaptive-trajectories')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT / 'outputs/relaxation-study')
    sub = parser.add_subparsers(dest='action', required=True)
    init = sub.add_parser('init')
    init.add_argument('--flow-cli', type=Path, default=ROOT / 'build/cartmesh2d_flow_cli')
    init.add_argument('--mesh-cli', type=Path, default=ROOT / 'build/cartmesh2d_cli')
    run = sub.add_parser('run')
    run.add_argument('--batch', required=True)
    run.add_argument('--cases', nargs='+', default=['training'])
    run.add_argument('--methods', nargs='+', choices=METHODS+['coupled'], default=METHODS)
    run.add_argument('--pairs', nargs='+', help='u,p pairs; omitted means manifest grid')
    run.add_argument('--iterations', type=int, default=3000)
    run.add_argument('--seconds', type=float, default=60)
    run.add_argument('--repeats', type=int, default=1)
    run.add_argument('--workers', type=int, choices=[1, 2], default=2)
    run.add_argument('--reference', action='store_true')
    run.add_argument('--linear-policy', choices=['adaptive', 'strict'], help='explicit reference/diagnostic override')
    run.add_argument('--reference-tolerance', type=float, help='explicit tighter-than-candidate reference override; failures retained')
    selection = sub.add_parser('select')
    selection.add_argument('--batches', nargs='+', required=True)
    selection.add_argument('--output', type=Path, required=True)
    reporting = sub.add_parser('report')
    reporting.add_argument('--destination', type=Path, default=ROOT / 'artifacts/current/relaxation-study')
    bootstrap = sub.add_parser('bootstrap-reference')
    bootstrap.add_argument('--case', required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.action == 'init':
        manifest = initialise(root, args.flow_cli.resolve(), args.mesh_cli.resolve())
        print(json.dumps({k: v['cells'] for k, v in manifest['cases'].items()}))
        return
    manifest = load(root)
    if args.action == 'run':
        cases = [case for value in args.cases for case in manifest.get(value, [value])]
        pairs = [list(map(float, pair.split(','))) for pair in args.pairs] if args.pairs else [
            [u, p] for u in manifest['alphaU'] for p in manifest['alphaP']]
        if any(len(pair) != 2 or any(not 0 < v <= 1 for v in pair) for pair in pairs):
            parser.error('Each pair must have two factors in (0,1]')
        if args.reference and args.methods != ['coupled']:
            parser.error('References must use coupled')
        run_batch(root, manifest, args.batch, cases, args.methods, pairs,
                  dict(iterations=args.iterations, seconds=args.seconds), args.repeats, args.workers,
                  args.reference, args.linear_policy, args.reference_tolerance)
    elif args.action == 'select':
        result = select(root, manifest, args.batches, args.output)
        print(json.dumps({m: rows[:3] for m, rows in result['rankings'].items()}, indent=2))
    elif args.action == 'bootstrap-reference':
        bootstrap_reference(root, manifest, args.case)
    else:
        report(root, manifest, args.destination)


if __name__ == '__main__':
    main()
