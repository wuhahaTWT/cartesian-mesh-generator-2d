#!/usr/bin/env python3
"""Drive native mesher/solver and check output/lifecycle contracts, not PDEs."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import queue
import signal
import subprocess
import tempfile
import threading

parser = argparse.ArgumentParser()
parser.add_argument('--cli', required=True)
parser.add_argument('--mesh-cli', required=True)
parser.add_argument('--keep', help='Fresh evidence directory; otherwise clean only after success.')
args = parser.parse_args()
cli = str(Path(args.cli).resolve())
mesher = str(Path(args.mesh_cli).resolve())
root = Path(args.keep) if args.keep else Path(tempfile.mkdtemp(prefix='cartmesh-compatible-cli-'))
if args.keep:
    root.mkdir(parents=True, exist_ok=False)


def read(prefix, suffix='.summary.json'):
    return json.loads(Path(str(prefix) + suffix).read_text())


def snapshot(prefix):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.glob(prefix.name + '.*') if p.is_file()}


def maximum_state_difference(left, right):
    return max(abs(a - b) for key in ('cells', 'faces')
               for lrow, rrow in zip(left[key], right[key])
               for a, b in zip(lrow, rrow))


def run(label, case='cavity', extra=(), code=0, compatible=True, nu='.1'):
    prefix = root / label
    cmd = [cli, '--mesh', str(mesh), '--output', str(prefix), '--case', case, '--nu', nu]
    if compatible:
        cmd += ['--discretization', 'compatible']
    cmd += list(extra)
    result = subprocess.run(cmd, text=True, capture_output=True, timeout=35)
    Path(str(prefix) + '.stdout').write_text(result.stdout)
    Path(str(prefix) + '.stderr').write_text(result.stderr)
    assert result.returncode == code, (label, result.returncode, result.stderr, result.stdout[-1000:])
    return prefix


def exported(prefix, expected):
    summary = read(prefix)
    assert summary['status'] == expected and summary['exportsComplete']
    assert summary['converged'] == (expected == 'converged')
    assert not summary['physicalAccuracyQualified'] and summary['contextBoundCheckpoint']
    assert summary['checkpointAvailable'] and Path(str(prefix) + '.checkpoint').is_file()
    assert summary['lastAcceptedAvailable']
    accepted, field = read(prefix, '.accepted.json'), read(prefix, '.fields.json')
    assert accepted['kind'] == 'accepted-iterate'
    assert accepted['cells'] == field['cells'] and accepted['faces'] == field['faces']
    assert len(accepted['cells']) == summary['cells'] and len(accepted['faces']) == summary['faces']
    assert all(len(c) == 9 for c in accepted['cells'])
    assert all(len(f) == 4 for f in accepted['faces'])
    if expected == 'cancelled':
        assert not summary['boundaryLoadsAvailable'] and not Path(str(prefix) + '.loads.json').exists()
    else:
        assert summary['boundaryLoadsAvailable']
        loads = read(prefix, '.loads.json')
        assert loads['sourceState'] == 'last-accepted'
        assert loads['acceptedIterations'] == summary['totalAcceptedIterations']
        assert loads['numericallyConverged'] == summary['numericallyConverged']
        assert not loads['physicalAccuracyQualified']
        assert loads['units']['forceAndTractionMoments'] == 'm^3/s^2'
        assert loads['boundaries'] and all(len(b['tractionMoments']) == 2 for b in loads['boundaries'])
        assert not Path(str(prefix) + '.loads.json.tmp').exists()
    if summary['resumed']:
        assert not Path(str(prefix) + '.seed.json').exists()
    else:
        assert read(prefix, '.seed.json')['kind'] == 'seed'
    with Path(str(prefix) + '.cells.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == summary['cells']
    for c, v in zip(rows, accepted['cells']):
        assert [float(c[k]) for k in ('u', 'uX', 'uY', 'v', 'vX', 'vY', 'p', 'pX', 'pY')] == v
    with Path(str(prefix) + '.faces.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == summary['faces']
    for f, v in zip(rows, accepted['faces']):
        assert [float(f[k]) for k in ('u0', 'us', 'v0', 'vs')] == v
    vtk = Path(str(prefix) + '.vtk').read_text()
    assert vtk.count('CELL_DATA ') == 1
    assert f'CELL_DATA {summary["cells"]}\n' in vtk
    assert vtk.count('VECTORS velocity double') == 1
    return summary


try:
    outline = root / 'square.xy'
    outline.write_text('0 0\n1 0\n1 1\n0 1\n')
    mesh_prefix = root / 'square'
    generated = subprocess.run([mesher, str(outline), str(mesh_prefix), '3', '.25', '.1',
                                'interior', str(root / 'foam'), '3'],
                               capture_output=True, text=True, timeout=20)
    (root / 'mesh.stdout').write_text(generated.stdout)
    (root / 'mesh.stderr').write_text(generated.stderr)
    assert generated.returncode == 0, generated.stderr
    mesh = root / 'square.solver.cm2d'

    cavity = run('cavity')
    solved = exported(cavity, 'converged')
    assert solved['controls']['linearInitialGuess'] == 'zero'
    assert solved['controls']['pressureInverse'] == 'viscous-mass'
    assert solved['controls']['initialPseudoStep'] == .1 and solved['controls']['maximumLinearRestarts'] == 50
    with Path(str(cavity) + '.residuals.csv').open() as stream:
        cavity_history = list(csv.DictReader(stream))
    assert cavity_history and all('linearInitialRelativeResidual' in row for row in cavity_history)
    explicit_zero = run('explicit-zero', extra=('--linear-initial-guess', 'zero'))
    zero_summary = exported(explicit_zero, 'converged')
    assert zero_summary['controls']['linearInitialGuess'] == 'zero'
    assert Path(str(explicit_zero) + '.checkpoint').read_bytes() == Path(str(cavity) + '.checkpoint').read_bytes()
    assert read(explicit_zero, '.loads.json') == read(cavity, '.loads.json')
    current = run('current-state', extra=('--linear-initial-guess', 'current-state'))
    current_summary = exported(current, 'converged')
    assert current_summary['controls']['linearInitialGuess'] == 'current-state'
    # Algebraically equivalent initial guesses may change the final rounded
    # iterate. This is a numerical-regression allowance, not a CFD gate.
    assert maximum_state_difference(read(current, '.accepted.json'), read(cavity, '.accepted.json')) < 1e-10
    assert max(abs(a - b) for a, b in zip(read(current, '.loads.json')['boundaryTraction'],
                                         read(cavity, '.loads.json')['boundaryTraction'])) < 1e-10
    budget = run('budget', extra=('--max-iterations', '1'), code=2)
    limited = exported(budget, 'nonlinear-budget')
    assert limited['acceptedIterations'] == 1 and not limited['numericallyConverged']
    assert read(budget, '.accepted.json')['cells'] != read(budget, '.seed.json')['cells']
    linear = run('linear-budget', extra=('--linear-restarts', '1', '--krylov-directions', '1'), code=2)
    assert read(linear)['status'] == 'linear-budget' and not read(linear)['lastAcceptedAvailable']
    assert not Path(str(linear) + '.fields.json').exists() and not Path(str(linear) + '.accepted.json').exists()
    assert read(linear, '.seed.json')['kind'] == 'seed'
    assert not read(linear)['boundaryLoadsAvailable'] and not Path(str(linear) + '.loads.json').exists()

    # Same pressure-referenced channel, unchanged equation and default cold
    # controls. Switching only the inverse must preserve full fields and loads.
    schur = ('--compatible-pressure-inverse', 'schur')
    aggregation = ('--compatible-pressure-inverse', 'schur-aggregation')
    channel = run('channel-schur', case='channel', extra=schur)
    exported(channel, 'converged')
    channel_aggregation = run('channel-aggregation', case='channel', extra=aggregation)
    channel_budget = run('channel-budget', case='channel', extra=schur + ('--max-iterations', '1'), code=2)
    channel_limited = exported(channel_budget, 'nonlinear-budget')
    assert channel_limited['totalAcceptedIterations'] == 1
    checkpoint_before = Path(str(channel_budget) + '.checkpoint').read_bytes()
    aggregation_resume = run('aggregation-resume', case='channel',
                             extra=aggregation + ('--restart', str(channel_budget) + '.checkpoint'))
    for candidate in (channel_aggregation, aggregation_resume):
        report = exported(candidate, 'converged')
        assert report['controls']['pressureInverse'] == 'diagonal-schur-aggregation'
        assert report['controls']['linearInitialGuess'] == 'zero'
        assert report['controls']['initialPseudoStep'] == .1 and report['controls']['maximumLinearRestarts'] == 50
        # Reuse the initial-guess roundoff allowance for another algebraically
        # equivalent solve; no new spatial/physical acceptance threshold.
        assert maximum_state_difference(read(candidate, '.accepted.json'), read(channel, '.accepted.json')) < 1e-10
        left, right = read(candidate, '.loads.json'), read(channel, '.loads.json')
        assert left['absolutePressureReference'] and right['absolutePressureReference']
        assert max(abs(a-b) for a,b in zip(left['boundaryTraction'], right['boundaryTraction'])) < 1e-10
        assert len(left['boundaries']) == len(right['boundaries'])
        for a,b in zip(left['boundaries'], right['boundaries']):
            assert a['face'] == b['face']
            for key in ('tractionMoments', 'pressureMoments'):
                assert max(abs(x-y) for am,bm in zip(a[key],b[key]) for x,y in zip(am,bm)) < 1e-10
    assert read(aggregation_resume)['resumed'] and read(aggregation_resume)['acceptedIterationsBefore'] == 1
    aggregation_failed = run('aggregation-linear-budget', case='channel',
                             extra=aggregation + ('--linear-restarts', '1', '--krylov-directions', '1'), code=2)
    assert read(aggregation_failed)['status'] == 'linear-budget'
    assert not read(aggregation_failed)['checkpointAvailable'] and not read(aggregation_failed)['lastAcceptedAvailable']
    assert read(aggregation_failed, '.seed.json')['kind'] == 'seed'
    for suffix in ('.checkpoint', '.accepted.json', '.fields.json', '.loads.json'):
        assert not Path(str(aggregation_failed) + suffix).exists()
    aggregation_failed_resume = run('aggregation-failed-resume', case='channel',
                                    extra=aggregation + ('--restart', str(channel_budget) + '.checkpoint',
                                                         '--linear-restarts', '1', '--krylov-directions', '1'), code=2)
    retained = exported(aggregation_failed_resume, 'linear-budget')
    assert retained['acceptedIterations'] == 0 and retained['totalAcceptedIterations'] == 1
    assert Path(str(aggregation_failed_resume) + '.checkpoint').read_bytes() == checkpoint_before
    assert Path(str(channel_budget) + '.checkpoint').read_bytes() == checkpoint_before
    closed = run('closed-aggregation', extra=aggregation, code=1)
    assert read(closed)['status'] == 'failed' and not read(closed)['checkpointAvailable']
    assert 'traction pressure reference' in read(closed)['reason']

    custom = run('custom', case='custom', extra=('--boundary', str(cavity) + '.boundaries'))
    exported(custom, 'converged')
    assert read(custom, '.accepted.json')['cells'] == read(cavity, '.accepted.json')['cells']
    assert read(custom, '.accepted.json')['faces'] == read(cavity, '.accepted.json')['faces']

    continued = run('continued', extra=('--restart', str(budget) + '.checkpoint'))
    resumed_summary = exported(continued, 'converged')
    assert resumed_summary['resumed'] and resumed_summary['acceptedIterationsBefore'] == 1
    assert resumed_summary['totalAcceptedIterations'] == solved['acceptedIterations']
    assert Path(str(continued) + '.checkpoint').read_bytes() == Path(str(cavity) + '.checkpoint').read_bytes()
    assert read(continued, '.loads.json') == read(cavity, '.loads.json')
    assert read(continued, '.accepted.json')['cells'] == read(cavity, '.accepted.json')['cells']
    assert read(continued, '.accepted.json')['faces'] == read(cavity, '.accepted.json')['faces']
    continued_current = run('continued-current', extra=('--restart', str(budget) + '.checkpoint',
                                                       '--linear-initial-guess', 'current-state'))
    current_resumed_summary = exported(continued_current, 'converged')
    assert current_resumed_summary['resumed'] and current_resumed_summary['controls']['linearInitialGuess'] == 'current-state'
    assert maximum_state_difference(read(continued_current, '.accepted.json'), read(cavity, '.accepted.json')) < 1e-10
    failed_resume = run('failed-resume', extra=('--restart', str(budget) + '.checkpoint',
                                              '--linear-restarts', '1', '--krylov-directions', '1'), code=2)
    failed_summary = exported(failed_resume, 'linear-budget')
    assert failed_summary['acceptedIterations'] == 0 and failed_summary['totalAcceptedIterations'] == 1
    assert Path(str(failed_resume) + '.checkpoint').read_bytes() == Path(str(budget) + '.checkpoint').read_bytes()
    mismatched = run('mismatch', extra=('--restart', str(budget) + '.checkpoint'), nu='.2', code=1)
    assert read(mismatched)['status'] == 'failed' and not read(mismatched)['checkpointAvailable']
    assert not Path(str(mismatched) + '.accepted.json').exists()
    seed_as_restart = run('seed-as-restart', extra=('--restart', str(budget) + '.seed.json'), code=1)
    assert not Path(str(seed_as_restart) + '.summary.json').exists()
    corrupt = root / 'corrupt.checkpoint'
    text = Path(str(budget) + '.checkpoint').read_text()
    index = text.index('STATE\n') + len('STATE\n') + 15
    corrupt.write_text(text[:index] + ('1' if text[index] == '0' else '0') + text[index+1:])
    damaged = run('corrupt-read', extra=('--restart', str(corrupt)), code=1)
    assert not Path(str(damaged) + '.summary.json').exists()

    # New opt-in must not alter the original default path.
    default = run('default', compatible=False, extra=('--max-iterations', '1500',))
    collocated = run('collocated', compatible=False, extra=('--discretization', 'collocated', '--max-iterations', '1500'))
    for suffix in ('.cells.csv', '.faces.csv'):
        assert Path(str(default) + suffix).read_bytes() == Path(str(collocated) + suffix).read_bytes()

    for label, extra in [('unsupported', ('--time-step', '.1')),
                         ('quadrature', ('--quadrature-order', '2')),
                         ('duplicate', ('--nu', '.2')),
                         ('bad-linear-initial', ('--linear-initial-guess', 'previous')),
                         ('bad-pressure-inverse', ('--compatible-pressure-inverse', 'auto')),
                         ('bad-selector', ('--discretization', 'collocated'))]:
        prefix = run(label, extra=extra, code=1)
        assert not Path(str(prefix) + '.summary.json').exists()
    missing = run('missing-boundary', case='custom', extra=('--boundary', str(root / 'missing')), code=1)
    assert not Path(str(missing) + '.summary.json').exists()
    legacy = root / 'legacy.json'
    legacy_text = json.dumps({'converged': False, 'status': 'failed'})
    legacy.write_text(legacy_text)
    run('legacy', code=1)
    assert legacy.read_text() == legacy_text
    assert not (root / 'legacy.summary.json').exists()
    (root / 'protected-load.loads.json').write_text('preserve existing load')
    run('protected-load', code=1)
    assert (root / 'protected-load.loads.json').read_text() == 'preserve existing load'
    before = snapshot(cavity)
    repeat = subprocess.run([cli, '--discretization', 'compatible', '--mesh', str(mesh),
                             '--output', str(cavity), '--case', 'cavity', '--nu', '.1'],
                            capture_output=True, text=True, timeout=10)
    assert repeat.returncode == 1 and 'already exists' in repeat.stderr
    assert snapshot(cavity) == before

    # POSIX signals are the CLI's native cooperative cancellation contract.
    if os.name == 'posix':
        prefix = root / 'cancel'
        cmd = [cli, '--discretization', 'compatible', '--mesh', str(mesh),
               '--output', str(prefix), '--case', 'cavity', '--nu', '.1']
        child = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        lines, events = [], queue.Queue()
        def collect():
            for line in child.stdout:
                lines.append(line)
                if '"type":"compatible-flow-progress"' in line:
                    events.put(line)
        reader = threading.Thread(target=collect, daemon=True)
        reader.start()
        try:
            events.get(timeout=20)
            child.send_signal(signal.SIGINT)
            code = child.wait(timeout=20)
            reader.join(timeout=2)
        finally:
            if child.poll() is None:
                child.kill()
                child.wait()
        Path(str(prefix) + '.stdout').write_text(''.join(lines))
        Path(str(prefix) + '.stderr').write_text(child.stderr.read())
        assert code == 130, (code, lines)
        cancelled = exported(prefix, 'cancelled')
        assert 1 <= cancelled['acceptedIterations'] < solved['acceptedIterations']
        assert not cancelled['converged']
        cancel_resume = run('cancel-resume', extra=('--restart', str(prefix) + '.checkpoint'))
        exported(cancel_resume, 'converged')
        assert Path(str(cancel_resume) + '.checkpoint').read_bytes() == Path(str(cavity) + '.checkpoint').read_bytes()

        # Uncatchable termination must leave the most recently published file.
        killed_prefix = root / 'killed'
        cmd = [cli, '--discretization', 'compatible', '--mesh', str(mesh),
               '--output', str(killed_prefix), '--case', 'cavity', '--nu', '.1']
        child = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        lines, events = [], queue.Queue()
        reader = threading.Thread(target=collect, daemon=True)
        reader.start()
        try:
            events.get(timeout=20)
            child.kill()
            code = child.wait(timeout=20)
            reader.join(timeout=2)
        finally:
            if child.poll() is None:
                child.kill()
                child.wait()
        Path(str(killed_prefix) + '.stdout').write_text(''.join(lines))
        Path(str(killed_prefix) + '.stderr').write_text(child.stderr.read())
        assert code == -signal.SIGKILL
        assert read(killed_prefix)['status'] == 'running' and not read(killed_prefix)['converged']
        recovered = run('kill-resume', extra=('--restart', str(killed_prefix) + '.checkpoint'))
        exported(recovered, 'converged')
        assert Path(str(recovered) + '.checkpoint').read_bytes() == Path(str(cavity) + '.checkpoint').read_bytes()

    print(json.dumps({'passed': True, 'cavityIterations': solved['acceptedIterations'],
                      'cavityCells': solved['cells'], 'signalTested': os.name == 'posix',
                      'aggregationColdAndResume': True, 'aggregationFailureRetention': True,
                      'evidenceDirectory': str(root)}))
except BaseException:
    print(f'Preserved failed CLI evidence: {root}')
    raise
else:
    if not args.keep:
        import shutil
        shutil.rmtree(root)
