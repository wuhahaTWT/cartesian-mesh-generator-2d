#!/usr/bin/env python3
"""Native coupled Picard solves with an explicit native block linear backend.

Python only orchestrates existing native assembly/recovery/residual evaluation.
The optional bounded NumPy backend solves the exported matrix as a reference;
it is not a second PDE discretization or a sparse-LU performance benchmark.
U=L=1 and pressure scale U^2: the existing 1e-9 research iteration criterion
requires both full uncondensed residual and state change. Not a product gate.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mesh')
    parser.add_argument('n', type=int)
    parser.add_argument('problem')
    parser.add_argument('viscosity', type=float)
    parser.add_argument('equation', choices=['ns', 'stokes'])
    parser.add_argument('boundary', choices=['closed', 'open', 'pressure', 'traction', 'pseudo-traction', 'normal-stress'])
    parser.add_argument('output', type=Path)
    parser.add_argument('--backend', choices=['ilu0', 'ic0', 'jacobi', 'dense-reference'], default='ilu0')
    parser.add_argument('--transport', type=Path, default=root / 'build/native-laminar-p1-transport')
    parser.add_argument('--block', type=Path, default=root / 'build/native-laminar-block-precondition')
    parser.add_argument('--fields', type=Path, default=root / 'build/native-laminar-state-compare')
    parser.add_argument('--initial-state', default='zero')
    parser.add_argument('--iterations', type=int, default=40)
    parser.add_argument('--tolerance', type=float, default=1e-9)
    parser.add_argument('--linear-tolerance', type=float, default=1e-13)
    parser.add_argument('--restarts', type=int, default=50)
    parser.add_argument('--quadrature', type=int, choices=range(4, 13), default=6)
    args = parser.parse_args()
    for value in (args.viscosity, args.tolerance, args.linear_tolerance):
        if not math.isfinite(value) or value <= 0:
            parser.error('viscosity and tolerances must be positive and finite')
    if args.iterations < 1 or args.restarts < 1 or args.linear_tolerance > .01:
        parser.error('invalid iteration control')
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    history = []
    controls = [args.mesh, str(args.n), args.problem, str(args.viscosity),
                args.equation, args.boundary, str(args.quadrature)]
    previous = args.initial_state
    record = {'controls': controls, 'backend': args.backend, 'completed': False,
              'initialState': previous, 'iterations': history, 'pid': os.getpid(),
              'iterationTolerance': args.tolerance, 'linearTolerance': args.linear_tolerance,
              'maximumIterations': args.iterations, 'maximumRestarts': args.restarts,
              'lastCompletedIterate': None, 'nativeBinaries': {}, 'sourceSha256': {}}

    def save():
        record['totalSeconds'] = time.monotonic() - start
        tmp = args.output / 'result.json.tmp'
        tmp.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
        tmp.replace(args.output / 'result.json')

    def run_json(command, target, stderr):
        t = time.monotonic()
        result = subprocess.run([str(x) for x in command], text=True, capture_output=True)
        stderr.write(result.stderr)
        stderr.flush()
        target.with_suffix(target.suffix + '.stdout').write_text(result.stdout)
        if result.returncode != 0:
            raise RuntimeError(f'{command[0]} exited {result.returncode}; see {target.name}.stdout and stderr')
        value = json.loads(result.stdout)
        target.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
        return value, time.monotonic() - t

    try:
        binaries = {'transport': args.transport, 'fields': args.fields}
        if args.backend != 'dense-reference':
            binaries['block'] = args.block
        record['nativeBinaries'] = {k: {'path': str(v), 'sha256': sha(v)} for k, v in binaries.items()}
        sources = ['native-laminar-oseen-block.py', 'native-laminar-block-precondition.cpp',
                   'native-laminar-state-compare.cpp', 'native-laminar-open-boundary.cpp', 'native-laminar-p1-transport.cpp',
                   'native-laminar-p1-oseen.cpp', 'native-laminar-p1-stress.cpp',
                   'native-laminar-hybrid-stokes-p1.cpp']
        record['sourceSha256'] = {f: sha(root / 'artifacts/current' / f) for f in sources}
        if previous != 'zero':
            record['initialStateSha256'] = sha(previous)
        if Path(args.mesh).is_file():
            record['meshSha256'] = sha(args.mesh)
        geometry = args.output / 'geometry'
        with (args.output / 'geometry.stderr').open('w') as err:
            record['geometry'], record['geometryProcessSeconds'] = run_json(
                [args.fields, 'geometry', args.mesh, args.n, geometry],
                args.output / 'geometry.json', err)
        save()
        for i in range(args.iterations):
            prefix = args.output / f'iteration-{i:03d}'
            step = {'iteration': i, 'previous': previous}
            history.append(step)
            save()
            native_args = [*controls, previous, str(prefix)]
            with Path(str(prefix) + '.stderr').open('w') as err:
                step['assembly'], step['assemblyProcessSeconds'] = run_json(
                    [args.transport, 'assemble', *native_args], Path(str(prefix) + '.assemble.json'), err)
                save()
                t = time.monotonic()
                if args.backend == 'dense-reference':
                    import numpy as np
                    size = step['assembly']['unknowns']
                    # Hard resource cap for the optional small reference only.
                    # It does not reject geometry or change the native backend.
                    if size > 2175:
                        raise RuntimeError('dense reference limited to 2175 unknowns; use a native backend')
                    z = np.fromfile(str(prefix) + '.entries', dtype=[('i', '<i8'), ('j', '<i8'), ('v', '<f8')])
                    matrix = np.zeros((size, size))
                    np.add.at(matrix, (z['i'], z['j']), z['v'])
                    rhs = np.fromfile(str(prefix) + '.rhs', dtype='<f8')
                    x = np.linalg.solve(matrix, rhs)
                    for _ in range(2):
                        x += np.linalg.solve(matrix, rhs - matrix @ x)
                    residual = rhs - matrix @ x
                    rhs_norm = float(np.linalg.norm(rhs))
                    relative = float(np.linalg.norm(residual)) / rhs_norm if rhs_norm else float(np.linalg.norm(residual))
                    ok = bool(np.all(np.isfinite(x)) and math.isfinite(relative) and relative <= args.linear_tolerance)
                    step['linear'] = {'converged': ok, 'unknowns': size, 'relative_residual': relative,
                                      'maximum_residual': float(np.max(abs(residual))), 'method': 'dense-reference'}
                    x.astype('<f8').tofile(str(prefix) + ('.solution' if ok else '.candidate'))
                    del matrix, z, rhs, x, residual
                    code = 0 if ok else 2
                else:
                    command = [args.block, prefix, str(geometry) + '.cells.csv', args.backend,
                               'gauge' if args.boundary == 'closed' else 'outlet', prefix,
                               str(args.linear_tolerance), str(args.restarts), str(args.viscosity)]
                    result = subprocess.run([str(x) for x in command], text=True, capture_output=True)
                    err.write(result.stderr)
                    err.flush()
                    Path(str(prefix) + '.linear.stdout').write_text(result.stdout)
                    code = result.returncode
                    linear_record = Path(str(prefix) + '.json')
                    if linear_record.exists():
                        step['linear'] = json.loads(linear_record.read_text())
                step['linearProcessSeconds'] = time.monotonic() - t
                step['linearExitCode'] = code
                save()
                if code != 0:
                    record['failure'] = f'linear solve exited {code}; candidate was not recovered or accepted'
                    save()
                    return 2 if code == 2 else 1
                if not step['linear']['converged'] or not Path(str(prefix) + '.solution').is_file():
                    raise RuntimeError('linear success has no confirmed solution')
                step['recovery'], step['recoveryProcessSeconds'] = run_json(
                    [args.transport, 'recover', *native_args], Path(str(prefix) + '.recover.json'), err)
                candidate = str(prefix) + '.recover.state'
                step['nonlinearCheck'], step['checkProcessSeconds'] = run_json(
                    [args.transport, 'check', *controls, candidate, prefix], Path(str(prefix) + '.check.json'), err)
                check = step['nonlinearCheck']
                metrics = [check['internalResidual'], check['freeFaceResidual'], check['maxDivergence'],
                           step['recovery']['stateCoefficientChange']]
                if not all(math.isfinite(v) and v >= 0 for v in metrics):
                    raise RuntimeError('invalid physical iteration metric')
                step['equationResidualMax'] = max(metrics[:3])
                step['stateSha256'] = sha(candidate)
                previous = candidate
                record['lastCompletedIterate'] = candidate
                save()
                print(json.dumps({'iteration': i, 'residual': step['equationResidualMax'], 'change': metrics[3],
                                  'linearSeconds': step['linearProcessSeconds'],
                                  'krylovProducts': step['linear'].get('krylov_products')}), flush=True)
                if step['equationResidualMax'] <= args.tolerance and (args.equation == 'stokes' or metrics[3] <= args.tolerance):
                    record['completed'] = True
                    record['finalState'] = candidate
                    save()
                    return 0
        record['failure'] = 'nonlinear iteration budget exhausted; no final acceptance'
        save()
        return 2
    except KeyboardInterrupt:
        record['failure'] = 'interrupted; no final acceptance'
        save()
        return 130
    except Exception as exc:
        record['failure'] = str(exc)
        save()
        raise


if __name__ == '__main__':
    raise SystemExit(main())
