#!/usr/bin/env python3
"""Serial, resource-bounded physical-time runs on one existing final mesh.

Every command, binary/mesh hash and failed invocation is retained. This driver
checks discrete balance; a completed run is not an accuracy qualification.
Use separate directories to compare time steps without overwriting evidence.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import time

import verify_transient_flow as audit


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cli', type=Path, default=Path('build/cartmesh2d_flow_cli'))
    p.add_argument('--mesh', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--case', choices=['taylor-green', 'channel', 'cavity', 'external'], default='taylor-green')
    p.add_argument('--dt', type=float, nargs='+', default=[.04, .02, .01, .005])
    p.add_argument('--end-time', type=float, default=.2)
    p.add_argument('--nu', type=float, default=.1)
    p.add_argument('--speed', type=float, default=1)
    p.add_argument('--tolerance', type=float, default=1e-9)
    p.add_argument('--max-iterations', type=int, default=1000)
    p.add_argument('--timeout', type=float, default=180)
    a = p.parse_args()
    for value in [a.end_time, a.nu, a.speed, a.tolerance, a.timeout, *a.dt]:
        if not math.isfinite(value) or value <= 0:
            p.error('physical parameters and resource budgets must be positive and finite')
    if a.max_iterations < 10 or len(set(a.dt)) != len(a.dt):
        p.error('need at least 10 inner iterations and distinct time steps')
    for dt in a.dt:
        if not math.isclose(a.end_time/dt, round(a.end_time/dt), abs_tol=1e-10, rel_tol=0) or round(a.end_time/dt) < 1:
            p.error('end-time must be a positive integer multiple of each dt')
    a.output.mkdir(parents=True, exist_ok=False)
    report = {'scope': 'fixed time, discrete-balance verification; no engineering qualification',
              'mesh': str(a.mesh.resolve()), 'meshSha256': audit.native.sha256_file(a.mesh),
              'binarySha256': audit.native.sha256_file(a.cli), 'runs': []}
    for dt in a.dt:
        prefix = a.output / ('dt-' + str(dt)) / 'result'
        prefix.parent.mkdir()
        command = [str(a.cli.resolve()), '--mesh', str(a.mesh.resolve()), '--output', str(prefix.resolve()),
                   '--case', a.case, '--nu', str(a.nu), '--speed', str(a.speed), '--convection', 'limited-linear',
                   '--time-step', str(dt), '--steps', str(round(a.end_time/dt)), '--tolerance', str(a.tolerance),
                   '--max-iterations', str(a.max_iterations), '--profile']
        run = {'command': command, 'dt': dt, 'endTime': a.end_time, 'timeoutSeconds': a.timeout}
        start = time.monotonic()
        with Path(str(prefix)+'.stdout.log').open('w') as stdout, Path(str(prefix)+'.stderr.log').open('w') as stderr:
            try:
                completed = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=a.timeout)
                run['returnCode'] = completed.returncode
            except subprocess.TimeoutExpired:
                run['returnCode'] = None
                run['issue'] = 'timeout; process terminated; accepted checkpoint retained if available'
        run['elapsedSeconds'] = time.monotonic()-start
        run['valid'] = False
        if run['returnCode'] == 0:
            try:
                result = audit.verify(a.mesh, prefix, Path(str(prefix)+'.audit.json'))
                run.update(valid=True, counts=result['counts'], analyticErrors=result['analyticErrors'],
                           computedMaxCourant=result['computedMaxCourant'],
                           totalInnerIterations=sum(r['innerIterations'] for r in result['history']))
            except (ValueError, OSError, audit.native.VerificationError) as exc:
                run['issue'] = str(exc)
                audit.native.write_json(Path(str(prefix)+'.audit.json'), {'valid': False, 'issues': [str(exc)]})
        report['runs'].append(run)
        audit.native.write_json(a.output/'runs.json', report)
        print(json.dumps(run), flush=True)
        if not run['valid']:
            return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
