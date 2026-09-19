#!/usr/bin/env python3
"""Independent, deliberately restricted SST-2003m homogeneous-decay audit.

This validates a prescribed F1=1, zero-strain, zero-gradient ODE limit of the
transport block on actual CM2D cells. It is NOT a wall-bounded RANS benchmark.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import time
import verify_native_flow as native


def require(ok, message):
    if not ok:
        raise ValueError(message)


def finite(value):
    x = float(value)
    require(math.isfinite(x), 'nonfinite evidence')
    return x


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(mesh_path, prefix, dt, steps):
    require(math.isfinite(dt) and dt > 0 and isinstance(steps, int) and 0 < steps <= 100,
            'invalid time controls')
    mesh = native.read_cm2d(mesh_path)
    measured = native.measure(mesh, 1e-10, 1e-9)
    require(not measured.issues, str(measured.issues))
    def rows(suffix):
        with Path(str(prefix)+suffix).open() as f:
            return list(csv.DictReader(f))
    history, cells, faces = rows('.history.csv'), rows('.cells.csv'), rows('.faces.csv')
    require(len(history) == steps and len(cells) == len(mesh.cells) and len(faces) == len(mesh.edges),
            'incomplete output')
    k, w = .02, 4.
    max_history_error = 0.
    for i, row in enumerate(history):
        previous_k, previous_w = k, w
        w = 2*w/(1+math.sqrt(1+4*dt*.075*w))
        k /= 1+dt*.09*w
        require(int(row['step']) == i+1 and 1 <= int(row['outerIterations']) <= 80, 'invalid accepted step')
        require(math.isclose(finite(row['time']), (i+1)*dt, abs_tol=1e-14) and
                math.isclose(finite(row['dt']), dt, abs_tol=1e-14), 'mixed physical times')
        err = max(abs(finite(row['k'])-k), abs(finite(row['omega'])-w))
        max_history_error = max(max_history_error, err)
        require(err < 1e-8 and 0 <= finite(row['maxNonlinearRateResidual']) < 1e-9,
                'nonlinear accepted-step mismatch')
    rk, rw = [], []
    final_k, final_w = [], []
    for i, row in enumerate(cells):
        require(int(row['cell']) == i, 'cell order mismatch')
        a = measured.areas[i]
        for key, exact in zip(('area', 'x', 'y'), (a, *measured.centroids[i])):
            require(math.isclose(finite(row[key]), exact, rel_tol=1e-9, abs_tol=1e-11), 'geometry mismatch')
        kv, wv, pk, pw = map(finite, (row['k'], row['omega'], row['previousK'], row['previousOmega']))
        require(kv >= 0 and wv > 0, 'inadmissible turbulence state')
        require(abs(kv-k) < 1e-8 and abs(wv-w) < 1e-8 and abs(pk-previous_k) < 1e-8 and
                abs(pw-previous_w) < 1e-8, 'field or previous state differs from discrete ODE root')
        rk.append(a*((kv-pk)/dt+.09*wv*kv))
        rw.append(a*((wv-pw)/dt+.075*wv*wv))
        final_k.append(kv); final_w.append(wv)
    max_flux = 0.
    for i, row in enumerate(faces):
        require(int(row['face']) == i, 'face order mismatch')
        qk, qw = finite(row['kFlux']), finite(row['omegaFlux'])
        max_flux = max(max_flux, abs(qk), abs(qw))
        require(abs(qk) < 1e-10 and abs(qw) < 1e-10, 'homogeneous zero-gradient face flux mismatch')
        edge = mesh.edges[i]
        rk[edge.owner] += qk; rw[edge.owner] += qw
        if edge.neighbour >= 0:
            rk[edge.neighbour] -= qk; rw[edge.neighbour] -= qw
    max_rate = max(max(abs(x)/a for x, a in zip(rk, measured.areas)),
                   max(abs(x)/a for x, a in zip(rw, measured.areas)))
    require(max_rate < 2e-9, 'original nonlinear PDE residual failed')
    final_time = dt*steps
    exact_w = 4/(1+.075*4*final_time)
    exact_k = .02*(1+.075*4*final_time)**(-.09/.075)
    return dict(valid=True, scope='homogeneous-decay-only', cells=len(cells), faces=len(faces),
                mesh=str(mesh_path), meshSha256=sha(mesh_path), prefix=str(prefix), timeStep=dt, steps=steps,
                time=final_time, maxAcceptedHistoryError=max_history_error, maxNonlinearRateResidual=max_rate,
                globalKBalance=math.fsum(rk), globalOmegaBalance=math.fsum(rw), maxFaceFlux=max_flux,
                kContinuousRms=native.weighted_l2((x-exact_k for x in final_k), measured.areas),
                omegaContinuousRms=native.weighted_l2((x-exact_w for x in final_w), measured.areas),
                exactK=exact_k, exactOmega=exact_w,
                files={suffix: sha(str(prefix)+suffix) for suffix in ('.cells.csv', '.faces.csv', '.history.csv')})


def study(probe, mesh, output):
    require(not output.exists(), 'output exists; preserve previous evidence')
    output.mkdir(parents=True)
    results = []
    for steps in (4, 8, 16):
        require(shutil.disk_usage(output).free >= 1.5*1024**3, 'disk floor 1.5 GiB reached')
        dt = .8/steps
        prefix = output/f'n{steps}'
        command = [str(probe), str(mesh), str(prefix), str(dt), str(steps)]
        start = time.monotonic()
        record = dict(command=command, timeoutSeconds=90, status='started')
        record_path = Path(str(prefix)+'.run.json')
        record_path.write_text(json.dumps(record, indent=2)+'\n')
        try:
            run = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
        except subprocess.TimeoutExpired as error:
            output_text = error.stdout or ''
            if isinstance(output_text, bytes):
                output_text = output_text.decode(errors='replace')
            Path(str(prefix)+'.log').write_text(output_text)
            record.update(status='timeout', wallSeconds=time.monotonic()-start)
            record_path.write_text(json.dumps(record, indent=2)+'\n')
            raise
        Path(str(prefix)+'.log').write_text(run.stdout)
        record.update(status='completed' if run.returncode == 0 else 'failed',
                      returncode=run.returncode, wallSeconds=time.monotonic()-start)
        record_path.write_text(json.dumps(record, indent=2)+'\n')
        require(run.returncode == 0, f'probe failed, retained {prefix}.log')
        result = audit(mesh, prefix, dt, steps)
        result.update(command=command, wallSeconds=time.monotonic()-start)
        Path(str(prefix)+'.audit.json').write_text(json.dumps(result, indent=2)+'\n')
        results.append(result)
    orders = {}
    for field in ('kContinuousRms', 'omegaContinuousRms'):
        errors = [r[field] for r in results]
        order = [math.log(a/b, 2) for a, b in zip(errors, errors[1:])]
        # Analytic BE first-order ODE limit: this bounded asymptotic window is
        # fixed before the study, not a general SST engineering quality gate.
        require(all(.8 < x < 1.2 for x in order), 'first-order decay time trend failed')
        orders[field] = order
    result = dict(valid=True, definition='SST-2003m rho=1, SI; F1=1 homogeneous zero-shear decay',
                  notValidated=['wall distance', 'wall conditions', 'velocity-pressure coupling',
                                'spatial SST gradients', 'engineering turbulence'],
                  source='https://tmbwg.github.io/turbmodels/sst.html#sst-2003',
                  probeSha256=sha(probe), runs=results, temporalOrders=orders)
    (output/'study.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(study(args.probe.resolve(), args.mesh.resolve(), args.output.resolve()), indent=2))
