#!/usr/bin/env python3
"""Euler scheme validation; --output retains actual fields and compact evidence."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools/verification'))
import verify_euler as euler


def run_suite(cli, root):
    records = []
    def run(mesh, name, scheme, order, extra=(), expected=0):
        prefix = root / name
        command = [str(cli), '--mesh', str(mesh), '--output', str(prefix), '--gas-r', '1',
                   '--flux', scheme, '--order', str(order), *map(str, extra)]
        started = time.perf_counter()
        result = subprocess.run(command, capture_output=True, text=True, timeout=90)
        (root / (name+'.log')).write_text(result.stdout+result.stderr)
        assert result.returncode == expected, (command, result.returncode, result.stdout, result.stderr)
        summary = json.loads(Path(str(prefix)+'.json').read_text())
        assert summary['fluxScheme'] == scheme and summary['order'] == order
        audit = euler.audit(mesh, prefix)
        record = dict(case=name, command=command, binarySha256=hashlib.sha256(cli.read_bytes()).hexdigest(), wallSeconds=time.perf_counter()-started,
                      summary=summary, audit=audit, meshSha256=hashlib.sha256(mesh.read_bytes()).hexdigest())
        records.append(record)
        return prefix, record
    sod_errors = {}
    for n in (50, 100, 200):
        mesh = root / f'sod{n}.solver.cm2d'; euler.rectangle(mesh, n, 4)
        for scheme, order in [('rusanov', 1), ('hllc', 1), ('hllc', 2)]:
            name = f'sod{n}-{scheme}-o{order}'
            prefix, record = run(mesh, name, scheme, order)
            rows = list(csv.DictReader(Path(str(prefix)+'.cells.csv').open()))
            errors = {key: math.fsum(float(r['area'])*abs(float(r[key])-euler.sod(float(r['x']), .2)[k]) for r in rows)/.1
                      for key, k in [('rho', 0), ('u', 1), ('p', 2)]}
            record['areaWeightedL1'] = errors; sod_errors[n, scheme, order] = errors
            print(name, errors, flush=True)
        # Same physical checkpoint, same method and target => identical accepted trajectory.
        if n == 50:
            first, _ = run(mesh, 'hllc2-interrupted', 'hllc', 2, ['--max-steps', '20'], 2)
            resumed, _ = run(mesh, 'hllc2-resumed', 'hllc', 2, ['--restart', str(first)+'.checkpoint'])
            assert Path(str(resumed)+'.checkpoint').read_bytes() == (root/'sod50-hllc-o2.checkpoint').read_bytes()
            # Switching numerical controls is allowed; the final state remains physically bound.
            run(mesh, 'hllc2-switched', 'hllc', 2, ['--restart', str(root/'sod50-rusanov-o1.checkpoint'), '--end-time', '.21'])
    unit_scaling = []
    rho_scale, pressure_scale = 1.225, 101325.
    velocity_scale = math.sqrt(pressure_scale/rho_scale)
    for scheme, order in [('rusanov',1),('hllc',2)]:
        prefix, record = run(root/'sod50.solver.cm2d', f'sod50-SI-{scheme}-o{order}', scheme, order,
                             ['--density',rho_scale,'--pressure',pressure_scale,'--end-time',.2/velocity_scale])
        rows = list(csv.DictReader(Path(str(prefix)+'.cells.csv').open()))
        base = list(csv.DictReader((root/f'sod50-{scheme}-o{order}.cells.csv').open()))
        scales = dict(rho=rho_scale,u=velocity_scale,v=velocity_scale,p=pressure_scale)
        error = max(abs(float(a[k])/scale-float(b[k]))/(1+abs(float(b[k])))
                    for a,b in zip(rows,base) for k,scale in scales.items())
        # Unit changes must only introduce accumulated roundoff, not a different
        # physical trajectory. Per-step binary64 allowance; all flux equations
        # are separately audited with dimensionally covariant normalization.
        assert error < 256*math.ulp(1.)*record['summary']['steps'], ('SI scaling',error)
        record['maximumUnitScalingDifference'] = error
        unit_scaling.append(dict(scheme=scheme,order=order,maximumNormalizedDifference=error))
        print(record['case'], 'unit scaling', error, flush=True)
    for key in ('rho', 'u', 'p'):
        assert sod_errors[200, 'hllc', 2][key] < sod_errors[100, 'hllc', 2][key] < sod_errors[50, 'hllc', 2][key]
        assert sod_errors[200, 'hllc', 2][key] < .75*sod_errors[200, 'rusanov', 1][key], 'same-grid Sod improvement'
    vortex = []
    end, gamma, beta = .5, 1.4, 5.
    for n in (32, 64, 128):
        mesh = root/f'vortex{n}.solver.cm2d'; euler.rectangle(mesh, n, n, 20, 20)
        prefix, record = run(mesh, f'vortex{n}-hllc-o2', 'hllc', 2, ['--case', 'vortex', '--end-time', end])
        rows = list(csv.DictReader(Path(str(prefix)+'.cells.csv').open()))
        errors = [0., 0., 0., 0.]
        for row in rows:
            x, y = float(row['x'])-10-end, float(row['y'])-10-end
            radius2 = x*x+y*y
            factor = beta/(2*math.pi)*math.exp((1-radius2)/2)
            temperature = 1-(gamma-1)*beta*beta/(8*gamma*math.pi**2)*math.exp(1-radius2)
            rho = temperature**(1/(gamma-1))
            exact = [rho, 1-factor*y, 1+factor*x, rho**gamma]
            for k, key in enumerate(['rho', 'u', 'v', 'p']):
                errors[k] += float(row['area'])*abs(float(row[key])-exact[k])/400
        record['areaWeightedL1'] = dict(zip(['rho', 'u', 'v', 'p'], errors))
        record['reference'] = 'analytic isentropic vortex sampled at cell centroids; same centre-sampled initialization, domain 20x20, beta=5, t=0.5'
        vortex.append(errors)
        print(record['case'], record['areaWeightedL1'], flush=True)
    rates = [math.log2(a/b) for a,b in zip(vortex[-2],vortex[-1])]
    # Smooth fields test order, unlike the discontinuous Sod wave. This bounded
    # 64->128 regression requires >1.4 observed L1 order, not a universal error.
    assert all(rate > 1.4 for rate in rates), ('vortex observed orders', rates)
    # A changed reported scheme/stage must be rejected, even if summary says success.
    prefix = root/'sod100-hllc-o2'; mesh = root/'sod100.solver.cm2d'
    for suffix, field in [('.faces.csv', 'energy'), ('.faces.csv', 'hllcFallbackStages'), ('.cells.csv', 'previousRhoE'), ('.history.csv', 'mass')]:
        path = Path(str(prefix)+suffix); original = path.read_text(); rows = list(csv.DictReader(original.splitlines()))
        index = -1 if suffix == '.history.csv' else 10
        rows[index][field] = str(float(rows[index][field])+.1) if field != 'hllcFallbackStages' else '1'
        with path.open('w') as out:
            writer = csv.DictWriter(out, fieldnames=rows[0].keys()); writer.writeheader(); writer.writerows(rows)
        try:
            try: euler.audit(mesh, prefix)
            except ValueError: pass
            else: raise AssertionError('independent audit accepted mutation: '+field)
        finally: path.write_text(original)
    evidence = dict(scope='native ideal-gas inviscid Euler development validation', cases=records,
                    vortexObservedOrders=dict(zip(['rho','u','v','p'],rates)),
                    restartIdentical=True, numericalMethodSwitchChecked=True, independentMutationChecks=4, unitScaling=unit_scaling,
                    limitations=['no viscous stress, heat conduction or turbulence', 'no general external-flow or engineering qualification',
                                 'audit reconstructs the last accepted step and checks complete exported time history, not full-run replay'])
    (root/'validation.json').write_text(json.dumps(evidence, indent=2)+'\n')
    return evidence


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--cli', required=True, type=Path); parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output:
        args.output.mkdir(parents=True, exist_ok=True); run_suite(args.cli.resolve(), args.output.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='cm2d-euler-accuracy-') as d: run_suite(args.cli.resolve(), Path(d))
