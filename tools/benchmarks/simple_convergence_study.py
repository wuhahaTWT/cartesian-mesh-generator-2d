#!/usr/bin/env python3
"""Run retained native SIMPLE/coupled cases at explicit stopping controls.

Python only launches existing native binaries and compares their exported fields.
It does not reconstruct PDEs or grant physical accuracy from a residual value.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fields(prefix):
    path = Path(str(prefix) + '.cells.csv')
    if not path.exists():
        return None
    with path.open(newline='') as stream:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]


def differences(left, right, speed, pressure_scale):
    if not left or not right or len(left) != len(right):
        return None
    for a, b in zip(left, right):
        if any(a[k] != b[k] for k in ('cell', 'x', 'y', 'area')):
            raise RuntimeError('Field comparison requires identical mesh/order')
    area = sum(a['area'] for a in right)
    return {
        'maxVelocityComponentOverUref': max(abs(a[k]-b[k])/speed for a, b in zip(left, right) for k in ('u', 'v')),
        'areaL2VelocityOverUref': (sum(a['area']*((a['u']-b['u'])**2+(a['v']-b['v'])**2) for a, b in zip(left, right))/area)**.5/speed,
        'maxPressureOverReferenceScale': max(abs(a['p']-b['p']) for a, b in zip(left, right))/pressure_scale,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--variants', type=Path, default=ROOT/'outputs/simple-pc-variants')
    parser.add_argument('--binary', type=Path, default=ROOT/'build/Release/cartmesh2d_flow_cli.exe')
    parser.add_argument('--suite', choices=['simple', 'gpu'], default='simple')
    parser.add_argument('--seconds', type=float, default=120)
    args = parser.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    binary = args.binary.resolve(); variants = args.variants.resolve()
    cases = {
        'annulus-small': ('annulus-l5', .1, .5, 'face-limited-linear'),
        'annulus-fine': ('annulus-l7', .1, .5, 'face-limited-linear'),
        'dfg': ('dfg20-l7', .001, .3, 'limited-linear'),
    }
    runs = []
    if args.suite == 'simple':
        # The source ablations use exactly the original upwind study controls.
        for variant in ('baseline', 'velocity', 'pressure', 'both'):
            runs.append((f'dfg-ablation-{variant}', 'dfg', binary if variant == 'baseline' else variants/f'flow-{variant}.exe',
                         ['--convection', 'upwind'], 1e-8, 'simple'))
        for tol in (.05, .01, .001, 1e-6, 1e-8):
            runs.append((f'dfg-research-tol-{tol:g}', 'dfg', variants/'flow-both.exe', [], tol, 'simple'))
        for case in ('annulus-small', 'annulus-fine', 'dfg'):
            runs.append((f'{case}-product-simple', case, binary, [], 1e-8, 'simple'))
            runs.append((f'{case}-product-coupled', case, binary, ['--linear-policy', 'adaptive'], 1e-8, 'coupled'))
        runs += [
            ('dfg-research-anderson', 'dfg', variants/'flow-both.exe', ['--steady-acceleration', 'anderson'], 1e-8, 'simple'),
            ('dfg-research-anderson-adaptive', 'dfg', variants/'flow-both.exe', ['--steady-acceleration', 'anderson', '--linear-policy', 'adaptive'], 1e-8, 'simple'),
            ('dfg-research-coupled', 'dfg', variants/'flow-both.exe', ['--linear-policy', 'adaptive'], 1e-8, 'coupled'),
        ]
    else:
        for case in ('annulus-small',):
            for pc in ('jacobi', 'ic0', 'aggregation'):
                runs.append((f'{case}-cpu-{pc}', case, binary, ['--pressure-preconditioner', pc], 1e-6, 'simple'))
                if pc != 'aggregation':
                    runs.append((f'{case}-gpu-{pc}', case, binary,
                                 ['--pressure-preconditioner', pc, '--pressure-backend', 'torch-cuda'], 1e-6, 'simple'))
            for backend in ('cpu', 'torch-cuda'):
                runs.append((f'{case}-coupled-{backend}', case, binary,
                             ['--pressure-backend', backend, '--linear-policy', 'adaptive'], 1e-8, 'coupled'))
    result = {
        'format': 'cartmesh2d-simple-convergence-study-v1',
        'head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'dirtyStatusAtStart': subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True),
        'scope': 'native stopping/iteration comparison on identical retained inputs; no new spatial accuracy qualification',
        'serialRuns': True, 'budgetSecondsPerRun': args.seconds, 'runs': [],
        'startupPolicy': 'one --help warmup per executable; every timed CFD run is a fresh process including GPU initialization',
        'sourceSha256': {str(p): sha(ROOT/p) for p in (
            Path('src/fv/Incompressible2D.cpp'), Path('src/fv/FlowCoupled2D.cpp'),
            Path('src/fv/TorchPressureModule2D.cpp'), Path('include/cartmesh2d/fv/detail/FlowConvergence2D.hpp'),
            Path('tools/benchmarks/simple_convergence_study.py'))},
        'gpuModuleSha256': sha(ROOT/'build/Release/cartmesh2d_torch_pressure.dll') if (ROOT/'build/Release/cartmesh2d_torch_pressure.dll').exists() else None,
        'torchRootEnvironment': os.environ.get('CARTMESH2D_TORCH_ROOT'),
    }
    record = out/'summary.json'

    def save():
        record.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')

    for executable in dict.fromkeys(run[2] for run in runs):
        subprocess.run([str(executable), '--help'], cwd=ROOT, stdout=subprocess.DEVNULL,
                       stderr=subprocess.PIPE, check=True, timeout=max(120, args.seconds))
    for name, case, executable, extra, tolerance, coupling in runs:
        stem, nu, speed, scheme = cases[case]
        mesh = ROOT/'outputs/laminar-inputs'/f'{stem}.solver.cm2d'
        boundary = mesh.with_name(f'{stem}.boundaries')
        prefix = out/name
        command = [str(executable), '--mesh', str(mesh), '--boundary', str(boundary), '--case', 'custom',
                   '--nu', str(nu), '--speed', str(speed), '--convection', scheme, '--coupling', coupling,
                   '--velocity-relaxation', '.6', '--pressure-relaxation', '.25', '--pressure-corrections', '4',
                   '--linear-policy', 'strict', '--tolerance', str(tolerance), '--max-iterations', '5000',
                   '--profile', '--output', str(prefix), *extra]
        if args.suite == 'gpu':
            command += ['--pressure-snapshot', str(prefix)+'.pressure.bin']
        run = {'name': name, 'case': case, 'command': command, 'binarySha256': sha(executable),
               'meshSha256': sha(mesh), 'boundarySha256': sha(boundary), 'tolerance': tolerance,
               'changesSpatialDiscretization': executable.name in ('flow-pressure.exe', 'flow-both.exe')}
        result['runs'].append(run); save()
        start = time.perf_counter()
        with Path(str(prefix)+'.log').open('w', encoding='utf-8') as stream:
            try:
                run['exitCode'] = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                                                 timeout=args.seconds).returncode
            except subprocess.TimeoutExpired:
                run['exitCode'] = 'timeout'
        run['wallSeconds'] = time.perf_counter()-start
        for key, suffix in [('summary', '.json'), ('performance', '.performance.json')]:
            path = Path(str(prefix)+suffix)
            if path.exists():
                run[key] = json.loads(path.read_text(encoding='utf-8'))
        run['accepted'] = run['exitCode'] == 0 and run.get('summary', {}).get('converged') is True
        run['prefix'] = str(prefix)
        save(); print(name, run['exitCode'], run.get('summary', {}).get('iterations'), f"{run['wallSeconds']:.3f}s", flush=True)
    lookup = {r['name']: r for r in result['runs']}
    for run in result['runs']:
        if run['accepted'] and run['case'] == 'dfg':
            load = next((x for x in run['summary'].get('namedWallLoads', []) if x['name'] == 'cylinder'), None)
            if load is not None:
                # Retained DFG Re20 normalization: Umean=.2, diameter=.1.
                run['cylinderLoad'] = {'Cd': load['forceX']/.002, 'Cl': load['forceY']/.002}
    if args.suite == 'simple':
        reference = lookup['dfg-research-tol-1e-08']
        if reference['accepted']:
            for run in result['runs']:
                if run['name'].startswith('dfg-research-') and run['accepted']:
                    s = run['summary']; scale = .3**2+.001*.3/s['domainHeight']
                    run['differenceFromSameResearchStrictSolution'] = differences(fields(run['prefix']), fields(reference['prefix']), .3, scale)
                    run['wallForceDifference'] = {k: s[k]-reference['summary'][k] for k in ('wallForceX', 'wallForceY')}
        for case in ('annulus-small', 'annulus-fine', 'dfg'):
            left, right = lookup[f'{case}-product-simple'], lookup[f'{case}-product-coupled']
            if left['accepted'] and right['accepted']:
                _, nu, speed, _ = cases[case]
                right['differenceFromSameProductSIMPLE'] = differences(fields(right['prefix']), fields(left['prefix']),
                    speed, speed**2+nu*speed/right['summary']['domainHeight'])
    else:
        for run in result['runs']:
            if '-gpu-' in run['name'] and run['accepted']:
                ref = lookup[run['name'].replace('-gpu-', '-cpu-')]
                if ref['accepted']:
                    run['differenceFromMatchedCPU'] = differences(fields(run['prefix']), fields(ref['prefix']), .5, .25+.05/run['summary']['domainHeight'])
                    run['fullCostSpeedupCPUoverGPU'] = ref['wallSeconds']/run['wallSeconds']
    save()


if __name__ == '__main__':
    main()
