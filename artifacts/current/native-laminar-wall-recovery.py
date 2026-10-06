#!/usr/bin/env python3
"""Native wall recovery study. Python orchestrates and solves exported matrices only."""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, help='new output directory or new case prefix within it')
    parser.add_argument('--case', help='run one named control; default runs all 40')
    parser.add_argument('--stress-binary', type=Path, default=ROOT/'build/native-laminar-p1-stress')
    parser.add_argument('--wall-binary', type=Path, default=ROOT/'build/native-laminar-wall-recovery')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    cases = []
    for mesh, problem in [('square','noslip'), ('sheared','noslip-sheared')]:
        for n in [2,4,8]:
            for load, stress in [('cell','laplace'),('cell','symmetric'),('lift','laplace'),('lift','symmetric')]:
                cases.append((f'{mesh}-{n}-{load}-{stress}', mesh, str(n), problem, '0', load, stress, '8'))
        for load in ['cell','lift']:
            cases.append((f'{mesh}-8-gradient-{load}', mesh, '8', problem, '1', load, 'symmetric', '8'))
        cases.append((f'{mesh}-16-lift-symmetric', mesh, '16', problem, '0', 'lift', 'symmetric', '8'))
    for radial, angular in [(1,8),(2,16),(4,32)]:
        for load in ['cell','lift']:
            cases.append((f'ring-{radial}-{angular}-{load}', str(out/f'ring-{radial}-{angular}.cm2d'), '0', 'hydrostatic', '1', load, 'symmetric', '8'))
    for load in ['cell','lift']:
        cases.append((f'cut-4-hydro-{load}', 'cut', '4', 'hydrostatic', '1', load, 'symmetric', '8'))
    cases += [('square-4-strong-gradient','square','4','noslip','10000','lift','symmetric','8'),
              ('square-4-lift-symmetric-q12','square','4','noslip','0','lift','symmetric','12')]
    if args.case:
        cases = [case for case in cases if case[0] == args.case]
        if not cases:
            raise SystemExit('unknown case')
    for label, *_ in cases:
        if any(out.glob(label+'.*')):
            raise SystemExit('case output exists: '+label)
    rings = set()
    records = []
    for label, *control in cases:
        if label.startswith('ring-'):
            _, radial, angular, _ = label.split('-')
            mesh = Path(control[0])
            if mesh not in rings:
                subprocess.run([str(args.wall_binary), 'write-ring', radial, angular, str(mesh)], check=True)
                rings.add(mesh)
        prefix = out/label
        command = [*control, str(prefix)]
        started = time.monotonic()
        with Path(str(prefix)+'.stderr').open('w') as err:
            assembly = json.loads(subprocess.check_output([str(args.stress_binary),'assemble',*command], stderr=err, text=True))
            # Local reproduction is intentionally bounded. Do not silently switch
            # to a new sparse method or an approximate solve for larger inputs.
            if assembly['unknowns'] > 2300:
                raise RuntimeError('local dense linear-algebra budget exceeded')
            entries = np.fromfile(str(prefix)+'.entries', dtype=[('i','<i8'),('j','<i8'),('v','<f8')])
            matrix = np.zeros((assembly['unknowns'],)*2)
            np.add.at(matrix, (entries['i'],entries['j']), entries['v'])
            rhs = np.fromfile(str(prefix)+'.rhs', dtype='<f8')
            solution = np.linalg.solve(matrix,rhs)
            for _ in range(2):
                solution += np.linalg.solve(matrix,rhs-matrix@solution)
            if not np.all(np.isfinite(solution)):
                raise RuntimeError('nonfinite solution')
            solution.tofile(str(prefix)+'.solution')
            linear = {'maxResidual':float(np.max(abs(matrix@solution-rhs))),
                      'relativeResidual':float(np.linalg.norm(matrix@solution-rhs)/np.linalg.norm(rhs))}
            recovery = json.loads(subprocess.check_output([str(args.stress_binary),'recover',*command], stderr=err, text=True))
            wall = json.loads(subprocess.check_output([str(args.wall_binary),*command], stderr=err, text=True))
        record = {'label':label, 'args':command, 'assembly':assembly, 'linear':linear, 'recovery':recovery, 'wall':wall, 'seconds':time.monotonic()-started,
                  'binariesSha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.stress_binary,args.wall_binary]}}
        Path(str(prefix)+'.result.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
        records.append(record)
        print(label,wall['reactionPressureRms'],flush=True)
    manifest = out/(('case-'+args.case) if args.case else 'complete')
    manifest.with_suffix('.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')

if __name__ == '__main__':
    main()
