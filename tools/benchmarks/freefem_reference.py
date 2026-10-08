#!/usr/bin/env python3
"""Run the native FreeFEM 2D annulus reference on unchanged CM2D polygons.

This is an explicit external solver adapter, not a reconstruction of the native
FV equations. PDE assembly and sparse solves execute inside FreeFEM/UMFPACK.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = 'b1e524c8ca6a9b44cb9eff780459bacfc12a3ad7'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--boundary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='Output prefix')
    parser.add_argument('--freefem', type=Path, required=True)
    parser.add_argument('--converter', type=Path, default=ROOT / 'build/cartmesh2d_freefem_mesh')
    parser.add_argument('--space', choices=['cr', 'th'], default='cr')
    parser.add_argument('--trace', choices=['facets', 'exact'], default='facets')
    parser.add_argument('--stress', choices=['laplacian', 'symmetric'], default='laplacian')
    parser.add_argument('--inertia', type=int, choices=[0, 1], default=1)
    args = parser.parse_args()
    if args.space == 'th' and args.trace == 'facets':
        parser.error('P2 shared corner DOFs cannot satisfy both discontinuous facet velocities; use a declared continuous control.')
    if args.space == 'cr' and args.stress == 'symmetric':
        parser.error('CR symmetric-gradient coercivity is not established here; this reference uses the standard broken Laplacian.')
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = Path(str(output) + '.run.json')
    for suffix in ['.run.json', '.summary.json', '.cells.csv', '.iterations.csv', '.fe-solution.txt']:
        if Path(str(output) + suffix).exists():
            parser.error('Output prefix already contains a run; choose a new prefix to preserve prior fields and failures.')
    converted = Path(str(output) + '-input')
    script = ROOT / 'tools/benchmarks' / ('freefem_' + args.space + '.edp')
    files = [args.mesh.resolve(), args.boundary.resolve(), args.converter.resolve(), args.freefem.resolve(),
             Path(__file__), ROOT / 'tools/benchmarks/freefem_mesh.cpp',
             ROOT / 'tools/benchmarks/freefem_core.idp', script]
    for path in files:
        if not path.is_file():
            parser.error('Required input or executable is missing: ' + str(path))
    record = dict(status='preparing', sourceHead=subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        sourceAndInputSha256={str(p): sha(p) for p in files},
        sourceHasUncommittedChanges=bool(subprocess.check_output(
            ['git', 'status', '--porcelain', '--untracked-files=normal'], cwd=ROOT, text=True).strip()),
        reference=dict(package='FreeFEM 4.15', sourceCommit=REFERENCE,
                       sourceUrl='https://github.com/FreeFem/FreeFem-sources/tree/' + REFERENCE,
                       license='FreeFEM LGPL-2.1; packaged components have their own notices'),
        controls=dict(space=args.space, trace=args.trace, stress=args.stress,
                      inertia=args.inertia, nu=.1, speed=.5, tolerance=1e-8), runs=[])
    def save():
        metadata.write_text(json.dumps(record, indent=2) + '\n')
    def run(name, command):
        started = time.perf_counter()
        log = Path(str(output) + '.' + name + '.log')
        with log.open('w') as stream:
            result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        record['runs'].append(dict(name=name, command=command, returnCode=result.returncode,
                                   wallSeconds=time.perf_counter() - started, log=str(log), sha256=sha(log)))
        if result.returncode:
            record['status'] = 'failed-' + name
        save()
        return result.returncode
    save()
    code = run('export', [str(args.converter.resolve()), str(args.mesh.resolve()),
                          str(args.boundary.resolve()), str(converted)])
    if code:
        return code
    record['meshMapping'] = json.loads(Path(str(converted) + '.mapping.json').read_text())
    record['status'] = 'solving'
    save()
    code = run('solve', [str(args.freefem.resolve()), '-nw', '-v', '1', str(script),
                        '--input', str(converted), '--output', str(output), '--inertia', str(args.inertia),
                        '--boundary', args.trace, '--symmetric', str(int(args.stress == 'symmetric'))])
    summary = Path(str(output) + '.summary.json')
    if summary.exists():
        record['solution'] = json.loads(summary.read_text())
        record['status'] = 'converged' if code == 0 and record['solution']['converged'] else 'failed-solve'
        if not record['solution']['converged'] and code == 0:
            code = 2
        record['fields'] = {suffix: sha(str(output) + suffix)
                            for suffix in ['.cells.csv', '.fe-solution.txt', '.iterations.csv']
                            if Path(str(output) + suffix).exists()}
    elif code == 0:
        record['status'] = 'missing-summary'
        code = 3
    record['scope'] = dict(
        topology='Every original polygon/edge/boundary retained; positive fan triangulation adds original FV centres.',
        boundary='CR/facets preserves original face velocity means; exact is a different boundary control, imposed by FE interpolation at Dirichlet degrees of freedom.',
        fieldMapping='FreeFEM quadrature computes true triangle means, then aggregates them to unchanged native polygons; native FV values are centroid estimates.',
        pressure='Exact pressure point constraint and mean subtraction, without a pressure mass penalty; all original continuity residuals, including the gauge row, are checked.',
        accuracy='Circle error is a physical approximation comparison for facet data; exact trace restricts the analytic solution to the polygon domain, with FE boundary interpolation error.',
        convergence='Newton increment and native free-DOF weak residual normalized by declared velocity/pressure/equation scales; not native FV acceptance or a spatial-accuracy certificate.')
    save()
    print(json.dumps({'result': str(metadata), 'status': record['status'], 'solution': record.get('solution')}, indent=2))
    return code


if __name__ == '__main__':
    raise SystemExit(main())
