#!/usr/bin/env python3
"""Bounded serial curved-duct grid study, with independent conservation audits.

The selected outline must have vertical openings at its x-extrema. Grid changes
are reported separately from convergence/conservation; no reference accuracy is
claimed. Final meshes, raw fields and failed runs remain under the output root.
"""
import argparse
import json
import math
from pathlib import Path
import shutil
import time

import verify_native_flow as native


def generate(args):
    root = args.output.resolve()
    if root.exists() and not args.verify_only:
        raise ValueError('output must be a fresh directory')
    root.mkdir(parents=True, exist_ok=args.verify_only)
    previous = json.loads((root/'summary.json').read_text()) if args.verify_only else None
    geometry = args.geometry.resolve(strict=True)
    mesh_cli, flow_cli = args.mesh_cli.resolve(strict=True), args.flow_cli.resolve(strict=True)
    report = dict(valid=False, cases=[], runs=[], issues=[], meshSensitivity=[],
                  scope='Steady curved duct geometry/conservation audit and grid sensitivity; no reference accuracy qualification.',
                  geometry=dict(path=str(geometry), sha256=native.sha256_file(geometry)),
                  executables={str(p):native.sha256_file(p) for p in (mesh_cli, flow_cli)},
                  controls=dict(nu=args.nu, speed=args.speed, tolerance=1e-9, maxIterations=4000,
                                convection=args.convection, timeoutSeconds=args.timeout))

    def save():
        native.write_json(root/'summary.json', report)

    def execute(command, label):
        command = list(map(str, command))
        if previous is not None:
            matches = [r for r in previous['runs'] if r['command'] == command]
            if len(matches) != 1 or matches[0]['status'] != 'passed':
                raise ValueError(f'{label}: no unique successful recorded invocation matching these controls')
            report['runs'].append(matches[0])
            return matches[0]
        if shutil.disk_usage(root).free < 5*1024**3:
            raise ValueError('resource stop: less than 5 GiB free')
        start = time.monotonic()
        record = native.run(command, root/'logs'/label, args.timeout)
        record['elapsedSeconds'] = time.monotonic()-start
        report['runs'].append(record)
        save()
        if record['status'] != 'passed':
            raise ValueError(f'{label} failed; recorded artifacts retained')
        return record

    try:
        for level in args.levels:
            label = f'nozzle-l{level}'
            directory = root/label
            directory.mkdir(exist_ok=args.verify_only)
            prefix = directory/'mesh'
            execute([mesh_cli, geometry, prefix, level, 1/30, .1, 'interior',
                     directory/'openfoam', level, 0], label+'-mesh')
            mesh_path = Path(str(prefix)+'.solver.cm2d')
            output = directory/'flow'
            stage = execute([flow_cli, '--mesh', mesh_path, '--output', output,
                             '--case', 'duct', '--nu', args.nu, '--speed', args.speed,
                             '--tolerance', '1e-9', '--max-iterations', '4000',
                             '--pressure-preconditioner', 'aggregation', '--convection', args.convection], label+'-flow')
            item = native.verify_case(mesh_path, output, 'duct', args.nu, args.speed,
                                      native.argument_parser().parse_args(['--max-iterations', '4000']))
            item.update(label=label, mesh=str(mesh_path), prefix=str(output), stage=stage)
            report['cases'].append(item)
            save()
            if not item['valid']:
                raise ValueError(f'{label} independent audit failed')
            mesh = native.read_cm2d(mesh_path)
            measured = native.measure(mesh, 1e-11, 1e-9)
            bc = native.flow_boundaries(mesh, measured, 'duct', args.speed)
            faces, _ = native.read_faces(Path(str(output)+'.faces.csv'), mesh)
            cells = native.read_cells(Path(str(output)+'.cells.csv'), mesh, measured, 'duct')
            def pressure(role):
                selected = [e for e in mesh.edges if bc['roles'][e.id] == role]
                def length(e):
                    a, z = mesh.vertices[e.v0], mesh.vertices[e.v1]
                    return math.hypot(a[0]-z[0], a[1]-z[1])
                return sum(faces[e.id]['pressure']*length(e) for e in selected)/sum(map(length, selected))
            metrics = dict(pressureDrop=pressure('inlet')-pressure('outlet'),
                           maxSpeed=max(c['speed'] for c in cells),
                           wallForceX=item['native']['wallForceX'],
                           wallForceY=item['native']['wallForceY'])
            sensitivity = dict(label=label, cells=len(mesh.cells), h=measured.characteristic_h,
                               metrics=metrics, relativeChange=None)
            if report['meshSensitivity']:
                old = report['meshSensitivity'][-1]['metrics']
                sensitivity['relativeChange'] = {
                    key:abs(metrics[key]-old[key])/max(abs(metrics[key]), 1e-30)
                    for key in ('pressureDrop', 'maxSpeed', 'wallForceX')}
            report['meshSensitivity'].append(sensitivity)
            save()
        report['valid'] = True
    except (ValueError, OSError, KeyError) as exc:
        report['issues'].append(str(exc))
    save()
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verify-only', action='store_true', help='reaudit retained fields with identical recorded run controls')
    parser.add_argument('--geometry', type=Path, default=native.REPO/'examples/complex/nozzle_profile.xy')
    parser.add_argument('--mesh-cli', type=Path, default=native.REPO/'build/cartmesh2d_cli')
    parser.add_argument('--flow-cli', type=Path, default=native.REPO/'build/cartmesh2d_flow_cli')
    parser.add_argument('--levels', type=int, nargs='+', default=[5, 6, 7])
    parser.add_argument('--nu', type=float, default=.1)
    parser.add_argument('--speed', type=float, default=1.)
    parser.add_argument('--convection', choices=['upwind', 'limited-linear'], default='limited-linear')
    parser.add_argument('--timeout', type=int, default=180)
    args = parser.parse_args()
    if (len(args.levels) < 3 or args.levels != sorted(set(args.levels)) or
            not all(3 <= level <= 9 for level in args.levels) or
            not all(math.isfinite(x) and x > 0 for x in (args.nu, args.speed)) or
            not 1 <= args.timeout <= 600):
        parser.error('use at least three increasing levels in [3,9], positive finite nu/speed, timeout 1..600 s')
    try:
        result = generate(args)
    except (ValueError, OSError) as exc:
        parser.exit(1, str(exc)+'\n')
    print(json.dumps({key:result[key] for key in ('valid', 'issues', 'meshSensitivity')}, indent=2))
    raise SystemExit(0 if result['valid'] else 1)
