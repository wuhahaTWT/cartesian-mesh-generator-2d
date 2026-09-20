#!/usr/bin/env python3
"""Real curved Cut-cell duct: steady audit, accepted transient steps and restart."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools/verification'))
import verify_native_flow as steady
import verify_transient_flow as transient
import verify_scalar_transport as scalar


def run(command, success=True):
    result = subprocess.run(list(map(str, command)), capture_output=True, text=True, timeout=25)
    assert (result.returncode == 0) == success, result.stdout + result.stderr
    return result


def main(cli, mesh_cli, transport_cli):
    with tempfile.TemporaryDirectory(prefix='cartmesh-duct-') as name:
        root = Path(name)
        prefix = root / 'nozzle'
        mesh_path = root / 'nozzle.solver.cm2d'
        run([mesh_cli, ROOT / 'examples/complex/nozzle_profile.xy', prefix,
             4, 1/30, .1, 'interior', root / 'foam', 4, 0])
        mesh = steady.read_cm2d(mesh_path)
        measured = steady.measure(mesh, 1e-11, 1e-9)
        assert not measured.issues, measured.issues
        bc = steady.flow_boundaries(mesh, measured, 'duct', 1.)
        # A genuine curved domain, not the old filled rectangle shortcut.
        area = measured.total_area
        xmin, ymin, xmax, ymax = measured.bounds
        assert area < .9 * (xmax-xmin)*(ymax-ymin)
        assert all(role in bc['roles'] for role in ('inlet', 'outlet', 'wall'))
        assert any(bc['roles'][e.id] == 'wall' and
                   abs(mesh.vertices[e.v1][1]-mesh.vertices[e.v0][1]) > 1e-6
                   and abs(mesh.vertices[e.v1][0]-mesh.vertices[e.v0][0]) > 1e-6
                   for e in mesh.edges)

        def flow(output, extra=(), case='duct', success=True):
            return run([cli, '--mesh', mesh_path, '--output', output, '--case', case,
                        '--nu', .1, '--speed', 1, '--tolerance', 1e-8,
                        '--max-iterations', 400, '--pressure-preconditioner', 'aggregation',
                        *extra], success)

        output = root / 'steady'
        flow(output)
        report = steady.verify_case(mesh_path, output, 'duct', .1, 1., steady.argument_parser().parse_args([]))
        assert report['valid'], report['issues']
        assert report['benchmark']['status'] == 'not-qualified'  # no invented physical reference
        with Path(str(output)+'.faces.csv').open() as stream:
            faces = list(csv.DictReader(stream))
        assert all(float(faces[i]['flux']) == 0 for i, role in enumerate(bc['roles']) if role == 'wall')
        inlet = sum(float(faces[i]['flux']) for i, role in enumerate(bc['roles']) if role == 'inlet')
        assert abs(inlet + ymax-ymin) < 1e-12
        rejection = flow(root/'wrong', case='channel', success=False).stderr
        assert 'rectangular outer boundary' in rejection or 'not axis aligned' in rejection, rejection

        whole, first, resumed = (root/k for k in ('whole', 'first', 'resumed'))
        flow(whole, ['--time-step', .02, '--steps', 2])
        flow(first, ['--time-step', .02, '--steps', 1])
        flow(resumed, ['--time-step', .02, '--steps', 1,
                       '--restart', str(first)+'.checkpoint'])
        for output in (whole, resumed):
            report = transient.verify(mesh_path, output, root/(output.name+'-audit.json'))
            assert report['valid'], report
        assert Path(str(whole)+'.checkpoint').read_bytes() == Path(str(resumed)+'.checkpoint').read_bytes()

        # The App exposes synchronized passive heat transport for this same
        # flow. An initially uniform temperature with insulated walls stays
        # uniform while its nonuniform carrier accelerates through the throat.
        boundary = root / 'temperature.csv'
        with boundary.open('w') as stream:
            writer = csv.writer(stream)
            writer.writerow(['face', 'type', 'value', 'inflowValue'])
            for e in mesh.edges:
                if e.neighbour < 0:
                    writer.writerow([e.id, 'flux', 0, 300])
        heat = root / 'thermal'
        run([transport_cli, '--mesh', mesh_path, '--output', heat, '--boundary', boundary,
             '--evolve-flow', 'duct', '--flow-nu', .1, '--flow-speed', 1,
             '--pressure-preconditioner', 'aggregation', '--diffusivity', .1,
             '--initial', 300, '--dt', .02, '--steps', 2])
        thermal = json.loads(Path(str(heat)+'.json').read_text())
        assert thermal['flowCase'] == 'duct' and thermal['evolvingFlow']
        assert abs(thermal['minValue']-300) < 1e-9 and abs(thermal['maxValue']-300) < 1e-9
        assert scalar.verify(heat)['valid']

        # Rotate the actual mesh. Horizontal openings must fail explicitly;
        # selecting duct must not quietly turn a curved tip into an inlet.
        lines = mesh_path.read_text().splitlines()
        n = int(lines[1].split()[1])
        for i in range(2, 2+n):
            ident, x, y = lines[i].split()
            lines[i] = f'{ident} {-float(y):.17g} {float(x):.17g}'
        mesh_path.write_text('\n'.join(lines)+'\n')
        failure = flow(root/'rotated', success=False)
        assert 'left inlet and right pressure outlet' in failure.stderr, failure.stderr
        print(json.dumps({'valid': True, 'cells': len(mesh.cells), 'steadyAndTransient': True,
                          'restartIdentical': True, 'unsupportedRotatedPortsRejected': True}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cli', type=Path, required=True)
    parser.add_argument('--mesh-cli', type=Path, required=True)
    parser.add_argument('--transport-cli', type=Path, required=True)
    args = parser.parse_args()
    main(args.cli.resolve(), args.mesh_cli.resolve(), args.transport_cli.resolve())
