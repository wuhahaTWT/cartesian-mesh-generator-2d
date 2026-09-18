#!/usr/bin/env python3
"""Small independent CLI regressions for flow physics and explicit failure paths."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile


def rectangle(path, nx, ny, width=4.0, separated=False):
    vertices = [(i * width / nx, j / ny)
                for j in range(ny + 1) for i in range(nx + 1)]
    edges, cells, by_pair = [], [], {}
    for j in range(ny):
        for i in range(nx):
            if separated and i in (2, 3):
                continue
            ids = [j * (nx + 1) + i, j * (nx + 1) + i + 1,
                   (j + 1) * (nx + 1) + i + 1, (j + 1) * (nx + 1) + i]
            incident = []
            for a, b in zip(ids, ids[1:] + ids[:1]):
                pair = tuple(sorted((a, b)))
                if pair not in by_pair:
                    by_pair[pair] = len(edges)
                    edges.append([len(edges), a, b, len(cells), -1, 2])
                else:
                    edges[by_pair[pair]][4:] = [len(cells), 0]
                incident.append(by_pair[pair])
            cells.append([len(cells), 0, 0, width / nx / ny, 4, *ids, 4, *incident])
    if separated:
        for edge in edges:
            a, b = vertices[edge[1]], vertices[edge[2]]
            if edge[4] == -1 and a[0] == b[0] and 0 < a[0] < width:
                edge[5] = 1  # Stationary embedded divider faces.
    text = ['CM2D 1', f'VERTICES {len(vertices)}']
    text += [f'{i} {x:.17g} {y:.17g}' for i, (x, y) in enumerate(vertices)]
    text += [f'EDGES {len(edges)}', *(' '.join(map(str, e)) for e in edges),
             f'CELLS {len(cells)}', *(' '.join(map(str, c)) for c in cells),
             'AUDIT 0 0 0 0 0 0 0', 'END']
    path.write_text('\n'.join(text) + '\n')


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


parser = argparse.ArgumentParser()
parser.add_argument('--cli', required=True)
parser.add_argument('--mesh-cli', required=True)
args = parser.parse_args()
cli = str(Path(args.cli).resolve())
with tempfile.TemporaryDirectory(prefix='cartmesh-flow-') as name:
    root = Path(name)

    def run(label, mesh, case='channel', extra=(), code=0, error_contains=None):
        prefix = root / label
        result = subprocess.run([cli, '--mesh', str(mesh), '--output', str(prefix),
                                 '--case', case, '--nu', '.01', '--speed', '1',
                                 '--max-iterations', '700', *extra],
                                text=True, capture_output=True, timeout=40)
        assert result.returncode == code, (label, result.stdout[-1000:], result.stderr)
        if code == 1:
            if error_contains:
                assert error_contains in result.stderr, result.stderr
            assert not prefix.with_suffix('.json').exists()
            return None
        data = json.loads(prefix.with_suffix('.json').read_text())
        field = rows(prefix.with_suffix('.cells.csv'))
        flux = rows(prefix.with_suffix('.faces.csv'))
        imbalance = [0.0] * len(field)
        for face in flux:
            q = float(face['flux'])
            assert math.isfinite(q)
            imbalance[int(face['owner'])] += q
            if int(face['neighbour']) >= 0:
                imbalance[int(face['neighbour'])] -= q
        assert max(map(abs, imbalance)) < 1e-8
        for cell in field:
            assert all(math.isfinite(float(cell[k])) for k in ('u', 'v', 'p', 'speed'))
            assert abs(float(cell['speed']) - math.hypot(float(cell['u']), float(cell['v']))) < 1e-12
        assert data['converged'] is (code == 0)
        assert len(json.loads(prefix.with_suffix('.fields.json').read_text())['cells']) == len(field)
        assert prefix.with_suffix('.vtk').read_text().count('CELL_DATA ') == 1
        return data, field

    errors, pressure_errors = [], []
    for ny in (8, 16):
        mesh = root / f'channel{ny}.solver.cm2d'
        rectangle(mesh, ny * 4, ny)
        data, field = run(f'channel{ny}', mesh)
        if ny == 8:
            # Profiling must observe, never alter the physical solve or exports.
            run('profiled', mesh, extra=('--profile',))
            for suffix in ('.json', '.fields.json', '.cells.csv', '.faces.csv',
                           '.residuals.csv', '.vtk'):
                assert (root / ('profiled' + suffix)).read_bytes() == (
                    root / (f'channel{ny}' + suffix)).read_bytes(), suffix
            profile = json.loads((root / 'profiled.performance.json').read_text())
            assert profile['simpleIterations'] == data['iterations']
            assert profile['converged'] is data['converged']
            assert profile['momentumSolves'] == 2 * data['iterations']
            assert profile['pressureSolves'] == 4 * data['iterations']
            assert profile['momentumIterations'] >= profile['maxMomentumIterations'] > 0
            assert profile['pressureIterations'] >= profile['maxPressureIterations'] > 0
            for key in ('readAndMeshSeconds', 'solveSeconds', 'momentumLinearSolveSeconds', 'pressureLinearSolveSeconds'):
                assert math.isfinite(profile[key]) and profile[key] >= 0, (key, profile[key])
            assert profile['solveSeconds'] >= profile['momentumLinearSolveSeconds'] + profile['pressureLinearSolveSeconds']
            assert profile['pressurePreconditioner'] == 'ic0'
            jacobi, jacobi_field = run('jacobi', mesh, extra=('--pressure-preconditioner', 'jacobi'))
            assert jacobi['pressurePreconditioner'] == 'jacobi'
            # Different Krylov paths must solve the same physical equations.
            for first, second in zip(field, jacobi_field):
                assert max(abs(float(first[k]) - float(second[k])) for k in ('u', 'v', 'p')) < 1e-8
        error = math.sqrt(sum((float(c['u']) - 4 * float(c['y']) * (1 - float(c['y']))) ** 2
                              for c in field) / len(field))
        errors.append(error)
        assert error < .025
        # Fit p(x) away from the pressure outlet; expected dp/dx=-8 nu Umax/H^2.
        interior = [c for c in field if .5 < float(c['x']) < 3.5]
        mean_x = sum(float(c['x']) for c in interior) / len(interior)
        mean_p = sum(float(c['p']) for c in interior) / len(interior)
        slope = sum((float(c['x']) - mean_x) * (float(c['p']) - mean_p) for c in interior)
        slope /= sum((float(c['x']) - mean_x) ** 2 for c in interior)
        pressure_errors.append(abs(slope / -.08 - 1))
    assert errors[1] < errors[0] / 2.5, errors
    assert pressure_errors[1] < .02 and pressure_errors[1] < pressure_errors[0] / 2, pressure_errors
    run('limited', mesh, extra=('--max-iterations', '1', '--profile'), code=2)
    limited = json.loads((root / 'limited.performance.json').read_text())
    assert limited['converged'] is False and limited['simpleIterations'] == 1
    for label, options in [('negative', ('--nu', '-1')), ('nan', ('--speed', 'nan')),
                           ('overflow', ('--speed', '1e200')), ('count', ('--max-iterations', '1.5')),
                           ('preconditioner', ('--pressure-preconditioner', 'unknown'))]:
        run(label, mesh, extra=options, code=1)
    run('no-solid', mesh, case='external', code=1)
    separated = root / 'separated.solver.cm2d'
    rectangle(separated, 6, 2, separated=True)
    run('disconnected', separated, case='external', code=1,
        error_contains='one connected fluid region')
    cavity = root / 'cavity.solver.cm2d'
    rectangle(cavity, 12, 12, 1.0)
    data, field = run('cavity', cavity, case='cavity')
    assert abs(float(field[0]['p'])) < 1e-12  # Closed-domain pressure gauge.
    centre = min(field, key=lambda c: (float(c['x'])-.5)**2 + (float(c['y'])-.5)**2)
    assert -.3 < float(centre['u']) < -.05  # Clockwise recirculation, not zero flow.
    # Real regression: a 3596-cell generated channel previously exhausted the
    # restarted pressure Krylov solver before its first SIMPLE step.
    outline = root / 'channel.xy'
    outline.write_text('0 0\n4 0\n4 1\n0 1\n')
    generated_prefix = root / 'generated-channel'
    generated = subprocess.run([str(Path(args.mesh_cli).resolve()), str(outline),
                                str(generated_prefix), '6', str(1 / 62), '.1',
                                'interior', str(root / 'openfoam'), '6', '0'],
                               text=True, capture_output=True, timeout=30)
    assert generated.returncode == 0, generated.stderr
    data, field = run('pressure-long-mode', root / 'generated-channel.solver.cm2d')
    assert len(field) > 3000 and data['globalRelativeImbalance'] < 1e-8
    print('Flow CLI: Poiseuille refinement, pressure gradient, closed-cavity gauge, '
          'independent flux balance and explicit failures verified.', errors)
