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
        assert data['pressureDiscretization'] == 'shared-face-gauss'
        assert data['convection'] in ('upwind', 'limited-linear')
        for face in flux:
            assert all(math.isfinite(float(face[k])) for k in
                       ('pressure','advectionX','advectionY','diffusionX','diffusionY'))
        assert len(json.loads(prefix.with_suffix('.fields.json').read_text())['cells']) == len(field)
        assert prefix.with_suffix('.vtk').read_text().count('CELL_DATA ') == 1
        return data, field

    errors, pressure_errors, wall_force_errors = [], [], []
    for ny in (8, 16):
        mesh = root / f'channel{ny}.solver.cm2d'
        rectangle(mesh, ny * 4, ny)
        data, field = run(f'channel{ny}', mesh)
        assert data['viscousStress'] == 'symmetric'
        assert data['forceDefinition'] == 'shared-face-newtonian-traction'
        # Poiseuille: both walls receive +x shear, total 8*nu*Umax*L/H=.32.
        wall_force_errors.append(abs(data['wallViscousForceX'] - .32))
        assert abs(data['wallViscousForceY']) < 1e-8
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
            # Zero-residual pressure passes may skip preconditioning. Require
            # bounded accounting without assuming every pressure solve builds
            # or consumes a factor; direct cache tests cover actual reuse.
            assert 0 <= profile['pressureFactorizations'] <= profile['pressureSolves']
            assert 0 <= profile['pressureFactorReuses'] <= profile['pressureSolves']
            assert (profile['pressureFactorizations'] +
                    profile['pressureFactorReuses'] <= profile['pressureSolves'])
            assert profile['momentumIterations'] >= profile['maxMomentumIterations'] > 0
            assert profile['pressureIterations'] >= profile['maxPressureIterations'] > 0
            for key in ('readAndMeshSeconds', 'solveSeconds', 'momentumLinearSolveSeconds', 'pressureLinearSolveSeconds'):
                assert math.isfinite(profile[key]) and profile[key] >= 0, (key, profile[key])
            assert profile['solveSeconds'] >= profile['momentumLinearSolveSeconds'] + profile['pressureLinearSolveSeconds']
            assert profile['pressurePreconditioner'] == 'ic0'
            aggregation, aggregation_field = run(
                'aggregation', mesh,
                extra=('--pressure-preconditioner', 'aggregation', '--profile'))
            aggregation_profile = json.loads(
                (root / 'aggregation.performance.json').read_text())
            assert aggregation_profile['pressurePreconditioner'] == 'aggregation'
            assert 0 < aggregation_profile['pressureHierarchyBuilds'] <= aggregation_profile['simpleIterations']
            assert 0 <= aggregation_profile['pressureHierarchyReuses'] <= (
                aggregation_profile['pressureSolves'] - aggregation_profile['pressureHierarchyBuilds'])
            assert aggregation_profile['maxPressureHierarchyLevels'] > 1
            assert aggregation_profile['maxPressureCoarseCells'] <= 32
            for first, second in zip(field, aggregation_field):
                assert max(abs(float(first[k]) - float(second[k]))
                           for k in ('u', 'v', 'p')) < 1e-8
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
    assert wall_force_errors[1] < wall_force_errors[0], wall_force_errors
    assert pressure_errors[1] < .02 and pressure_errors[1] < pressure_errors[0] / 2, pressure_errors
    run('limited', mesh, extra=('--max-iterations', '1', '--profile'), code=2)
    limited = json.loads((root / 'limited.performance.json').read_text())
    assert limited['converged'] is False and limited['simpleIterations'] == 1
    for label, options in [('negative', ('--nu', '-1')), ('nan', ('--speed', 'nan')),
                           ('overflow', ('--speed', '1e200')), ('count', ('--max-iterations', '1.5')),
                           ('preconditioner', ('--pressure-preconditioner', 'unknown')),
                           ('convection', ('--convection', 'unknown')),
                           ('stress', ('--viscous-stress', 'unknown'))]:
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
    high, high_field = run('cavity-high', cavity, case='cavity',
                           extra=('--convection','limited-linear'))
    assert high['convection'] == 'limited-linear'
    laplacian, laplacian_field = run('cavity-laplacian', cavity, case='cavity',
        extra=('--convection','limited-linear','--viscous-stress','laplacian'))
    assert laplacian['viscousStress']=='laplacian'
    assert laplacian['forceDefinition']=='reconstructed-newtonian-traction'
    assert max(abs(float(a['u'])-float(b['u'])) for a,b in zip(high_field,laplacian_field)) > 1e-5
    assert abs(float(high_field[0]['p'])) < 1e-12
    # The selection must actually change nonlinear transport, not only metadata.
    assert max(abs(float(a['u'])-float(b['u'])) for a,b in zip(field,high_field)) > 1e-3
    run('high-limit', cavity, case='cavity',
        extra=('--convection','limited-linear','--max-iterations','1'), code=2)
    # A single pressure boundary must retain both incoming and outgoing flux.
    # Driven parallel counterflow has an analytic solution; do not clip its
    # negative normal velocity or confuse this with a usual positive inlet.
    reverse_errors=[]
    for n in (8,16):
        reverse_mesh=root/f'counterflow-{n}.solver.cm2d'
        rectangle(reverse_mesh,n,n,1.)
        if n==8:
            run('reverse-rejected',reverse_mesh,case='counterflow',code=1,
                error_contains='backflow unsupported')
            run('bad-reverse-mode',reverse_mesh,case='counterflow',
                extra=('--outlet-backflow','unknown'),code=1)
        summary,field=run(f'reverse-{n}',reverse_mesh,case='counterflow',
            extra=('--nu','.1','--outlet-backflow','normal-inlet',
                   '--convection','limited-linear','--tolerance','1e-8'))
        assert summary['outletBackflow']=='normal-inlet'
        face_rows=rows(root/f'reverse-{n}.faces.csv')
        # On this Cartesian fixture right-boundary owners have x=1-h/2.
        outlet=[f for f in face_rows if int(f['neighbour'])<0 and
                float(field[int(f['owner'])]['x'])>1-.51/n and abs(float(f['flux']))>1e-10]
        incoming=[f for f in outlet if float(f['flux'])<0]
        assert incoming and any(float(f['flux'])>0 for f in outlet)
        assert summary['outletBackflowFaces']==len(incoming)
        assert abs(summary['outletInflow']+sum(float(f['flux']) for f in incoming))<1e-12
        for face in incoming:
            cell=field[int(face['owner'])]
            assert abs(float(face['advectionX'])-float(face['flux'])*float(cell['u']))<1e-12
            assert float(face['advectionY'])==0 and float(face['pressure'])==0
        reverse_errors.append(math.sqrt(sum(float(c['area'])*
            ((float(c['u'])-1-2*math.cos(2*math.pi*float(c['y'])))**2+float(c['v'])**2)
            for c in field)))
    assert reverse_errors[1]<reverse_errors[0]/3,reverse_errors
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
