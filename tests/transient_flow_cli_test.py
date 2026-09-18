#!/usr/bin/env python3
"""Small transient Taylor-Green, restart, and checkpoint CLI regressions."""
import argparse
import csv
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile


def rectangle(path, n=14):
    vertices = [(i / n, j / n) for j in range(n + 1) for i in range(n + 1)]
    edges, cells, by_pair = [], [], {}
    for j in range(n):
        for i in range(n):
            ids = [j * (n + 1) + i, j * (n + 1) + i + 1,
                   (j + 1) * (n + 1) + i + 1, (j + 1) * (n + 1) + i]
            incident = []
            for a, b in zip(ids, ids[1:] + ids[:1]):
                key = tuple(sorted((a, b)))
                if key not in by_pair:
                    by_pair[key] = len(edges)
                    edges.append([len(edges), a, b, len(cells), -1, 2])
                else:
                    edges[by_pair[key]][4:] = [len(cells), 0]
                incident.append(by_pair[key])
            cells.append([len(cells), 0, 0, 1 / n / n, 4, *ids, 4, *incident])
    text = ['CM2D 1', f'VERTICES {len(vertices)}']
    text += [f'{i} {x:.17g} {y:.17g}' for i, (x, y) in enumerate(vertices)]
    text += [f'EDGES {len(edges)}'] + [' '.join(map(str, edge)) for edge in edges]
    text += [f'CELLS {len(cells)}'] + [' '.join(map(str, cell)) for cell in cells]
    text += ['AUDIT 0 0 0 0 0 0 0', 'END']
    path.write_text('\n'.join(text) + '\n')


def invoke(cli, mesh, output, dt=None, steps=None, extra=(), restart=None):
    command = [str(cli), '--mesh', str(mesh), '--output', str(output),
               '--case', 'taylor-green', '--nu', '.01', '--speed', '1',
               '--tolerance', '1e-9', '--max-iterations', '150']
    if dt is not None:
        command += ['--time-step', str(dt)]
    if steps is not None:
        command += ['--steps', str(steps)]
    if restart is not None:
        command += ['--restart', str(restart)]
    return subprocess.run(command + list(extra), text=True, capture_output=True, timeout=20)


def read_rows(path):
    return list(csv.DictReader(path.open()))


def state_error(path, time=.2):
    rows = read_rows(path.with_suffix('.cells.csv'))
    amplitude = math.exp(-2 * .01 * math.pi * math.pi * time)
    x0, y0 = float(rows[0]['x']), float(rows[0]['y'])
    p0 = .25 * amplitude * amplitude * (math.cos(2 * math.pi * x0) + math.cos(2 * math.pi * y0))
    error = 0.0
    dot = norm = exact_norm = 0.0
    for row in rows:
        x, y = float(row['x']), float(row['y'])
        u = amplitude * math.sin(math.pi * x) * math.cos(math.pi * y)
        v = -amplitude * math.cos(math.pi * x) * math.sin(math.pi * y)
        p = .25 * amplitude * amplitude * (math.cos(2 * math.pi * x) + math.cos(2 * math.pi * y)) - p0
        error += (float(row['u']) - u) ** 2 + (float(row['v']) - v) ** 2
        observed = float(row['p'])
        dot += observed * p
        norm += observed * observed
        exact_norm += p * p
    return math.sqrt(error / len(rows)), dot / math.sqrt(norm * exact_norm)


def solution_distance(first, second):
    left, right = read_rows(first.with_suffix('.cells.csv')), read_rows(second.with_suffix('.cells.csv'))
    values = []
    for a, b in zip(left, right):
        values.extend(float(a[key]) - float(b[key]) for key in ('u', 'v', 'p'))
    return math.sqrt(sum(value * value for value in values) / len(values))


def main(cli):
    with tempfile.TemporaryDirectory(prefix='cartmesh-transient-') as name:
        root = Path(name)
        mesh = root / 'taylor.solver.cm2d'
        rectangle(mesh)

        # Fixed-grid runs check physical decay and self-convergence of the
        # accepted solution; they do not claim spatial convergence.
        runs = {}
        for dt in (.04, .02, .01):
            prefix = root / ('dt' + str(dt).replace('.', '_'))
            result = invoke(cli, mesh, prefix, dt, int(round(.2 / dt)))
            assert result.returncode == 0, (dt, result.stderr)
            summary = json.loads(prefix.with_suffix('.json').read_text())
            assert summary['completedSteps'] == int(round(.2 / dt))
            assert summary['acceptedTime'] > .19 and summary['time'] > .19
            history = list(csv.DictReader(prefix.with_suffix('.time-history.csv').open()))
            energies = [float(row['kineticEnergy']) for row in history]
            assert all(math.isfinite(value) and value > 0 for value in energies)
            assert all(a > b for a, b in zip(energies, energies[1:]))
            runs[dt] = (prefix, state_error(prefix))
        assert runs[.01][1][1] > .98  # analytic pressure field has positive correlation
        coarse_change = solution_distance(runs[.04][0], runs[.02][0])
        fine_change = solution_distance(runs[.02][0], runs[.01][0])
        assert fine_change < coarse_change

        first = root / 'split_first'
        assert invoke(cli, mesh, first, .04, 2).returncode == 0
        split = root / 'split_final'
        assert invoke(cli, mesh, split, .04, 3, restart=first.with_suffix('.checkpoint')).returncode == 0
        continuous = root / 'continuous'
        assert invoke(cli, mesh, continuous, .04, 5).returncode == 0
        for suffix in ('.cells.csv', '.faces.csv', '.fields.json', '.checkpoint'):
            assert split.with_suffix(suffix).read_bytes() == continuous.with_suffix(suffix).read_bytes(), suffix

        failed = root / 'failed'
        result = invoke(cli, mesh, failed, .04, 1, ('--max-iterations', '1'))
        assert result.returncode == 2, result.stderr
        failed_summary = json.loads(failed.with_suffix('.json').read_text())
        assert failed_summary['completedSteps'] == 0 and failed_summary['acceptedTime'] == 0
        assert re.search(r'^TIME 0(?:\.0+)?$', failed.with_suffix('.checkpoint').read_text(), re.MULTILINE)

        tampered_config = root / 'tampered_config.checkpoint'
        config = first.with_suffix('.checkpoint').read_text().replace(
            '"taylor-green" 0.01 ', '"taylor-green" 0.011 ', 1)
        tampered_config.write_text(config)
        assert invoke(cli, mesh, root / 'bad_config', .04, 1, restart=tampered_config).returncode == 1
        tampered_mesh = root / 'tampered_mesh.checkpoint'
        mesh_text = re.sub(r'^(CELL 0 \S+ \S+ )\S+', r'\g<1>0.006',
                           first.with_suffix('.checkpoint').read_text(), count=1, flags=re.MULTILINE)
        tampered_mesh.write_text(mesh_text)
        assert invoke(cli, mesh, root / 'bad_mesh', .04, 1, restart=tampered_mesh).returncode == 1

        # An exception after starting a new run must invalidate an
        # earlier successful summary under that same prefix.
        assert invoke(cli, mesh, continuous, .04, 1, ('--manufactured-pressure-slope', '1')).returncode == 1
        failed_again = json.loads(continuous.with_suffix('.json').read_text())
        assert failed_again['status'] == 'failed' and failed_again['converged'] is False
        assert failed_again['acceptedTime'] == 0

        assert invoke(cli, mesh, root / 'missing_steps', .04).returncode == 1
        assert invoke(cli, mesh, root / 'missing_dt', None, 1).returncode == 1
        assert invoke(cli, mesh, root / 'zero_dt', 0, 1).returncode == 1
        assert invoke(cli, mesh, root / 'steady_taylor').returncode == 1

    print('Transient flow CLI: decay, fixed-grid self-convergence, restart identity, and checkpoint rejection paths verified.')
    print('raw dt=.04/.02/.01 velocity errors:', [runs[dt][1][0] for dt in (.04, .02, .01)])
    print('raw self-distance coarse/fine:', coarse_change, fine_change,
          'pressure correlation dt=.01:', runs[.01][1][1])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cli', required=True)
    parser.add_argument('--mesh-cli')  # accepted for parity with other CLI regressions
    main(parser.parse_args().cli)
