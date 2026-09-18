#!/usr/bin/env python3
"""Plot actual manufactured-flow polygons, fields and independently audited errors."""
import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'verification'))
from verify_native_flow import read_cm2d, manufactured_sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path, action='append', required=True,
                        help='Repeat for each mesh kind and convection scheme')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--field', choices=('speed', 'pressure'), default='speed')
    parser.add_argument('--baseline', type=Path,
                        help='Optional independently audited Cartesian series for comparison')
    args = parser.parse_args()
    groups = []
    for path in args.summary:
        report = json.loads(path.read_text())
        if not report['valid'] or not report['cases']:
            raise ValueError(f'Expected a passed manufactured-flow series: {path}')
        if report.get('format') == 'cartmesh2d-manufactured-flow-runner-v1':
            grouped = {}
            for item in report['cases']:
                case = dict(item['verification'], label=item['label'])
                grouped.setdefault((item['meshKind'], item['scheme']), []).append(case)
            series = grouped.values()
        else:
            series = [report['cases']]
        for cases in series:
            if len(cases) < 2 or any(not c['valid'] or c['case'] != 'manufactured' for c in cases):
                raise ValueError(f'At least two passed manufactured meshes are required per series: {path}')
            cases.sort(key=lambda c: c['counts']['cells'])
            scheme = cases[0]['momentumAudit']['convection']
            kind = 'Warped' if cases[0]['label'].startswith('warped-') else 'Cartesian'
            groups.append((kind, scheme, cases))
    candidates = [cases[-1] for kind, scheme, cases in groups
                  if kind == 'Warped' and scheme == 'limited-linear']
    if not candidates:
        raise ValueError('A warped limited-linear series is required for the field view')
    case = candidates[0]
    mesh = read_cm2d(Path(case['mesh']))
    with Path(case['prefix'] + '.cells.csv').open() as stream:
        rows = sorted(csv.DictReader(stream), key=lambda row: int(row['cell']))
    polygons = [[mesh.vertices[i] for i in c.vertices] for c in mesh.cells]
    speed = np.array([float(row['speed']) for row in rows])
    pressure = np.array([float(row['p']) for row in rows])
    errors = []
    pressure_errors = []
    slope = case.get('manufacturedPressureSlope', 0.0)
    gauge = manufactured_sample(float(rows[0]['x']), float(rows[0]['y']),
                                case['speed'], case['nu'], slope)['p']
    for row in rows:
        exact = manufactured_sample(float(row['x']), float(row['y']), case['speed'], case['nu'], slope)
        errors.append(np.hypot(float(row['u']) - exact['u'], float(row['v']) - exact['v']))
        pressure_errors.append(abs(float(row['p']) - (exact['p'] - gauge)))
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), layout='constrained')
    views = [
        (axes[0, 0], speed, 'viridis', 'Computed forced vortex', 'Speed (m/s)'),
        (axes[0, 1], np.array(errors), 'magma', 'Error against the analytic velocity', '|U - U exact| (m/s)'),
    ] if args.field == 'speed' else [
        (axes[0, 0], pressure, 'coolwarm', 'Computed kinematic pressure', 'p (m2/s2)'),
        (axes[0, 1], np.array(pressure_errors), 'magma', 'Error against the analytic pressure', '|p - p exact| (m2/s2)'),
    ]
    for ax, values, cmap, title, label in views:
        collection = PolyCollection(polygons, array=values, cmap=cmap,
                                    edgecolors=(.1, .15, .2, .22), linewidths=.14)
        ax.add_collection(collection)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect('equal')
        ax.set_title(f'{title}\n{len(rows):,} actual warped cells, limited reconstruction')
        ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
        fig.colorbar(collection, ax=ax, label=label, shrink=.8)
    sample = rows[::max(1, len(rows) // 160)]
    axes[0, 0].quiver([float(r['x']) for r in sample], [float(r['y']) for r in sample],
                      [float(r['u']) for r in sample], [float(r['v']) for r in sample],
                      color='white', scale=13, width=.002)
    for ax, metric, label in [(axes[1, 0], 'velocityL2Relative', 'Velocity L2 / Uref'),
                               (axes[1, 1], 'pressureL2Relative', 'Pressure L2 / Uref squared')]:
        for kind, scheme, cases in groups:
            h = [c['meshMeasurement']['characteristicH'] for c in cases]
            error = [c['benchmark'][metric] for c in cases]
            order = np.log(error[-2] / error[-1]) / np.log(h[-2] / h[-1])
            color = '#b46027' if scheme == 'limited-linear' else '#416986'
            ax.loglog(h, error, 'o-' if kind == 'Cartesian' else 's--', color=color,
                      label=f'{kind}, {scheme}; last order {order:.2f}')
        if args.baseline:
            baseline = json.loads(args.baseline.read_text())
            if not baseline['valid'] or not baseline['cases']:
                raise ValueError('Baseline must have passed independent verification')
            cases = sorted(baseline['cases'], key=lambda c: c['counts']['cells'])
            ax.loglog([c['meshMeasurement']['characteristicH'] for c in cases],
                      [c['benchmark'][metric] for c in cases], 'x:', color='#416986',
                      label='Previous pressure boundary, Cartesian')
        ax.set_xlabel('h = sqrt(total area / cells) (m)'); ax.set_ylabel(label)
        ax.set_title('Area-weighted error, three mesh sizes')
        ax.grid(alpha=.2, which='both'); ax.legend(fontsize=8)
    fig.suptitle('Native 2D Navier-Stokes: manufactured-solution verification', fontsize=15)
    fig.supxlabel(f'Stationary walls; analytic volume forcing; pressure slope = {slope:g}. This smooth case does not qualify cavity or turbulent-flow accuracy.', fontsize=9)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    main()
