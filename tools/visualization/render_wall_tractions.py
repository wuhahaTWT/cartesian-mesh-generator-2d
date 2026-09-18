#!/usr/bin/env python3
"""Plot real native polygons and independently integrated wall-traction errors."""
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
from verify_native_flow import read_cm2d


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manufactured', type=Path, required=True)
    parser.add_argument('--tractions', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    mms = json.loads(args.manufactured.read_text())
    audit = json.loads(args.tractions.read_text())
    if not mms['valid'] or not audit['valid']:
        raise ValueError('Expected verified MMS and valid traction readback')
    case = max((c for c in mms['cases'] if c['meshKind'] == 'warped'),
               key=lambda c: c['verification']['counts']['cells'])['verification']
    mesh = read_cm2d(Path(case['mesh']))
    with Path(case['prefix'] + '.cells.csv').open() as stream:
        cells = sorted(csv.DictReader(stream), key=lambda x: int(x['cell']))
    polygons = [[mesh.vertices[i] for i in c.vertices] for c in mesh.cells]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
    collection = PolyCollection(polygons, array=np.array([float(c['speed']) for c in cells]),
                                cmap='viridis', edgecolor='#e4e8e980', linewidth=.12)
    axes[0, 0].add_collection(collection)
    axes[0, 0].set(xlim=(0, 1), ylim=(0, 1), xlabel='x', ylabel='y',
                   title=f"Actual warped mesh / {len(cells):,} cells")
    axes[0, 0].set_aspect('equal')
    fig.colorbar(collection, ax=axes[0, 0], label='Computed speed (m/s)')
    for kind, color in [('cartesian', '#2377a4'), ('warped', '#b96535')]:
        series = sorted((c for c in audit['cases'] if c.get('meshKind') == kind),
                        key=lambda c: c['counts']['cells'])
        axes[0, 1].loglog([c['counts']['cells'] for c in series],
                         [c['tractionRmsPerFace']['total'] for c in series],
                         'o-', color=color, label=kind)
        fine = series[-1]
        fm = read_cm2d(Path(fine['mesh']))
        rows = []
        for row in fine['wallFaces']:
            e = fm.edges[row['face']]
            a, b = fm.vertices[e.v0], fm.vertices[e.v1]
            if abs(a[1]) < 1e-10 and abs(b[1]) < 1e-10:
                rows.append(((a[0]+b[0])/2, row['viscousActual'][0]))
        rows.sort()
        axes[1, 1].plot(*zip(*rows), '.', color=color, label=kind)
    axes[0, 1].set(xlabel='Final cells', ylabel='Wall traction RMS (m²/s²)',
                   title='Analytic forced flow: total wall traction error')
    for scheme, marker in [('upwind', 's'), ('limited-linear', 'o')]:
        series = sorted((c for c in audit['cases'] if c['case'] == 'channel' and c['scheme'] == scheme),
                        key=lambda c: c['counts']['cells'])
        axes[1, 0].loglog([c['counts']['cells'] for c in series],
                         [abs(c['forceTotals']['viscous']['error'][0] / c['channelViscousExpectedForceX'])
                          for c in series], marker+'-', label=scheme)
    axes[1, 0].set(xlabel='Final cells', ylabel='Relative viscous wall force error',
                   title='Poiseuille channel: analytic force = 0.32')
    x = np.linspace(0, 1, 401)
    axes[1, 1].plot(x, 2*np.pi*case['nu']*case['speed']*np.sin(np.pi*x)**2,
                    color='#333333', linewidth=1, label='analytic')
    axes[1, 1].set(xlabel='x along bottom wall', ylabel='Viscous traction x (m²/s²)',
                   title='Finest forced flow: bottom-wall shear')
    for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.grid(alpha=.18, which='both'); ax.legend(frameon=False)
    fig.suptitle('Conservative symmetric stress | actual native results', fontsize=16)
    fig.supxlabel('Independent face audit passed. Original cavity Ghia refinement gate remains failed.', fontsize=10)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=165)


if __name__ == '__main__':
    main()
