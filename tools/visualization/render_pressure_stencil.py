#!/usr/bin/env python3
"""Render an accepted CM2D flow field and its actual boundary-pressure stencil."""
import argparse
import csv
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'verification'))
import verify_native_flow as native

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--mesh', type=Path, required=True)
parser.add_argument('--prefix', type=Path, required=True)
parser.add_argument('--cell', type=int, required=True)
parser.add_argument('--time', type=float)
parser.add_argument('--view-radius', type=float)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
summary=json.loads(Path(str(args.prefix)+'.json').read_text())
if not summary.get('converged') or summary.get('status')!='converged':
    raise ValueError('Require an accepted flow result')
mode=summary.get('pressureBoundaryReconstruction')
if mode not in ('one-sided-linear','one-sided-linear-2ring','one-sided-linear-adaptive'):
    raise ValueError('This plot requires an extrapolated boundary-pressure operator')
actual_time=summary.get('acceptedTime',0)
if args.time is not None and not math.isclose(args.time,actual_time,rel_tol=1e-12,abs_tol=1e-12):
    raise ValueError('Requested time differs from the accepted field')
mesh = native.read_cm2d(args.mesh)
measured = native.measure(mesh, 1e-10, 1e-9)
if measured.issues:
    raise ValueError(measured.issues)
with Path(str(args.prefix) + '.cells.csv').open() as stream:
    rows = list(csv.DictReader(stream))
if len(rows) != len(mesh.cells):
    raise ValueError('Field and mesh dimensions differ')
polygons = [[mesh.vertices[v] for v in cell.vertices] for cell in mesh.cells]

def adjacent(i):
    return {e.neighbour if e.owner == i else e.owner
            for e in (mesh.edges[f] for f in mesh.cells[i].edges)
            if e.neighbour >= 0}

if not 0 <= args.cell < len(mesh.cells):
    raise ValueError('Unknown target cell')
# The displayed fixture must have unknown pressure on its boundary; a known
# outlet pressure would add a different row to the least-squares matrix.
if summary['case']!='external' or not any(mesh.edges[f].patch==1 for f in mesh.cells[args.cell].edges):
    raise ValueError('Select an external-flow embedded wall cell')
direct = adjacent(args.cell)
extended = set().union(*(adjacent(i) for i in direct)) - direct - {args.cell}
centre = measured.centroids[args.cell]

def condition(ids):
    xx = xy = yy = 0.
    for i in ids:
        dx, dy = (measured.centroids[i][k] - centre[k] for k in (0, 1))
        length = math.hypot(dx, dy)
        dx, dy = dx / length, dy / length
        xx += dx * dx; xy += dx * dy; yy += dy * dy
    root = math.hypot(xx - yy, 2 * xy)
    maximum=.5*(xx+yy+root)
    determinant=xx*yy-xy*xy
    return maximum*maximum/determinant if determinant>0 else math.inf

second = extended.copy()
extra=set()
if mode=='one-sided-linear-adaptive':
    visited={args.cell}|direct|second; frontier=second.copy()
    for ring in range(3,7):
        if condition(direct|second|extra)<=16:break
        fresh=set().union(*(adjacent(i) for i in frontier))-visited
        if not fresh:break
        extra.update(fresh);visited.update(fresh);frontier=fresh
elif mode=='one-sided-linear' and math.isfinite(condition(direct)) and 1/(condition(direct)+2+1/condition(direct))>64*sys.float_info.epsilon:
    second=set()
extended=second|extra

fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), layout='constrained')
ax = axes[0]
field = PolyCollection(polygons, array=[float(r['speed']) for r in rows],
                       cmap='viridis', edgecolors='#73818c', linewidths=.10)
ax.add_collection(field); ax.autoscale_view()
ax.set(aspect='equal', xlabel='x [m]', ylabel='y [m]',
       title=f'Actual accepted flow: {len(rows):,} cells | '+('steady' if not actual_time else f't = {actual_time:g} s'))
if args.view_radius is not None:
    if args.view_radius<=0:raise ValueError('view-radius must be positive')
    ax.set_xlim(centre[0]-args.view_radius,centre[0]+args.view_radius)
    ax.set_ylim(centre[1]-args.view_radius,centre[1]+args.view_radius)
ax.plot(*centre, 'o', color='#e65c28', markeredgecolor='white', markersize=6)
fig.colorbar(field, ax=ax, label='Speed [m/s]', shrink=.85)

ax = axes[1]
for ids, color, label in ((direct, '#377eb8', 'Direct neighbours'),
                         (second, '#e68632', 'Second ring'),
                         (extra, '#4b9b54', 'Adaptive extra samples'),
                         ({args.cell}, '#ce4c60', 'Boundary cell')):
    for i in sorted(ids):
        ax.add_collection(PolyCollection([polygons[i]], facecolors=color,
                          edgecolors='#3c4855', linewidths=1, alpha=.24))
        c = measured.centroids[i]
        ax.plot(*c, 'o', color=color, markersize=6)
        if i != args.cell:
            ax.plot([centre[0], c[0]], [centre[1], c[1]], color=color,
                    linewidth=1.3, linestyle='--' if i in extended else '-')
        if i == args.cell:
            ax.annotate(f'cell {i}', c, xytext=(6, 6), textcoords='offset points', fontsize=9)
    ax.plot([], [], 'o-', color=color, label=label)
ax.autoscale_view(); ax.margins(.15)
ax.set(aspect='equal', xlabel='x [m]', ylabel='y [m]',
       title=f'Pressure samples: {len(direct)} → {len(direct | extended)}\n'
             f'LS condition: direct {condition(direct):.1f}, 2-ring {condition(direct|second):.1f}, final {condition(direct | extended):.2f}')
ax.legend(loc='lower left', fontsize=9, framealpha=.95)
fig.suptitle('Boundary pressure reconstruction on a real cut-cell mesh', fontsize=15)
fig.supxlabel('Same geometry and solver tolerance; no cell removal or flux clipping.\n'
              'Local stability regression, not cylinder-flow accuracy certification.', fontsize=10)
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, dpi=160)
plt.close(fig)
