#!/usr/bin/env python3
"""Render actual rectilinear meshes and audited flow, including failed refinements."""
import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'verification'))
import verify_native_flow as native

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--study', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
study = json.loads(args.study.read_text())
cases = study['cases']
if not cases or not all(c['valid'] for c in cases):
    raise ValueError('Each displayed flow needs its independent audit')
selected = next(c for c in cases if c['stretch'] > 0)
mesh = native.read_cm2d(Path(selected['mesh']))
polygons = [[mesh.vertices[v] for v in c.vertices] for c in mesh.cells]

def rows(case, suffix):
    with Path(case['prefix'] + suffix).open() as stream:
        return list(csv.DictReader(stream))

def label(case):
    return f"{case['nx']} x {case['ny']}, " + ('graded' if case['stretch'] else 'uniform')

fig, ax = plt.subplots(2, 3, figsize=(16, 9), layout='constrained')
ax[0, 0].add_collection(PolyCollection(polygons, facecolors='#f4f6f8',
    edgecolors='#5e7387', linewidths=.35))
ax[0, 0].set(xlim=(0, 1), ylim=(0, .3), xlabel='x [m]', ylabel='y [m]',
    title=f"Actual near-wall mesh: {selected['cells']} cells")
ax[0, 0].set_aspect('equal')
ax[0, 0].plot([0, .5], [0, 0], '--', color='#cd9c45', lw=3)
ax[0, 0].plot([.5, 1], [0, 0], color='#20252a', lw=4)
field = rows(selected, '.cells.csv')
collection = PolyCollection(polygons, array=[float(r['speed']) for r in field],
    cmap='viridis', edgecolors='none', rasterized=True)
ax[0, 1].add_collection(collection)
ax[0, 1].set(xlim=(0, 1), ylim=(0, 1), xlabel='x [m]', ylabel='y [m]', title='Returned speed [m/s]')
ax[0, 1].set_aspect('equal')
fig.colorbar(collection, ax=ax[0, 1], shrink=.8)
for case in cases:
    field = rows(case, '.cells.csv')
    # An actual cell-centre column, without invented near-wall interpolation.
    column = min({float(r['x']) for r in field}, key=lambda x: (abs(x-.75), x))
    profile = sorted((r for r in field if float(r['x'])==column), key=lambda r: float(r['y']))
    ax[0, 2].plot([float(r['u']) for r in profile], [float(r['y']) for r in profile], '.-', label=label(case))
    samples = case['plateWallSamples']
    for axis, key in [(ax[1, 0], 'Cf'), (ax[1, 1], 'yPlus')]:
        axis.plot([r['xFromLeadingEdge'] for r in samples], [r[key] for r in samples], '.-', label=label(case))
    history = rows(case, '.history.csv')
    ax[1, 2].semilogy([int(r['iteration']) for r in history],
        [float(r['momentumResidual']) for r in history], label=label(case))
ax[0, 2].set(xlabel='u / U_infinity', ylabel='y [m]', ylim=(0, .3), title=f'Cell-centre profile, x = {column:.6f} m')
ax[1, 0].set(xlabel='Distance from leading edge [m]', ylabel='Cf', title='Discrete wall friction')
ax[1, 1].set(xlabel='Distance from leading edge [m]', ylabel='y+', title='Owner-centre wall distance')
ax[1, 2].set(xlabel='SIMPLE iteration', ylabel='Momentum residual', title='Audited completed solves')
for axis in [ax[0, 2], *ax[1]]:
    axis.grid(alpha=.2)
    axis.legend(fontsize=8)
failed = ', '.join(f"{f['grid'][0]} x {f['grid'][1]} graded" for f in study['failures'])
fig.suptitle('Near-wall refinement | steady incompressible SST-2003m | Re_plate = 500', fontsize=17)
fig.supxlabel(f"Completed: {len(cases)} | Rejected solves: {failed or 'none'}\n"
    'Same physical inputs; no high-Re reference or mesh-independence claim. Failed cases have no accepted flow fields.', fontsize=10)
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, dpi=150)
