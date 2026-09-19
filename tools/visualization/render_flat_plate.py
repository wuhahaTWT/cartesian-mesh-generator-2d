#!/usr/bin/env python3
"""Plot audited mixed-boundary SST diagnostics, not a TMR benchmark claim."""
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
if not study['valid'] or not all(c['valid'] for c in study['cases']):
    raise ValueError('independent final-state audits required')
cases = study['cases']
finest = max(c['cells'] for c in cases)
selected = [next(c for c in cases if c['cells']==finest and c['flatPlateTop']==mode)
            for mode in ('pressure-farfield', 'symmetry')]

def rows(case, suffix):
    with Path(case['prefix']+suffix).open() as stream:
        return list(csv.DictReader(stream))

fields = [rows(c, '.cells.csv') for c in selected]
speed_max = max(float(r['speed']) for field in fields for r in field)
fig, axes = plt.subplots(2, 3, figsize=(16, 9), layout='constrained')
for ax, case, field, key, title in zip(axes[0],
        [selected[0], selected[1], selected[0]], [fields[0], fields[1], fields[0]],
        ['speed', 'speed', 'p'], ['Open pressure top: speed', 'Symmetry top: speed', 'Open pressure top: p / rho']):
    mesh = native.read_cm2d(Path(case['mesh']))
    polygons = [[mesh.vertices[v] for v in c.vertices] for c in mesh.cells]
    collection = PolyCollection(polygons, array=[float(r[key]) for r in field],
        cmap='viridis' if key=='speed' else 'coolwarm', edgecolors='#627280', linewidths=.2, rasterized=True)
    if key=='speed': collection.set_clim(0, speed_max)
    ax.add_collection(collection); ax.autoscale(); ax.set_aspect('equal')
    ax.plot([0, .5], [0, 0], '--', color='#cd9c45', linewidth=2)
    ax.plot([.5, 1], [0, 0], color='#20252a', linewidth=3)
    ax.set(title=title, xlabel='x [m]', ylabel='y [m]')
    fig.colorbar(collection, ax=ax, shrink=.85, label='m/s' if key=='speed' else 'm2/s2')

for case in cases:
    sample = case['plateWallSamples']
    label = f"{case['cells']} cells, "+('open' if case['flatPlateTop']=='pressure-farfield' else 'symmetry')
    style = '-' if case['flatPlateTop']=='pressure-farfield' else '--'
    for ax, key in zip(axes[1, :2], ['Cf', 'yPlus']):
        ax.plot([r['xFromLeadingEdge'] for r in sample], [r[key] for r in sample],
                style, marker='.', label=label)
for ax, title, ylabel in zip(axes[1, :2], ['Discrete wall friction', 'Owner-centre wall distance'], ['Cf', 'y+']):
    ax.set(title=title, xlabel='x - leading edge [m]', ylabel=ylabel)
    ax.grid(alpha=.2); ax.legend(fontsize=8)

history = rows(selected[0], '.history.csv')
ax = axes[1, 2]
for key, label in [('momentumResidual','Momentum'), ('kCellResidual','k'), ('omegaCellResidual','omega')]:
    ax.semilogy([int(r['iteration']) for r in history],
               [max(float(r[key]),1e-18) for r in history], label=label)
ax.set(title='Open top: current equation residuals', xlabel='SIMPLE iteration',
       ylabel='Diagonal-scaled residual (equation-specific units)')
ax.legend(); ax.grid(alpha=.2)
fig.suptitle(f'Experimental SST-2003m flat-plate boundaries | {finest:,} actual cells in field panels', fontsize=16)
fig.supxlabel('U=1 m/s, nu=0.001 m2/s, inlet k=0.001 m2/s2, omega=2 1/s; plate x=0.5 to 1 (Re_length=500).\n'
    'Dashed bottom: symmetry; solid bottom: no slip. Discrete-equation audits passed; no high-Re, TMR or physical-accuracy qualification.', fontsize=9)
fig.savefig(args.output, dpi=150); plt.close(fig)
