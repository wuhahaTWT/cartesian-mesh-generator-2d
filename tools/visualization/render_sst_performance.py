#!/usr/bin/env python3
"""Plot measured SST phases alongside the actual independently audited mesh/field."""
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
if not study['valid']:
    raise ValueError('audited, numerically identical fields required')
runs = study['runs']
selected = runs['medium-after']
mesh = native.read_cm2d(Path(selected['mesh']))
polygons = [[mesh.vertices[i] for i in cell.vertices] for cell in mesh.cells]
with Path(selected['prefix'] + '.cells.csv').open() as stream:
    fields = list(csv.DictReader(stream))
if len(fields) != len(polygons):
    raise ValueError('field/mesh size mismatch')
fig, axes = plt.subplots(2, 2, figsize=(13, 8), layout='constrained')
axes[0, 0].add_collection(PolyCollection(polygons, facecolors='#f4f6f8', edgecolors='#536b80', linewidths=.3))
speed = PolyCollection(polygons, array=[float(r['speed']) for r in fields], cmap='viridis', edgecolors='none')
axes[0, 1].add_collection(speed)
fig.colorbar(speed, ax=axes[0, 1], label='Speed [m/s]', shrink=.8)
for ax in axes[0]:
    ax.set(xlim=(-.5, 2), ylim=(0, .025), xlabel='x [m]', ylabel='y [m]')
axes[0, 0].set_title('Actual 5,120-cell mesh (vertical zoom)')
axes[0, 1].set_title('Returned field unchanged byte for byte')

pairs = [('profiled', 'final-small'), ('medium-before', 'medium-after')]
for offset, index, label, color in [(-.18, 0, 'Before', '#a1abb5'), (.18, 1, 'After', '#327da2')]:
    values = [runs[p[index]]['native']['solveSeconds'] for p in pairs]
    bars = axes[1, 0].bar([x+offset for x in range(2)], values, .36, label=label, color=color)
    axes[1, 0].bar_label(bars, fmt='%.2f s', padding=3)
axes[1, 0].set(xticks=[0, 1], xticklabels=['2,560 cells', '5,120 cells'], ylabel='Native solve time [s]', title='Same inputs, iterations and convergence gates')
axes[1, 0].legend()

def phases(run):
    total = run['native']['solveSeconds']
    p = run['native']['performance']
    solve, evaluate = p['scalarSolves'], p['scalarEvaluations']
    setup = solve['setupSeconds'] + evaluate['setupSeconds']
    flux = solve['faceFluxSeconds'] + evaluate['faceFluxSeconds']
    linear = solve['linearSeconds'] + evaluate['linearSeconds']
    flow_linear = p['pressureLinearSeconds'] + p['momentumLinearSeconds']
    values = [setup, linear, flux, p['sstUpdateSeconds']-setup-linear-flux,
              flow_linear, total-p['sstUpdateSeconds']-flow_linear]
    if min(values) < 0:
        raise ValueError('overlapping or invalid timing decomposition')
    return values

before, after = [phases(runs[n]) for n in pairs[1]]
labels = ['Scalar setup', 'Scalar linear solves', 'Scalar face fluxes', 'Other SST work', 'Flow linear solves', 'Other flow work']
colors = ['#d29d53', '#4479a6', '#629f90', '#98b0bb', '#7f78a8', '#c9ccd0']
bottom = [0., 0.]
for i, (label, color) in enumerate(zip(labels, colors)):
    values = [before[i], after[i]]
    axes[1, 1].bar([0, 1], values, bottom=bottom, color=color, label=label, width=.5)
    bottom = [a+b for a, b in zip(bottom, values)]
axes[1, 1].set(xticks=[0, 1], xticklabels=['Before', 'After'], ylabel='Seconds', title='5,120 cells: disjoint measured phases')
axes[1, 1].legend(fontsize=8, ncols=2, loc='upper right')
for ax in axes[1]:
    ax.set_ylim(0, runs['medium-before']['native']['solveSeconds']*1.18)
    ax.grid(axis='y', alpha=.2)
    ax.set_axisbelow(True)
fig.suptitle('SST-2003m | Re_plate = 10 million | lower overhead, identical numerical results', fontsize=15)
fig.supxlabel('Single serial observations on one macOS host; independent field audits passed.\nNo change in physical accuracy; previous 12,800-cell timeout remains unqualified.', fontsize=10)
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, dpi=150)
