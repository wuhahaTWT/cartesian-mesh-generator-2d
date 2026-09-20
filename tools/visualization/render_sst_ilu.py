#!/usr/bin/env python3
"""Plot recorded SST preconditioner timings and an audited actual near-wall field."""
import argparse
import csv
import hashlib
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
run = study['acceptedRuns']['cells5120']['ilu0']
if not run['audit']['valid'] or not run['converged']:
    raise ValueError('audited converged medium field required')
mesh_path = Path(run['mesh'])
field_path = Path(run['files']['.cells.csv']['path'])
for path, expected in [(mesh_path, run['meshSha256']),
                       (field_path, run['files']['.cells.csv']['sha256'])]:
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError('evidence hash mismatch: ' + str(path))
mesh = native.read_cm2d(mesh_path)
polygons = [[mesh.vertices[i] for i in cell.vertices] for cell in mesh.cells]
with field_path.open() as stream:
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
axes[0, 1].set_title('Audited ILU(0) field (vertical zoom)')
for offset, method, label, color in [(-.18, 'jacobi', 'Jacobi', '#a1abb5'), (.18, 'ilu0', 'ILU(0)', '#327da2')]:
    values = [study['acceptedRuns'][key][method]['solveSeconds'] for key in ('cells2560', 'cells5120')]
    bars = axes[1, 0].bar([x+offset for x in range(2)], values, .36, label=label, color=color)
    axes[1, 0].bar_label(bars, fmt='%.2f s', padding=3)
axes[1, 0].set(xticks=[0, 1], xticklabels=['2,560 cells', '5,120 cells'], ylabel='Native solve time [s]',
               title='Complete solves: original gates passed', ylim=(0, 46))
axes[1, 0].legend()
fixed = study['fixed12800Throughput']
values = [fixed[key]['solveSeconds'] for key in ('jacobi', 'ilu0')]
bars = axes[1, 1].bar([0, 1], values, color=['#a1abb5', '#327da2'], width=.5)
axes[1, 1].bar_label(bars, fmt='%.2f s', padding=3)
axes[1, 1].set(xticks=[0, 1], xticklabels=['Jacobi', 'ILU(0)'], ylabel='Native solve time [s]',
               title='12,800 cells / 100 iterations: NOT converged', ylim=(0, 58))
for ax in axes[1]:
    ax.grid(axis='y', alpha=.2)
    ax.set_axisbelow(True)
fig.suptitle('SST-2003m | Re_plate = 10 million | optional scalar ILU(0)', fontsize=15)
fig.supxlabel('Same host, serial observations. Field agreement checked; not a physical accuracy qualification.\nFull 12,800-cell run: 90-second timeout, no accepted final field.', fontsize=10)
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, dpi=150)
