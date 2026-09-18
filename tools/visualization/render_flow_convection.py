#!/usr/bin/env python3
"""Plot real native fields and independently verified convection comparisons."""
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


def load(path):
    report = json.loads(path.read_text())
    if not report['cases'] or any(not case['valid'] for case in report['cases']):
        raise ValueError(f'Refusing fields with failed independent case checks: {path}')
    return {case['label']: case for case in report['cases']}


def polygons_and_fields(case):
    mesh = read_cm2d(Path(case['mesh']))
    prefix = Path(case['prefix'])
    with Path(str(prefix)+'.cells.csv').open() as f:
        rows = sorted(csv.DictReader(f), key=lambda row: int(row['cell']))
    polygons = [[mesh.vertices[i] for i in cell.vertices] for cell in mesh.cells]
    return polygons, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upwind', type=Path, required=True)
    parser.add_argument('--limited', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    base, high = load(args.upwind), load(args.limited)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), layout='constrained')
    for ax, label, title in [(axes[0, 0], 'circle-re20', 'Re20 circle: limited reconstruction'),
                             (axes[0, 1], 'cavity-l6', 'Re100 cavity: limited reconstruction')]:
        polys, rows = polygons_and_fields(high[label])
        speed = np.array([float(row['speed']) for row in rows])
        collection = PolyCollection(polys, array=speed, cmap='viridis',
                                    edgecolors=(.12,.16,.21,.22), linewidths=.18)
        ax.add_collection(collection)
        ax.autoscale_view(); ax.set_aspect('equal')
        if label == 'circle-re20':
            ax.set_xlim(-3, 5); ax.set_ylim(-3, 3)
            ax.text(.02,.02,'Near-field crop of the full 42 x 42 domain',
                    transform=ax.transAxes,fontsize=8,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
        else:
            stride = max(1, len(rows)//180)
            sample = rows[::stride]
            ax.quiver([float(r['x']) for r in sample], [float(r['y']) for r in sample],
                      [float(r['u']) for r in sample], [float(r['v']) for r in sample],
                      color='white', scale=10, width=.0022)
        ax.set_title(f'{title}\n{len(rows):,} actual polygon cells')
        ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
        fig.colorbar(collection, ax=ax, label='Speed (m/s)', shrink=.85)
    ax = axes[1, 0]
    for label, color in [('cavity-l4','#6e99b5'),('cavity-l5','#378079'),('cavity-l6','#b46027')]:
        samples = [s for s in high[label]['benchmark']['samples'] if s['field']=='u']
        samples.sort(key=lambda s:s['coordinate'])
        ax.plot([s['actual'] for s in samples], [s['coordinate'] for s in samples],
                color=color,label=f"Limited, {high[label]['counts']['cells']:,} cells")
    samples = [s for s in base['cavity-l6']['benchmark']['samples'] if s['field']=='u']
    samples.sort(key=lambda s:s['coordinate'])
    ax.plot([s['actual'] for s in samples], [s['coordinate'] for s in samples],
            '--',color='#6a6a6a',label='Upwind, fine')
    ax.scatter([s['reference'] for s in samples], [s['coordinate'] for s in samples],
               s=20,facecolors='none',edgecolors='black',label='Ghia et al. (1982)')
    ax.set_xlabel('u / lid speed'); ax.set_ylabel('y / H'); ax.set_title('Cavity vertical centreline')
    ax.grid(alpha=.18); ax.legend(fontsize=8)
    ax = axes[1, 1]
    labels = ['cavity-l4','cavity-l5','cavity-l6']
    for cases, name, color in [(base,'First-order upwind','#64778b'),
                               (high,'Limited linear reconstruction','#b46027')]:
        ax.loglog([cases[l]['counts']['cells'] for l in labels],
                  [cases[l]['benchmark']['centrelineRmse'] for l in labels],
                  'o-',color=color,label=name)
    ax.set_xlabel('Final cell count'); ax.set_ylabel('Combined u/v centreline RMSE / lid speed')
    ax.set_title('Same geometry and stopping criteria')
    ax.grid(alpha=.18,which='both'); ax.legend(fontsize=8)
    report=json.loads(args.limited.read_text())
    if not report['sequenceChecks']['valid']:
        ax.text(.04,.06,'Limited reconstruction: finest-mesh\nrefinement trend NOT qualified',transform=ax.transAxes,
                fontsize=9,color='#9a3f24',bbox={'facecolor':'white','alpha':.9,'edgecolor':'#d8b3a6'})
    fig.suptitle('Shared-face pressure and bounded convection: native 2D FVM',fontsize=16)
    fig.supxlabel('Real CM2D / CSV fields. Steady laminar checks; no turbulence or grid-independence certification.',fontsize=9)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    main()
