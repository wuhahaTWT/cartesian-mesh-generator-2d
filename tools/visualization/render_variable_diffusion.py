#!/usr/bin/env python3
"""Plot audited variable-diffusion fields on actual native final meshes."""
import argparse
import csv
import json
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools/verification'))
import verify_native_flow as native


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    study = json.loads(args.study.read_text())
    if not study['valid']:
        raise ValueError('A fully audited study is required')
    case = study['shear']['spatial']['limited-linear'][-1]
    prefix = case['prefix']
    info = json.loads(Path(prefix+'.json').read_text())
    mesh = native.read_cm2d(Path(info['mesh']))
    with Path(prefix+'.cells.csv').open() as stream:
        cells = list(csv.DictReader(stream))
    with Path(prefix+'.faces.csv').open() as stream:
        faces = list(csv.DictReader(stream))
    if len(cells) != len(mesh.cells) or len(faces) != len(mesh.edges):
        raise ValueError('Field dimensions do not match actual mesh')
    for row, cell in zip(cells, mesh.cells):
        if int(row['cell']) != cell.id:
            raise ValueError('Field ordering mismatch')
    polygons = [[mesh.vertices[v] for v in c.vertices] for c in mesh.cells]
    values = [float(r['value']) for r in cells]
    errors = [abs(float(r['value'])-float(r['exact'])) for r in cells]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
    for ax, data, cmap, title in (
        (axes[0,0], values, 'viridis', f'Computed scalar on {len(cells):,} actual cells'),
        (axes[0,1], errors, 'magma', 'Absolute error against manufactured solution')):
        collection = PolyCollection(polygons, array=data, cmap=cmap, edgecolors='none', rasterized=True)
        ax.add_collection(collection)
        ax.autoscale(); ax.set_aspect('equal')
        ax.set(title=title, xlabel='x', ylabel='y')
        fig.colorbar(collection, ax=ax, shrink=.75)
    segments = [[mesh.vertices[e.v0], mesh.vertices[e.v1]] for e in mesh.edges]
    collection = LineCollection(segments, array=[float(r['diffusivity']) for r in faces], cmap='cividis', linewidths=.7)
    ax = axes[1,0]; ax.add_collection(collection)
    ax.set(xlim=(-.015,.29), ylim=(-.015,.29), xlabel='x', ylabel='y', title='Actual cut-cell detail; color = prescribed face D')
    ax.set_aspect('equal'); fig.colorbar(collection, ax=ax, shrink=.75)
    ax = axes[1,1]
    for domain, color in [('square','#377b90'), ('shear','#aa6442')]:
        for scheme, style in [('upwind','--'), ('limited-linear','-')]:
            rows = study[domain]['spatial'][scheme]
            ax.loglog([r['cells']**-.5 for r in rows], [r['l2Error'] for r in rows], 'o'+style,
                      color=color, label=f'{domain}, {scheme}')
    ax.set(xlabel='Effective spacing sqrt(area / cells); area = 1', ylabel='Area-weighted RMS error', title='Same PDE and stopping criteria on all grids')
    ax.grid(True, which='both', alpha=.18); ax.legend(fontsize=9)
    fig.suptitle('Conservative transport with spatial diffusivity', fontsize=16)
    fig.supxlabel('D(x) = 0.08(1+x), U = (0.35, 0), exact scalar = sin(pi x) sin(pi y)\n'
                  'Independent geometry, face constitutive flux and cell balance audits passed.\n'
                  'Prescribed coefficients on frozen carrier; not a turbulence model.', fontsize=10)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=150); plt.close(fig)


if __name__ == '__main__':
    main()
