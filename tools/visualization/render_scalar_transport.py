#!/usr/bin/env python3
"""Render actual CM2D scalar fields and archived independent refinement evidence."""
import argparse
import csv
import json
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--circle-prefix',type=Path,required=True)
parser.add_argument('--refinement',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
summary=json.loads(Path(str(args.circle_prefix)+'.json').read_text())
mesh=native.read_cm2d(Path(summary['mesh']))
with Path(str(args.circle_prefix)+'.cells.csv').open() as f:
    rows=list(csv.DictReader(f))
assert len(rows)==len(mesh.cells)
polygons=[[mesh.vertices[v] for v in cell.vertices] for cell in mesh.cells]
values=[float(row['value']) for row in rows]
assert min(values)>=0 and max(values)<=1  # this figure's fixed colour range must not conceal overshoots
ref=json.loads(args.refinement.read_text())
fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
for ax,extent,title in ((axes[0,0],(-2,6,-2,2),'Actual cut-cell temperature field'),
                        (axes[0,1],(.65,1.9,-.65,.65),'Near-wall / downstream detail')):
    collection=PolyCollection(polygons,array=values,cmap='YlOrRd',edgecolors='#8a929b',linewidths=.16)
    collection.set_clim(0,1);ax.add_collection(collection)
    ax.set(xlim=extent[:2],ylim=extent[2:],xlabel='x',ylabel='y',title=title,aspect='equal')
    fig.colorbar(collection,ax=ax,label=r'$\theta$')
ax=axes[1,0]
for scheme,series in ref['spatial'].items():
    ax.loglog([1/r['cells']**.5 for r in series],[r['l2Error'] for r in series],'o-',label=scheme)
ax.set(xlabel='h = sqrt(area / cells)',ylabel='Volume-weighted scalar L2 error',title='Spatial MMS: three independently read meshes')
ax.grid(True,which='both',alpha=.18);ax.legend()
ax=axes[1,1];series=ref['temporal']
ax.loglog([r['dt'] for r in series],[r['l2Error'] for r in series],'o-',color='#2166ac',label='Backward Euler')
ax.set(xlabel='Time step [s]',ylabel='Scalar L2 error against analytic decay',title='Temporal refinement: common final time 0.5 s')
ax.grid(True,which='both',alpha=.18);ax.legend()
fig.suptitle('Native passive thermal transport\n4,880 cells | frozen carrier at t = 0.5 s | D = 0.1 | wall 1, inflow 0',fontsize=14)
args.output.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(args.output,dpi=160);plt.close(fig)
