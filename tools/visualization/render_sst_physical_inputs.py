#!/usr/bin/env python3
"""Plot audited high-Re diagnostic fields; retain independent-audit failures."""
import argparse
import csv
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--study',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
study=json.loads(args.study.read_text())
cases=[c for c in study['cases'] if c['plateReynolds']==1e7]
if not cases or not all(c['valid'] for c in cases):raise ValueError('independent accepted fields required')
selected=max(cases,key=lambda c:c['cells'])
mesh=native.read_cm2d(Path(selected['mesh']))
polygons=[[mesh.vertices[v] for v in cell.vertices] for cell in mesh.cells]
with Path(selected['prefix']+'.cells.csv').open() as stream:field=list(csv.DictReader(stream))
fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
axes[0,0].add_collection(PolyCollection(polygons,facecolors='#f4f6f8',edgecolors='#536b80',linewidths=.35))
axes[0,0].set(xlim=(-.5,2),ylim=(0,.025),xlabel='x [m]',ylabel='y [m]',title=f"Actual {selected['cells']}-cell mesh (vertical zoom)")
axes[0,0].plot([-.5,0],[0,0],'--',color='#d2a23f',lw=3)
axes[0,0].plot([0,2],[0,0],color='#222',lw=3)
collection=PolyCollection(polygons,array=[float(r['speed']) for r in field],cmap='viridis',edgecolors='none')
axes[0,1].add_collection(collection)
axes[0,1].set(xlim=(-.5,2),ylim=(0,.025),xlabel='x [m]',ylabel='y [m]',title='Returned speed [m/s] (vertical zoom)')
fig.colorbar(collection,ax=axes[0,1],shrink=.8)
for c in cases:
 samples=c['plateWallSamples'];label=f"{c['nx']} x {c['ny']}, stretch={c['stretch']}"
 axes[1,0].plot([r['xFromLeadingEdge'] for r in samples],[r['Cf'] for r in samples],'.-',label=label)
 axes[1,1].plot([r['xFromLeadingEdge'] for r in samples],[r['yPlus'] for r in samples],'.-',label=label)
axes[1,0].set(xlabel='Distance from leading edge [m]',ylabel='Cf',title='Wall friction remains grid-sensitive')
axes[1,1].set(xlabel='Distance from leading edge [m]',ylabel='y+',yscale='log',title='First-cell wall resolution remains insufficient')
axes[1,1].axhline(1,color='#777',ls='--',lw=1,label='y+ = 1 reference')
for ax in axes[1]:ax.grid(alpha=.2);ax.legend(fontsize=8)
failed=', '.join(c['name'] for c in study['failures']) or 'none'
fig.suptitle('Incompressible SST-2003m | Re_plate = 10 million | diagnostic only',fontsize=16)
fig.supxlabel(f'Independent audit rejected: {failed}\nNo TMR reference match or high-Re accuracy qualification; rejected fields are not plotted.',fontsize=10)
args.output.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(args.output,dpi=150)
