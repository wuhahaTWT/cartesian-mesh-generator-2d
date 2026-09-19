#!/usr/bin/env python3
"""Actual signed outlet flux and independently measured counterflow errors."""
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
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--circle-audit',type=Path,required=True)
p.add_argument('--counterflow-audits',type=Path,nargs=3,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
circle=json.loads(a.circle_audit.read_text())
series=[json.loads(path.read_text()) for path in a.counterflow_audits]
if not all(r['valid'] for r in [circle,*series]):
    raise ValueError('Only verified outputs may be rendered')
def rows(prefix,suffix):
    with Path(str(prefix)+suffix).open() as stream:
        return list(csv.DictReader(stream))
mesh=native.read_cm2d(Path(circle['mesh']))
measured=native.measure(mesh,1e-10,1e-9)
geo=native.face_geometry(mesh,measured)
cells=rows(circle['prefix'],'.cells.csv');faces=rows(circle['prefix'],'.faces.csv')
summary=json.loads(Path(circle['prefix']+'.json').read_text())
fig,axes=plt.subplots(1,3,figsize=(15.5,5.3),layout='constrained',gridspec_kw={'width_ratios':[1.25,1,1]})
ax=axes[0]
polygons=[[mesh.vertices[v] for v in cell.vertices] for cell in mesh.cells]
pc=PolyCollection(polygons,array=[float(r['u']) for r in cells],cmap='viridis',
                  edgecolors='#777f86',linewidths=.08)
ax.add_collection(pc);ax.autoscale_view();ax.set_aspect('equal')
xmin,ymin,xmax,ymax=measured.bounds
for e,g,row in zip(mesh.edges,geo,faces):
    if e.neighbour<0 and abs(g.centre[0]-xmax)<1e-10:
        q=float(row['flux']);vx=q/math.hypot(*g.area_vector)
        ax.arrow(g.centre[0],g.centre[1],.14*vx,0,head_width=.035,
                 length_includes_head=True,color='#cd522b' if q<0 else '#2374ad',linewidth=.7)
ax.set(xlim=(xmin-.1,xmax+.4),ylim=(ymin-.1,ymax+.1),xlabel='x [m]',ylabel='y [m]',
       title=f'Compact cylinder: Re = 40, {len(cells):,} cells\n{summary["outletBackflowFaces"]} incoming pressure-outlet faces')
fig.colorbar(pc,ax=ax,label='Streamwise velocity [m/s]',shrink=.7,pad=.02)
ax=axes[1]
ys=[i/300 for i in range(301)]
ax.plot([1+2*math.cos(2*math.pi*y) for y in ys],ys,color='#30363d',label='Exact')
errors=[];sizes=[]
for color,r in zip(['#d18a24','#4f9b82','#4776ac'],series):
    data=rows(r['prefix'],'.cells.csv')
    x=min((float(c['x']) for c in data),key=lambda x:abs(x-.5))
    profile=sorted((c for c in data if abs(float(c['x'])-x)<1e-12),key=lambda c:float(c['y']))
    ax.plot([float(c['u']) for c in profile],[float(c['y']) for c in profile],'.',color=color,
            markersize=4,label=f'{len(data):,} cells')
    errors.append(r['benchmark']['velocityL2OverSpeed']);sizes.append(1/math.sqrt(len(data)))
ax.axvline(0,color='#9aa4ae',linewidth=.7)
ax.set(xlabel='Streamwise velocity / Uref',ylabel='y / H',title='Forced counterflow: sampled near x/H = 0.5')
ax.legend(fontsize=9);ax.grid(alpha=.15)
ax=axes[2];ax.loglog(sizes,errors,'o-',color='#4776ac')
orders=[math.log(errors[i]/errors[i+1])/math.log(sizes[i]/sizes[i+1]) for i in range(2)]
ax.set(xlabel='h / H',ylabel='Velocity vector L2 / Uref',
       title=f'Full-field spatial error\nObserved orders: {orders[0]:.2f} / {orders[1]:.2f}')
ax.grid(True,which='both',alpha=.2)
fig.suptitle('Pressure outlet with simultaneous inflow and outflow',fontsize=15)
fig.supxlabel('Opt-in normal-inlet model: prescribed pressure, zero incoming tangential velocity.\n'
              'Manufactured counterflow has an analytic body force; the compact cylinder is not a far-field accuracy benchmark.',fontsize=10)
a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=150);plt.close(fig)
