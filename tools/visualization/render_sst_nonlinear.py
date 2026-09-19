#!/usr/bin/env python3
"""Plot audited nonlinear SST fields on actual cells; no invented flow solution."""
import argparse,csv,json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();study=json.loads(a.study.read_text());r=study['cases']['shear'][-1]
if not study['valid'] or not r['valid']:raise ValueError('independent audit required')
mesh=native.read_cm2d(Path(r['mesh']))
def rows(suffix):
    with Path(r['prefix']+suffix).open() as f:return list(csv.DictReader(f))
cells=rows('.cells.csv');history=rows('.history.csv')
polygons=[[mesh.vertices[i] for i in c.vertices] for c in mesh.cells]
fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
for ax,key,title,cmap,label in zip(axes.flat,['k','omega','nuT'],
        ['Turbulent kinetic energy','Specific dissipation rate','Turbulent viscosity'],
        ['viridis','cividis','magma'],['k [m2/s2]','omega [1/s]','nu_t [m2/s]']):
    collection=PolyCollection(polygons,array=[float(row[key]) for row in cells],cmap=cmap,
                              edgecolors='#697681',linewidths=.12,rasterized=True)
    ax.add_collection(collection);ax.autoscale();ax.set_aspect('equal')
    ax.set(title=title,xlabel='x [m]',ylabel='y [m]');fig.colorbar(collection,ax=ax,shrink=.8,label=label)
ax=axes[1,1]
for key,label in [('kCellResidual','k'),('omegaCellResidual','omega')]:
    ax.semilogy([int(row['iteration']) for row in history],[float(row[key]) for row in history],label=label)
ax.axhline(1e-10,ls='--',color='#777777',lw=.8,label='Probe stopping target')
ax.grid(alpha=.2);ax.legend();ax.set(title='Residual at the updated nonlinear state',
    xlabel='Nonlinear iteration',ylabel='Max. diagonal-scaled cell residual')
fig.suptitle(f'Nonlinear SST transport on {len(cells):,} actual Cut-cells',fontsize=16)
fig.supxlabel('Prescribed U=(y,0), one backward-Euler step dt=0.02 s; resolved bottom-wall k/omega boundary.\n'
    'Fixed carrier only: no momentum feedback, no boundary-layer accuracy or engineering RANS qualification.',fontsize=10)
fig.savefig(a.output,dpi=150);plt.close(fig)
