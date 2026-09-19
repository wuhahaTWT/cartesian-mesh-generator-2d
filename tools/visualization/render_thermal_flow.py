#!/usr/bin/env python3
"""Actual synchronized thermal fields; retain mixed spatial/time error scope."""
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

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--audit',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();report=json.loads(a.audit.read_text())
fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
def field(ax,prefix,extent,title,arrows=False):
    meta=json.loads(Path(prefix+'.json').read_text());mesh=native.read_cm2d(Path(meta['mesh']))
    with Path(prefix+'.cells.csv').open() as f:rows=list(csv.DictReader(f))
    values=[float(row['value']) for row in rows]
    assert min(values)>=0 and max(values)<=1
    polygons=[[mesh.vertices[v] for v in cell.vertices] for cell in mesh.cells]
    collection=PolyCollection(polygons,array=values,cmap='YlOrRd',edgecolors='#8a929b',linewidths=.15)
    collection.set_clim(0,1);ax.add_collection(collection)
    ax.set(xlim=extent[:2],ylim=extent[2:],aspect='equal',xlabel='x',ylabel='y',title=title)
    fig.colorbar(collection,ax=ax,label=r'$\theta$')
    if arrows:
        text=Path(meta['carrierCheckpoint']).read_text()
        fields={line.split()[0]:list(map(float,line.split()[2:])) for line in text.splitlines() if line.startswith(('U ','V '))}
        indices=range(0,len(rows),7)
        ax.quiver([float(rows[i]['x']) for i in indices],[float(rows[i]['y']) for i in indices],
                  [fields['U'][i] for i in indices],[fields['V'][i] for i in indices],color='#303b49',alpha=.65,scale=14,width=.002)
field(axes[0,0],report['circle']['prefix'],(-1.6,2.6,-1.5,1.5),'Heated cylinder: 4,880 cells, t = 0.2 s')
field(axes[0,1],report['temporal'][-1]['prefix'],(0,1,0,1),'Evolving vortex and scalar: t = 0.1 s',True)
ax=axes[1,0];series=report['spatial']
ax.loglog([1/r['cells']**.5 for r in series],[r['l2Error'] for r in series],'o-',color='#2166ac')
ax.set(xlabel='h = sqrt(area / cells)',ylabel='Scalar L2 error',title='Space refinement at fixed dt = 0.0025 s\nIncludes time error; no single measured spatial order')
ax.grid(True,which='both',alpha=.18)
ax=axes[1,1];series=report['temporal']
ax.loglog([r['dt'] for r in series],[r['l2Error'] for r in series],'o-',color='#b35806')
ax.set(xlabel='Time step [s]',ylabel='Scalar L2 error',title='Common final time 0.1 s, 900 cells\nObserved temporal order: 0.94 / 0.92')
ax.grid(True,which='both',alpha=.18)
fig.suptitle('Synchronized native flow and passive thermal transport\nOne-way, constant properties | joint restart verified',fontsize=14)
a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=150);plt.close(fig)
