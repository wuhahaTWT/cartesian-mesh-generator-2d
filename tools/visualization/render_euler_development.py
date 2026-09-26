#!/usr/bin/env python3
"""Plot only retained numerical fields and analytic reference solutions."""
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
from matplotlib.ticker import NullFormatter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_euler as euler
import verify_native_flow as geometry


def render(root, destination):
    evidence=json.loads((root/'validation.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(13.4,9.6),layout='constrained')
    xx=[i/1600 for i in range(1601)];exact=[euler.sod(x,.2) for x in xx]
    styles=[('rusanov',1,'#d77b31','Rusanov / first order'),('hllc',1,'#8d99a6','HLLC-HLLE / first order'),('hllc',2,'#007f88','HLLC-HLLE / second order')]
    for ax,key,k in [(axes[0,0],'rho',0),(axes[0,1],'p',2)]:
        ax.plot(xx,[q[k] for q in exact],color='#202833',lw=1.2,label='Exact Riemann solution')
        for scheme,order,color,label in styles:
            prefix=root/f'sod200-{scheme}-o{order}';euler.audit(root/'sod200.solver.cm2d',prefix)
            rows=list(csv.DictReader(Path(str(prefix)+'.cells.csv').open()))[:200]
            ax.plot([float(r['x']) for r in rows],[float(r[key]) for r in rows],color=color,lw=1.6 if order==2 else 1.2,label=label)
        ax.set(xlabel='x',ylabel='Density' if key=='rho' else 'Pressure',title=f'Sod shock tube | 200 x 4 actual cells | t = 0.2',xlim=(0,1))
        ax.grid(alpha=.18);ax.legend(fontsize=8)
    ax=axes[1,0];meshpath=root/'vortex128.solver.cm2d';prefix=root/'vortex128-hllc-o2';euler.audit(meshpath,prefix)
    mesh=geometry.read_cm2d(meshpath);rows=list(csv.DictReader(Path(str(prefix)+'.cells.csv').open()))
    polygons=[[mesh.vertices[i] for i in c.vertices] for c in mesh.cells]
    collection=PolyCollection(polygons,array=[float(r['rho']) for r in rows],cmap='viridis',edgecolors='none')
    ax.add_collection(collection);ax.set(xlim=(7,14),ylim=(7,14),aspect='equal',xlabel='x',ylabel='y',title='Isentropic vortex | actual density | 128 x 128 | t = 0.5')
    fig.colorbar(collection,ax=ax,shrink=.8,label='Density')
    ax=axes[1,1];runs=[r for r in evidence['cases'] if r['case'].startswith('vortex')]
    nn=[int(r['case'].split('-')[0][6:]) for r in runs]
    for key,color in [('rho','#007f88'),('u','#d77b31'),('p','#7869a9')]:
        errors=[r['areaWeightedL1'][key] for r in runs]
        order=evidence['vortexObservedOrders'][key]
        ax.loglog(nn,errors,'o-',color=color,lw=1.6,label=f'{key}: observed order {order:.2f} (64 to 128)')
    anchor=runs[1]['areaWeightedL1']['rho']
    ax.loglog(nn,[anchor*(64/n)**2 for n in nn],'--',color='#9da6ac',lw=1,label='Second-order reference slope')
    ax.set(xlabel='Cells per side',ylabel='Area-weighted L1 error',title='Smooth-vortex grid refinement | second-order method',xticks=nn)
    ax.set_xticklabels([str(n) for n in nn]);ax.xaxis.set_minor_formatter(NullFormatter());ax.grid(alpha=.18,which='both');ax.legend(fontsize=8)
    fig.suptitle('Native 2D compressible Euler: conservation, shock resolution and smooth-flow accuracy',fontsize=15,color='#20333a')
    fig.supxlabel('Ideal gas / inviscid development branch. Real solver fields; analytic references shown explicitly. No general engineering qualification.',fontsize=9)
    destination.parent.mkdir(parents=True,exist_ok=True);fig.savefig(destination,dpi=160);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--study',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();render(args.study,args.output)
