#!/usr/bin/env python3
"""Plot actual circular Couette fields, mesh and independently audited errors."""
import argparse
import json
import math
from pathlib import Path
import sys
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_native_flow as native

def main(args):
    report=json.loads((args.study/'summary.json').read_text())
    if not report['valid']:raise ValueError('Study has unresolved acceptance issues')
    mesh=native.read_cm2d(args.mesh);measured=native.measure(mesh,1e-11,1e-9)
    cells=native.read_cells(Path(str(args.cutcell_prefix)+'.cells.csv'),mesh,measured,'custom')
    audit=json.loads(Path(str(args.cutcell_prefix)+'-audit.json').read_text())
    if not audit['valid']:raise ValueError('Cut-cell field has not passed independent audit')
    fig,axes=plt.subplots(1,3,figsize=(15.8,5.3),layout='constrained')
    polygons=[[mesh.vertices[i] for i in c.vertices] for c in mesh.cells]
    coll=PolyCollection(polygons,array=[math.hypot(c['u'],c['v']) for c in cells],cmap='viridis',edgecolors='#182830',linewidths=.18)
    axes[0].add_collection(coll);axes[0].autoscale();axes[0].set_aspect('equal');fig.colorbar(coll,ax=axes[0],label='Speed [m/s]',shrink=.77)
    axes[0].set(xlabel='x [m]',ylabel='y [m]',title=f'Actual Cut-cell solution: {len(cells):,} cells')
    fine=report['cases'][-1];mf=native.read_cm2d(Path(fine['mesh']));measure=native.measure(mf,1e-11,1e-9)
    cf=native.read_cells(Path(fine['prefix']+'.cells.csv'),mf,measure,'custom')
    radial={}
    for c in cf:
        r=math.hypot(c['x'],c['y']);key=round(r,10)
        radial.setdefault(key,[]).append((-c['u']*c['y']+c['v']*c['x'])/r)
    radius=sorted(radial)
    axes[1].plot(radius,[math.fsum(radial[r])/len(radial[r]) for r in radius],'o',ms=3,label=f'Fine ring grid ({len(cf):,})')
    exact=[.5+.5*i/200 for i in range(201)]
    axes[1].plot(exact,[-r/3+1/(3*r) for r in exact],color='black',lw=1,label='Circular analytic solution')
    axes[1].set(xlabel='Radius [m]',ylabel='Azimuthal velocity [m/s]',title='Solver reference on ring quadrilaterals');axes[1].legend(fontsize=8);axes[1].grid(alpha=.2)
    count=[c['metrics']['cells'] for c in report['cases']]
    axes[2].loglog(count,[c['metrics']['velocityL2'] for c in report['cases']],'o-',label='Smooth wall samples')
    previous=[json.loads((args.baseline/f'l{k}'/'metrics.json').read_text()) for k in (5,6,7)]
    axes[2].loglog([c['cells'] for c in previous],[c['velocityL2'] for c in previous],'s--',label='Piecewise constant wall')
    axes[2].set(xlabel='Ring quadrilateral cells',ylabel='Velocity L2 / inner-wall speed',title='Refined geometry and mesh');axes[2].legend(fontsize=8);axes[2].grid(which='both',alpha=.2)
    fig.suptitle('Incompressible rotating annulus: inner wall moves, outer wall fixed',fontsize=16)
    fig.supxlabel('Ri = 0.5 m; Ro = 1 m; inner speed = 0.5 m/s; nu = 0.1 m²/s.\nRing-grid reference accuracy does not qualify Cut-cell pressure accuracy.',fontsize=10)
    args.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(args.output,dpi=150);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--mesh',type=Path,required=True)
    p.add_argument('--cutcell-prefix',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
