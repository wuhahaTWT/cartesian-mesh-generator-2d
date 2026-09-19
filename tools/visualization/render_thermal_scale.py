#!/usr/bin/env python3
"""Real coupled fields and continuous/discrete-time-reference errors."""
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
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--allow-incomplete',action='store_true',help='Show completed audited cases with a prominent failure notice')
    a=p.parse_args();study=json.loads(a.study.read_text())
    incomplete=not study.get('valid')
    if (incomplete and not a.allow_incomplete) or not study['cases'] or not all(c.get('valid') for c in study['cases']):
        raise ValueError('Requires a completed independently audited study')
    rows=[dict(c['reference']) for c in study['cases']]
    for coarse,fine in zip(rows,rows[1:]):
        if fine['h']>=coarse['h']:raise ValueError('Non-refining completed cases')
        for field in ('scalar','velocity'):
            e1,e2=coarse[field+'BackwardEulerL2'],fine[field+'BackwardEulerL2']
            fine[field+'BackwardEulerL2Order']=math.log(e1/e2)/math.log(coarse['h']/fine['h']) if e1>0 and e2>0 else None
    case=study['cases'][-1];prefix=case['prefix']
    info=json.loads(Path(prefix+'.json').read_text());mesh=native.read_cm2d(Path(info['mesh']))
    with Path(prefix+'.cells.csv').open() as f:cells=list(csv.DictReader(f))
    if len(cells)!=len(mesh.cells):raise ValueError('actual field and mesh mismatch')
    polys=[[mesh.vertices[v] for v in c.vertices] for c in mesh.cells]
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    ax=axes[0,0]
    pc=PolyCollection(polys,array=[float(c['value']) for c in cells],cmap='inferno',edgecolors='none',rasterized=True)
    ax.add_collection(pc);ax.autoscale_view();ax.set_aspect('equal')
    ax.set(xlabel='x [m]',ylabel='y [m]',title=f'Computed scalar: {len(cells):,} cells, t = {info["time"]:g} s')
    fig.colorbar(pc,ax=ax,shrink=.7,label='Scalar perturbation (normalized)')
    ax=axes[1,0];bounds=(.475,.525,.475,.525)
    zoom=[poly for poly in polys if max(x for x,y in poly)>=bounds[0] and min(x for x,y in poly)<=bounds[1]
          and max(y for x,y in poly)>=bounds[2] and min(y for x,y in poly)<=bounds[3]]
    ax.add_collection(PolyCollection(zoom,facecolors='white',edgecolors='#6f7d89',linewidths=.5))
    ax.set(xlim=bounds[:2],ylim=bounds[2:],xlabel='x [m]',ylabel='y [m]',title='Actual final CM2D mesh: local detail');ax.set_aspect('equal')
    for ax,field,title in [(axes[0,1],'scalar','Scalar error'),(axes[1,1],'velocity','Velocity vector error [m/s]')]:
        hs=[r['h'] for r in rows]
        ax.loglog(hs,[r[field+'ContinuousL2'] for r in rows],'o-',color='#d08332',label='Continuous exact solution')
        ax.loglog(hs,[r[field+'BackwardEulerL2'] for r in rows],'s-',color='#277f87',label='Space-continuous backward Euler')
        ax.loglog(hs,[r[field+'TimeDiscretizationL2'] for r in rows],':',color='#777777',label='Reference time-discretization gap')
        orders=[r[field+'BackwardEulerL2Order'] for r in rows[1:]]
        order_label=', '.join(f'{x:.2f}' if x is not None else 'undefined' for x in orders) or 'not available'
        ax.set(xlabel='h / H',ylabel='Area-weighted RMS',title=title+'\nCompleted-data error slope: '+order_label)
        ax.grid(True,which='both',alpha=.18);ax.legend(fontsize=8)
    fig.suptitle('Thermal vortex: '+('INCOMPLETE SCALE STUDY' if incomplete else 'completed three-grid study'),fontsize=15,
                 color='#9d3434' if incomplete else '#222222')
    notice=''
    if incomplete:
        completed={r['cells'] for r in rows}
        missing=[m['cellsAcross']**2 for m in study.get('meshes',[]) if m['cellsAcross']**2 not in completed]
        if missing and any(r.get('timedOut') for r in study['runs']):
            notice='Incomplete: '+', '.join(f'{n:,}' for n in missing)+f' cells exceeded the {study["timeoutSeconds"]:g} s wrapper budget.\n'
        else:notice='Study NOT passed: '+ '; '.join(study.get('issues',[]))+'\n'
    fig.supxlabel(f'Fixed dt = {info["timeStep"]:g} s, {info["steps"]} accepted steps; same final carrier as separately solved flow.\n'
                  +notice+
                  'Short-time, constant-property analytic case. No buoyancy, turbulence or long-time accuracy qualification.',fontsize=10)
    a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=155);plt.close(fig)
    print(a.output)

if __name__=='__main__':main()
