#!/usr/bin/env python3
"""Render real thermal fields and fixed-grid continuous-time reference errors."""
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
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_native_flow as native


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();study=json.loads(a.study.read_text())
    if not study['valid'] or not all(c['valid'] for c in study['cases']):raise ValueError('Audited time study required')
    rows=study['temporal']['series'];case=study['cases'][-1]
    info=json.loads(Path(case['prefix']+'.json').read_text());mesh=native.read_cm2d(Path(info['mesh']))
    with Path(case['prefix']+'.cells.csv').open() as f:scalar=list(csv.DictReader(f))
    with Path(case['flowPrefix']+'.cells.csv').open() as f:flow=list(csv.DictReader(f))
    if len(scalar)!=len(flow) or len(flow)!=len(mesh.cells):raise ValueError('Actual field/mesh mismatch')
    polygons=[];values=[];errors=[]
    amp=rows[-1]['speed']*math.exp(-2*math.pi**2*rows[-1]['nu']*rows[-1]['time'])
    for cell,s,f in zip(mesh.cells,scalar,flow):
        if cell.id!=int(s['cell']) or cell.id!=int(f['cell']):raise ValueError('Cell ordering mismatch')
        x,y=float(f['x']),float(f['y']);u=amp*math.sin(math.pi*x)*math.cos(math.pi*y);v=-amp*math.cos(math.pi*x)*math.sin(math.pi*y)
        polygons.append([mesh.vertices[k] for k in cell.vertices]);values.append(float(s['value']))
        errors.append(math.hypot(float(f['u'])-u,float(f['v'])-v))
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    for ax,data,cmap,title,label in [(axes[0,0],values,'inferno',f'Computed scalar, {len(mesh.cells):,} cells','Normalized scalar'),(axes[1,0],errors,'viridis','Actual velocity error vs continuous solution','Vector error [m/s]')]:
        artist=PolyCollection(polygons,array=data,cmap=cmap,edgecolors='none',rasterized=True)
        ax.add_collection(artist);ax.set(xlim=(0,1),ylim=(0,1),xlabel='x [m]',ylabel='y [m]',title=title);ax.set_aspect('equal');fig.colorbar(artist,ax=ax,shrink=.7,label=label)
    for ax,key,label in [(axes[0,1],'scalar','Scalar'),(axes[1,1],'velocity','Velocity vector')]:
        dt=[r['dt'] for r in rows]
        for suffix,text,color in [('ContinuousL2','Continuous exact solution','#b97833'),('BackwardEulerL2','Space-continuous backward Euler','#277f87')]:
            ax.loglog(dt,[r[key+suffix] for r in rows],'o-',color=color,label=text)
        orders=[r[key+'ContinuousL2Order'] for r in rows[1:]]
        ax.set(xlabel='Time step [s]',ylabel='Area-weighted RMS',title=label+' error; observed slopes '+', '.join(f'{o:.2f}' if o is not None else 'undefined' for o in orders));ax.grid(True,which='both',alpha=.18);ax.legend(fontsize=9)
    fig.suptitle('Fixed mesh and final time: time-step refinement',fontsize=15)
    fig.supxlabel(f'Same {len(mesh.cells):,}-cell mesh, t = {rows[-1]["time"]:g} s, flow tolerance {rows[-1]["flowTolerance"]:g}; dt = '+ ' / '.join(f'{r["dt"]:g}' for r in rows)+' s.\n'+
                  f'Fields shown at dt = {rows[-1]["dt"]:g} s. All coupled carriers match separately solved flow.\n'
                  'Short-time, constant-property analytic case; slopes include spatial and iterative error.',fontsize=10)
    a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=155);plt.close(fig)


if __name__=='__main__':main()
