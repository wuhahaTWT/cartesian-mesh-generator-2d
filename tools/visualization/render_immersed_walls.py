#!/usr/bin/env python3
"""Compare actual dense original-wall samples and final Cartesian flow fields."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/cartmesh2d-immersed-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


def table(path):
    with path.open() as source:
        return list(csv.DictReader(source))


def render(pairs, destination):
    fig,axes=plt.subplots(2,len(pairs),figsize=(5*len(pairs),8.2),squeeze=False,layout='constrained')
    provenance={'pairs':[],'description':'Computed fields and dense boundary samples; no synthesized flow'}
    for column,(label,base,current) in enumerate(pairs):
        base,current=Path(base),Path(current)
        summary=json.loads((current/'summary.json').read_text());g=summary['grid'];m=summary['metrics']
        ref=summary['normalization']['reference_velocity'];nx,ny=g['nx'],g['ny']
        previous=json.loads((base/'summary.json').read_text())
        if previous['normalization']['reference_velocity']!=ref or any(previous['grid'][key]!=g[key] for key in ('nx','ny','length','height')) or any(previous['controls'][key]!=summary['controls'][key] for key in ('viscosity','drive')) or (base/'boundary.xy').read_bytes()!=(current/'boundary.xy').read_bytes():
            raise ValueError('Each pair must use the same reference speed, grid and original boundary')
        rows=table(current/'cells.csv')
        field=lambda key:np.array([float(row[key]) for row in rows]).reshape(ny,nx)
        speed=np.hypot(field('u'),field('v'))/ref
        speed=np.ma.masked_where(field('classification')==1,speed)
        x=np.linspace(0,g['length'],nx+1);y=np.linspace(0,g['height'],ny+1)
        top,bottom=axes[:,column]
        color=top.pcolormesh(x,y,speed,cmap='viridis',vmin=0,vmax=1,shading='flat')
        top.set_facecolor('#c6cdd2')
        top.add_collection(LineCollection([[(q,0),(q,g['height'])] for q in x]+[[(0,q),(g['length'],q)] for q in y],colors='#243847',linewidths=.3,alpha=.45))
        loops=[];points=[]
        for line in (current/'boundary.xy').read_text().splitlines()+['']:
            if line.strip():points.append(tuple(map(float,line.split())))
            elif points:loops.append(np.array(points));points=[]
        if len(loops)!=1:
            raise ValueError('This perimeter-comparison view requires one loop per case')
        all_points=np.vstack(loops);minimum=all_points.min(axis=0);maximum=all_points.max(axis=0)
        for loop in loops:
            ring=np.vstack([loop,loop[0]]);top.plot(ring[:,0],ring[:,1],color='#f7bc39',lw=1.8)
        margin=3*max(g['dx'],g['dy'])
        top.set(xlim=(minimum[0]-margin,maximum[0]+margin),ylim=(minimum[1]-margin,maximum[1]+margin),aspect='equal',xlabel='x',ylabel='y',title=f'{label} | {nx} x {ny}: computed speed / Uref')
        fig.colorbar(color,ax=top,shrink=.8)
        hashes={}
        for source,name,color in ((base,'Brinkman','#d68029'),(current,'Coupled wall','#236ea5')):
            walls=table(source/'walls.csv')
            point=np.array([[float(row[k]) for k in ('x','y')] for row in walls])
            arc=np.concatenate(([0.],np.cumsum(np.linalg.norm(np.diff(point,axis=0),axis=1))))
            values=100*np.hypot([float(row['u']) for row in walls],[float(row['v']) for row in walls])/ref
            bottom.semilogy(arc/arc[-1],np.ma.masked_less_equal(values,0),color=color,lw=1.5,label=f'{name}: max {max(values):.4g}%')
            hashes[name]={'directory':str(source.resolve()),'sha256':{p:hashlib.sha256((source/p).read_bytes()).hexdigest() for p in ('cells.csv','walls.csv','boundary.xy','summary.json')}}
        bottom.set(xlabel='fraction of original wall perimeter',ylabel='sampled wall speed / Uref (%)',title=f"Normal max: {100*m['wall_normal_speed_max']/ref:.4g}% | wall time {summary['controls']['wall_penalty_time']:g}")
        bottom.legend(loc='best',fontsize=9);bottom.grid(alpha=.22,which='both')
        provenance['pairs'].append({'label':label,'sources':hashes})
    fig.suptitle('Pure Cartesian CFD: pressure and wall force solved together\nFull square cells retained; curves use independent dense wall samples, including corners',fontsize=15)
    fig.supxlabel('Gray cells hold solid auxiliary values. Each pair has the same grid, reference speed and original geometry.\nSmall wall velocity and algebraic residuals do not establish field accuracy or grid independence.',fontsize=10)
    destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(destination,dpi=150);plt.close(fig)
    destination.with_suffix(destination.suffix+'.json').write_text(json.dumps(provenance,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pair',nargs=3,action='append',required=True,metavar=('LABEL','BASELINE','CURRENT'))
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();render(args.pair,args.output)
