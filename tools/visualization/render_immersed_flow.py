#!/usr/bin/env python3
"""Plot actual Cartesian immersed fields; never synthesize flow samples."""
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


def render(source, destination):
    source, destination=Path(source),Path(destination)
    summary=json.loads((source/'summary.json').read_text())
    g,m=summary['grid'],summary['metrics']
    with (source/'cells.csv').open() as f:
        rows=list(csv.DictReader(f))
    nx,ny=g['nx'],g['ny']
    field=lambda name:np.array([float(row[name]) for row in rows]).reshape(ny,nx)
    u,v,classification=field('u'),field('v'),field('classification')
    x=np.linspace(0,g['length'],nx+1);y=np.linspace(0,g['height'],ny+1)
    xc=(x[1:]+x[:-1])/2;yc=(y[1:]+y[:-1])/2
    loops=[];points=[]
    for line in (source/'boundary.xy').read_text().splitlines()+['']:
        if line.strip():points.append(tuple(map(float,line.split())))
        elif points:loops.append(np.array(points));points=[]
    with (source/'history.csv').open() as f:
        history=list(csv.DictReader(f))
    fig=plt.figure(figsize=(13.5,8.5),layout='constrained')
    grid=fig.add_gridspec(2,2,height_ratios=[1,1.15])
    top=fig.add_subplot(grid[0,:]);mesh=fig.add_subplot(grid[1,0]);res=fig.add_subplot(grid[1,1])
    speed=np.hypot(u,v)
    plot=top.pcolormesh(x,y,speed,cmap='viridis',shading='flat',rasterized=True)
    flow_u=np.ma.masked_where(classification==1,u);flow_v=np.ma.masked_where(classification==1,v)
    top.streamplot(xc,yc,flow_u,flow_v,density=[2.5,.75],color='white',linewidth=.6,arrowsize=.7)
    top.set(xlabel='x',ylabel='y',aspect='equal',title=f"Actual velocity | {nx} x {ny} uncut Cartesian cells | {summary['stop_reason']}")
    fig.colorbar(plot,ax=top,label='speed')
    mesh.pcolormesh(x,y,classification,cmap=matplotlib.colors.ListedColormap(['#edf5fa','#a9b7c6','#f2b366']),vmin=0,vmax=2,shading='flat')
    segments=[[(q,0),(q,g['height'])] for q in x]+[[(0,q),(g['length'],q)] for q in y]
    mesh.add_collection(LineCollection(segments,colors='#617487',linewidths=.4))
    for loop in loops:
        ring=np.vstack([loop,loop[0]])
        top.plot(ring[:,0],ring[:,1],color='#ffcc5c',lw=1.3)
        mesh.plot(ring[:,0],ring[:,1],color='#182936',lw=1.4)
    if loops:
        points=np.vstack(loops);centre=points.mean(axis=0);radius=max(np.ptp(points[:,0]),np.ptp(points[:,1]))
        mesh.set_xlim(max(0,centre[0]-1.2*radius),min(g['length'],centre[0]+1.2*radius))
        mesh.set_ylim(max(0,centre[1]-1.2*radius),min(g['height'],centre[1]+1.2*radius))
    mesh.set(xlabel='x',ylabel='y',aspect='equal',title='Full lattice near original wall (no cut cells)')
    step=np.array([float(row['step']) for row in history])
    if len(step):
        for name,color in [('momentum','#236fa0'),('continuity','#bd5149'),('field_change','#668c36')]:
            values=np.array([float(row[name]) for row in history])
            label=name+(' (zeros at 1e-17)' if np.any(values==0) else '')
            res.semilogy(step,np.maximum(values,1e-17),label=label,color=color)
        res.axhline(summary['controls']['steady_tolerance'],color='#236fa0',ls=':',lw=.8)
        res.legend(fontsize=9)
    res.set(xlabel='accepted iteration',ylabel='normalized residual',title='Accepted iteration residuals')
    res.grid(alpha=.2)
    ref=summary['normalization']['reference_velocity']
    fig.suptitle(f"Cartesian immersed-boundary prototype: {summary['case']}\n"
                 f"Periodic x; {summary['controls'].get('wall_method','brinkman')} | sampled wall speed / Uref = {100*m['wall_speed_max']/ref:.3f}%",fontsize=14)
    fig.supxlabel(f"Solid auxiliary values retained. Grid {summary['timing']['grid_seconds']:.4f}s + boundary {summary['timing']['boundary_seconds']:.3f}s + solve {summary['timing']['solve_seconds']:.2f}s.\n"
                  'Development result: finite wall penalty and spatial error remain; physical accuracy is not qualified.',fontsize=10)
    destination.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(destination,dpi=160)
    plt.close(fig)
    evidence={'source':str(source.resolve()),'image':str(destination.resolve()),'inputs':{}}
    for name in ['cells.csv','u.csv','v.csv','walls.csv','history.csv','boundary.xy','summary.json']:
        evidence['inputs'][name]=hashlib.sha256((source/name).read_bytes()).hexdigest()
    if (source/'wall-markers.csv').exists():
        evidence['inputs']['wall-markers.csv']=hashlib.sha256((source/'wall-markers.csv').read_bytes()).hexdigest()
    destination.with_suffix(destination.suffix+'.json').write_text(json.dumps(evidence,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    render(args.directory,args.output)
