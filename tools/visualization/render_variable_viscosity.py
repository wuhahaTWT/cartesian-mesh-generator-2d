#!/usr/bin/env python3
"""Render actual prescribed-viscosity MMS fields and audited refinement errors."""
import argparse
import csv
import json
import math
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_native_flow as native


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();study=json.loads(args.study.read_text())
    if not study['valid']:raise ValueError('Audited study required')
    series=study['cases']['symmetric'];fine=series[-1];prefix=fine['prefix']
    mesh=native.read_cm2d(Path(fine['mesh']))
    with Path(prefix+'.cells.csv').open() as f:cells=list(csv.DictReader(f))
    with Path(prefix+'.faces.csv').open() as f:faces=list(csv.DictReader(f))
    if len(cells)!=len(mesh.cells) or len(faces)!=len(mesh.edges):raise ValueError('Mesh/field mismatch')
    polygons=[[mesh.vertices[i] for i in c.vertices] for c in mesh.cells]
    speed=[math.hypot(float(c['u']),float(c['v'])) for c in cells]
    error=[math.hypot(float(c['u'])-float(c['exactU']),float(c['v'])-float(c['exactV'])) for c in cells]
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    for ax,values,cmap,title in [(axes[0,0],speed,'viridis',f'Computed velocity magnitude, {len(cells):,} cells'),(axes[0,1],error,'magma','Actual vector velocity error')]:
        artist=PolyCollection(polygons,array=values,cmap=cmap,edgecolors='none',rasterized=True)
        ax.add_collection(artist);ax.autoscale();ax.set_aspect('equal');ax.set(title=title,xlabel='x [m]',ylabel='y [m]')
        fig.colorbar(artist,ax=ax,shrink=.75,label='m/s')
    segments=[[mesh.vertices[e.v0],mesh.vertices[e.v1]] for e in mesh.edges]
    artist=LineCollection(segments,array=[float(f['viscosity']) for f in faces],cmap='cividis',linewidths=.6)
    ax=axes[1,0];ax.add_collection(artist);ax.set(xlim=(0,.3),ylim=(0,.3),xlabel='x [m]',ylabel='y [m]',title='Final mesh detail; same viscosity in stress and matrix');ax.set_aspect('equal')
    fig.colorbar(artist,ax=ax,shrink=.75,label='Kinematic viscosity [m2/s]')
    ax=axes[1,1];h=[r['meshMeasurement']['characteristicH'] for r in series]
    for name,key,color in [('Velocity / Uref','velocityL2Relative','#327a91'),('Pressure / Uref^2','pressureL2Relative','#b46a3e'),('Wall viscous traction / Uref^2',None,'#6c618e')]:
        values=[r['benchmark'][key] if key else r['wallTraction']['l2']/r['speed']**2 for r in series]
        ax.loglog(h,values,'o-',color=color,label=name)
    ax.set(xlabel='Effective spacing [m]',ylabel='Normalized RMS error',title='Full symmetric stress: refinement of fields and wall force');ax.grid(True,which='both',alpha=.18);ax.legend(fontsize=9)
    fig.suptitle('Spatial viscosity in conservative incompressible flow',fontsize=16)
    fig.supxlabel('Prescribed nu(x)=0.1(1+x), Uref=1; unit-square forced vortex, same 1e-9 solve tolerance.\n'
                  'Independent geometry, source, face stress, wall force and momentum balance checks passed.\n'
                  'A variable-coefficient laminar foundation; not SST or engineering flow qualification.',fontsize=10)
    args.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(args.output,dpi=150);plt.close(fig)


if __name__=='__main__':main()
