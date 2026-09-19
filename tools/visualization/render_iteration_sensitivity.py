#!/usr/bin/env python3
"""Plot the measured flow-tolerance study, including its qualification limit."""
import argparse
import csv
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_native_flow as native


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();evidence=json.loads(args.evidence.read_text())
    rows=evidence['finalTimeComparison']['cases']
    if len(rows)!=3 or not all(c['valid'] for c in evidence['independentAudits']):
        raise ValueError('Three-grid comparison and independent audits required')
    final=rows[-1];prefix=Path(final['tight']['prefix'])
    summary=json.loads(Path(str(prefix)+'.json').read_text())
    if summary['time']!=final['time']:raise ValueError('Field time mismatch')
    mesh=native.read_cm2d(Path(final['mesh']))
    with Path(str(prefix)+'.cells.csv').open() as f:fields=list(csv.DictReader(f))
    if len(fields)!=len(mesh.cells):raise ValueError('Field/mesh count mismatch')
    amplitude=final['expectedBEAmplitude'];polys=[];errors=[]
    for cell,f in zip(mesh.cells,fields):
        if int(f['cell'])!=cell.id:raise ValueError('Field ordering mismatch')
        x,y=float(f['x']),float(f['y'])
        u=amplitude*math.sin(math.pi*x)*math.cos(math.pi*y)
        v=-amplitude*math.cos(math.pi*x)*math.sin(math.pi*y)
        errors.append(math.hypot(float(f['u'])-u,float(f['v'])-v))
        polys.append([mesh.vertices[k] for k in cell.vertices])
    fig,(ax,chart)=plt.subplots(1,2,figsize=(12,5.8),layout='constrained')
    artist=PolyCollection(polys,array=errors,cmap='magma',edgecolors='none',rasterized=True)
    ax.add_collection(artist);ax.set(xlim=(0,1),ylim=(0,1),xlabel='x [m]',ylabel='y [m]',
                                   title='102,400 cells: actual velocity error\nTighter iteration tolerance, t = 0.01 s')
    ax.set_aspect('equal');fig.colorbar(artist,ax=ax,shrink=.75,label='Vector error vs backward-Euler reference [m/s]')
    xs=[r['cells'] for r in rows]
    for key,color,label in [('baseline','#bd7035','Iteration tolerance 1e-8'),('tight','#257b83','Iteration tolerance 1e-10')]:
        chart.loglog(xs,[r[key]['velocityBackwardEulerL2'] for r in rows],'o-',color=color,label=label)
    chart.set(xticks=xs,xticklabels=['10,000','40,000','102,400'],xlabel='Actual cell count',ylabel='Area-weighted velocity RMS error [m/s]',
              title='Same physical time, grid, time step and scheme')
    chart.minorticks_off();chart.grid(True,alpha=.18);chart.legend(fontsize=9)
    fig.suptitle('Iteration tolerance can obscure a grid-refinement comparison',fontsize=14)
    fig.supxlabel('Fixed dt = 0.005 s, t = 0.01 s. Reference is continuous in space and backward Euler in time.\n'
                  'Flow-only diagnosis; the original synchronized thermal qualification remains NOT PASSED.\n'
                  'The tight 102,400-cell run uses two checkpoint-linked steps; total runtime is reported, not a 300 s pass.',fontsize=10)
    args.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(args.output,dpi=155);plt.close(fig)


if __name__=='__main__':main()
