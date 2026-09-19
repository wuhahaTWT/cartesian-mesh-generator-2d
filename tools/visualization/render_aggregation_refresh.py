#!/usr/bin/env python3
"""Render actual mesh fields and the measured coefficient-refresh comparison."""
import argparse
import csv
import json
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_native_flow as native


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    evidence=json.loads(args.evidence.read_text())
    cases=evidence['cases']
    if not all(c['newAudit']['valid'] for c in cases):
        raise ValueError('Independent audit required')
    case=next(c for c in cases if c['name']=='n320')
    mesh=native.read_cm2d(Path(case['run']['command'][case['run']['command'].index('--mesh')+1]))
    with Path(case['newPrefix']+'.cells.csv').open() as f:rows=list(csv.DictReader(f))
    if len(rows)!=len(mesh.cells):raise ValueError('Cell count mismatch')
    polygons=[];speeds=[]
    for cell,row in zip(mesh.cells,rows):
        if cell.id!=int(row['cell']):raise ValueError('Cell ordering mismatch')
        polygons.append([mesh.vertices[i] for i in cell.vertices])
        speeds.append(np.hypot(float(row['u']),float(row['v'])))
    fig,(ax,chart)=plt.subplots(1,2,figsize=(12,5.8),layout='constrained')
    artist=PolyCollection(polygons,array=np.asarray(speeds),cmap='viridis',edgecolors='none',rasterized=True)
    ax.add_collection(artist);ax.set(xlim=(0,1),ylim=(0,1),xlabel='x [m]',ylabel='y [m]',title='102,400 actual cells\nVelocity magnitude, t = 0.005 s');ax.set_aspect('equal')
    fig.colorbar(artist,ax=ax,shrink=.75,label='Speed [m/s]')
    x=np.arange(len(cases));width=.36
    for shift,key,label,color in [(-width/2,'oldPerformance','Full hierarchy setup','#a6aeb8'),(width/2,'newPerformance','Current coefficients, reused grouping','#367f91')]:
        values=[c[key]['pressureLinearSolveSeconds'] for c in cases]
        bars=chart.bar(x+shift,values,width,label=label,color=color)
        chart.bar_label(bars,fmt='%.2f',fontsize=9)
    chart.set(xticks=x,xticklabels=['Circle\n4,880','Vortex\n10,000','Vortex\n40,000','Vortex\n102,400'],ylabel='Pressure linear solve + setup [s]',title='Same equations and stopping tolerances')
    chart.legend(fontsize=8);chart.set_ylim(top=max(c[k]['pressureLinearSolveSeconds'] for c in cases for k in ['oldPerformance','newPerformance'])*1.3)
    fig.suptitle('Pressure multigrid: measure the total benefit of reusing grouping',fontsize=14)
    fig.supxlabel('Circle: fresh paired runs. Vortex baseline: earlier same-machine runs; one timing each.\n10k / 40k: two steps. 102,400: one step. No general speed or thermal-qualification claim.',fontsize=10)
    args.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(args.output,dpi=155);plt.close(fig)


if __name__=='__main__':main()
