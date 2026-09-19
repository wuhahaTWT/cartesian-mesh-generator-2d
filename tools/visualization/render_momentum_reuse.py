#!/usr/bin/env python3
"""Plot actual flow cells and paired timings from the momentum-reuse study."""
import argparse
import csv
import json
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'tools/verification'))
import verify_native_flow as native


def path(value):
    p=Path(value)
    return p if p.is_absolute() else REPO/p


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    evidence=json.loads(args.evidence.read_text())
    comparisons=evidence['comparisons']
    if not comparisons or not all(c['valid'] and not c['mismatches'] for c in comparisons.values()):
        raise ValueError('Requires independently audited byte-identical paired cases')
    run=comparisons['circle']['runs']['new']
    command=run['command'];prefix=path(command[command.index('--output')+1])
    mesh=native.read_cm2d(path(run['mesh']))
    with Path(str(prefix)+'.cells.csv').open() as f:cells=list(csv.DictReader(f))
    if len(cells)!=len(mesh.cells):raise ValueError('Field/mesh size mismatch')
    polygons=[];speed=[]
    for cell,field in zip(mesh.cells,cells):
        if int(field['cell'])!=cell.id:raise ValueError('Cell ordering mismatch')
        poly=[mesh.vertices[v] for v in cell.vertices]
        if max(x for x,y in poly)<-3 or min(x for x,y in poly)>3 or max(y for x,y in poly)<-3 or min(y for x,y in poly)>3:continue
        polygons.append(poly);speed.append((float(field['u'])**2+float(field['v'])**2)**.5)
    fig,(ax,chart)=plt.subplots(1,2,figsize=(12,5.4),layout='constrained')
    artist=PolyCollection(polygons,array=speed,cmap='viridis',edgecolors='#617283',linewidths=.25,rasterized=True)
    ax.add_collection(artist);ax.set(xlim=(-3,3),ylim=(-3,3),xlabel='x [m]',ylabel='y [m]',
        title=f'Actual Cut-cell flow, {len(cells):,} cells\nLocal mesh and velocity at t = 0.5 s')
    ax.set_aspect('equal');fig.colorbar(artist,ax=ax,shrink=.8,label='Speed [m/s]')
    labels=[];old=[];new=[]
    for label,c in comparisons.items():
        r=c['runs']['new'];cmd=r['command'];p=path(cmd[cmd.index('--output')+1])
        summary=json.loads(Path(str(p)+'.json').read_text())
        labels.append(f'{label}\n{summary["cells"]:,} cells')
        old.append(c['runs']['old']['elapsedSeconds']);new.append(r['elapsedSeconds'])
    xs=list(range(len(labels)))
    a=chart.bar([x-.19 for x in xs],old,width=.38,label='Before',color='#a9b0b7')
    b=chart.bar([x+.19 for x in xs],new,width=.38,label='Reuse current-state assembly',color='#287f87')
    chart.bar_label(a,fmt='%.2f',fontsize=8,padding=3);chart.bar_label(b,fmt='%.2f',fontsize=8,padding=3)
    chart.set(xticks=xs,xticklabels=labels,ylabel='Wall time [s], including I/O',
              title='Serial paired runs, same controls and solver stops')
    chart.set_ylim(0,max(old+new)*1.2);chart.grid(axis='y',alpha=.15);chart.legend(fontsize=8)
    fig.suptitle('Reuse numerical work without changing the solved fields',fontsize=14)
    fig.supxlabel('Compared fields, face fluxes and residual histories are byte-identical; independent equation audits pass.\n'
                  'One paired run per case. These cases do not certify general engineering CFD accuracy.',fontsize=10)
    args.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(args.output,dpi=155);plt.close(fig)


if __name__=='__main__':main()
