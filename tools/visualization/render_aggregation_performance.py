#!/usr/bin/env python3
"""Render real fields alongside serial baseline/candidate benchmark measurements."""
import argparse
import json
from pathlib import Path
import statistics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from render_native_flow import read_cm2d, read_cells


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--comparisons',nargs='+',type=Path,required=True)
    p.add_argument('--audit',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    audit=json.loads(args.audit.read_text())
    if audit.get('valid') is not True: raise ValueError('Independent field audit required')
    pairs=[json.loads(path.read_text()) for path in args.comparisons]
    if not all(d['pairs'] and all(all(p['byteIdentical'].values()) for p in d['pairs']) for d in pairs):
        raise ValueError('Requires byte-identical verified comparison outputs')
    if not any(r['prefix']==audit['prefix'] for d in pairs for r in d['runs'] if r['label']=='candidate'):
        raise ValueError('Field audit must belong to a plotted candidate')
    mesh=read_cm2d(Path(audit['mesh']))
    cells=read_cells(Path(audit['prefix']+'.cells.csv'),len(mesh.cells))
    polygons=[[mesh.vertices[v] for v in c.vertices] for c in mesh.cells]
    fig,axes=plt.subplots(2,2,figsize=(12,7.3),layout='constrained')
    ax=axes[0,0]
    pc=PolyCollection(polygons,array=[c['u'] for c in cells],cmap='viridis',edgecolors='none',rasterized=True)
    ax.add_collection(pc);ax.autoscale_view();ax.set_aspect('equal')
    ax.set(title=f'Actual converged channel: {len(cells):,} cells',xlabel='x [m]',ylabel='y [m]')
    fig.colorbar(pc,ax=ax,shrink=.65,label='Streamwise velocity [m/s]')
    ax=axes[1,0]
    xmin,ymin,xmax,ymax=audit['meshMeasurement']['bounds']
    cx=(xmin+xmax)*.5; width=(ymax-ymin)*.025
    # A narrow wall-adjacent view exposes the actual fine Cartesian cells.
    zoom=(cx-width,cx+width,ymin,ymin+2*width)
    visible=[poly for poly in polygons if max(x for x,y in poly)>=zoom[0] and min(x for x,y in poly)<=zoom[1]
             and max(y for x,y in poly)>=zoom[2] and min(y for x,y in poly)<=zoom[3]]
    ax.add_collection(PolyCollection(visible,facecolors='white',edgecolors='#73818e',linewidths=.45))
    ax.set(xlim=zoom[:2],ylim=zoom[2:],xlabel='x [m]',ylabel='y [m]',title='Actual final mesh near the wall')
    ax.set_aspect('equal')
    for ax,key,title,factor in [(axes[0,1],'elapsedSeconds','End-to-end elapsed time [s]',1),
                                (axes[1,1],'peakRssBytes','Peak resident memory [MiB]',1/1048576)]:
        labels=[]
        for i,d in enumerate(pairs):
            labels.append(f"{d['runs'][0]['summary']['cells']:,}\n"+d['runs'][0]['summary']['case'])
            for label,offset,color in [('baseline',-.18,'#8c9ba8'),('candidate',.18,'#318781')]:
                values=[r[key]*factor for r in d['runs'] if r['label']==label]
                mean=statistics.mean(values)
                ax.bar(i+offset,mean,width=.32,color=color,label=label if i==0 else None)
                ax.scatter([i+offset]*len(values),values,color='#273543',s=10,zorder=3)
                ax.text(i+offset,mean,f'{mean:.2f}',ha='center',va='bottom',fontsize=9)
        ax.set(xticks=range(len(labels)),xticklabels=labels,title=title)
        ax.margins(y=.16);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True);ax.legend(fontsize=9)
    fig.suptitle('Pressure multigrid: same computed fields, less temporary assembly storage',fontsize=14)
    fig.supxlabel('Same meshes, equations, tolerances and iterations. Dots are individual serial runs.\n'
                  'Small timing differences do not establish a general speedup; 500k is one pair.',fontsize=10)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=155,facecolor='white');plt.close(fig)
    print(args.output)

if __name__=='__main__':main()
