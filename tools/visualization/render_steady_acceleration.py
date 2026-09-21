#!/usr/bin/env python3
"""Render actual mesh/fields, all residual rows and measured solve costs."""
import argparse,csv,json,math,sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_native_flow as native

def main(args):
    study=json.loads((args.study/'report.json').read_text())
    if study.get('complete') is not True:raise ValueError('Study is incomplete')
    for name,digest in study['artifactSha256'].items():
        path=args.study/name
        if not path.resolve().is_relative_to(args.study.resolve()) or native.sha256_file(path)!=digest:raise ValueError('Changed artifact: '+name)
    case=study['cases'][0];command=case['runs'][0]['execution']['command'];meshpath=Path(command[command.index('--mesh')+1])
    if native.sha256_file(meshpath)!=case['meshSha256']:raise ValueError('Changed mesh')
    mesh=native.read_cm2d(meshpath);measured=native.measure(mesh,1e-11,1e-9)
    values=[];histories=[]
    for mode in ['none','anderson']:
        prefix=args.study/('half-'+mode)
        audit=native.verify_case(meshpath,prefix,'custom',.1,1,native.argument_parser().parse_args(['--max-iterations','20000']))
        if not audit['valid']:raise ValueError(audit['issues'])
        values.append(native.read_cells(Path(str(prefix)+'.cells.csv'),mesh,measured,'custom'))
        histories.append(list(csv.DictReader(Path(str(prefix)+'.residuals.csv').open())))
    polygons=[[mesh.vertices[i] for i in cell.vertices] for cell in mesh.cells]
    speed=[math.hypot(r['u'],r['v']) for r in values[1]]
    delta=[math.hypot(a['u']-b['u'],a['v']-b['v']) for a,b in zip(*values)]
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    for ax,data,title in [(axes[0,0],speed,'Accelerated half-channel: actual Cut-cell field'),(axes[0,1],delta,'Velocity difference from ordinary SIMPLE')]:
        collection=PolyCollection(polygons,array=data,cmap='viridis',edgecolors='#536467',linewidths=.14)
        ax.add_collection(collection);ax.set(xlim=(0,2),ylim=(0,1),aspect='equal',xlabel='x [m]',ylabel='y [m]',title=title)
        fig.colorbar(collection,ax=ax,label='m/s',shrink=.8)
    for history,label in zip(histories,['SIMPLE','Safeguarded Anderson']):
        axes[1,0].semilogy([int(r['iteration']) for r in history],[float(r['momentumResidual']) for r in history],label=label,lw=1)
    axes[1,0].axhline(1e-9,color='black',linestyle=':',label='Unchanged nonlinear target')
    axes[1,0].set(xlabel='Iteration',ylabel='Original momentum residual',title='Every recorded iterate; no smoothing');axes[1,0].legend(fontsize=8);axes[1,0].grid(alpha=.2)
    cases=study['cases'][:3];x=list(range(len(cases)))
    for shift,label,mode in [(-.18,'SIMPLE','none'),(.18,'Safeguarded Anderson','anderson')]:
        seconds=[next(r['performance']['solveSeconds'] for r in c['runs'] if r['mode']==mode) for c in cases]
        bars=axes[1,1].bar([i+shift for i in x],seconds,.36,label=label);axes[1,1].bar_label(bars,fmt='%.2f s',fontsize=8)
    axes[1,1].set(xticks=x,xticklabels=['Half-channel','Nozzle','Cylinder Re20'],ylabel='Measured solve wall time [s]',title='One serial run per mode on this Mac');axes[1,1].legend(fontsize=8)
    fig.suptitle('Steady acceleration with unchanged equations and acceptance gates',fontsize=15)
    fig.supxlabel('Timing is case/hardware-specific. Final ordinary SIMPLE step certifies convergence. Field differences include iterative stopping error.',fontsize=9)
    args.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(args.output,dpi=145);plt.close(fig)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--study',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);main(parser.parse_args())
