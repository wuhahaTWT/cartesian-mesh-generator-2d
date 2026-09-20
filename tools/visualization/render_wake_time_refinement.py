#!/usr/bin/env python3
"""Plot retained fields and temporal sensitivity without implying qualification."""
import argparse,csv,json,math,sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_transient_flow as transient
n=transient.native


def main(args):
    report=json.loads((args.study/'study.json').read_text())
    if report.get('valid') is not True or report.get('accuracyQualification')!='not-qualified':
        raise ValueError('Expected a completed observation study without a physical qualification claim')
    meshpath=Path(report['mesh']);mesh=n.read_cm2d(meshpath);m=n.measure(mesh,1e-11,1e-9)
    if n.sha256_file(meshpath)!=report['meshSha256']:raise ValueError('Mesh digest changed')
    fields=[];histories=[]
    for run in report['runs']:
        folder=args.study/f"n{run['steps']}"
        for name,digest in run['artifactSha256'].items():
            if Path(name).name!=name or n.sha256_file(folder/name)!=digest:raise ValueError('Retained run artifact differs: '+name)
        a=transient.verify(meshpath,folder/'flow',folder/'figure-audit.json')
        if a['controls']!=report['controls']:raise ValueError('Physical controls differ')
        fields.append(n.read_cells(folder/'flow.cells.csv',mesh,m,'external'))
        histories.append(list(csv.DictReader((folder/'flow.time-history.csv').open())))
    polygons=[[mesh.vertices[i] for i in c.vertices] for c in mesh.cells]
    fine=fields[-1];coarse=fields[0];speed=[math.hypot(c['u'],c['v']) for c in fine]
    delta=[math.hypot(a['u']-b['u'],a['v']-b['v']) for a,b in zip(coarse,fine)]
    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    for ax,values,title in [(axes[0,0],speed,'Final speed: smallest time step'),(axes[0,1],delta,'Velocity difference: largest vs smallest step')]:
        collection=PolyCollection(polygons,array=values,cmap='viridis',edgecolors='#667477',linewidths=.12)
        ax.add_collection(collection);ax.set(xlim=(-2,12),ylim=(-4,4),aspect='equal',title=title,xlabel='x [m]',ylabel='y [m]')
        fig.colorbar(collection,ax=ax,label='m/s',shrink=.85)
    for run,history in zip(report['runs'],histories):
        axes[1,0].plot([float(r['time']) for r in history],[float(r['forceY']) for r in history],label=f"dt={run['dt']:g} s")
    axes[1,0].set(xlabel='Physical time [s]',ylabel='Fy / density / depth [m³/s²]',title='Recorded transverse force at every accepted step')
    axes[1,0].legend();axes[1,0].grid(alpha=.25)
    steps=[r['dt'] for r in report['runs'][:-1]];u=report['controls']['speed']
    for key,scale,label in [('velocityL2',u,'Velocity L2 / Uref'),('pressureL2',u*u,'Pressure L2 / Uref²')]:
        order=report['observedDifferenceOrders'][key]
        axes[1,1].loglog(steps,[d[key]/scale for d in report['successiveDifferences']],'o-',label=f'{label}; difference order={order:.3f}' if order is not None else label)
    axes[1,1].set(xlabel='Larger dt of each successive pair [s]',ylabel='Normalized field difference',title='Two successive differences: observed sensitivity')
    axes[1,1].legend(fontsize=8);axes[1,1].grid(alpha=.25,which='both')
    fig.suptitle(f"Same audited wake state at t={report['startTime']:g} s; {report['duration']:g} s window on {len(mesh.cells):,} cells",fontsize=15)
    fig.supxlabel('Earlier accumulated temporal error and mesh error remain. These comparisons do not qualify long-time shedding frequency or spatial accuracy.',fontsize=10)
    args.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(str(args.output)+'.png',dpi=145);plt.close(fig)
    evidence={**report,'studySha256':n.sha256_file(args.study/'study.json'),'plotSourceSha256':n.sha256_file(Path(__file__))}
    n.write_json(Path(str(args.output)+'.json'),evidence)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
