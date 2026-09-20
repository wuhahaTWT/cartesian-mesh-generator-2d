#!/usr/bin/env python3
"""Render the actual final mesh/field and all retained physical-time monitors."""
import argparse,csv,json,math,sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_transient_flow as transient
n=transient.native


def main(args):
    checked=json.loads(args.verified.read_text());source=Path(checked['source']);history=Path(checked['historyFile'])
    if checked.get('valid') is not True or n.sha256_file(source)!=checked['sourceSha256'] or n.sha256_file(history)!=checked['historySha256']:
        raise ValueError('Source trajectory/history differs from the verified record')
    original=json.loads(source.read_text())
    # A budget-stopped study may contain a later, unaudited attempt. Render the
    # final verified segment, never that candidate field.
    label=checked['segments'][-1]['label']
    candidates=[r for r in original['runs'] if r.get('valid') is True and r.get('time')==checked['endTime'] and
        Path(r['command'][r['command'].index('--output')+1]).parent.name==label]
    if len(candidates)!=1:raise ValueError('Cannot identify the final verified segment')
    run=candidates[0];command=run['command']
    prefix=Path(command[command.index('--output')+1]);mesh_path=Path(original['mesh'])
    if n.sha256_file(mesh_path)!=checked['meshSha256']:raise ValueError('Final mesh changed')
    final=transient.verify(mesh_path,prefix,Path(str(prefix)+'.plot-audit.json'))
    if final['time']!=checked['endTime']:raise ValueError('Final field time differs')
    mesh=n.read_cm2d(mesh_path);m=n.measure(mesh,1e-11,1e-9);g=n.face_geometry(mesh,m)
    cells=n.read_cells(Path(str(prefix)+'.cells.csv'),mesh,m,'external')
    summary=json.loads(Path(str(prefix)+'.json').read_text())
    if summary['case']!='external' or summary['nu']!=.02 or summary['speed']!=1:
        raise ValueError('This cylinder observation figure requires the recorded Re100 controls')
    boundary=n.flow_boundaries(mesh,m,'external',1)
    u=[c['u'] for c in cells];v=[c['v'] for c in cells]
    gu=n.reconstruct_gradient(mesh,m,g,u,boundary['u'],boundary['fixedU'])
    gv=n.reconstruct_gradient(mesh,m,g,v,boundary['v'],boundary['fixedV'])
    omega=[b[0]-a[1] for a,b in zip(gu,gv)]
    polygons=[[mesh.vertices[i] for i in c.vertices] for c in mesh.cells]
    rows=list(csv.DictReader(history.open()));times=[float(r['time']) for r in rows]
    fx=[float(r['forceX']) for r in rows];fy=[float(r['forceY']) for r in rows]
    if len(rows)!=checked['acceptedSteps'] or times[-1]!=checked['endTime']:raise ValueError('Monitor count/time differs')
    fig,axes=plt.subplots(2,2,figsize=(15,9),layout='constrained')
    for ax,values,title,cmap in [(axes[0,0],[math.hypot(a,b) for a,b in zip(u,v)],'Actual final velocity magnitude','viridis'),
        (axes[0,1],omega,'Vorticity reconstructed from final velocity','RdBu_r')]:
        collection=PolyCollection(polygons,array=values,cmap=cmap,edgecolors='#647276',linewidths=.12)
        if cmap=='RdBu_r':collection.set_clim(-max(map(abs,omega)),max(map(abs,omega)))
        ax.add_collection(collection);ax.set(xlim=(-2,12),ylim=(-4,4),aspect='equal',title=title,xlabel='x [m]',ylabel='y [m]')
        fig.colorbar(collection,ax=ax,label='Speed [m/s]' if cmap=='viridis' else 'Vorticity [1/s]',shrink=.8)
    axes[1,0].plot(times,fx);axes[1,0].set(title='Streamwise force'+(', including startup impulse' if checked['startTime']==0 else ''),xlabel='Physical time [s]',ylabel='Fx / density / depth [m³/s²]',yscale='log')
    axes[1,1].plot(times,fy);axes[1,1].set(title='Transverse force history',xlabel='Physical time [s]',ylabel='Fy / density / depth [m³/s²]')
    for ax in axes[1]:ax.grid(alpha=.25)
    title=f"Native Re100 cylinder: {len(mesh.cells):,} Cut-cell cells, {checked['acceptedSteps']:,} accepted steps to t={checked['endTime']:g} s"
    if checked.get('complete') is False:title+=f"\nStopped before requested t={checked['requestedEndTime']:g} s; trailing attempt not audited"
    fig.suptitle(title,fontsize=16)
    fig.supxlabel('Every retained segment-final state and restart chain audited. Monitor curves are recorded data. No mesh/time independence or saturated shedding qualification.',fontsize=10)
    args.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(str(args.output)+'.png',dpi=145);plt.close(fig)
    peaks=[dict(time=times[i],forceY=fy[i]) for i in range(1,len(fy)-1) if times[i]>5 and fy[i]>0 and fy[i]>fy[i-1] and fy[i]>=fy[i+1]]
    evidence={**checked,'vorticityDefinition':'Least-squares gradient of exported final cell velocity; dv/dx-du/dy',
        'vorticityRange':[min(omega),max(omega)],'positiveLocalMonitorPeaks':peaks,
        'plotSourceSha256':n.sha256_file(Path(__file__)),'verifiedRecordSha256':n.sha256_file(args.verified),
        'finalFieldAuditSha256':n.sha256_file(Path(str(prefix)+'.plot-audit.json'))}
    Path(str(args.output)+'.json').write_text(json.dumps(evidence,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verified',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
