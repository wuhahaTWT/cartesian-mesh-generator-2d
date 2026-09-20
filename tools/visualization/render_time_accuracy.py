#!/usr/bin/env python3
"""Re-audit native Taylor-Green time refinement and render measured fields."""
import argparse,csv,json,math,sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_transient_flow as audit
import compare_transient_steps as compare

def main(root,output,binary_snapshot=None):
    report={'valid':True,'scope':'One Taylor-Green decay, fixed grids; temporal self-convergence, not general CFD qualification','groups':{}}
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for name in ('cartesian','warped'):
        runset=root/name/'runs.json';source=json.loads(runset.read_text());mesh=Path(source['mesh'])
        comparison=compare.compare(runset,root/name/'comparison.json',binary_snapshot)
        runs=[]
        for entry in source['runs']:
            command=entry['command'];prefix=Path(command[command.index('--output')+1])
            result=audit.verify(mesh,prefix,Path(str(prefix)+'.audit.json'));assert result['valid']
            runs.append({'dt':entry['dt'],'elapsedSeconds':entry['elapsedSeconds'],'command':command,
                         'analyticErrors':result['analyticErrors'],'computedMaxCourant':result['computedMaxCourant'],
                         'independentAuditSha256':audit.native.sha256_file(Path(str(prefix)+'.audit.json'))})
        report['groups'][name]={'comparison':comparison,'runs':runs,'runsetSha256':audit.native.sha256_file(runset)}
        axes[0,1].loglog([r['dt'] for r in runs],[r['analyticErrors']['relativeVelocityL2'] for r in runs],'o-',label=name)
        distances=comparison['adjacentDistances']
        axes[1,0].loglog([d['coarseDt'] for d in distances],[d['vectorL2'] for d in distances],'o-',label=name)
    mesh_path=root/'warped.solver.cm2d';mesh=audit.native.read_cm2d(mesh_path);measured=audit.native.measure(mesh,1e-11,1e-9)
    fine=root/'warped/dt-0.005/result';tight=root/'tight/dt-0.005/result'
    original=audit.native.read_cells(Path(str(fine)+'.cells.csv'),mesh,measured,'taylor-green')
    tightened=audit.native.read_cells(Path(str(tight)+'.cells.csv'),mesh,measured,'taylor-green')
    checked=audit.verify(mesh_path,tight,root/'tight-audit.json');assert checked['valid']
    baseline=json.loads(Path(str(fine)+'.json').read_text());tighter=json.loads(Path(str(tight)+'.json').read_text())
    def same_physics(other):
        for key in ('case','cells','nu','speed','convection','viscousStress','temporalDiscretization','temporalFaceInterpolation','velocityRelaxation','time'):
            a,b=baseline[key],other[key]
            if isinstance(a,(int,float)):
                if not audit.native.close(a,b,1e-12,1e-10):raise ValueError('Different comparison '+key)
            elif a!=b:raise ValueError('Different comparison '+key)
    same_physics(tighter)
    if baseline['case']!='taylor-green' or not tighter['tolerance']<baseline['tolerance'] or tighter['dt']!=baseline['dt']:
        raise ValueError('Tightening must retain dt and decrease the stopping tolerance')
    report['tightening']={'tolerances':[baseline['tolerance'],tighter['tolerance']],
      'maxVelocityDifference':max(math.hypot(a['u']-b['u'],a['v']-b['v']) for a,b in zip(original,tightened)),
      'maxPressureDifference':max(abs(a['p']-b['p']) for a,b in zip(original,tightened)),
      'auditSha256':audit.native.sha256_file(root/'tight-audit.json')}
    polygons=[[mesh.vertices[i] for i in cell.vertices] for cell in mesh.cells]
    field=PolyCollection(polygons,array=[math.hypot(c['u'],c['v']) for c in original],cmap='viridis',edgecolors='#24313b',linewidths=.15)
    axes[0,0].add_collection(field);axes[0,0].autoscale();axes[0,0].set_aspect('equal')
    axes[0,0].set(title='Actual warped grid and final speed',xlabel='x [m]',ylabel='y [m]');fig.colorbar(field,ax=axes[0,0],label='Speed [m/s]')
    adaptive=audit.verify(mesh_path,root/'adaptive',root/'adaptive-audit.json');assert adaptive['valid']
    summary=json.loads((root/'adaptive.json').read_text())
    same_physics(summary)
    report['adaptive']={'summary':summary,'analyticErrors':adaptive['analyticErrors'],'auditSha256':audit.native.sha256_file(root/'adaptive-audit.json')}
    for label,prefix in [('fixed dt=0.005',fine),(f"adaptive CFL <= {summary['targetCourant']:g}",root/'adaptive')]:
        result=audit.verify(mesh_path,prefix,Path(str(prefix)+'.audit.json'))
        history=result['analyticDecay']['history']
        axes[1,1].plot([r['time'] for r in history],[100*r['relativeEnergyError'] for r in history],label=label)
    axes[0,1].set(xlabel='Fixed time step [s]',ylabel='Velocity L2 / exact final velocity L2',title='Analytic error includes space and time')
    axes[1,0].set(xlabel='Coarser dt of adjacent pair [s]',ylabel='Velocity difference L2 [m/s]',title='Temporal self-convergence on each fixed grid')
    axes[1,1].set(xlabel='Physical time [s]',ylabel='Energy error relative to exact [%]',title='CFL control is not a temporal error estimator')
    for ax in (axes[0,1],axes[1,0],axes[1,1]):ax.grid(alpha=.2,which='both');ax.legend(fontsize=8)
    fig.suptitle(f"Native incompressible Taylor-Green decay | backward Euler\nnu = {baseline['nu']:g} m²/s, initial speed scale = {baseline['speed']:g} m/s, target time = {baseline['time']:g} s",fontsize=14)
    output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(str(output)+'.png',dpi=145);plt.close(fig)
    report['sourceSha256']=audit.native.sha256_file(Path(__file__))
    Path(str(output)+'.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'valid':True,'tightening':report['tightening'],'adaptive':adaptive['analyticErrors']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--binary-snapshot',type=Path)
    a=p.parse_args();main(a.root,a.output,a.binary_snapshot)
