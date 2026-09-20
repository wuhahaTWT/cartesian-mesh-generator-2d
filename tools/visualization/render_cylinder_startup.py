#!/usr/bin/env python3
"""Plot independently re-audited startup observations, without qualification claims."""
import argparse,json,math,sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_transient_flow as audit

def main(args):
    records=[]
    fig,axes=plt.subplots(1,3,figsize=(16,4.9),layout='constrained')
    for mesh,prefix in [(args.coarse_mesh,args.coarse),(args.fine_mesh,args.fine)]:
        q=json.loads(Path(str(prefix)+'.json').read_text());a=audit.verify(mesh,prefix,Path(str(prefix)+'-plot-audit.json'));assert a['valid']
        if q['case']!='external' or q['nu']!=.02 or q['speed']!=1 or not math.isclose(q['time'],5):
            raise ValueError('This observation plot requires the recorded Re100 startup controls and t=5')
        records.append({'mesh':str(mesh),'prefix':str(prefix),'summary':q,'auditSha256':audit.native.sha256_file(Path(str(prefix)+'-plot-audit.json'))})
        for ax,key in zip(axes[1:],('forceX','forceY')):
            ax.plot([r['time'] for r in a['history']],[r[key] for r in a['history']],label=f"{q['cells']:,} cells")
    mesh=audit.native.read_cm2d(args.fine_mesh);m=audit.native.measure(mesh,1e-11,1e-9)
    cells=audit.native.read_cells(Path(str(args.fine)+'.cells.csv'),mesh,m,'external')
    collection=PolyCollection([[mesh.vertices[i] for i in c.vertices] for c in mesh.cells],
        array=[math.hypot(c['u'],c['v']) for c in cells],cmap='viridis',edgecolors='#37464a',linewidths=.18)
    axes[0].add_collection(collection);axes[0].set(xlim=(-2,7),ylim=(-3,3),aspect='equal',xlabel='x [m]',ylabel='y [m]',title=f"Actual final fine-grid field, t={records[-1]['summary']['time']:g} s")
    fig.colorbar(collection,ax=axes[0],label='Speed [m/s]',shrink=.7)
    axes[1].set(title='Streamwise force during startup',xlabel='Physical time [s]',ylabel='Fx / density / depth [m³/s²]',yscale='log')
    axes[2].set(title='Strong coarse-grid transverse bias',xlabel='Physical time [s]',ylabel='Fy / density / depth [m³/s²]')
    for ax in axes[1:]:ax.grid(alpha=.25);ax.legend()
    fig.suptitle('Circular cylinder startup: D=2 m, U=1 m/s, nu=0.02 m²/s (Re=100)',fontsize=15)
    fig.supxlabel('Different meshes and adaptive histories. Discrete balances pass; this is not a spatial/time convergence or vortex-shedding qualification.',fontsize=10)
    args.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(str(args.output)+'.png',dpi=145);plt.close(fig)
    Path(str(args.output)+'.json').write_text(json.dumps({'valid':True,'scope':'Startup observations only; no force or shedding accuracy qualification','runs':records,'sourceSha256':audit.native.sha256_file(Path(__file__))},indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ('coarse-mesh','fine-mesh','coarse','fine','output'):p.add_argument('--'+k,type=Path,required=True)
    main(p.parse_args())
