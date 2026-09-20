#!/usr/bin/env python3
"""Actual initial state and cylinder response; no shedding-accuracy qualification."""
import argparse,json,math,sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_transient_flow as audit


def main(args):
    result=audit.verify(args.mesh,args.prefix,Path(str(args.prefix)+'.plot-audit.json'))
    if not result['initialVortex']:raise ValueError('A retained initial vortex is required')
    if result['case']!='external':raise ValueError('This plot requires an external-flow startup')
    mesh=audit.native.read_cm2d(args.mesh);m=audit.native.measure(mesh,1e-11,1e-9)
    cells=audit.native.read_cells(Path(str(args.prefix)+'.cells.csv'),mesh,m,'external')
    state={line.split()[0]:list(map(float,line.split()[2:])) for line in Path(str(args.prefix)+'.initial.checkpoint').read_text().splitlines() if line.split()[0] in ('U','V')}
    polygons=[[mesh.vertices[i] for i in c.vertices] for c in mesh.cells]
    fig,axes=plt.subplots(1,3,figsize=(16,4.8),layout='constrained')
    for ax,speeds,title in [(axes[0],[math.hypot(u,v) for u,v in zip(state['U'],state['V'])],'Retained initial speed, t=0'),
        (axes[1],[math.hypot(c['u'],c['v']) for c in cells],f"Computed speed, t={result['time']:g} s")]:
        collection=PolyCollection(polygons,array=speeds,cmap='viridis',edgecolors='#49616a',linewidths=.13)
        ax.add_collection(collection);ax.set(xlim=(-1.5,7),ylim=(-2.5,2.5),aspect='equal',title=title,xlabel='x [m]',ylabel='y [m]')
        fig.colorbar(collection,ax=ax,label='Speed [m/s]',shrink=.67)
    rows=result['history'];amplitude=100*abs(result['initialVortex']['definition']['peakSpeed'])/result['controls']['speed']
    axes[2].plot([r['time'] for r in rows],[r['forceY'] for r in rows],label=f'Explicit {amplitude:g}% initial vortex')
    baseline=None
    if args.baseline:
        baseline=audit.verify(args.mesh,args.baseline,Path(str(args.baseline)+'.initial-comparison-audit.json'))
        if baseline['initialVortex'] or baseline['time']!=result['time'] or baseline['controls']!=result['controls']:
            raise ValueError('Baseline physical controls or target time differ')
        if len(baseline['history'])!=len(rows) or any(not math.isclose(a[k],b[k],rel_tol=1e-12,abs_tol=1e-14)
            for a,b in zip(rows,baseline['history']) for k in ('time','dt')):raise ValueError('Baseline accepted time steps differ')
        axes[2].plot([r['time'] for r in baseline['history']],[r['forceY'] for r in baseline['history']],label='Unperturbed symmetric start')
    axes[2].set(title='Transverse force response',xlabel='Time [s]',ylabel='Fy / density / depth [m³/s²]');axes[2].grid(alpha=.2);axes[2].legend(fontsize=9)
    fig.suptitle(f"Explicit compact initial condition on {len(mesh.cells):,} actual Cut-cell cells",fontsize=15)
    fig.supxlabel('One startup trajectory. Initial/final discrete audits pass; force accuracy, long-time shedding and mesh/time independence remain unqualified.',fontsize=10)
    args.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(str(args.output)+'.png',dpi=145);plt.close(fig)
    evidence=dict(valid=True,scope='Initial state and short transient response only; not shedding qualification',
        initial=result['initialVortex'],time=result['time'],cells=len(mesh.cells),finalForceY=rows[-1]['forceY'],
        maxAbsoluteForceY=max(abs(r['forceY']) for r in rows),meshSha256=audit.native.sha256_file(args.mesh),
        auditSha256=audit.native.sha256_file(Path(str(args.prefix)+'.plot-audit.json')),
        sourceSha256=audit.native.sha256_file(Path(__file__)))
    if baseline:
        evidence['unperturbedFinalForceY']=baseline['history'][-1]['forceY']
        evidence['baselineAuditSha256']=audit.native.sha256_file(Path(str(args.baseline)+'.initial-comparison-audit.json'))
    Path(str(args.output)+'.json').write_text(json.dumps(evidence,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ('mesh','prefix','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--baseline',type=Path);main(p.parse_args())
