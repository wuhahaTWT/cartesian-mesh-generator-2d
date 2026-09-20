#!/usr/bin/env python3
"""Plot actual audited mesh/fields, physical profiles, and descriptive references."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--mesh-case',default='fine',help='case whose actual mesh is displayed')
    args=parser.parse_args()
    study=json.loads(args.study.read_text())
    cases=study['cases']
    if not cases or not all(c['equationAuditValid'] and not c['physicalAccuracyQualified'] for c in cases):
        raise ValueError('audited descriptive study required')
    # Display only the exact files whose strict audits produced the study.
    for c in cases:
        for record in c['files']:
            if hashlib.sha256(Path(record['path']).read_bytes()).hexdigest()!=record['sha256']:
                raise ValueError('field provenance differs: '+record['path'])
    colors=['#8f9eae','#b68d52','#247c88','#9a5e82']
    fig,axes=plt.subplots(2,3,figsize=(15,8),layout='constrained')
    selected=next(c for c in cases if c['name']==args.mesh_case)
    mesh=native.read_cm2d(Path(selected['files'][0]['path']))
    cells_path=next(Path(r['path']) for r in selected['files'] if r['path'].endswith('.cells.csv'))
    with cells_path.open() as stream:rows=list(csv.DictReader(stream))
    polys=[];values=[]
    for cell,row in zip(mesh.cells,rows):
        vertices=[mesh.vertices[v] for v in cell.vertices]
        if max(v[0] for v in vertices)>-.03 and min(v[0] for v in vertices)<.2 and min(v[1] for v in vertices)<.003:
            polys.append(vertices);values.append(float(row['u'])/selected['physicalInputs']['speed'])
    collection=PolyCollection(polys,array=values,cmap='viridis',edgecolors='#8495a1',linewidths=.3)
    axes[0,0].add_collection(collection)
    axes[0,0].set(xlim=(-.03,.2),ylim=(0,.003),xlabel='x [m]',ylabel='y [m]',title=f"Actual {selected['name']} mesh and u/U (vertical zoom)")
    fig.colorbar(collection,ax=axes[0,0],shrink=.8,label='u / U')
    reference=study.get('reference')
    for c,color in zip(cases,colors):
        label=f"{c['name']}: {c['cells']:,} cells"
        samples=c['wallSamples']
        axes[0,1].plot([s['x'] for s in samples],[s['Cf'] for s in samples],color=color,label=label)
        axes[1,0].plot([s['x'] for s in samples],[s['yPlus'] for s in samples],color=color,label=label)
        for i,p in enumerate(c['profiles']):
            axes[1,i+1].plot([v['uOverU'] for v in p['samples']],[v['y'] for v in p['samples']],'.-',
                             color=color,markersize=3,label=label)
    if reference:
        cf=reference['tables']['cf_plate_sstv.dat']['CFL3D']
        axes[0,1].plot([v[0] for v in cf],[v[1] for v in cf],'k--',lw=1,label='TMR SST-Vm / M=0.2 (context)')
        for i,p in enumerate(cases[0]['profiles']):
            table=next(v for k,v in reference['tables']['flatplate_u_sstv.dat'].items() if float(k.split('=')[1])==p['x'])
            axes[1,i+1].plot([v[0] for v in table],[v[1] for v in table],'k--',lw=1,label='TMR SST-Vm / M=0.2')
    axes[0,1].set(xlim=(0,2),ylim=(.002,.006),xlabel='x [m]',ylabel='Cf',title='Skin friction (leading-edge peak cropped)')
    axes[1,0].set(xlim=(0,2),xlabel='x [m]',ylabel='y+',yscale='log',title='First-cell wall resolution')
    axes[1,0].axhline(1,color='#333',ls=':',lw=1)
    for i,p in enumerate(cases[0]['profiles']):
        axes[1,i+1].set(xlim=(0,1.04),ylim=(0,.04),xlabel='u / U',ylabel='y [m]',title=f"Velocity profile at x={p['x']}")
    axes[0,2].axis('off')
    lines=['SAME-MODEL SENSITIVITY','']
    for comparison in study['comparisons']:
        lines += [comparison['from']+' → '+comparison['to']]
        for station in comparison['stations']:
            lines.append(f"  x={station['x']:.2f}: Cf {100*station['signedRelativeCfChange']:+.2f}%")
        lines.append(f"  Integrated Cd {100*comparison['signedRelativeCdChange']:+.2f}%")
        lines.append('')
    lines += ['Changes are observations, not error estimates.','Grid independence is not established.']
    axes[0,2].text(0,1,'\n'.join(lines),va='top',fontsize=10,linespacing=1.35)
    for axis in [axes[0,1],*axes[1]]:
        axis.grid(alpha=.18)
        axis.legend(fontsize=7,loc='best')
    fig.suptitle('Native SST-2003m | Re based on plate length = 10 million | Physical checks',fontsize=15)
    fig.supxlabel('All native fields passed independent equation audits. TMR uses a different SST variant; reference differences are not accuracy certification.',fontsize=10)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=150)


if __name__=='__main__':main()
