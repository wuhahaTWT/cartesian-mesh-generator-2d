#!/usr/bin/env python3
"""Plot recorded SST pressure timings; never treat fixed-work fields as accepted."""
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

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--study',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
a=parser.parse_args();study=json.loads(a.study.read_text())
def checked(record):
    path=Path(record['path'])
    if hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:
        raise ValueError('evidence changed: '+str(path))
    return path
mesh=native.read_cm2d(checked(study['mesh']))
for run in study['runs'].values():
    for item in run['files']:checked(item)
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
polygons=[]
for cell in mesh.cells:
    points=[mesh.vertices[i] for i in cell.vertices]
    if min(x for x,y in points)<1 and max(x for x,y in points)>.95 and min(y for x,y in points)<.003:
        polygons.append(points)
axes[0,0].add_collection(PolyCollection(polygons,facecolors='#f1f4f6',edgecolors='#627b89',linewidths=.65))
axes[0,0].set(xlim=(.95,1),ylim=(0,.003),xlabel='x [m]',ylabel='y [m]',title=f'Same actual {len(mesh.cells):,}-cell mesh (vertical zoom)')
colors={'ic0':'#778b98','aggregation':'#a17b47','jacobi':'#955665'}
means=study['fixedWorkMeans'];methods=list(means)
for i,method in enumerate(methods):
    row=means[method];p=row['pressureSeconds'];total=row['nativeSeconds']
    axes[0,1].bar(i,p,color=colors[method],label='Pressure' if i==0 else None)
    axes[0,1].bar(i,total-p,bottom=p,color='#d8dfe2',label='Other work' if i==0 else None)
    axes[0,1].text(i,total+.35,f'{total:.2f} s',ha='center')
    values=[r['diagnostics']['solveSeconds'] for n,r in study['runs'].items() if n.startswith(method+'-') and r['completedFixed100']]
    axes[0,1].scatter([i]*len(values),values,color='#252c31',s=12,zorder=4)
axes[0,1].set(xticks=range(len(methods)),xticklabels=[m.upper() for m in methods],ylabel='Native solve time [s]',title='100 SIMPLE steps: repeated fixed-work timings')
axes[0,1].set_ylim(0,max(r['nativeSeconds'] for r in means.values())*1.18)
axes[0,1].legend(fontsize=9);axes[0,1].grid(axis='y',alpha=.15)
full=study['runs'].get('aggregation-full')
if full and 'diagnostics' in full:
    history=next((r for r in full['files'] if r['path'].endswith('.history.csv')),None)
    if history:
        with checked(history).open() as f:rows=list(csv.DictReader(f))
        for key,label,color in [('momentumResidual','Momentum','#357e91'),('kCellResidual','k local','#ab8142'),('omegaCellResidual','omega local','#876397')]:
            axes[1,0].semilogy([int(r['iteration']) for r in rows],[float(r[key]) for r in rows],label=label,color=color)
        axes[1,0].axhline(1e-7,color='#357e91',ls=':',lw=1)
        axes[1,0].axhline(1e-9,color='#666',ls=':',lw=1)
        axes[1,0].legend(fontsize=9)
    axes[1,0].set(xlabel='SIMPLE iteration',ylabel='Residual (different original definitions)',title='Aggregation: complete-run history')
    axes[1,0].grid(alpha=.15)
else:
    axes[1,0].axis('off');axes[1,0].text(.05,.5,'No completed full-run history available.')
axes[1,1].axis('off')
lines=['Same physics, initial fields and stopping limits.','IC0 remains the default.','', 'Jacobi: pressure iteration limit; no 100-step result.', 'Fixed-work timings are not convergence proof.','']
if full and 'diagnostics' in full:
    d=full['diagnostics'];accepted=full.get('audit',{}).get('valid') is True
    lines += [f"Full aggregation run: {d['stopReason']}",f"{d['iterations']:,} steps; {d['solveSeconds']:.2f} native seconds",'Independent original-equation audit: '+('PASS' if accepted else 'NOT ACCEPTED')]
    if 'fullComparison' in study:
        b=study['fullComparison']['baselineDiagnostics'];lines += ['',f"Earlier IC0: {b['iterations']:,} steps; {b['solveSeconds']:.2f} s",'Full-run times: one observation per method.']
lines += ['','Physical accuracy is not qualified by this test.']
axes[1,1].text(.03,.98,'\n'.join(lines),va='top',fontsize=10.5,linespacing=1.5)
fig.suptitle('SST-2003m | Re_plate = 10 million | Pressure-solver comparison',fontsize=15)
fig.supxlabel('Local serial macOS measurements. No tolerance relaxation, cell removal or hidden fallback.',fontsize=10)
a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=150)
