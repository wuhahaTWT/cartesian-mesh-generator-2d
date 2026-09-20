#!/usr/bin/env python3
"""Render real SST mesh data and full convergence histories from a hashed study."""
import argparse,csv,hashlib,json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();study=json.loads(a.study.read_text())
def checked(item):
    path=Path(item['path'])
    if hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:raise ValueError('study data hash mismatch')
    return path
reference=study['runs']['original'];candidate=study['runs']['completion']
chosen=candidate if candidate.get('fields') else reference
mesh=native.read_cm2d(checked(study['mesh']))
with checked(chosen['fields']).open() as stream:cells=list(csv.DictReader(stream))
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
polys=[];values=[]
for cell,row in zip(mesh.cells,cells):
    polygon=[mesh.vertices[i] for i in cell.vertices]
    if max(v[0] for v in polygon)>-.006 and min(v[0] for v in polygon)<.016 and min(v[1] for v in polygon)<.0004:
        polys.append(polygon);values.append(float(row['omega']))
collection=PolyCollection(polys,array=values,cmap='viridis',norm=LogNorm(),edgecolors='#728292',linewidths=.4)
axes[0,0].add_collection(collection);axes[0,0].set(xlim=(-.006,.016),ylim=(0,.0004),xlabel='x [m]',ylabel='y [m]',title='Actual final mesh and omega near the leading edge')
axes[0,0].xaxis.set_major_locator(MaxNLocator(5))
fig.colorbar(collection,ax=axes[0,0],label='Omega [1/s]')
for name,run,color in [('Original schedule',reference,'#7d8793'),('Completion sweep',candidate,'#347b78')]:
    if not run.get('history'):continue
    with checked(run['history']).open() as stream:rows=list(csv.DictReader(stream))
    x=[int(r['iteration']) for r in rows]
    axes[0,1].semilogy(x,[float(r['momentumResidual']) for r in rows],label=name,color=color)
    axes[1,0].semilogy(x,[float(r['omegaCellResidual']) for r in rows],label=name,color=color)
for axis,title,limit in [(axes[0,1],'Original momentum residual',1e-7),(axes[1,0],'Original omega local residual',1e-9)]:
    axis.axhline(limit,color='black',ls='--',lw=.8,label='Unchanged stopping limit')
    axis.set(xlabel='SIMPLE iteration',ylabel='Residual',title=title);axis.legend(fontsize=8);axis.grid(alpha=.2)
axes[1,1].axis('off')
lines=['Same mesh, physics and original acceptance gates','']
for name,run in [('Original',reference),('Completion',candidate)]:
    lines += [f"{name}: {'equation audit passed' if run['auditValid'] else 'not accepted'}",f"  SIMPLE iterations: {run['iterations']:,}",f"  Native solve: {run['solveSeconds']:.2f} s"]
    if 'completionPasses' in run:lines.append(f"  Completion sweeps: {run['completionPasses']}")
    lines.append('')
lines+=['Wall times are single-run observations.','They are not a controlled speedup benchmark.','Physical accuracy and mesh independence remain unqualified.']
axes[1,1].text(.02,.96,'\n'.join(lines),va='top',fontsize=11,linespacing=1.5)
fig.suptitle('12,800 cells | SST-2003m | Bounded completion of turbulence updates',fontsize=14)
fig.supxlabel('Final field is accepted only if the independent original-equation audit passes. Vertical mesh zoom.',fontsize=10)
a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=150)
