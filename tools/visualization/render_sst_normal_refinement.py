#!/usr/bin/env python3
"""Show actual nested meshes and bounded-solver evidence, including failure."""
import argparse,csv,hashlib,json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();study=json.loads(a.study.read_text())
def checked(record):
    path=Path(record['path'])
    if hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('evidence hash differs')
    return path
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for axis,key in zip(axes[0],('original','refined')):
    record=study['meshes'][key];mesh=native.read_cm2d(checked(record))
    polygons=[]
    for c in mesh.cells:
        points=[mesh.vertices[v] for v in c.vertices]
        if max(v[0] for v in points)>.95 and min(v[0] for v in points)<1.:
            if min(v[1] for v in points)<.003:polygons.append(points)
    axis.add_collection(PolyCollection(polygons,facecolors='#f2f5f7',edgecolors='#557381',linewidths=.6))
    axis.axhspan(0,study['geometry']['preservedBandTop'],color='#c28e43',alpha=.35,label='Unchanged first four rows')
    axis.set(xlim=(.95,1.),ylim=(0,.003),xlabel='x [m]',ylabel='y [m]',title=f"{key.capitalize()}: {len(mesh.cells):,} actual cells (vertical zoom)")
    axis.legend(fontsize=8,loc='upper left')
run=study['run']
if run.get('history'):
    with checked(run['history']).open() as stream:rows=list(csv.DictReader(stream))
    for key,color,label in [('momentumResidual','#317c88','Momentum'),('kCellResidual','#a47749','k local'),('omegaCellResidual','#975c86','omega local')]:
        axes[1,0].semilogy([int(r['iteration']) for r in rows],[float(r[key]) for r in rows],color=color,label=label)
    axes[1,0].axhline(1e-7,color='#317c88',ls=':',lw=1)
    axes[1,0].axhline(1e-9,color='#555',ls=':',lw=1)
    axes[1,0].set(xlabel='SIMPLE iteration',ylabel='Residual (different original definitions)',title='Refined-grid solver history')
    axes[1,0].legend(fontsize=9);axes[1,0].grid(alpha=.2)
else:axes[1,0].text(.1,.5,'No returned full history; see retained raw log.')
axes[1,1].axis('off')
d=run['diagnostics'];status='Original-equation audit PASSED' if run['audit']['valid'] else 'NOT ACCEPTED'
lines=[status,'',f"Stop: {d['stopReason']}",f"Iterations: {d['iterations']:,}",f"Native solve: {d['solveSeconds']:.2f} s",'',
       'Same domain, streamwise mesh, first four rows,',
       'physics, initial values and stopping limits.',
       'Each remaining normal interval bisected once.','']
if 'physicalComparison' in study:
    lines.append('Accepted-field friction changes:')
    for station in study['physicalComparison']['stations']:
        lines.append(f"  x={station['x']:.2f}: {100*station['signedRelativeCfChange']:+.2f}%")
else:
    lines+=['No accepted physical comparison is reported.',
            'An unfinished field cannot establish grid sensitivity.']
axes[1,1].text(.03,.97,'\n'.join(lines),va='top',fontsize=11,linespacing=1.5)
fig.suptitle('Controlled wall-normal refinement | Re_plate = 10 million | SST-2003m',fontsize=15)
fig.supxlabel('Geometry, original-equation convergence and physical accuracy are separate checks. No thresholds were relaxed.',fontsize=10)
a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=150)
