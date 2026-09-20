from pathlib import Path
import json,math,sys
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import NullFormatter
sys.path.insert(0,'tools/verification')
import verify_native_flow as n
import verify_rotating_annulus as ring
root=Path('outputs/native-flow/annulus-facet-study');study=json.loads((root/'study.json').read_text());assert study['valid']
fig,axes=plt.subplots(2,2,figsize=(13,10),layout='constrained');samples=[]
for case in study['cases']:
 folder=root/f"f{case['facets']}";meshpath=folder/'mesh.solver.cm2d';prefix=folder/'flow'
 assert n.sha256_file(meshpath)==case['meshSha256']
 assert n.sha256_file(folder/'audit.json')==case['auditSha256']
 audit=n.verify_case(meshpath,prefix,'custom',.1,.5,n.argument_parser().parse_args(['--max-iterations','6000']));assert audit['valid'],audit['issues'];n.write_json(folder/'figure-audit.json',audit)
 metrics=ring.metrics(meshpath,prefix)
 for k,v in metrics.items():assert v==case['metrics'][k],(k,v,case['metrics'][k])
 mesh=n.read_cm2d(meshpath);m=n.measure(mesh,1e-11,1e-9);cells=n.read_cells(Path(str(prefix)+'.cells.csv'),mesh,m,'custom')
 pressure=lambda r:r*r/18-2*math.log(r)/9-1/(18*r*r)
 gauge=pressure(math.hypot(cells[0]['x'],cells[0]['y']));error=[(c['p']-pressure(math.hypot(c['x'],c['y']))+gauge)/.25 for c in cells]
 samples.append(dict(facets=case['facets'],metrics=metrics,maximumAbsoluteNormalizedPressureError=max(map(abs,error)),reauditSha256=n.sha256_file(folder/'figure-audit.json')))
 if case['facets'] in (64,256):
  ax=axes[0,0 if case['facets']==64 else 1]
  polygons=[[mesh.vertices[i] for i in c.vertices] for c in mesh.cells]
  coll=PolyCollection(polygons,array=error,cmap='RdBu_r',norm=SymLogNorm(linthresh=.005,vmin=-2,vmax=2),edgecolors='#586567',linewidths=.12)
  # The full range is shown: warn rather than clipping an unexpected larger value.
  assert max(map(abs,error))<=2
  ax.add_collection(coll);ax.set(xlim=(-1.04,1.04),ylim=(-1.04,1.04),aspect='equal',xlabel='x [m]',ylabel='y [m]',title=f"{case['facets']} facets per wall / {len(cells):,} actual Cut-cell cells")
  fig.colorbar(coll,ax=ax,label='Signed pressure error / Uwall² (symmetric log)',shrink=.85,ticks=[-1,-.1,-.01,0,.01,.1,1])
x=[c['facets'] for c in study['cases']]
for key,label in [('velocityL2','Velocity L2 / Uwall'),('pressureL2','Pressure L2 / Uwall²'),('exactCircularNormalVelocityMismatchL2','Circle vs polygon wall: normal velocity / Uwall')]:
 axes[1,0].loglog(x,[c['metrics'][key] for c in study['cases']],'o-',label=label)
axes[1,0].set(xlabel='Facets per circular wall',ylabel='Normalized error',title='Fixed Cartesian background; geometry changes');axes[1,0].legend(fontsize=8);axes[1,0].grid(alpha=.25,which='both')
for key,label in [('rotorTorqueError','Inner-wall torque'),('housingTorqueError','Outer-wall torque')]:axes[1,1].semilogy(x,[c['metrics'][key] for c in study['cases']],'o-',label=label)
axes[1,1].set(xlabel='Facets per circular wall',ylabel='Relative torque error',title='Torque errors are not monotonic');axes[1,1].legend();axes[1,1].grid(alpha=.25,which='both')
for ax in axes[1]:
 ax.set_xticks(x,labels=list(map(str,x)));ax.xaxis.set_minor_formatter(NullFormatter())
fig.suptitle('Circular Couette reference vs polygon Cut-cell geometry',fontsize=16)
fig.supxlabel('Ri=0.5 m, Ro=1 m, inner speed=0.5 m/s, nu=0.1 m²/s. Pressure gauge: cell 0.\nGeometry and wall-sampling sensitivity; not a three-grid flow accuracy qualification.',fontsize=10)
prefix=Path('artifacts/current/native-annulus-facets');fig.savefig(str(prefix)+'.png',dpi=145);plt.close(fig)
n.write_json(Path(str(prefix)+'.json'),dict(valid=True,scope=study['scope'],controls=study['controls'],cases=samples,source=str(root/'study.json'),sourceSha256=n.sha256_file(root/'study.json'),figureDriverSha256=n.sha256_file(Path(__file__)),reference=ring.REFERENCE,notRun=['checkMesh']))
print(json.dumps(samples))
