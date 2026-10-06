#!/usr/bin/env python3
"""Summarize native RT1 transport experiments; no PDE is assembled in Python.
Run after native-laminar-p1-oseen.py and the retained readback scripts.
"""
from pathlib import Path
import hashlib,json,math
root=Path(__file__).resolve().parents[2];base=root/'outputs/cloud-laminar/p1-convection'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def localize(x):
 if isinstance(x,str):return x.replace(str(root)+'/', '')
 if isinstance(x,list):return [localize(v) for v in x]
 if isinstance(x,dict):return {k:localize(v) for k,v in x.items()}
 return x
final=read(base/'final-state-evidence.json');affine=read(base/'affine-final-readback.json')
for k,r in final.items():
 if k.startswith('split-rt1-'):
  prior=read(base/k.removeprefix('split-')/'result.json')
  r['predecessorIterations']=len(prior['iterations']);r['predecessorSeconds']=prior['totalSeconds']
  r['totalFromZeroSecondsIncludingUnsplitPredecessor']=r['totalSeconds']+prior['totalSeconds']
 else:r['totalFromZeroSecondsIncludingUnsplitPredecessor']=r['totalSeconds']
negative={}
for f in sorted(base.glob('poly-*/result.json')):
 if 'rotation' in f.parent.name:
  r=read(f)
  if r.get('completed'):negative[f.parent.name]=r['iterations'][-1]['nonlinearCheck']
poly={}
for f in sorted(base.glob('split-poly-*/result.json')):
 r=read(f)
 if r.get('completed'):poly[f.parent.name]=r['iterations'][-1]['nonlinearCheck']
inventory=[]
for f in sorted(base.glob('*/result.json')):
 r=read(f)
 if 'iterations' not in r:continue
 inventory.append({'case':f.parent.name,'completedAggregate':r.get('completed',False),'iterationsInAggregate':len(r['iterations']),'secondsInAggregate':r.get('totalSeconds'),'sha256':sha(f),'failure':r.get('failure')})
closedStokes={}
for f in sorted((root/'outputs/cloud-laminar/p1-stress').glob('closed-*.result.json')):
 r=read(f);closedStokes[f.stem]={'args':r['args'],'recovery':r['recovery'],'totalSeconds':r['totalSeconds']}
orders={}
for geometry in ('square','sheared'):
 for nu in (1,.1):
  candidates=[v['metrics'] for k,v in final.items() if geometry in k and 'cylinder' not in k and (('nu01' in k)==(nu==.1))]
  candidates.sort(key=lambda m:m['cells']);a,b=candidates[-2:]
  orders[f'{geometry}-nu{nu}']={key:math.log(a[key]/b[key],2) for key in ['velocityRms','pressureRms','rt1VelocityRms','manufacturedTractionMomentRms']}
sources=['native-laminar-p1-oseen.cpp','native-laminar-p1-transport.cpp','native-laminar-p1-oseen.py','native-laminar-p1-stress.cpp','native-laminar-hybrid-stokes-p1.cpp']
result={'schema':'native-laminar-compatible-convection-v1','status':'research only; production equations/defaults/quality gates unchanged','method':{
 'unknowns':'P1 cell/face velocity, P1 kinematic static pressure; symmetric weak-gradient viscosity; RT1 compatible load and advector/transported/test velocity; same B/-B^T coupling; local 8x8 elimination',
 'convection':'Conservative RT1 fan-volume derivative plus radial upwind jumps and original-face hybrid flux. One normal flux per interface. Natural outlet adds beta.n uF.vF.',
 'quadratureFix':'Split each edge exactly at the P1 beta.n sign-change root before Gauss quadrature; no fitted epsilon or parameter sweep. Same rule for assembly and energy identity.',
 'energy':'Pressure work cancels with the same weak divergence. Convection identity includes fan/face upwind dissipation, divergence work and boundary kinetic-energy flux. Raw extension reactions and physical boundary traction work are distinct.',
 'boundary':'True embedded wall has both velocity trace modes zero. Inlet (1,0); horizontal sides symmetry v=0 and natural tangential traction; right outlet natural (2 nu sym G - p I)n=0. Static pressure, not Bernoulli.',
 'boundaryLimitation':'Natural stress outlet differs from the original product pressure-Dirichlet/zero-normal-gradient-velocity stencil. These are coupled Re20 controlled subproblems on original meshes, not a repair certification of the exact original branch.',
 'iterationControl':'Picard with no relaxation. Both complete uncondensed residual and coefficient change <=1e-9, research normalization U=L=1, p scale U^2; not a product or physical acceptance threshold.',
 'maxima':'Pressure range includes all P1 cell vertices; velocity maxima are sampled P2 potential/RT1 values, not certified continuous maxima.'},
 'finalCases':final,'lastRefinementOrders':orders,'affineControls':{'count':len(affine),'maxVelocityRms':max(x['velocityRms'] for x in affine.values()),'maxPressureRms':max(x['pressureRms'] for x in affine.values()),'nativeReadbacks':affine},
 'negativeCellTestRotation':negative,'quadraticIrregularControls':poly,'closedStokesFactorControls':closedStokes,
 'unsplitQuadratureComparison':read(base/'quadrature-comparison.json'),'splitQuadratureComparison':read(base/'split-quadrature-comparison.json'),
 'runInventory':inventory,'inventoryCaution':'Mutable snapshot indexes of several early runs are incomplete although final state/logs exist; retained as partial, not counted as completed. Three affected small cases rerun from zero with immutable markers. Final selected states natively read back, exact state-byte roundtrip and complete cell/face row counts verified.',
 'cost':'Serial nonexclusive wall times include native assembly, sparse LU, two algebraic residual corrections, native recovery/check and hashing/export. Exclude compilation and extra readback. Split continuation cost includes unsplit predecessor explicitly; no claim that continuation-only cost is from zero. LU memory scaling remains unresolved.',
 'limitations':['RT1 solenoidal reconstructed transport space is affine exact but does not contain general quadratic velocity; irregular quadratic Poiseuille controls are not exact.','P2 potential velocity and RT1 transport velocity have different error orders; neither may be substituted for the other.','Closed no-slip wall traction moments only about order1.3; not qualified.','Only two original fixed-polygon external space levels, no outer-domain or tolerance qualification yet.','No product integration, physical-time or backflow qualification, no completed original-BC Re20 comparison.','No new full101 native/frontend198/Electron/Mac/package validation.'],
 'validation':{'nativeWerror':True,'nativeCtest':'flow_face, flow_boundary, flow_branch_certificate: 3/3; log retained','geometry':'original makeFvMesh2D gate, areas/atomic faces unchanged; no deleted cells, repaired contour or threshold change'},
 'sourceSha256':{s:sha(root/'artifacts/current'/s) for s in sources},'binarySha256':sha(root/'build/native-laminar-p1-transport-atomic'),
 'meshSha256':{str(p.relative_to(root)):sha(p) for p in [root/'outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid0.solver.cm2d',root/'outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1.solver.cm2d']}}
probe=base/'block-coarse-probe/verified-result.json'
if probe.exists():result['independentStokesBlockProbe']=read(probe)
(root/'artifacts/current/native-laminar-p1-convection.json').write_text(json.dumps(localize(result),indent=2)+'\n')
print('summary cases',len(final),'affine',len(affine),'run indexes',len(inventory),flush=True)
