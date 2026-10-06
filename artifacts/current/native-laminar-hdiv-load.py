#!/usr/bin/env python3
"""Native RT0 load reconstruction experiments; Python only orchestrates and compares fields."""
import hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/cloud-laminar/hdiv-load'
BINARY=ROOT/'build/native-laminar-hdiv-load'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(p):return np.genfromtxt(p,delimiter=',',names=True)
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 build=['g++','-std=c++20','-O2','-Wall','-Wextra','-Wpedantic','-Wno-unused-parameter','-Werror','-Iinclude','artifacts/current/native-laminar-hdiv-load.cpp','build/libcartmesh2d_fv.a','build/libcartmesh2d.a','-o',str(BINARY)]
 audits={'cylinder':'outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1.solver.cm2d','annulus':'outputs/cloud-laminar/annulus-joint-7.solver.cm2d'}
 if '--summarize' not in sys.argv:
  subprocess.run(build,cwd=ROOT,check=True)
  with (OUT/'runs.jsonl').open('w') as f:subprocess.run([str(BINARY),str(OUT)],cwd=ROOT,stdout=f,check=True)
  for label,path in audits.items():
   with (OUT/(label+'-audit.json')).open('w') as f:subprocess.run([str(BINARY),'--audit',path],cwd=ROOT,stdout=f,check=True)
 runs=[json.loads(line) for line in (OUT/'runs.jsonl').read_text().splitlines()]
 good=[r for r in runs if 'rejected' not in r]
 groups=[]
 for r in good:
  if not r['label'].endswith('-0'):continue
  label=r['label'][:-2];base=table(OUT/(r['label']+'.cells.csv'));bf=table(OUT/(r['label']+'.faces.csv'))
  item={'label':label,'cells':r['cells'],'baseVelocityRms':r['velocityRms'],'basePressureRms':r['pressureRms'],'gradientControls':[]}
  for value in (1,10000):
   target=next(x for x in good if x['label']==label+'-'+str(value));c=table(OUT/(target['label']+'.cells.csv'));f=table(OUT/(target['label']+'.faces.csv'))
   assert np.array_equal(c['area'],base['area']) and np.array_equal(c['cell'],base['cell'])
   a=c['area'];du=np.hypot(c['u']-base['u'],c['v']-base['v']);df=np.hypot(f['u']-bf['u'],f['v']-bf['v'])
   dp=c['p']-base['p']-c['p_exact_plus_gauge']+base['p_exact_plus_gauge'];dp-=np.dot(a,dp)/sum(a)
   item['gradientControls'].append({'amplitude':value,'velocityRms':float(np.sqrt(np.dot(a,du*du)/sum(a))),'velocityMax':float(max(du)),'faceVelocityMax':float(max(df)),'pressureShiftRms':float(np.sqrt(np.dot(a,dp*dp)/sum(a))),'pressureShiftMax':float(max(abs(dp)))})
  groups.append(item)
 lifted=[r for r in good if '-lift-' in r['label']];lg=[c for g in groups if '-lift' in g['label'] for c in g['gradientControls']]
 record={'schema':1,'scope':'Native research-only k=0 hybrid Laplacian Stokes with RT0 H(div) body-force reconstruction. No product/default/geometry/quality-gate change. Not old collocated solver repair or full Navier-Stokes qualification.',
 'method':'Centroid fan, each triangle R(v)=a_j+b(x-c_T), b=D_T(v)/2. Exact outward face normal moments and internal radial normal continuity, with remaining circulation minimizing integral |R(v)-v_T|^2. Integrate arbitrary f(x).R(v) by 6x6 Duffy Gauss quadrature. Same native k=0 viscous stiffness and same B/-B^T for lifted and cell-load controls.',
 'identity':'For every local test basis, integral grad(phi).R(v) = sum_F mean_F(phi)*(v_F.S_TF) - mean_T(phi)*B_T(v). Assembly cancels shared-face terms; gradient forces change only projected pressure for fixed Dirichlet velocity.',
 'units':'Lref=1 m, Uref=1 m/s, nu=1 m^2/s. phi=Uref^2*((x/Lref)^3+(x/Lref)*(y/Lref)^2); lambda dimensionless=0,1,10000. f=-nu Laplacian(u)+lambda grad(phi) in m/s^2. u error m/s, p error m^2/s^2, div s^-1, pressure work m^4/s^3. Local load identity coefficient per unit velocity basis m^3/s^2; trace residual per unit velocity m. No amplitude is an acceptance threshold.',
 'manufacturedProblems':{'hydrostatic':'u=0,p=lambda phi','polynomial':'u=(x^2,-2xy) in unit reference coordinates,p=lambda phi','vortex':'u=(sin(pi x)cos(pi y),-cos(pi x)sin(pi y)),p=lambda phi'},
 'boundary':'Exact manufactured velocity face averages on all boundary faces. Zero velocity hydrostatic control is impermeable; nonzero manufactured traces are not physical cylinder or annulus boundary conditions.',
 'reference':'https://arxiv.org/abs/2203.07180v3 Sections 2.4 and 3.2, read via PDF. Implementation is a direct two-dimensional k=0 constrained lift, not full paper code/high-order or convection implementation.',
 'buildCommand':build,'solveCount':len(good),'rejectedMeshes':[r for r in runs if 'rejected' in r],'nativeSolveSeconds':sum(r['seconds'] for r in good),'maximumLiftedGradientVelocityRms':max(c['velocityRms'] for c in lg),'maximumLiftedGradientFaceVelocity':max(c['faceVelocityMax'] for c in lg),'maximumLiftedPressureShiftRms':max(c['pressureShiftRms'] for c in lg),'maximumLiftedDivergence':max(r['maxDivergence'] for r in lifted),'maximumLiftedLinearResidual':max(r['linearResidual'] for r in lifted),'groups':groups,'realMeshReconstructionChecks':[],
 'limits':['Only low-order Stokes solves <=144 cells; dense <=1500 unknowns. Full original 17260-cell and 6606-cell checks cover local reconstruction and integration identity, not global flow solve.','Lowest-order viscosity still gives poor sheared vortex pressure accuracy; gradient invariance does not qualify general pressure accuracy.','Centroid fan must have positive triangles and closed radial endpoints within construction roundoff. No polygon alteration or welding. Not proved for arbitrary non-star-shaped cells.','No convection, traction boundary, symmetric-stress qualification, physical time, old nonlinear branch selection, Electron/macOS/packaged app validation.']}
 for label,path in audits.items():
  r=json.loads((OUT/(label+'-audit.json')).read_text());r['meshSha256']=sha(ROOT/path);record['realMeshReconstructionChecks'].append(r)
 record['sourceHashes']={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'artifacts/current/native-laminar-hdiv-load.cpp',ROOT/'artifacts/current/native-laminar-hybrid-stokes.cpp']}
 record['binarySha256']=sha(BINARY);record['libraryHashes']={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'build/libcartmesh2d_fv.a',ROOT/'build/libcartmesh2d.a']}
 record['rawOutputHashes']={p.name:sha(p) for p in sorted(OUT.glob('*')) if p.is_file() and p.suffix in ('.csv','.json','.jsonl')}
 (ROOT/'artifacts/current/native-laminar-hdiv-load.json').write_text(json.dumps(record,indent=2)+'\n')
 print(json.dumps({k:record[k] for k in ['solveCount','nativeSolveSeconds','maximumLiftedGradientVelocityRms','maximumLiftedGradientFaceVelocity','maximumLiftedPressureShiftRms','maximumLiftedDivergence']}))
if __name__=='__main__':main()
