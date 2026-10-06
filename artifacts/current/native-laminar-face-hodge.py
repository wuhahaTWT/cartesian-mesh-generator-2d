#!/usr/bin/env python3
"""Run native face-Hodge assembly; only linear algebra and field comparison here.
Manufactured pressure diffusion is a separate structural test with Dirichlet
pressure on all boundaries. Closed-wall projection uses zero flux on ALL
boundaries. Neither run solves Re20 nor chooses a nonlinear physical branch.
"""
import hashlib,json,os,subprocess,time,sys
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/cloud-laminar/face-hodge'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(p):return np.genfromtxt(p,delimiter=',',names=True)
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 binary=ROOT/'build/native-laminar-face-hodge'
 command=['g++','-std=c++20','-O2','-Wall','-Wextra','-Wpedantic','-Werror','-Iinclude','artifacts/current/native-laminar-face-hodge.cpp','build/libcartmesh2d_fv.a','build/libcartmesh2d.a','-o',str(binary)]
 subprocess.run(command,cwd=ROOT,check=True)
 cases=[('skew8','skew:8'),('skew16','skew:16'),('skew32','skew:32'),('cylinder4716','outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid0.solver.cm2d'),('cylinder17260','outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1.solver.cm2d'),('cylinder66912','outputs/cloud-laminar/cylinder-joint/far-20.solver.cm2d'),('annulus576','outputs/cloud-laminar/annulus-joint-5.solver.cm2d'),('annulus6606','outputs/cloud-laminar/annulus-joint-7.solver.cm2d')]
 result={'schema':1,'scope':'Research-only joint face metric/pressure/divergence prototype. Not Re20 or a product default.','formula':'H_i=R_i R_i^T/V_i + tr(R_i R_i^T/V_i)/2 * (I-N_i(N_i^T N_i)^-1 N_i^T); N_i^T R_i=V_i I; B=signed incidence.','units':'q integral velocity normal: m^2/s; H dimensionless in 2D; 0.5 q^T H q: m^4/s^2. Projection step is unit 1 s; pressure multiplier m^2/s^2. Manufactured diffusion uses unit response 1 s, pressure amplitude 1 m^2/s^2.','limits':'Projection norm is face-Hodge kinetic energy, not proof for the existing collocated momentum norm. No convection, no-slip tangential stress or physical-time integration is included. Pressure manufactured Dirichlet traces do not represent the cylinder flow boundary conditions.','buildCommand':command,'binarySha256':sha(binary),'sourceSha256':sha(ROOT/'artifacts/current/native-laminar-face-hodge.cpp'),'runs':[]}
 if '--resume' in sys.argv:
  prior=json.loads((OUT/'runs.json').read_text())
  if prior['sourceSha256']!=result['sourceSha256'] or prior['binarySha256']!=result['binarySha256']:raise RuntimeError('changed native prototype; cannot resume')
  result['runs']=prior['runs']
 for label,mesh in cases:
  if any(r['label']==label for r in result['runs']):continue
  prefix=OUT/label;started=time.monotonic()
  run=[str(binary),mesh,str(prefix)]
  native=subprocess.run(run,cwd=ROOT,check=True,capture_output=True,text=True)
  (OUT/(label+'.native.json')).write_text(native.stdout)
  record=json.loads(native.stdout);record.update(label=label,command=run)
  if not mesh.startswith('skew:'):record['meshSha256']=sha(ROOT/mesh)
  c=table(prefix.with_suffix('.cells.csv'));f=table(prefix.with_suffix('.faces.csv'));t=table(prefix.with_suffix('.hodge.csv'))
  n=len(c);nf=len(f);H=sp.coo_matrix((t['value'],(t['row'].astype(int),t['col'].astype(int))),shape=(nf,nf)).tocsc()
  interior=f['neighbour']>=0;boundary=~interior;fid=np.arange(nf);owner=f['owner'].astype(int);neighbour=f['neighbour'][interior].astype(int)
  B=sp.coo_matrix((np.r_[np.ones(nf),-np.ones(sum(interior))],(np.r_[owner,neighbour],np.r_[fid,fid[interior]])),shape=(n,nf)).tocsc()
  # Same H and exact transpose B^T give a symmetric mixed Poisson system.
  K=sp.bmat([[H,-B.T],[-B,None]],format='csc')
  start=time.monotonic();lu=splu(K);factor_seconds=time.monotonic()-start
  record['poissonUnknowns']=nf+n;record['poissonFactorSeconds']=factor_seconds;record['factorNonzeros']=lu.L.nnz+lu.U.nnz
  outputs={}
  for kind in ('affine','nonlinear'):
   trace=np.where(boundary,f[kind+'Pressure'],0)
   source=np.zeros(n) if kind=='affine' else c['sourceIntegral']
   rhs=np.r_[-trace,-source];sol=lu.solve(rhs)
   for _ in range(2):sol+=lu.solve(rhs-K@sol)
   q,p=sol[:nf],sol[nf:];exactp=c[kind+'Pressure'];exactq=f[kind+'Flux'];dp=p-exactp;dq=q-exactq
   record[kind]={'pressureRms':float(np.sqrt(np.dot(c['area'],dp*dp)/sum(c['area']))),'pressureMax':float(max(abs(dp))),'fluxHodgeRms':float(np.sqrt(max(0,dq@(H@dq))/sum(c['area']))),'algebraicMaxResidual':float(max(abs(K@sol-rhs))),'maximumCellDivergenceResidual':float(max(abs(B@q-source)/c['area']))}
   outputs[kind+'Flux']=q;outputs[kind+'Pressure']=p
  del lu,K
  # Homogeneous normal perturbations on every boundary; no geometry or cell omitted.
  Hi=H[interior,:][:,interior];Bi=B[:,interior];Bg=Bi[:-1,:]
  K=sp.bmat([[Hi,Bg.T],[Bg,None]],format='csc')
  q0=f['trialFlux'][interior];rhs=np.r_[Hi@q0,np.zeros(n-1)]
  start=time.monotonic();lu=splu(K);sol=lu.solve(rhs)
  for _ in range(2):sol+=lu.solve(rhs-K@sol)
  q=sol[:sum(interior)];p=np.r_[sol[sum(interior):],0];dq=q0-q
  e0=.5*q0@(Hi@q0);e1=.5*q@(Hi@q);loss=.5*dq@(Hi@dq)
  record['projection']={'seconds':time.monotonic()-start,'energyBefore':float(e0),'energyAfter':float(e1),'removedEnergy':float(loss),'energyIdentityError':float(e0-e1-loss),'pressureWork':float(p@(Bi@q)),'maximumDivergence':float(max(abs(Bi@q)/c['area'])),'algebraicMaxResidual':float(max(abs(K@sol-rhs))),'boundaryNormalFluxMax':0.0}
  outputs['projectedFlux']=np.zeros(nf);outputs['projectedFlux'][interior]=q;outputs['projectionPressure']=p
  np.savez_compressed(str(prefix)+'.fields.npz',**outputs)
  record['seconds']=time.monotonic()-started
  record['outputHashes']={p.name:sha(p) for p in OUT.glob(label+'.*')}
  result['runs'].append(record)
  (OUT/'runs.json').write_text(json.dumps(result,indent=2)+'\n')
  print(label,json.dumps({'affine':record['affine'],'nonlinear':record['nonlinear'],'projection':record['projection'],'seconds':record['seconds']}),flush=True)
 accepted_branch_control(result)
 (ROOT/'artifacts/current/native-laminar-face-hodge.json').write_text(json.dumps(result,indent=2)+'\n')
def accepted_branch_control(result):
 # An already conservative flux must be invariant under projection, regardless
 # of whether its collocated state belongs to an anomalous or bounded branch.
 prefix=OUT/'cylinder17260';c=table(prefix.with_suffix('.cells.csv'));f=table(prefix.with_suffix('.faces.csv'));t=table(prefix.with_suffix('.hodge.csv'))
 n=len(c);nf=len(f);inside=f['neighbour']>=0;fid=np.arange(nf)
 H=sp.coo_matrix((t['value'],(t['row'].astype(int),t['col'].astype(int))),shape=(nf,nf)).tocsc()
 B=sp.coo_matrix((np.r_[np.ones(nf),-np.ones(sum(inside))],(np.r_[f['owner'].astype(int),f['neighbour'][inside].astype(int)],np.r_[fid,fid[inside]])),shape=(n,nf)).tocsc()
 Hi=H[inside,:][:,inside];Bg=B[:-1,inside];K=sp.bmat([[Hi,Bg.T],[Bg,None]],format='csc');lu=splu(K)
 records=[]
 for name in ('far-20-fixed128-grid1-face-limited-linear','far-20-fixed128-grid1-short-nu01-from-nu1'):
  path=ROOT/'outputs/cloud-laminar/cylinder-joint'/str(name+'.faces.csv');data=table(path)
  if not np.array_equal(data['face'].astype(int),fid):raise RuntimeError('face ordering mismatch')
  q=data['flux'];rhs=np.r_[np.zeros(sum(inside)),-(B@q)[:-1]];sol=lu.solve(rhs)
  for _ in range(2):sol+=lu.solve(rhs-K@sol)
  dq=np.zeros(nf);dq[inside]=sol[:sum(inside)]
  np.savez_compressed(OUT/(name+'.projection.npz'),inputFlux=q,correction=dq)
  records.append({'input':str(path.relative_to(ROOT)),'sha256':sha(path),'originalMaxDivergence':float(max(abs(B@q)/c['area'])),'fluxCorrectionHodgeRms':float(np.sqrt(max(0,dq@(H@dq))/sum(c['area']))),'maximumFluxCorrection':float(max(abs(dq))),'fixedBoundaryFluxCorrectionMax':float(max(abs(dq[~inside]))),'projectedMaxDivergence':float(max(abs(B@(q+dq))/c['area'])),'netBoundaryFlux':float(sum(B@q)),'gaugeCellArea':float(c['area'][-1]),'note':'Fixed boundary net-flux roundoff is left in the gauge cell, not discarded; all-cell divergence includes it.','outputSha256':sha(OUT/(name+'.projection.npz'))})
 result['driverSha256']=sha(Path(__file__))
 result['productLibraryHashes']={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'build/libcartmesh2d_fv.a',ROOT/'build/libcartmesh2d.a')}
 result['acceptedBranchControl']={'conclusion':'Both conservative accepted face fluxes remain unchanged to their residual scale. Projection alone cannot repair or select the collocated nonlinear branch; a compatible momentum/viscous/convection discretization is still required.','runs':records}
if __name__=='__main__':main()
