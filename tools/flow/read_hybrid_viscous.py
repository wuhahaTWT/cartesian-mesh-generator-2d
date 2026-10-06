import pathlib,csv,gzip,json,time,numpy as np,scipy.linalg as la
import argparse
p=argparse.ArgumentParser();p.add_argument('--output',type=pathlib.Path,required=True);args=p.parse_args();root=args.output;out={}
for name in ['original264','phase-L4','phase-L5','phase-L6']:
 out[name]={}
 for method in ['hybrid','linear']:
  prefix=root/name/method;report=json.loads(prefix.with_suffix('.json').read_text());faces=list(csv.DictReader(prefix.with_suffix('.faces.csv').open()));walls=[x for x in faces if int(x['neighbour'])==-1];metrics={}
  for state in ['analytic-centroid','steady']:
   data=[x for x in walls if x['state']==state];err=den=torqueErr=torqueDen=shearErr=shearDen=0.;groups={}
   for x in data:
    vals={k:float(v) for k,v in x.items() if k!='state'};f=np.array([vals['fluxX'],vals['fluxY']]);r=np.array([vals['referenceX'],vals['referenceY']]);err+=np.linalg.norm(f-r);den+=np.linalg.norm(r);torque=vals['x']*f[1]-vals['y']*f[0];tangent=np.array([-vals['sy'],vals['sx']]);tangent/=np.linalg.norm(tangent);shearErr+=abs(np.dot(f-r,tangent));shearDen+=abs(np.dot(r,tangent));torqueErr+=abs(torque-vals['referenceTorque']);torqueDen+=abs(vals['referenceTorque']);side='inner' if vals['x']*vals['sx']+vals['y']*vals['sy']<0 else 'outer';g=groups.setdefault(side,{'torque':0.,'referenceTorque':0.,'work':0.,'referenceWork':0.,'referenceWorkGivenFaceValue':0.,'forceX':0.,'forceY':0.});g['torque']+=torque;g['referenceTorque']+=vals['referenceTorque'];g['work']+=vals['work'];g['referenceWork']+=vals['referenceWork'];vr=.3+.2/(vals['x']**2+vals['y']**2);g['referenceWorkGivenFaceValue']+=np.dot(np.array([-vr*vals['y'],vr*vals['x']]),r);g['forceX']+=f[0];g['forceY']+=f[1]
   metrics[state]={'wallTractionRelativeL1':err/den,'wallShearRelativeL1':shearErr/shearDen,'wallTorqueRelativeL1':torqueErr/torqueDen,'walls':groups}
  report['wallMetrics']=metrics
  if method=='hybrid' or name=='original264':
   t=time.monotonic();n=2*report['cells'];a=np.zeros((n,n));mass=np.zeros(n)
   with gzip.open(prefix.with_suffix('.cell-matrix.csv.gz'),'rt') as f:
    for x in csv.DictReader(f):i=int(x['row']);a[i,int(x['column'])]=float(x['value']);mass[i]=float(x['area'])
   s=a/np.sqrt(mass[:,None]*mass[None,:]);sym=(s+s.T)/2;v=la.eigh(sym,subset_by_index=[0,11],eigvals_only=True);spec={'massScaledRelativeAsymmetry':float(abs(s-s.T).max()/abs(s).max()),'lowestSymmetricPartEigenvaluesPerSecond':v.tolist(),'matrixReadbackAndSpectrumSeconds':time.monotonic()-t}
   if name=='original264':
    e=la.eigvals(s);spec.update(minimumRealEigenvaluePerSecond=float(e.real.min()),negativeRealEigenvalues=int((e.real<0).sum()),maximumImaginaryEigenvaluePerSecond=float(abs(e.imag).max()))
    ev,vec=la.eigh(sym,subset_by_index=[0,0]);np.savetxt(prefix.with_suffix('.lowest-mode.csv'),np.column_stack((np.arange(n),vec[:,0]/np.sqrt(mass))),delimiter=',',header='component,velocity',comments='')
   report['cellSpectrum']=spec
   if name=='original264':
    native=json.loads(prefix.with_suffix('.json').read_text());report['cellSpectrum']['operatorScaleRelativeMinimum']=float(v[0]/np.max(np.abs(s)))
  out[name][method]=report;print(name,method,'vel',report['velocityCentroidRelativeL1'],'traction',metrics['steady']['wallTractionRelativeL1'],'spec',report.get('cellSpectrum',{}).get('lowestSymmetricPartEigenvaluesPerSecond',[])[:2],flush=True)
 (root/'summary.json').write_text(json.dumps(out,indent=2))
