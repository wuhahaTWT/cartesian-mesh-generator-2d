#!/usr/bin/env python3
"""Native Oseen assembly, native RT1 convection and full nonlinear residual.
Python orchestrates and solves the exported sparse linear systems only.
No product acceptance rule: 1e-9 is an explicit dimensionless research iteration
control with U=L=1; pressure scale U^2, time scale L/U. It is tightened separately.
All iterates and failures retained; no automatic branch selection or relaxation.
"""
import sys,subprocess,json,time,hashlib,os
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
root=Path(__file__).resolve().parents[2]
exe=Path(os.environ.get('CARTMESH_P1_OSEEN_BINARY',str(root/'build/native-laminar-p1-transport')))
name,n,problem,nu,equation,boundary,label,*options=sys.argv[1:]
limit=int(options[0]) if options else 40
tol=float(options[1]) if len(options)>1 else 1e-9
folder=root/'outputs/cloud-laminar/p1-convection'/label
folder.mkdir(parents=True,exist_ok=False)
base=[name,n,problem,nu,equation,boundary,os.environ.get('CARTMESH_P1_QUADRATURE','6')];previous=os.environ.get('CARTMESH_P1_INITIAL_STATE','zero');history=[];start=time.monotonic()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
record={'args':base,'pid':os.getpid(),'maxIterations':limit,'initialState':previous,'researchIterationTolerance':tol,'nativeSourceSha256':{f:sha(root/'artifacts/current'/f) for f in ['native-laminar-p1-transport.cpp','native-laminar-p1-oseen.cpp','native-laminar-p1-stress.cpp','native-laminar-hybrid-stokes-p1.cpp']},'binaryPath':str(exe),'binarySha256':sha(exe),'iterations':history,'completed':False}
def save():
 record['totalSeconds']=time.monotonic()-start
 tmp=folder/'result.json.tmp'
 tmp.write_text(json.dumps(record,indent=2)+'\n')
 tmp.replace(folder/'result.json')
def native(mode,args,err):
 value=json.loads(subprocess.check_output([str(exe),mode,*args],stderr=err,text=True))
 Path(args[-1]+'.'+mode+'.json').write_text(json.dumps(value,indent=2)+'\n')
 return value
try:
 for i in range(limit):
  prefix=folder/f'iteration-{i:03d}';args=[*base,previous,str(prefix)]
  step={'iteration':i,'previous':previous};history.append(step);save()
  with open(str(prefix)+'.stderr','w') as err:
   step['assembly']=native('assemble',args,err);save()
   ts=time.monotonic();z=np.fromfile(str(prefix)+'.entries',dtype=[('i','<i8'),('j','<i8'),('v','<f8')]);size=step['assembly']['unknowns']
   A=coo_matrix((z['v'],(z['i'],z['j'])),shape=(size,size)).tocsc();A.sum_duplicates();b=np.fromfile(str(prefix)+'.rhs');lu=splu(A);x=lu.solve(b)
   for _ in range(2):x+=lu.solve(b-A@x)
   if not np.all(np.isfinite(x)):raise RuntimeError('nonfinite sparse solution')
   x.tofile(str(prefix)+'.solution');step['linear']={'seconds':time.monotonic()-ts,'nnz':A.nnz,'factorNnz':lu.L.nnz+lu.U.nnz,'residualMax':float(np.max(abs(A@x-b)))}
   del lu,A,z,x,b
   step['recovery']=native('recover',args,err);previous=str(prefix)+'.recover.state';save()
   step['nonlinearCheck']=native('check',[*base,previous,str(prefix)],err)
   step['hashes']={p.name:sha(p) for p in sorted(folder.glob(prefix.name+'.*')) if p.is_file() and p.suffix!='.stderr'}
   check=step['nonlinearCheck'];step['equationResidualMax']=max(check['internalResidual'],check['freeFaceResidual'],check['maxDivergence'])
   Path(str(prefix)+'.complete.json').write_text(json.dumps(step,indent=2)+'\n')
   save();print(json.dumps({'label':label,'iteration':i,'linearSeconds':step['linear']['seconds'],'change':step['recovery']['stateCoefficientChange'],'residual':step['equationResidualMax'],'maxSpeed':check['maxSpeed'],'pressureRange':check['pressureRange'],'wallFx':check['wallFx'],'energyError':check['energyBalance'],'convectionIdentityError':check['convectionIdentityError']}),flush=True)
   if step['equationResidualMax']<=tol and (equation=='stokes' or step['recovery']['stateCoefficientChange']<=tol):
    record['completed']=True;record['finalState']=previous;save();(folder/'completed.json').write_text(json.dumps(record,indent=2)+'\n');break
 if not record['completed']:record['failure']='iteration budget exhausted; no final acceptance';save();sys.exit(2)
except Exception as exc:
 record['failure']=str(exc);save();raise
