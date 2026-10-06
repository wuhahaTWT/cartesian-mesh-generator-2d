#!/usr/bin/env python3
"""Native P1 symmetric-stress assembly/recovery; Python only sparse linear algebra."""
import sys,subprocess,json,time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
root=Path(__file__).resolve().parents[2];bin=root/'build/native-laminar-p1-stress'
name,n,problem,lam,load,stress,label=sys.argv[1:];prefix=root/'outputs/cloud-laminar/p1-stress'/label
prefix.parent.mkdir(parents=True,exist_ok=True)
if any(prefix.parent.glob(prefix.name+'.*')):raise SystemExit('output prefix already exists; choose a new label to preserve raw evidence')
args=[name,n,problem,lam,load,stress,'6',str(prefix)]
t=time.monotonic()
with open(str(prefix)+'.stderr','w') as err:
 a=json.loads(subprocess.check_output([str(bin),'assemble',*args],stderr=err,text=True));Path(str(prefix)+'.assemble.json').write_text(json.dumps(a,indent=2)+'\n')
 print(label,'assembled',a['unknowns'],flush=True)
 ts=time.monotonic();z=np.fromfile(str(prefix)+'.entries',dtype=[('i','<i8'),('j','<i8'),('v','<f8')]);A=coo_matrix((z['v'],(z['i'],z['j'])),shape=(a['unknowns'],)*2).tocsc();A.sum_duplicates();b=np.fromfile(str(prefix)+'.rhs');lu=splu(A);x=lu.solve(b)
 for _ in range(2):x+=lu.solve(b-A@x)
 assert np.all(np.isfinite(x));x.tofile(str(prefix)+'.solution');linear={'seconds':time.monotonic()-ts,'nnz':A.nnz,'factorNnz':lu.L.nnz+lu.U.nnz,'residualMax':float(np.max(abs(A@x-b)))}
 r=json.loads(subprocess.check_output([str(bin),'recover',*args],stderr=err,text=True));record={'args':args,'assembly':a,'linear':linear,'recovery':r,'totalSeconds':time.monotonic()-t};Path(str(prefix)+'.result.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)
