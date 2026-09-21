#!/usr/bin/env python3
"""The optional system factor must retain full steady/transient acceptance."""
import argparse,csv,json,subprocess,sys,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
from verify_pressure_openings import rectangle,pressure_boundaries
import verify_native_flow as native
import verify_transient_flow as transient
p=argparse.ArgumentParser();p.add_argument('--cli',required=True);a=p.parse_args();cli=str(Path(a.cli).resolve())
with tempfile.TemporaryDirectory(prefix='cm2d-cholesky-') as directory:
 root=Path(directory);mesh=root/'channel.solver.cm2d';rectangle(mesh,8);template=root/'template';bc=root/'boundaries'
 r=subprocess.run([cli,'--mesh',str(mesh),'--case','duct','--export-boundaries',str(template)],capture_output=True,text=True,timeout=20);assert r.returncode==0,r.stderr
 pressure_boundaries(template,bc)
 def run(method,label,extra=()):
  prefix=root/label;command=[cli,'--mesh',str(mesh),'--case','custom','--boundary',str(bc),'--output',str(prefix),'--nu','.1','--speed','1','--tolerance','1e-10','--max-iterations','5000','--pressure-preconditioner',method,'--convection','face-limited-linear','--linear-policy','adaptive','--profile']
  return prefix,subprocess.run(command+list(extra),capture_output=True,text=True,timeout=30)
 if sys.platform!='darwin':
  prefix,r=run('cholesky','unsupported');assert r.returncode!=0 and 'only on macOS' in r.stderr;assert not prefix.with_suffix('.json').exists()
 else:
  for mode,extra in [('steady',()),('transient',('--time-step','.02','--steps','2'))]:
   prefixes=[]
   for method in ['aggregation','cholesky']:
    prefix,r=run(method,mode+'-'+method,extra);assert r.returncode==0,r.stderr;prefixes.append(prefix)
    summary=json.loads(prefix.with_suffix('.json').read_text());assert summary['converged'] and summary['strictLinearFinal'] and summary['pressurePreconditioner']==method
    if mode=='steady':audit=native.verify_case(mesh,prefix,'custom',.1,1,native.argument_parser().parse_args(['--max-iterations','5000']))
    else:
     assert summary['strictAcceptedSteps']==2
     audit=transient.verify(mesh,prefix,root/(mode+'-'+method+'-audit.json'))
    assert audit['valid'],audit
    profile=json.loads(prefix.with_suffix('.performance.json').read_text())
    if method=='cholesky':
     assert profile['pressureCholeskyBuilds']>0 and profile['pressureHierarchyBuilds']==0 and profile['maxPressureIterations']<=2
     repeated,rr=run(method,mode+'-repeated',extra);assert rr.returncode==0,rr.stderr
     for suffix in ['.cells.csv','.faces.csv','.fields.json','.residuals.csv','.vtk']:
      assert prefix.with_suffix(suffix).read_bytes()==repeated.with_suffix(suffix).read_bytes(),suffix
   x,y=[list(csv.DictReader(prefix.with_suffix('.cells.csv').open())) for prefix in prefixes]
   assert max(abs(float(i[k])-float(j[k])) for i,j in zip(x,y) for k in ['u','v','p'])<2e-7
  prefix,r=run('cholesky','partial',('--max-iterations','1'));assert r.returncode==2,r.stderr
  assert not json.loads(prefix.with_suffix('.json').read_text())['converged']
print('System pressure backend retains true equations, strict final gates, physical time and explicit platform limits')
