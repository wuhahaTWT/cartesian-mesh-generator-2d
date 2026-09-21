#!/usr/bin/env python3
"""Coarse initial iterates cannot bypass fine-grid equations or target binding."""
import argparse,csv,json,subprocess,sys,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
from verify_pressure_openings import rectangle,pressure_boundaries
from prolongate_regular_flow import interpolate,centres
import verify_native_flow as native
p=argparse.ArgumentParser();p.add_argument('--cli',required=True);args=p.parse_args();cli=str(Path(args.cli).resolve())
with tempfile.TemporaryDirectory(prefix='cm2d-initial-guess-') as directory:
 root=Path(directory);coarse=root/'coarse.solver.cm2d';fine=root/'fine.solver.cm2d';rectangle(coarse,4);rectangle(fine,8)
 meshes={x:native.read_cm2d(x) for x in [coarse,fine]}
 source_points=centres(meshes[coarse]);affine=[{'cell':str(i),'x':str(x),'y':str(y),'u':str(1+2*x-3*y),'v':str(-2+x+y),'p':str(5-4*x+2*y)} for i,(x,y) in enumerate(source_points)]
 points,values=interpolate(meshes[coarse],affine,meshes[fine])
 for (x,y),v in zip(points,values):assert max(abs(a-b) for a,b in zip(v,[1+2*x-3*y,-2+x+y,5-4*x+2*y]))<2e-14
 def run(mesh,label,guess=None,extra=()):
  bc=root/(mesh.stem+'-boundary');template=root/(mesh.stem+'-template')
  if not bc.exists():
   r=subprocess.run([cli,'--mesh',str(mesh),'--case','duct','--export-boundaries',str(template)],capture_output=True,text=True,timeout=20);assert r.returncode==0,r.stderr
   pressure_boundaries(template,bc)
  prefix=root/label
  command=[cli,'--mesh',str(mesh),'--case','custom','--boundary',str(bc),'--output',str(prefix),'--nu','.1','--speed','1','--tolerance','1e-10','--max-iterations','5000','--pressure-preconditioner','aggregation','--steady-acceleration','anderson','--linear-policy','adaptive']
  if guess:command+=['--initial-guess',str(guess)]
  r=subprocess.run(command+list(extra),capture_output=True,text=True,timeout=30)
  return prefix,r
 coarse_prefix,r=run(coarse,'coarse');assert r.returncode==0,r.stderr
 cold,r=run(fine,'cold');assert r.returncode==0,r.stderr
 rows=list(csv.DictReader(Path(str(coarse_prefix)+'.cells.csv').open()));points,values=interpolate(meshes[coarse],rows,meshes[fine])
 records=[[i,*point,*value] for i,(point,value) in enumerate(zip(points,values))]
 def save(path,rows):
  with path.open('w') as f:
   w=csv.writer(f);w.writerow(['cell','x','y','u','v','p']);w.writerows(rows)
 guess=root/'guess.csv';save(guess,records);warm,r=run(fine,'warm',guess);assert r.returncode==0,r.stderr
 for prefix in [cold,warm]:
  summary=json.loads(Path(str(prefix)+'.json').read_text());assert summary['converged'] and summary['strictLinearFinal']
  assert max(summary[k] for k in ['momentumResidual','velocityChange','pressureChange'])<1e-10
  audit=native.verify_case(fine,prefix,'custom',.1,1,native.argument_parser().parse_args(['--max-iterations','5000']));assert audit['valid'],audit['issues']
 a=list(csv.DictReader(Path(str(cold)+'.cells.csv').open()));b=list(csv.DictReader(Path(str(warm)+'.cells.csv').open()))
 assert max(abs(float(x[k])-float(y[k])) for x,y in zip(a,b) for k in ['u','v','p'])<2e-7
 for index,mutate in enumerate([lambda r:r.pop(),lambda r:r[0].__setitem__(1,99),lambda r:r[0].__setitem__(0,1),lambda r:r[0].__setitem__(3,'nan')]):
  bad=[row[:] for row in records];mutate(bad);file=root/f'bad-{index}.csv';save(file,bad);prefix,r=run(fine,f'rejected-{index}',file);assert r.returncode!=0;assert not Path(str(prefix)+'.json').exists()
 _,r=run(fine,'transient-rejected',guess,['--time-step','.01','--steps','1']);assert r.returncode!=0 and 'actual steady' in r.stderr
print('Initial iterate mapping, fine-grid equations, strict gates, target binding and transient refusal pass')
