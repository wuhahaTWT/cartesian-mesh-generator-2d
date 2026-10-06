#!/usr/bin/env python3
"""Prepare native inputs and call coupled Euler API harness; no independent PDE."""
import argparse,pathlib,subprocess,json,time,hashlib,shutil
p=argparse.ArgumentParser();p.add_argument('--common',action='store_true');p.add_argument('--restart-prefix',type=pathlib.Path);p.add_argument('--case',choices=['channel','cylinder'],required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--end',type=float,required=True);p.add_argument('--half',action='store_true');p.add_argument('--budget',type=float,default=240);p.add_argument('--schemes',nargs='+',default=['corrected','hybrid-heat','hybrid']);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
if a.case=='channel':mesh=pathlib.Path('outputs/recomputed-heated-smooth/mesh.solver.cm2d');boundary=pathlib.Path('outputs/recomputed-heated-smooth/full.boundaries');baseDt=5e-8;warmTime=2.5e-8;warmDt=1.25e-8;geom=mesh.parent/'channel.xy'
else:mesh=pathlib.Path('outputs/curved-cylinder/N48-W16/mesh.solver.cm2d');boundary=pathlib.Path('outputs/coupled-diffusion/cylinder-input/preparation.boundaries');baseDt=4e-9;warmTime=4e-9;warmDt=4e-9;geom=mesh.parent/'body.xy'
for f in [mesh,boundary,geom]:shutil.copyfile(f,a.output/f.name)
index={'codeCommit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'binarySha256':hashlib.sha256(pathlib.Path('build/cartmesh2d_euler_diffusion_benchmark').read_bytes()).hexdigest(),'case':a.case,'end':a.end,'half':a.half,'runs':[]}
if a.common:
 prefix=a.output/'comparison';cmd=['build/cartmesh2d_euler_diffusion_benchmark',str(a.output/mesh.name),str(a.output/boundary.name),str(prefix),'all',str(a.end),str(baseDt/(2 if a.half else 1)),str(a.budget)]
 if a.restart_prefix:cmd.append(str(a.restart_prefix))
 t=time.monotonic()
 with prefix.with_suffix('.log').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
 index['runs'].append({'command':cmd,'returncode':r.returncode,'fullProcessSeconds':time.monotonic()-t});(a.output/'run.json').write_text(json.dumps(index,indent=2));print(index,flush=True);raise SystemExit(r.returncode)
# Independent mode records stability diagnostics; actual accepted clocks may differ.
for scheme in a.schemes:
 checkpoint=None
 for phase,end,dt in [('warmup',warmTime,warmDt/(2 if a.half else 1)),('target',a.end,baseDt/(2 if a.half else 1))]:
  prefix=a.output/(scheme+'-'+phase);cmd=['build/cartmesh2d_euler_diffusion_benchmark',str(a.output/mesh.name),str(a.output/boundary.name),str(prefix),scheme,str(end),str(dt),str(a.budget)];
  if checkpoint:cmd.append(str(checkpoint))
  t=time.monotonic()
  with prefix.with_suffix('.log').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
  meta={'command':cmd,'returncode':r.returncode,'fullProcessSeconds':time.monotonic()-t}
  if prefix.with_suffix('.json').exists():meta['native']=json.loads(prefix.with_suffix('.json').read_text())
  prefix.with_suffix('.process.json').write_text(json.dumps(meta,indent=2));index['runs'].append(meta);(a.output/'run.json').write_text(json.dumps(index,indent=2));print(scheme,phase,meta,flush=True)
  if r.returncode:break
  checkpoint=prefix.with_suffix('.checkpoint')
