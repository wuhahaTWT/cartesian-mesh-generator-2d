#!/usr/bin/env python3
"""Call native viscous benchmark on retained real meshes; no PDE reconstruction."""
import argparse, pathlib, subprocess, json, hashlib, time, gzip
p=argparse.ArgumentParser();p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--mesh-root',type=pathlib.Path,default=pathlib.Path('outputs/curved-wall'));p.add_argument('--benchmark',type=pathlib.Path,default=pathlib.Path('build/cartmesh2d_hybrid_viscous_benchmark'));a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
cases={'original264':'ring-128-L4','phase-L4':'ring-128-L4-phase0.375','phase-L5':'ring-128-L5-phase0.375','phase-L6':'ring-128-L6-phase0.375'}
index={'commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'benchmarkSha256':hashlib.sha256(a.benchmark.read_bytes()).hexdigest(),'cases':{}}
for name,d in cases.items():
 root=a.output/name;root.mkdir(exist_ok=True);source=a.mesh_root/d
 for f in ['ring.xy','mesh.solver.cm2d','mesh.cm2d','mesh.construction-quality.json','mesh.sizing.json','mesh.resolution.json']:
  if (source/f).exists():(root/f).write_bytes((source/f).read_bytes())
 for method in ['hybrid','linear']:
  prefix=root/method;cmd=[str(a.benchmark),str(root/'mesh.solver.cm2d'),str(prefix),method];t=time.monotonic()
  with prefix.with_suffix('.log').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
  meta={'command':cmd,'returncode':r.returncode,'fullProcessSeconds':time.monotonic()-t};prefix.with_suffix('.process.json').write_text(json.dumps(meta,indent=2));index['cases'][name+'/'+method]=meta
  print(name,method,meta,flush=True)
  f=prefix.with_suffix('.cell-matrix.csv')
  if f.exists():
   raw=f.read_bytes();z=f.with_suffix(f.suffix+'.gz')
   with gzip.open(z,'wb') as out:out.write(raw)
   assert gzip.open(z,'rb').read()==raw;f.unlink()
  (a.output/'run.json').write_text(json.dumps(index,indent=2))
