#!/usr/bin/env python3
"""Call native HMM benchmark on preserved Solver meshes; read results only."""
import argparse,hashlib,json,pathlib,subprocess,time,gzip
p=argparse.ArgumentParser();p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--mesh-root',type=pathlib.Path,default=pathlib.Path('outputs/curved-wall'));p.add_argument('--benchmark',type=pathlib.Path,default=pathlib.Path('build/cartmesh2d_hybrid_heat_benchmark'));a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
meshes={'original264':'ring-128-L4','phase-L4':'ring-128-L4-phase0.375','phase-L5':'ring-128-L5-phase0.375','phase-L6':'ring-128-L6-phase0.375'}
report={'scope':'new cloud isolated native HMM scalar heat; NOT coupled Euler','sourceDiffSha256':hashlib.sha256(subprocess.check_output(['git','diff','HEAD','--'])).hexdigest(),'baseCommit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'benchmarkSha256':hashlib.sha256(a.benchmark.read_bytes()).hexdigest(),'cases':{}}
for name,d in meshes.items():
 m=a.mesh_root/d/'mesh.solver.cm2d';root=a.output/name;root.mkdir(exist_ok=True)
 for f in ['ring.xy','mesh.solver.cm2d','mesh.cm2d','mesh.construction-quality.json','mesh.sizing.json','mesh.resolution.json']:
  source=m.parent/f
  if source.exists():(root/f).write_bytes(source.read_bytes())
 for mode in ['trace','isothermal']:
  for label,dt in [('dt5ms',.005),('dt2p5ms',.0025)]+([('dt1p25ms',.00125)] if name=='phase-L6' else []):
   prefix=root/(mode+'-'+label);cmd=[str(a.benchmark),str(m),str(prefix),mode,'.05',str(dt)];
   if name=='original264' and mode=='trace' and label=='dt5ms':cmd+=['cell-matrix']
   t=time.monotonic()
   with prefix.with_suffix('.log').open('w') as log:result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
   metadata={'command':cmd,'returncode':result.returncode,'fullProcessSeconds':time.monotonic()-t,'meshSha256':hashlib.sha256(m.read_bytes()).hexdigest()};prefix.with_suffix('.process.json').write_text(json.dumps(metadata,indent=2))
   if prefix.with_suffix('.json').exists():metadata['native']=json.loads(prefix.with_suffix('.json').read_text())
   report['cases'][name+'/'+prefix.name]=metadata;(a.output/'run.json').write_text(json.dumps(report,indent=2));print(name,prefix.name,result.returncode,round(metadata['fullProcessSeconds'],3),flush=True)
# History compression only after every native writer has exited. Hash readback.
for f in a.output.rglob('*.history.csv'):
 raw=f.read_bytes();z=f.with_suffix(f.suffix+'.gz')
 with gzip.open(z,'wb') as out:out.write(raw)
 assert gzip.open(z,'rb').read()==raw;f.unlink()
