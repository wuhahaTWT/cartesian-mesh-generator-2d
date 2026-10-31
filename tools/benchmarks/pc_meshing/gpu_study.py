"""Serial CPU/GPU cost experiments using the real native solver and frozen inputs."""
import argparse,json,os,subprocess,threading,time
from pathlib import Path
from study import ROOT,BASE,FLOW,run,dump,read
OUT=BASE/'gpu-v2'
TORCH=r'D:\gas-turbine-envs\windows-runtime\Lib\site-packages\torch'
def native(mesh,coupling,preconditioner,backend,steps,budget,minimum=10,capture=False):
 name=f'{mesh}-{coupling}-{preconditioner}-{backend}-n{steps}-min{minimum}'
 prefix=OUT/name;prefix.parent.mkdir(parents=True,exist_ok=True)
 cmd=[FLOW,'--mesh',BASE/'meshes'/f'{mesh}.solver.cm2d','--case','custom','--boundary',BASE/'meshes'/f'{mesh}.boundaries','--nu','.1','--speed','.5','--convection','upwind','--viscous-stress','symmetric','--coupling',coupling,'--linear-policy','adaptive' if coupling=='coupled' else 'strict','--pressure-preconditioner',preconditioner,'--pressure-backend',backend,'--max-iterations',str(steps),'--minimum-iterations',str(minimum),'--tolerance','1e-8' if coupling=='coupled' else '1e-6','--profile','--output',prefix]
 if capture:cmd+=['--pressure-snapshot',str(prefix)+'.pressure.bin']
 return run('gpu-study-v2-'+name,cmd,budget,kind='gpu-study',mesh=mesh,coupling=coupling,preconditioner=preconditioner,backend=backend,maximumIterations=steps,minimumIterations=minimum,prefix=str(prefix.relative_to(ROOT)))
def main():
 p=argparse.ArgumentParser();p.add_argument('phase',choices=['native','replay']);a=p.parse_args()
 os.environ['CARTMESH2D_TORCH_ROOT']=TORCH
 OUT.mkdir(parents=True,exist_ok=True)
 if a.phase=='native':
  for floor in (10,1):native('shape160','coupled','aggregation','cpu',100,150,floor)
  for steps in (1,2):native('gmsh65k','coupled','aggregation','cpu',steps,180,1)
  for pc,backend in [('ic0','cpu'),('ic0','torch-cuda'),('aggregation','cpu')]:native('gmsh65k','simple',pc,backend,5,180,capture=backend=='cpu' and pc=='ic0')
  for backend in ('cpu','torch-cuda'):native('shape160','simple','ic0',backend,100,180)
  native('shape160','coupled','ic0','torch-cuda',100,180,1)
 else:
  snapshot=OUT/'gmsh65k-simple-ic0-cpu-n5-min10.pressure.bin'
  run('gpu-study-v2-pressure-replay65k',[ROOT/'build/Release/cartmesh2d_pressure_replay.exe','--snapshot',snapshot,'--gpu','--repeats','3','--output',OUT/'replay65k.json'],240,kind='gpu-study-replay')
if __name__=='__main__':main()
