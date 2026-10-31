"""Fluent staged timings and CPU/native-GPU comparison; no silent fallback.

2D native GPU is unsupported in this PyFluent/Fluent installation. The optional
one-layer 3D extrusion is a separate reference with symmetric end planes.
"""
from pathlib import Path
import argparse,json,math,re,time,sys
from study import ROOT,BASE,read,dump
sys.path.insert(0,str(ROOT/'tools/flow'))
import native_mesh as nm
OUT=BASE/'gpu-fluent'
EXE=r'E:\ANSYS Inc\v261\fluent\ntbin\win64\fluent.exe'
def extrude(name):
 out=OUT/'input3d';out.mkdir(parents=True,exist_ok=True);dest=out/(name+'.msh')
 if dest.exists():return dest
 m=nm.read_cm2d(BASE/'meshes'/f'{name}.solver.cm2d');n=len(m.vertices)
 groups={3:[],1001:[],1002:[]};names={3:'interior',1001:'front',1002:'back'}
 old=read(BASE/'fluent'/f'{name}-2cores-1e-8-coupled/input.interchange.json')
 keys={(round(u,10),round(v,10)):name for name,(u,v) in old['wallVelocities'].items() if name!='housing'}
 ids={name:int(name.split('_')[-1]) for name in keys.values()};names.update({v:k for k,v in ids.items()});names[1000]='housing'
 for e in m.edges:
  a,b=e.v0,e.v1
  # Native directed edge keeps owner on its left. 3D Fluent face normal
  # points INTO c0, so the side loop must point to that left-hand side.
  zone=3
  if e.neighbour<0:
   pa,pb=m.vertices[a],m.vertices[b];x=(pa[0]+pb[0])/2;y=(pa[1]+pb[1])/2
   if math.hypot(x,y)<.75:
    dx=pb[0]-pa[0];dy=pb[1]-pa[1];l=math.hypot(dx,dy)
    # owner-outward = (dy,-dx); moving inner velocity=(ny,-nx)*0.5
    u,v=-.5*dx/l,-.5*dy/l
    zone=ids[keys[(round(u,10),round(v,10))]]
   else:zone=1000
  groups.setdefault(zone,[]).append(([a+1,a+n+1,b+n+1,b+1],e.owner+1,e.neighbour+1))
 for c in m.cells:
  vs=list(c.vertices)
  area=sum(m.vertices[a][0]*m.vertices[b][1]-m.vertices[b][0]*m.vertices[a][1] for a,b in zip(vs,vs[1:]+vs[:1]))
  if area<0:vs.reverse()
  groups[1001].append(([v+1 for v in vs],c.id+1,0))
  groups[1002].append(([v+n+1 for v in reversed(vs)],c.id+1,0))
 H=lambda x:format(x,'x')
 with dest.open('x',encoding='ascii') as f:
  f.write(f'(0 "One layer symmetric extrusion; reference only")\n(2 3)\n(10 (0 1 {H(2*n)} 0))\n(10 (1 1 {H(2*n)} 1 3)(\n')
  for z in (0.,.1):
   for x,y in m.vertices:f.write(f'{x:.17g} {y:.17g} {z:.17g}\n')
  f.write(f'))\n(12 (0 1 {H(len(m.cells))} 0))\n(12 (2 1 {H(len(m.cells))} 1 7))\n(13 (0 1 {H(sum(map(len,groups.values())))} 0))\n')
  start=1
  for zone,rows in sorted(groups.items()):
   f.write(f'(13 ({H(zone)} {H(start)} {H(start+len(rows)-1)} {2 if zone==3 else 7 if zone in (1001,1002) else 3} 0)(\n')
   for nodes,c0,c1 in rows:f.write(' '.join(map(H,[len(nodes),*nodes,c0,c1]))+'\n')
   f.write('))\n');start+=len(rows)
  f.write('(45 (2 fluid fluid)())\n')
  for zone in sorted(groups):f.write(f'(45 ({zone} {"interior" if zone==3 else "symmetry" if zone in (1001,1002) else "wall"} {names[zone]})())\n')
 return dest

def execute(a):
 import ansys.fluent.core as f
 out=OUT/a.name;out.mkdir(parents=True,exist_ok=False);r={'name':a.name,'gpuRequested':a.gpu,'dimension':a.dimension,'cpuProcesses':a.cores,'precision':a.precision,'tolerance':1e-8,'stages':{},'scope':'same mesh, rest initialization; CPU/GPU native Fluent methods may differ'}
 start=time.monotonic();s=None
 def save():r['elapsedSeconds']=time.monotonic()-start;dump(out/'result.json',r)
 def mark(k,t):r['stages'][k]=time.monotonic()-t;save()
 try:
  t=time.monotonic()
  s=f.launch_fluent(fluent_path=EXE,dimension=a.dimension,precision=a.precision,gpu=[0] if a.gpu else False,processor_count=a.cores,ui_mode='no_gui',start_timeout=100,cwd=str(out),start_transcript=False,start_watchdog=True,env={'REMOTING_SERVER_ADDRESS':'127.0.0.1'})
  s.transcript.start(file_name=str(out/'transcript.log'));mark('launch',t)
  r['version']=str(s.get_fluent_version());x=s.settings
  t=time.monotonic()
  if a.prepare:
   mesh=extrude('gmsh65k');x.file.read_mesh(file_name=str(mesh));s.tui.mesh.check()
   text=(out/'transcript.log').read_text(encoding='utf-8',errors='replace')
   mins=re.findall(r'minimum volume.*?:\s*([-+0-9.eE]+)',text)
   assert mins and float(mins[-1])>0 and 'Mesh check failed' not in text,'Mesh check failed'
   r['minimumVolume']=float(mins[-1]);x.setup.models.viscous.model='laminar'
   material=x.setup.materials.fluid['air'];material.density.value=1.;material.viscosity.value=.1
   walls=read(BASE/'fluent/gmsh65k-2cores-1e-8-coupled/result.json')['configuredState']['walls']
   state={}
   for name,item in walls.items():
    if name=='housing':continue
    mom=item['momentum'].copy();mom['velocity_components']=[*mom['velocity_components'],{'option':'value','value':0.}]
    state[name]={'momentum':mom}
   x.setup.boundary_conditions.wall.set_state(state)
   x.solution.methods.p_v_coupling.flow_scheme='Coupled';x.solution.methods.discretization_scheme={'pressure':'standard','mom':'first-order-upwind'}
   for name in ('continuity','x-velocity','y-velocity','z-velocity'):x.solution.monitor.residual.equations[name].absolute_criteria=1e-8
  else:x.file.read_case(file_name=str(a.case.resolve()))
  mark('readAndSetup',t)
  r['methods']=x.solution.methods.get_state();r['controls']=x.solution.controls.get_state()
  r['runs']=[]
  for repeat in range(a.repeat):
   t=time.monotonic();x.solution.initialization.initialization_type='standard';x.solution.initialization.standard_initialize();mark('initialize'+str(repeat),t)
   t=time.monotonic();x.solution.run_calculation.iterate(iter_count=a.iterations);elapsed=time.monotonic()-t
   text=(out/'transcript.log').read_text(encoding='utf-8',errors='replace');hist=[]
   for row in re.findall(r'^\s*(\d+)\s+([0-9.eE+\-]+)\s+([0-9.eE+\-]+)\s+([0-9.eE+\-]+)',text,re.M):
    it,c,u,v=row;hist.append({'iteration':int(it),'continuity':float(c),'u':float(u),'v':float(v)})
   r['runs'].append({'solveSeconds':elapsed,'history':hist,'last':hist[-1] if hist else None,'converged':bool(hist and max(hist[-1][k] for k in ('continuity','u','v'))<=1e-8)})
   save()
  t=time.monotonic();x.file.write_case_data(file_name=str(out/'accepted.cas.h5'));mark('saveCaseData',t)
  r['status']='finished'
 except Exception as e:r.update(status='failed',error=repr(e));save();raise
 finally:
  if s:
   t=time.monotonic()
   try:s.exit(timeout=15,wait=15);r['sessionClosed']=True
   except Exception as e:r['exitError']=repr(e)
   mark('close',t)
  save()
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--name',required=True);p.add_argument('--case',type=Path);p.add_argument('--cores',type=int,default=1);p.add_argument('--gpu',action='store_true');p.add_argument('--dimension',type=int,default=2);p.add_argument('--precision',default='double');p.add_argument('--prepare',action='store_true');p.add_argument('--iterations',type=int,default=500);p.add_argument('--repeat',type=int,default=1)
 execute(p.parse_args())
