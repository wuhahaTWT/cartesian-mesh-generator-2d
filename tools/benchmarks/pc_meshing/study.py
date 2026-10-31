"""Local experiment orchestration and resource accounting; no FV equations.

All runs are immutable and resumable. Product defaults and solver kernels stay
unchanged. A timeout is retained as a timeout, never counted as convergence.
"""
from pathlib import Path
import argparse, datetime, hashlib, json, os, platform, subprocess, sys, time
import psutil

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'outputs/pc-mesh-study'
GEOMETRY=ROOT/'outputs/shape-replay/outputs/annulus-cross/original-physical-walls.xy'
ADAPTER=ROOT/'build/Release/cartmesh2d_annulus_mesh_cross.exe'
FLOW=ROOT/'build/Release/cartmesh2d_flow_cli.exe'
IMAGE='cartmesh2d-meshing-study:local'

def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temp.replace(path)

def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def run(name,command,budget,**metadata):
    receipt=BASE/'runs'/f'{name}.json';log=receipt.with_suffix('.log')
    if receipt.exists() and read(receipt).get('status')!='running':return read(receipt)
    if receipt.exists():raise RuntimeError(f'incomplete receipt {receipt}; inspect process before resuming')
    receipt.parent.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();record=dict(name=name,status='running',command=list(map(str,command)),
        startedUtc=datetime.datetime.now(datetime.timezone.utc).isoformat(),budgetSeconds=budget,**metadata)
    dump(receipt,record);print('START',name,flush=True)
    rss=0;cpus={};max_processes=0;gpu_samples=[];last_gpu=0
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONUTF8='1')
    with log.open('w',encoding='utf-8') as out:
        process=psutil.Popen(record['command'],cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,env=env)
        timeout=False
        while process.poll() is None:
            try:children=[process]+process.children(recursive=True)
            except psutil.Error:children=[process]
            current=0;max_processes=max(max_processes,len(children))
            for child in children:
                try:
                    c=child.cpu_times();cpus[child.pid]=max(cpus.get(child.pid,0),c.user+c.system)
                    current+=child.memory_info().rss
                except psutil.Error:pass
            rss=max(rss,current)
            if time.monotonic()-started>budget:
                timeout=True
                for child in reversed(children):
                    try:child.kill()
                    except psutil.Error:pass
                break
            time.sleep(.2)
        code=process.wait()
    record.update(status='timeout' if timeout or code==124 else 'finished',exitCode=code,
        wallSeconds=time.monotonic()-started,processTreeCpuSeconds=sum(cpus.values()),
        peakProcessTreeRssMiB=rss/1024**2,maxProcesses=max_processes,
        costScope='Windows process tree sampled every 0.2 s; Docker CLI excludes container CPU/RSS, see inner resource.json; no energy meter')
    dump(receipt,record);print('END',name,record['status'],code,round(record['wallSeconds'],3),flush=True)
    return record

def docker(arguments,cpus=4,memory='8g'):
    return ['docker','run','--rm','--network','none','--cpus',str(cpus),'--memory',memory,
        '-v',str(ROOT)+':/work','-w','/work',IMAGE,*arguments]

def generate():
    BASE.mkdir(parents=True,exist_ok=True)
    dump(BASE/'environment.json',dict(sourceCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        sourceDirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)),
        geometrySha256=sha(GEOMETRY),adapterSha256=sha(ADAPTER),flowSha256=sha(FLOW),
        os=platform.platform(),python=sys.version,logicalCpus=psutil.cpu_count(),physicalCpus=psutil.cpu_count(logical=False),
        ramGiB=psutil.virtual_memory().total/1024**3,
        physicalModel='Frozen 256-segment polygon annulus, facet tangential speed 0.5, nu=0.1, rho=1, rest start',
        energy='Not measured; process CPU-seconds and RSS are reported, no invented FLOPS or kWh'))
    specs=[]
    for n in (128,160,192,256,384,512,768):
        specs.append((f'shape{n}','shape',n,[]))
        specs.append((f'graded{n}','shape',n,['--wall-refine','1']))
    for n in (160,256,512):specs.append((f'legacy{n}','native',n,[]))
    for x,y in ((.25,0),(.25,.25),(.5,.5)):
        specs.append((f'phase128-{x}-{y}','shape',128,['--phase',str(x),str(y)]))
    for n in (160,256):specs.append((f'cap2-{n}','shape',n,['--area-ratio','2']))
    dump(BASE/'mesh-plan.json',specs)
    for name,mode,n,options in specs:
        prefix=BASE/'meshes'/name;prefix.parent.mkdir(exist_ok=True)
        run('mesh-'+name,[ADAPTER,mode,GEOMETRY,prefix,n,*options],180,
            kind='mesh',method=mode,n=n,options=options,prefix=str(prefix.relative_to(ROOT)))
    generate_gmsh()

def generate_gmsh():
    for name,mode,radial,sub,h in [('gmsh16k','mapped',32,2,.02),('gmsh65k','mapped',64,4,.01),
            ('gmsh196k','mapped',96,8,.006),('gmsh-tri-coarse','triangles',0,0,.035),('gmsh-tri-fine','triangles',0,0,.018)]:
        prefix=BASE/'meshes'/name
        args=['python3','tools/benchmarks/pc_meshing/gmsh_reference.py','--geometry','/work/'+str(GEOMETRY.relative_to(ROOT)).replace('\\','/'),
            '--output','/work/'+str(prefix.relative_to(ROOT)).replace('\\','/'),'--mode',mode,'--h',str(h)]
        if mode=='mapped':args+=['--radial',str(radial),'--subdivisions',str(sub)]
        r=run('gmsh-'+name,docker(args),180,kind='gmsh',prefix=str(prefix.relative_to(ROOT)))
        if r['exitCode']==0:
            run('mesh-'+name,[ADAPTER,'ready',GEOMETRY,prefix,256,str(prefix)+'.polygons'],180,
                kind='mesh',method='gmsh-'+mode,n=None,prefix=str(prefix.relative_to(ROOT)))

def generate_indexed():
    adapter=ADAPTER
    for n in (128,256,512,768):
        name=f'graded-indexed{n}';prefix=BASE/'meshes'/name
        run('mesh-'+name,[adapter,'shape',GEOMETRY,prefix,n,'--wall-refine','1','--distance-index','1'],180,
            kind='mesh',method='shape-refine-indexed',n=n,prefix=str(prefix.relative_to(ROOT)),adapterSha256=sha(adapter))

def solve_native(names,coupling='simple',policy='strict',tolerance='1e-6',budget=900,preconditioner='ic0',tag='',boundary_model='facet-wall',iterations=6000):
    for name in names:
        prefix=BASE/'meshes'/name;mesh=Path(str(prefix)+'.solver.cm2d')
        if not mesh.exists():print('SKIP rejected mesh',name,flush=True);continue
        boundary=Path(str(prefix)+'.boundaries')
        if not boundary.exists():
            r=run('export-boundary-'+name,[FLOW,'--mesh',mesh,'--case','annulus','--speed','.5','--export-boundaries',boundary],120,kind='input',mesh=name)
            if r['exitCode']!=0:raise RuntimeError('native annulus boundary export failed')
        if boundary_model=='analytic-trace':
            boundary=Path(str(prefix)+'.trace.boundaries')
            if not boundary.exists():run('analytic-trace-'+name,[sys.executable,ROOT/'tools/benchmarks/annulus_exact_trace.py',
                '--mesh',mesh,'--output',boundary],120,kind='input',mesh=name)
        suffix=('-'+tag if tag else '')+('-analytic-trace' if boundary_model=='analytic-trace' else '')
        output=BASE/'solutions'/f'{name}-{coupling}-{policy}-{preconditioner}-{tolerance}{suffix}'
        output.parent.mkdir(exist_ok=True)
        run('native-'+output.name,[FLOW,'--mesh',mesh,'--case','custom','--boundary',boundary,
            '--nu','.1','--speed','.5','--convection','upwind','--viscous-stress','symmetric','--coupling',coupling,
            '--linear-policy',policy,'--pressure-preconditioner',preconditioner,'--max-iterations',iterations,
            '--tolerance',tolerance,'--profile','--output',output],budget,kind='native',mesh=name,
            coupling=coupling,linearPolicy=policy,tolerance=float(tolerance),preconditioner=preconditioner,
            prefix=str(output.relative_to(ROOT)),threads=1,gpuUsed=False,boundaryModel=boundary_model,maximumIterations=iterations)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=['generate','gmsh','indexed','native'])
    p.add_argument('--names',nargs='+',default=['shape160','graded160','gmsh16k'])
    p.add_argument('--coupling',default='simple');p.add_argument('--linear-policy',default='strict')
    p.add_argument('--tolerance',default='1e-6');p.add_argument('--budget',type=int,default=900)
    p.add_argument('--preconditioner',default='ic0');p.add_argument('--tag',default='')
    p.add_argument('--boundary-model',choices=['facet-wall','analytic-trace'],default='facet-wall')
    p.add_argument('--iterations',type=int,default=6000);a=p.parse_args()
    if a.phase=='generate':generate()
    elif a.phase=='gmsh':generate_gmsh()
    elif a.phase=='indexed':generate_indexed()
    else:solve_native(a.names,a.coupling,a.linear_policy,a.tolerance,a.budget,a.preconditioner,a.tag,a.boundary_model,a.iterations)
