#!/usr/bin/env python3
"""Reproducible native CLI cases. Python only configures and reads output data;
no independent equation or topology reconstruction/audit chain.
"""
import argparse, csv, json, math, platform, subprocess, time, hashlib, shutil, tempfile
from pathlib import Path

def run(cmd, log):
    log.parent.mkdir(parents=True,exist_ok=True)
    log.with_suffix('.command.json').write_text(json.dumps([str(v) for v in cmd],indent=2))
    begin=time.monotonic()
    with log.open('w') as out:
        p=subprocess.run([str(v) for v in cmd],stdout=out,stderr=subprocess.STDOUT)
    return {'code':p.returncode,'seconds':time.monotonic()-begin,'log':str(log)}

def mesh_data(path):
    lines=path.read_text().splitlines();vertices={};edges=[];cells=0;orientation={}
    for i,line in enumerate(lines):
        if line.startswith('VERTICES '):
            for row in lines[i+1:i+1+int(line.split()[1])]:
                f=row.split();vertices[int(f[0])]=(float(f[1]),float(f[2]))
        if line.startswith('EDGES '):
            edges=[list(map(int,row.split())) for row in lines[i+1:i+1+int(line.split()[1])]]
        if line.startswith('CELLS '):
            cells=int(line.split()[1])
            for row in lines[i+1:i+1+cells]:
                f=row.split();owner=int(f[0]);n=int(f[4]);ids=list(map(int,f[5:5+n]))
                for a,b in zip(ids,ids[1:]+ids[:1]):orientation[owner,min(a,b),max(a,b)]=(a,b)
    return vertices,edges,cells,orientation

def configure(mesh,kind,root,events=True,speed=None):
    vertices,edges,cells,orientation=mesh_data(mesh)
    xmin=min(x for x,y in vertices.values());xmax=max(x for x,y in vertices.values())
    ymin=min(y for x,y in vertices.values());ymax=max(y for x,y in vertices.values())
    records=[];thermal=['face,type,value,inflowValue'];names={}
    for id,a,b,owner,neighbour,patch in edges:
        if neighbour>=0:continue
        a,b=orientation[owner,min(a,b),max(a,b)]
        ax,ay=vertices[a];bx,by=vertices[b];x=.5*(ax+bx);y=.5*(ay+by);sx=by-ay;sy=ax-bx
        left=abs(x-xmin)<1e-10;right=abs(x-xmax)<1e-10
        top=abs(y-ymax)<1e-10;bottom=abs(y-ymin)<1e-10
        name='wall';fk='wall';u=v=0;tk='flux';tv=0;ti=300
        if kind=='cavity':
            if top:name='lid';fk='moving-wall';u=float(speed) if speed else .2;tk='value';tv=301
            elif bottom:name='bottom';tk='value';tv=300
        else:
            if left:name='inlet';fk='velocity-inlet';u=float(speed) if speed else (1 if kind=='channel' else .2);tk='value';tv=300
            elif right:name='outlet';fk='pressure-opening'
            elif kind=='channel' or patch==1:tk='value';tv=301
            elif top or bottom:fk='symmetry';name='farfield'
        names[id]=name
        records.append(f'BOUNDARY {id} {owner} {x:.17g} {y:.17g} {sx:.17g} {sy:.17g} {fk} "{name}" {u} {v} 0')
        thermal.append(f'{id},{tk},{tv},{ti}')
    flow=root/'flow.boundaries';flow.write_text(f'CARTMESH2D_FLOW_BOUNDARIES 1\nCOUNTS {cells} {len(edges)} {len(records)}\n'+'\n'.join(records)+'\nEND\n')
    bc=root/'thermal.csv';bc.write_text('\n'.join(thermal)+'\n')
    event=root/'events.csv';changes={'channel':'1.37,inlet,value,300.5,300.5\n3.23,inlet,value,300,300\n',
      'cavity':'1.37,lid,value,300,300\n2.43,source,source,0.1,\n3.23,source,source,0,\n',
      'cylinder':'7.13,wall,value,300,300\n14.23,wall,value,301,300\n'}
    event.write_text('time,target,type,value,inflowValue\n'+changes[kind])
    (root/'face-names.json').write_text(json.dumps(names))
    return flow,bc,event,cells

def execute(args):
    root=Path(args.output).resolve();root.mkdir(parents=True,exist_ok=True)
    geom=root/'geometry.xy'
    if args.case=='cylinder':
        geom.write_text(''.join(f'{.15*math.cos(2*math.pi*i/64):.17g} {.15*math.sin(2*math.pi*i/64):.17g}\n' for i in range(64)))
        region='exterior';pad='3'
    else:geom.write_text('0 0\n'+('2 0\n2 1\n' if args.case=='channel' else '1 0\n1 1\n')+'0 1\n');region='interior';pad=str(1/14)
    mesh=root/'mesh.solver.cm2d'
    if not mesh.exists():
        r=run([args.mesh_cli,geom,root/'mesh',args.level,pad,'.1',region,root/'openfoam',args.level,'0'],root/'mesh.log')
        if r['code']:raise RuntimeError(r)
    flow,bc,event,cells=configure(mesh,args.case,root,speed=args.speed)
    prefix=root/args.label
    frozen=[]
    for input_path,suffix in [(flow,'.input.flow.boundaries'),(bc,'.input.thermal.csv'),(event,'.input.events.csv')]:
        target=Path(str(prefix)+suffix);shutil.copyfile(input_path,target);frozen.append(target)
    flow,bc,event=frozen
    cmd=[args.transport_cli,'--mesh',mesh,'--output',prefix,'--evolve-flow','custom','--flow-boundary',flow,
         '--boundary',bc,'--flow-nu',args.nu,'--flow-speed',args.speed or ('1' if args.case=='channel' else '.2'),'--flow-velocity-relaxation',args.relaxation,
         '--diffusivity',args.diffusivity,'--initial','300','--dt',args.dt,'--end-time',args.end,
         '--flow-max-iterations','1500','--flow-convection','limited-linear','--convection','limited-linear','--flux-correction',args.flux_correction,
         '--min-dt','0.000001','--max-courant',args.courant,'--max-step-retries',args.retries,'--time-error','on' if args.error else 'off',
         '--temperature-scale','1','--velocity-scale',args.speed or ('1' if args.case=='channel' else '.2'),'--time-rtol',args.rtol]
    if not args.no_events:cmd+=['--thermal-events',event]
    if args.restart:cmd+=['--restart',args.restart]
    provenance={'git':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        'sourceDiffSha256':hashlib.sha256(subprocess.check_output(['git','diff','--binary'])).hexdigest(),
        'binarySha256':hashlib.sha256(Path(args.transport_cli).read_bytes()).hexdigest()}
    r=run(cmd,root/(args.label+'.log'))
    meta={'case':args.case,'cells':cells,'run':r,'host':platform.uname()._asdict(),
          'provenance':provenance}
    if prefix.with_suffix('.json').exists():meta['summary']=json.loads(prefix.with_suffix('.json').read_text())
    if r['code']==0:
        field=list(csv.DictReader(prefix.with_suffix('.cells.csv').open()));faces=list(csv.DictReader(prefix.with_suffix('.faces.csv').open()))
        histories={suffix:list(csv.DictReader(Path(str(prefix)+suffix).open())) for suffix in
            ['.thermal-history.csv','.heat-history.csv','.attempt-history.csv']}
        accepted=histories['.thermal-history.csv'];heat=histories['.heat-history.csv'];attempts=histories['.attempt-history.csv']
        meta['outputIntegrity']=bool(len(field)==cells and len(accepted)==len(heat)==meta['summary']['completedSteps']
            and sum(row['reason']=='accepted' for row in attempts)==len(accepted)
            and accepted and float(accepted[-1]['time'])==float(heat[-1]['time'])==meta['summary']['acceptedTime'])
        if not meta['outputIntegrity']:
            meta['run']['code']=1;meta['run']['error']='Native output files are incomplete or inconsistent; retained for diagnosis'
        names=json.loads((root/'face-names.json').read_text());groups={}
        for f in faces:
            name=names.get(f['face']);
            if name is None:continue
            g=groups.setdefault(name,{'advectiveFlux':0.,'diffusiveFlux':0.,'volumeFlux':0.})
            for key in g:g[key]+=float(f[key])
        meta['boundaryGroups']=groups
        meta['heatContent']=sum(float(c['area'])*float(c['value']) for c in field)
        if 'outlet' in groups and groups['outlet']['volumeFlux']>0:
            meta['outletTemperature']=groups['outlet']['advectiveFlux']/groups['outlet']['volumeFlux']
    meta['outputSha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in root.glob(args.label+'.*')
                         if p.is_file() and not p.name.endswith('.metrics.json')}
    (root/(args.label+'.metrics.json')).write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta,indent=2));return r['code']

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mesh-cli',default='build/cartmesh2d_cli');p.add_argument('--transport-cli',default='build/cartmesh2d_transport_cli')
    p.add_argument('--case',choices=['channel','cavity','cylinder'],required=True);p.add_argument('--output',required=True)
    p.add_argument('--level',default='4');p.add_argument('--dt',default='.1');p.add_argument('--end',default='6');p.add_argument('--diffusivity',default='.1');p.add_argument('--nu',default='.1')
    p.add_argument('--flux-correction',choices=['unrestricted','bounded','bounded-spatial'],default='unrestricted')
    p.add_argument('--speed');p.add_argument('--relaxation',default='.6');p.add_argument('--courant',default='1');p.add_argument('--retries',default='18');p.add_argument('--error',action='store_true');p.add_argument('--rtol',default='.01');p.add_argument('--no-events',action='store_true');p.add_argument('--restart');p.add_argument('--label',default='result')
    p.add_argument('--isolate-live-output',action='store_true',help='Compute in temporary outputs/, then import closed files with hashes; avoids live workspace snapshot replacement')
    args=p.parse_args()
    if not args.isolate_live_output:raise SystemExit(execute(args))
    target=Path(args.output).resolve();stage=Path(tempfile.mkdtemp(prefix='cartmesh2d-thermal-'))/'outputs'
    args.output=str(stage)
    try:code=execute(args)
    finally:
        # Preserve failures as well as successes. The staging directory remains
        # available if importing files is interrupted; never copy active files.
        if stage.exists():shutil.copytree(stage,target,dirs_exist_ok=True)
    metrics=json.loads((target/(args.label+'.metrics.json')).read_text())
    for name,digest in metrics['outputSha256'].items():
        if hashlib.sha256((target/name).read_bytes()).hexdigest()!=digest:
            raise RuntimeError('Closed output import hash mismatch: '+name)
    raise SystemExit(code)
