#!/usr/bin/env python3
"""Static-pressure-driven Poiseuille verification on independent rectangular grids.

Kinematic pressure drop 4.8 over L=4, H=1, nu=.1 gives u=6*y*(1-y),
Q'=1. These manufactured rectangular meshes test the flow solver, not mesher
quality. Full independent momentum/continuity and checkpoint audits precede
reference comparisons; failures and commands remain in the output directory.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import shlex
import shutil
import time
import verify_native_flow as native
import verify_transient_flow as transient


def rectangle(path, n):
    nx,ny=4*n,n
    vertices=[(i/n,j/n) for j in range(ny+1) for i in range(nx+1)]
    edges,cells,pairs=[],[],{}
    for j in range(ny):
        for i in range(nx):
            a=j*(nx+1)+i;ids=[a,a+1,a+nx+2,a+nx+1];incident=[]
            for x,y in zip(ids,ids[1:]+ids[:1]):
                pair=tuple(sorted((x,y)))
                if pair not in pairs:
                    pairs[pair]=len(edges);edges.append([len(edges),x,y,len(cells),-1,2])
                else:edges[pairs[pair]][4:]=[len(cells),0]
                incident.append(pairs[pair])
            cells.append([len(cells),0,0,1/n**2,4,*ids,4,*incident])
    lines=['CM2D 1',f'VERTICES {len(vertices)}',*[f'{i} {x:.17g} {y:.17g}' for i,(x,y) in enumerate(vertices)],
           f'EDGES {len(edges)}',*[' '.join(map(str,e)) for e in edges],f'CELLS {len(cells)}',
           *[' '.join(map(str,c)) for c in cells],'AUDIT 0 0 0 0 0 0 0','END','']
    path.write_text('\n'.join(lines))


def pressure_boundaries(source,dest,left=4.8,right=0):
    lines=source.read_text().splitlines()
    for i,line in enumerate(lines):
        if not line.startswith('BOUNDARY '):continue
        t=shlex.split(line)
        if t[7] in ('velocity-inlet','pressure-outlet'):
            t[7]='pressure-opening';t[9:]=['0','0',str(left if float(t[5])<0 else right)]
        t[8]=json.dumps(t[8]);lines[i]=' '.join(t)
    dest.write_text('\n'.join(lines)+'\n')


def study(root,cli,levels=(8,16,32),convection='face-limited-linear',extended=True,tolerance=1e-10):
    if not levels or any(type(n) is not int or not 4<=n<=32 for n in levels) or any(a>=b for a,b in zip(levels,levels[1:])):
        raise ValueError('Use strictly increasing grid levels from 4 through 32')
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=False);cli=Path(cli).resolve(strict=True)
    report=dict(valid=False,accuracyQualification='case-specific rectangular Poiseuille only',
        reference=dict(nu=.1,length=4,height=1,pressureDrop=4.8,velocity='6*y*(1-y)',pressure='4.8*(1-x/4)',volumeFlowPerDepth=1),
        limits=dict(commandSeconds=180,totalSeconds=900,retainedBytes=50*1024**2,minimumFreeBytes=5*1024**3),
        targets=dict(fineVelocityL2=.002,fineFluxRelativeError=.003,pressureMaxError=1e-6,minimumObservedOrder=1.9),
        binarySha256=native.sha256_file(cli),runs=[],cases=[],issues=[])
    began=time.monotonic()
    def save():native.write_json(root/'study.json',report)
    def run(command,label,success=True):
        if shutil.disk_usage(root).free<5*1024**3:raise ValueError('below 5 GiB free space')
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file())>45*1024**2:raise ValueError('retained budget reached')
        remaining=900-(time.monotonic()-began)
        if remaining<=0:raise ValueError('total 900 second budget exhausted')
        start=time.monotonic();record=native.run(list(map(str,command)),root/'logs'/label,min(180,remaining))
        record['elapsedSeconds']=time.monotonic()-start;report['runs'].append(record);save()
        print(label,record['status'],round(record['elapsedSeconds'],3),flush=True)
        if (record['status']=='passed')!=success:raise ValueError(label+' unexpected status; see preserved logs')
    def solve(mesh,bc,prefix,extra=(),success=True):
        run([cli,'--mesh',mesh,'--case','custom','--boundary',bc,'--output',prefix,'--nu',.1,'--speed',1,
             '--convection',convection,'--pressure-preconditioner','aggregation','--tolerance',tolerance,
             '--max-iterations',10000,*extra],prefix.name,success)
        if not success:return
        if '--time-step' in extra:
            result=transient.verify(mesh,prefix,Path(str(prefix)+'.audit.json'))
        else:
            result=native.verify_case(mesh,prefix,'custom',.1,1,native.argument_parser().parse_args(['--max-iterations','10000']))
            native.write_json(Path(str(prefix)+'.audit.json'),result)
        if not result['valid']:raise ValueError(prefix.name+' independent audit failed: '+str(result.get('issues')))
    save()
    try:
        for n in levels:
            mesh=root/f'n{n}.solver.cm2d';rectangle(mesh,n)
            template=root/f'n{n}.template';bc=root/f'n{n}.boundaries';prefix=root/f'n{n}'
            run([cli,'--mesh',mesh,'--case','duct','--export-boundaries',template],f'n{n}-template')
            pressure_boundaries(template,bc);solve(mesh,bc,prefix)
            parsed=native.read_cm2d(mesh);measured=native.measure(parsed,1e-11,1e-9)
            rows=native.read_cells(Path(str(prefix)+'.cells.csv'),parsed,measured,'custom')
            faces,_=native.read_faces(Path(str(prefix)+'.faces.csv'),parsed)
            q=math.fsum(faces[e.id]['flux'] for e in parsed.edges if e.neighbour<0 and native.face_centre(parsed,e)[0]==4)
            err=native.weighted_l2([r['u']-6*r['y']*(1-r['y']) for r in rows],measured.areas)
            item=dict(n=n,cells=len(rows),velocityL2=err,flowPerDepth=q,flowRelativeError=abs(q-1),
                pressureMaxError=max(abs(r['p']-4.8*(1-r['x']/4)) for r in rows),
                transverseMax=max(abs(r['v']) for r in rows),summary=json.loads(Path(str(prefix)+'.json').read_text()))
            if item['pressureMaxError']>=1e-6:raise ValueError('linear pressure reference failed')
            if report['cases']:
                previous=report['cases'][-1]
                item['observedVelocityOrder']=math.log(previous['velocityL2']/err)/math.log(n/previous['n'])
                if item['observedVelocityOrder']<1.9:raise ValueError('velocity refinement order below 1.9')
            report['cases'].append(item);save()
            if extended and n==levels[0]:
                whole,first,resumed=(root/k for k in ('whole','first','resumed'))
                solve(mesh,bc,whole,('--time-step',.02,'--steps',2))
                solve(mesh,bc,first,('--time-step',.02,'--steps',1))
                solve(mesh,bc,resumed,('--time-step',.02,'--steps',1,'--restart',str(first)+'.checkpoint'))
                if Path(str(whole)+'.checkpoint').read_bytes()!=Path(str(resumed)+'.checkpoint').read_bytes():raise ValueError('checkpoint continuation differs')
                report['restartByteIdentical']=True
                equal=root/'equal.boundaries';pressure_boundaries(template,equal,7.25,7.25);solve(mesh,equal,root/'equal')
                with (root/'equal.cells.csv').open() as stream:
                    if any(math.hypot(float(r['u']),float(r['v']))>1e-10 or abs(float(r['p'])-7.25)>1e-9 for r in csv.DictReader(stream)):
                        raise ValueError('equal-pressure quiescent reference failed')
                # Changing physical type is incompatible with the checkpoint.
                solve(mesh,equal,root/'bad-restart',('--time-step',.02,'--steps',1,'--restart',str(first)+'.checkpoint'),False)
        if max(levels)>=32 and (report['cases'][-1]['velocityL2']>=.002 or report['cases'][-1]['flowRelativeError']>=.003):
            raise ValueError('fine Poiseuille reference target failed')
        report['artifactSha256']={str(p.relative_to(root)):native.sha256_file(p) for p in root.rglob('*') if p.is_file() and p.name!='study.json'}
        report['valid']=True;save()
    except Exception as error:report['issues'].append(str(error));save();raise
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cli',type=Path,default=Path('build/cartmesh2d_flow_cli'));parser.add_argument('--levels',type=int,nargs='+',default=[8,16,32])
    parser.add_argument('--convection',choices=['upwind','limited-linear','face-limited-linear'],default='face-limited-linear')
    args=parser.parse_args();print(json.dumps(study(args.output,args.cli,args.levels,args.convection)['cases'],indent=2))
