#!/usr/bin/env python3
"""Independent half-channel Poiseuille and frictionless acceleration checks."""
import argparse
import json
import math
from pathlib import Path
import shlex
import shutil
import time
import verify_native_flow as native
import verify_transient_flow as transient
from verify_pressure_openings import rectangle,pressure_boundaries


def symmetry_file(source,dest,both=False):
    lines=source.read_text().splitlines()
    for i,line in enumerate(lines):
        if not line.startswith('BOUNDARY '):continue
        t=shlex.split(line)
        if t[7]=='wall' and (both or float(t[6])>0):
            t[7]='symmetry';t[8]='symmetry-top' if float(t[6])>0 else 'symmetry-bottom'
        t[8]=json.dumps(t[8]);lines[i]=' '.join(t)
    dest.write_text('\n'.join(lines)+'\n')


def study(root,cli,levels=(4,8,16),convection='face-limited-linear'):
    if not levels or any(type(n) is not int or not 4<=n<=16 for n in levels) or any(a>=b for a,b in zip(levels,levels[1:])):
        raise ValueError('Use strictly increasing grid levels from 4 through 16')
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=False);cli=Path(cli).resolve(strict=True)
    report=dict(valid=False,binarySha256=native.sha256_file(cli),cases=[],runs=[],issues=[],
        reference=dict(length=4,height=1,nu=.1,pressureDrop=1.2,halfChannelVelocity='1.5*y*(2-y)',volumeFlow=1,plugAcceleration=.3),
        limits=dict(commandSeconds=180,totalSeconds=600,retainedBytes=30*1024**2,minimumFreeBytes=5*1024**3),
        targets=dict(fineVelocityL2=.002,fineFluxError=.003,minObservedOrder=1.9,pressureMaxError=1e-6,plugVelocityMaxError=1e-7),
        scope='Rectangular laminar half-channel and transient pressure-driven frictionless plug only; no general CFD accuracy qualification.')
    began=time.monotonic()
    def save():native.write_json(root/'study.json',report)
    def run(command,label,success=True):
        remaining=600-(time.monotonic()-began)
        if remaining<=0 or shutil.disk_usage(root).free<5*1024**3:raise ValueError('time/disk budget reached')
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file())>25*1024**2:raise ValueError('retained budget reached')
        start=time.monotonic();r=native.run(list(map(str,command)),root/'logs'/label,min(180,remaining));r['seconds']=time.monotonic()-start
        report['runs'].append(r);save();print(label,r['status'],round(r['seconds'],3),flush=True)
        if (r['status']=='passed')!=success:raise ValueError(label+' unexpected status')
    def solve(mesh,bc,prefix,extra=(),success=True):
        run([cli,'--mesh',mesh,'--case','custom','--boundary',bc,'--output',prefix,'--nu',.1,'--speed',1,
             '--max-iterations',12000,'--tolerance',1e-10,'--pressure-preconditioner','aggregation','--convection',convection,*extra],prefix.name,success)
        if not success:return
        if '--time-step' in extra:a=transient.verify(mesh,prefix,Path(str(prefix)+'.audit.json'))
        else:
            a=native.verify_case(mesh,prefix,'custom',.1,1,native.argument_parser().parse_args(['--max-iterations','12000']))
            native.write_json(Path(str(prefix)+'.audit.json'),a)
        if not a['valid']:raise ValueError(prefix.name+' failed independent audit: '+str(a.get('issues')))
    save()
    try:
        for n in levels:
            mesh=root/f'n{n}.solver.cm2d';rectangle(mesh,n);template=root/f'n{n}.template';pressure=root/f'n{n}.pressure';bc=root/f'n{n}.boundaries';prefix=root/f'n{n}'
            run([cli,'--mesh',mesh,'--case','duct','--export-boundaries',template],f'n{n}-template')
            pressure_boundaries(template,pressure,1.2,0);symmetry_file(pressure,bc);solve(mesh,bc,prefix)
            m=native.read_cm2d(mesh);measure=native.measure(m,1e-11,1e-9);cells=native.read_cells(Path(str(prefix)+'.cells.csv'),m,measure,'custom')
            faces,_=native.read_faces(Path(str(prefix)+'.faces.csv'),m);s=json.loads(Path(str(prefix)+'.json').read_text())
            q=math.fsum(faces[e.id]['flux'] for e in m.edges if e.neighbour<0 and native.face_centre(m,e)[0]==4)
            error=native.weighted_l2([r['u']-1.5*r['y']*(2-r['y']) for r in cells],measure.areas)
            p_error=max(abs(r['p']-1.2*(1-r['x']/4)) for r in cells)
            zero=max(abs(faces[b['face']]['flux']) for b in s['boundaryConditions'] if b['type']=='symmetry')
            if p_error>=1e-6 or zero!=0:raise ValueError('pressure reference or symmetry penetration failed')
            item=dict(n=n,cells=len(cells),velocityL2=error,flow=q,flowError=abs(q-1),pressureMaxError=p_error,
                      symmetryFluxMax=zero,iterations=s['iterations'])
            if report['cases']:
                prev=report['cases'][-1];item['observedOrder']=math.log(prev['velocityL2']/error)/math.log(n/prev['n'])
                if item['observedOrder']<1.9:raise ValueError('half-channel velocity order failed')
            report['cases'].append(item);save()
            if n==levels[0]:
                plug=root/'plug.boundaries';symmetry_file(pressure,plug,True)
                whole,first,resumed=(root/k for k in ('whole','first','resumed'))
                solve(mesh,plug,whole,('--time-step',.02,'--steps',3))
                solve(mesh,plug,first,('--time-step',.02,'--steps',1))
                solve(mesh,plug,resumed,('--time-step',.02,'--steps',2,'--restart',str(first)+'.checkpoint'))
                if Path(str(whole)+'.checkpoint').read_bytes()!=Path(str(resumed)+'.checkpoint').read_bytes():raise ValueError('plug restart differs')
                rows=native.read_cells(Path(str(whole)+'.cells.csv'),m,measure,'custom')
                err=max(math.hypot(r['u']-.018,r['v']) for r in rows)
                if err>=1e-7:raise ValueError('uniform acceleration reference failed')
                report['plug']=dict(time=.06,expectedVelocity=.018,maxVelocityError=err,restartByteIdentical=True)
                # A nonzero pressure value on a symmetry plane must fail.
                bad=root/'bad.boundaries';bad.write_text(bc.read_text().replace('symmetry "symmetry-top" 0 0 0','symmetry "symmetry-top" 0 0 1'))
                if bad.read_text()==bc.read_text():raise ValueError('invalid condition mutation failed')
                solve(mesh,bad,root/'invalid',success=False)
        fine=report['cases'][-1]
        if max(levels)>=16 and (fine['velocityL2']>=.002 or fine['flowError']>=.003):raise ValueError('fine half-channel targets failed')
        report['artifactSha256']={str(p.relative_to(root)):native.sha256_file(p) for p in root.rglob('*') if p.is_file() and p.name!='study.json'}
        report['valid']=True;save()
    except Exception as error:report['issues'].append(str(error));save();raise
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--cli',type=Path,default=Path('build/cartmesh2d_flow_cli'))
    a=p.parse_args();print(json.dumps(study(a.output,a.cli),indent=2))
