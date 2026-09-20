#!/usr/bin/env python3
"""Bounded named-boundary Poiseuille refinement and rotation study.

Reference is a 4 m by 1 m fully developed channel, Umax=1 m/s,
nu=.1 m2/s. Thus u=4*y*(1-y), p=pout+.8*(4-x), Q'=2/3.
The study checks an a-priori 2% fine-grid target for velocity, pressure
slope and total viscous wall force, not universal engineering qualification.
"""
import argparse
import json
import math
from pathlib import Path
import shlex
import shutil
import time

import verify_native_flow as native

REFERENCE = 'https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/1a114d602956fa0dd328155f9b45f93d_MIT2_25F13_Couet_and_Pois.pdf'


def transform_mesh(source, destination, angle):
    c, s = math.cos(angle), math.sin(angle)
    lines = source.read_text().splitlines()
    for i in range(2, 2+int(lines[1].split()[1])):
        ident, x, y = lines[i].split()
        x, y = float(x), float(y)
        lines[i] = f'{ident} {c*x-s*y:.17g} {s*x+c*y:.17g}'
    destination.write_text('\n'.join(lines)+'\n')


def transform_boundary(source, destination, angle, pressure):
    c, s = math.cos(angle), math.sin(angle)
    lines = source.read_text().splitlines()
    for i, line in enumerate(lines):
        if not line.startswith('BOUNDARY '): continue
        t = shlex.split(line)
        for j in (3,5,9):
            x,y = float(t[j]),float(t[j+1])
            t[j:j+2] = [f'{c*x-s*y:.17g}',f'{s*x+c*y:.17g}']
        if t[7] == 'pressure-outlet': t[11] = f'{pressure:.17g}'
        t[8] = json.dumps(t[8])
        lines[i] = ' '.join(t)
    destination.write_text('\n'.join(lines)+'\n')


def metrics(mesh_path, prefix, angle):
    mesh=native.read_cm2d(mesh_path);measured=native.measure(mesh,1e-11,1e-9)
    cells=native.read_cells(Path(str(prefix)+'.cells.csv'),mesh,measured,'custom')
    faces,_=native.read_faces(Path(str(prefix)+'.faces.csv'),mesh)
    summary=json.loads(Path(str(prefix)+'.json').read_text())
    bc=native.audit_explicit_boundaries(prefix,mesh,measured,summary)
    c,s=math.cos(angle),math.sin(angle)
    rows=[]
    for row in cells:
        x,y=c*row['x']+s*row['y'],-s*row['x']+c*row['y']
        u,v=c*row['u']+s*row['v'],-s*row['u']+c*row['v']
        rows.append(dict(x=x,y=y,u=u,v=v,p=row['p'],area=row['area']))
    weights=[r['area'] for r in rows];area=math.fsum(weights)
    xmean=math.fsum(r['x']*r['area'] for r in rows)/area
    pmean=math.fsum(r['p']*r['area'] for r in rows)/area
    slope=math.fsum(r['area']*(r['x']-xmean)*(r['p']-pmean) for r in rows)/math.fsum(r['area']*(r['x']-xmean)**2 for r in rows)
    q=-math.fsum(faces[b['face']]['flux'] for b in bc if b['type']=='velocity-inlet')
    force=c*summary['wallViscousForceX']+s*summary['wallViscousForceY']
    values=dict(velocityL2=native.weighted_l2([r['u']-4*r['y']*(1-r['y']) for r in rows],weights),
                transverseVelocityL2=native.weighted_l2([r['v'] for r in rows],weights),
                pressureL2Relative=native.weighted_l2([r['p']-(.25+.8*(4-r['x'])) for r in rows],weights)/3.2,
                fittedPressureGradient=slope,pressureGradientRelativeError=abs(slope+.8)/.8,
                flowRate=q,flowRateRelativeError=abs(q-2/3)/(2/3),
                viscousWallForce=force,viscousWallForceRelativeError=abs(force-3.2)/3.2)
    return values, rows


def generate(args):
    root=args.output.resolve()
    if root.exists(): raise ValueError('Choose a fresh output directory; failures are retained.')
    root.mkdir(parents=True)
    geometry=root/'channel.xy';geometry.write_text('0 0\n4 0\n4 1\n0 1\n')
    mesh_cli=args.mesh_cli.resolve(strict=True);flow_cli=args.flow_cli.resolve(strict=True)
    report=dict(valid=False,scope='Named-boundary fully developed laminar channel; reference accuracy and coordinate invariance for this case only.',
        reference=REFERENCE,definition=dict(length=4,height=1,peakSpeed=1,nu=.1,peakReynolds=10,pressureOutlet=.25,
        velocity='4*y*(1-y)',pressure='.25+.8*(4-x)',volumeFlux=2/3,wallViscousForce=3.2),
        predeclaredTargets=dict(fineVelocityL2=.02,finePressureGradientRelativeError=.02,fineWallForceRelativeError=.02,
        fineFlowRateRelativeError=.005,rotationMaxAbs=1e-6,iterativeVelocityChangeMaxAbs=1e-6),
        controls=dict(convection=args.convection,tolerance=1e-8,tightTolerance=1e-10,maxIterations=args.max_iterations,timeoutSeconds=args.timeout),
        sourceSha256={str(p):native.sha256_file(p) for p in (mesh_cli,flow_cli,Path(__file__))},runs=[],cases=[],issues=[])
    def save():native.write_json(root/'summary.json',report)
    def run(command,label):
        if shutil.disk_usage(root).free<5*1024**3:raise ValueError('resource guard: under 5 GiB free')
        start=time.monotonic();record=native.run(list(map(str,command)),root/'logs'/label,args.timeout)
        record['seconds']=time.monotonic()-start;report['runs'].append(record);save()
        print(label,record['status'],round(record['seconds'],2),flush=True)
        if record['status']!='passed':raise ValueError(label+' failed; see preserved logs')
    def solve(mesh,bc,prefix,label,tolerance):
        run([flow_cli,'--mesh',mesh,'--case','custom','--boundary',bc,'--output',prefix,'--nu',.1,'--speed',1,
            '--convection',args.convection,'--pressure-preconditioner','aggregation','--max-iterations',args.max_iterations,
            '--tolerance',tolerance],label)
        audit=native.verify_case(mesh,prefix,'custom',.1,1,native.argument_parser().parse_args(['--max-iterations',str(args.max_iterations)]))
        native.write_json(Path(str(prefix)+'.audit.json'),audit)
        if not audit['valid']:raise ValueError(label+': independent audit failed '+str(audit['issues']))
        return audit
    try:
        for level in args.levels:
            directory=root/f'l{level}';directory.mkdir()
            original=directory/'original.solver.cm2d';template=directory/'template.boundaries'
            run([mesh_cli,geometry,directory/'original',level,1/(2**level-2),.1,'interior',directory/'foam',level,0],f'l{level}-mesh')
            run([flow_cli,'--mesh',original,'--case','channel','--speed',1,'--export-boundaries',template],f'l{level}-boundaries')
            base_rows=None
            for label,angle in (('axis',0),('rotated',.63)):
                mesh=directory/(label+'.solver.cm2d');bc=directory/(label+'.boundaries');prefix=directory/label
                transform_mesh(original,mesh,angle);transform_boundary(template,bc,angle,.25)
                audit=solve(mesh,bc,prefix,f'l{level}-{label}',1e-8)
                values,rows=metrics(mesh,prefix,angle)
                item=dict(label=f'l{level}-{label}',level=level,angle=angle,mesh=str(mesh),prefix=str(prefix),
                          cells=audit['counts']['cells'],h=audit['meshMeasurement']['characteristicH'],metrics=values,
                          independentValid=audit['valid'],iterations=audit['native']['iterations'],
                          artifactSha256=audit['artifactSha256'],meshSha256=audit['meshSha256'])
                if base_rows is None:base_rows=rows
                else:
                    errors={k:max(abs(a[k]-b[k]) for a,b in zip(base_rows,rows)) for k in ('u','v','p')}
                    item['rotationMaxAbs']=errors
                    if max(errors.values())>1e-6:report['issues'].append(f'l{level}: rotation invariance failed')
                report['cases'].append(item);save()
                print(item['label'],values,flush=True)
            if level==args.levels[-1]:
                tight=directory/'tight'
                solve(directory/'axis.solver.cm2d',directory/'axis.boundaries',tight,f'l{level}-tight',1e-10)
                tight_metrics,tight_rows=metrics(directory/'axis.solver.cm2d',tight,0)
                report['iterationSensitivity']=dict(maxAbs={k:max(abs(a[k]-b[k]) for a,b in zip(base_rows,tight_rows)) for k in ('u','v','p')},
                    metrics=tight_metrics)
                if max(report['iterationSensitivity']['maxAbs'][k] for k in ('u','v'))>1e-6:
                    report['issues'].append('iterative velocity error exceeds predeclared target')
        axis=[c for c in report['cases'] if c['angle']==0]
        for key in ('velocityL2','pressureGradientRelativeError','viscousWallForceRelativeError','flowRateRelativeError'):
            sequence=[]
            for coarse,fine in zip(axis,axis[1:]):
                a,b=coarse['metrics'][key],fine['metrics'][key]
                sequence.append(dict(coarse=coarse['cells'],fine=fine['cells'],observedOrder=math.log(a/b)/math.log(coarse['h']/fine['h']) if min(a,b)>0 else None))
                if not b<a:report['issues'].append(key+' does not decrease under refinement')
            report.setdefault('refinement',{})[key]=sequence
        fine=axis[-1]['metrics']
        for key,limit in [('velocityL2',.02),('pressureGradientRelativeError',.02),('viscousWallForceRelativeError',.02),('flowRateRelativeError',.005)]:
            if fine[key]>limit:report['issues'].append(key+' misses predeclared fine-grid accuracy target')
        report['valid']=not report['issues'];save()
    except Exception as error:
        report['issues'].append(str(error));save();raise
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--levels',type=int,nargs='+',default=[5,6,7]);p.add_argument('--timeout',type=int,default=180)
    p.add_argument('--convection', choices=['limited-linear','face-limited-linear','upwind'],default='face-limited-linear')
    p.add_argument('--max-iterations',type=int,default=6000)
    p.add_argument('--mesh-cli',type=Path,default=native.REPO/'build/cartmesh2d_cli')
    p.add_argument('--flow-cli',type=Path,default=native.REPO/'build/cartmesh2d_flow_cli');args=p.parse_args()
    if len(args.levels)<3 or sorted(set(args.levels))!=args.levels or min(args.levels)<3 or max(args.levels)>9 or args.timeout<=0 or args.max_iterations<10:p.error('need at least three increasing levels in 3..9 and positive timeout')
    result=generate(args);print(json.dumps({'valid':result['valid'],'issues':result['issues']}));raise SystemExit(0 if result['valid'] else 1)
