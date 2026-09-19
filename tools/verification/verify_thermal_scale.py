#!/usr/bin/env python3
"""Bounded synchronized thermal scale study with two distinct analytic references.

The space-continuous backward-Euler reference removes the known temporal
truncation error for this unforced tangent vortex only. Its error still includes
spatial, iterative and carrier errors; it is not a continuum-accuracy shortcut.
"""
import argparse
import csv
import itertools
import json
import math
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import signal
import sys
import subprocess
import time
import verify_native_flow as native
import verify_scalar_transport as scalar
import verify_transient_flow as transient


def amplitude(rate, dt, steps):
    scalar.require(isinstance(steps,int) and not isinstance(steps,bool) and steps>0,
                   'reference requires a positive integer step count')
    scalar.require(math.isfinite(rate) and rate>0 and math.isfinite(dt) and dt>0,
                   'reference requires positive finite rate and dt')
    return math.exp(-steps*math.log1p(rate*dt))


def run_process_group(command, log, timeout):
    """Run a study command and terminate its entire descendant process group."""
    log.parent.mkdir(parents=True, exist_ok=True)
    record={'command':list(map(str,command)),
            'commandText':shlex.join(list(map(str,command))), 'returncode':None,
            'stdout':str(log.with_suffix('.stdout.log')),
            'stderr':str(log.with_suffix('.stderr.log')), 'timedOut':False}
    def text(value):
        return value.decode(errors='replace') if isinstance(value,bytes) else (value or '')
    try:
        process=subprocess.Popen(record['command'],cwd=native.REPO,text=True,
                                 stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                 start_new_session=os.name!='nt')
        try:
            stdout,stderr=process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            cleanup=[]
            try:
                if os.name=='nt':
                    try:
                        killed=subprocess.run(['taskkill','/PID',str(process.pid),'/T','/F'],
                                              text=True,capture_output=True,check=False,timeout=5)
                        if killed.returncode!=0:
                            cleanup.append(killed.stderr or killed.stdout or 'taskkill failed')
                            process.kill()
                    except subprocess.TimeoutExpired:
                        cleanup.append('taskkill timed out after 5s')
                        process.kill()
                else:
                    os.killpg(process.pid,signal.SIGKILL)
            except OSError as error:
                cleanup.append(str(error))
                try: process.kill()
                except OSError as fallback: cleanup.append(str(fallback))
            try:
                stdout_tail,stderr_tail=process.communicate(timeout=5)
            except subprocess.TimeoutExpired as tail_error:
                cleanup.append('post-cleanup wait timed out after 5s')
                try: process.kill()
                except OSError as error: cleanup.append(str(error))
                try:
                    stdout_tail,stderr_tail=process.communicate(timeout=1)
                except subprocess.TimeoutExpired as final_error:
                    cleanup.append('final wait timed out after 1s')
                    stdout_tail=getattr(final_error,'stdout',None) or getattr(tail_error,'stdout',None) or ''
                    stderr_tail=getattr(final_error,'stderr',None) or getattr(tail_error,'stderr',None) or ''
                    if process.stdout: process.stdout.close()
                    if process.stderr: process.stderr.close()
            stdout=text(stdout_tail)
            stderr=text(stderr_tail)
            record.update(returncode=process.returncode,status='failed',timedOut=True,
                          error=f'timed out after {timeout}s')
            if cleanup: record['cleanupErrors']=cleanup
        else:
            record.update(returncode=process.returncode,
                          status='passed' if process.returncode==0 else 'failed')
        Path(record['stdout']).write_text(text(stdout),encoding='utf-8')
        Path(record['stderr']).write_text(text(stderr),encoding='utf-8')
    except OSError as error:
        Path(record['stdout']).write_text('',encoding='utf-8')
        Path(record['stderr']).write_text(str(error)+'\n',encoding='utf-8')
        record.update(status='failed',error=str(error))
    return record


def read(prefix, suffix):
    return json.loads(Path(str(prefix)+suffix).read_text())


def reference_errors(prefix, flow_prefix):
    info=read(prefix,'.json'); flow=read(flow_prefix,'.json')
    scalar.require(info.get('verification')=='thermal-vortex' and info.get('evolvingFlow') is True,
                   'reference requires synchronized thermal vortex')
    scalar.require(not info.get('restart'),'full fixed-step trajectory from t=0 required')
    scalar.require(info.get('converged') is True and flow.get('converged') is True,
                   'reference requires converged outputs')
    scalar.require(flow.get('case')=='taylor-green' and flow.get('temporalDiscretization')=='backward-euler',
                   'reference carrier must use backward Euler Taylor-Green')
    steps=info['steps'];dt=scalar.finite(info['timeStep']);nu=scalar.finite(info['flowNu'])
    diffusivity=scalar.finite(info['diffusivity']);speed=scalar.finite(info['flowSpeed'])
    theta_be=amplitude(2*math.pi**2*diffusivity,dt,steps)
    velocity_be=speed*amplitude(2*math.pi**2*nu,dt,steps)
    scalar.require(speed>0 and flow['completedSteps']==steps,'carrier step count mismatch')
    for actual,expected in [(info['time'],dt*steps),(info['carrierTime'],dt*steps),
                            (flow['time'],dt*steps),(flow['dt'],dt),(flow['nu'],nu),(flow['speed'],speed)]:
        scalar.same(actual,expected,'mixed reference controls',absolute=1e-14,relative=1e-13)
    scalar.require(flow['convection']==info['flowConvection'],'carrier convection mismatch')
    relaxation=scalar.finite(info.get('flowVelocityRelaxation',.6))
    scalar.require(0<relaxation<=1,'invalid velocity relaxation')
    scalar.same(flow['velocityRelaxation'],relaxation,'carrier velocity relaxation mismatch')
    scalar.require(native.sha256_file(Path(str(prefix)+'.carrier.checkpoint'))==
                   native.sha256_file(Path(str(flow_prefix)+'.checkpoint')),
                   'coupled carrier differs from separately computed flow')
    theta_exact=math.exp(-2*math.pi**2*diffusivity*dt*steps)
    velocity_exact=speed*math.exp(-2*math.pi**2*nu*dt*steps)
    areas=[];terms={key:[] for key in ('scalarContinuousL2','scalarBackwardEulerL2',
                'velocityContinuousL2','velocityBackwardEulerL2','scalarTimeDiscretizationL2',
                'velocityTimeDiscretizationL2')}
    with Path(str(prefix)+'.cells.csv').open() as a,Path(str(flow_prefix)+'.cells.csv').open() as b:
        for i,(s,f) in enumerate(itertools.zip_longest(csv.DictReader(a),csv.DictReader(b))):
            scalar.require(s is not None and f is not None,'field cell count mismatch')
            scalar.require(int(s['cell'])==int(f['cell'])==i,'field cell ordering mismatch')
            for key in ('x','y','area'):scalar.same(s[key],f[key],'field geometry mismatch')
            x,y,area,value=map(scalar.finite,(s['x'],s['y'],s['area'],s['value']))
            scalar.require(0<x<1 and 0<y<1 and area>0,'reference requires unit-square interior cells')
            sx,sy=math.sin(math.pi*x),math.sin(math.pi*y)
            vx,vy=sx*math.cos(math.pi*y),-math.cos(math.pi*x)*sy
            u,v=scalar.finite(f['u']),scalar.finite(f['v'])
            areas.append(area)
            terms['scalarContinuousL2'].append(area*(value-theta_exact*sx*sy)**2)
            terms['scalarBackwardEulerL2'].append(area*(value-theta_be*sx*sy)**2)
            terms['velocityContinuousL2'].append(area*((u-velocity_exact*vx)**2+(v-velocity_exact*vy)**2))
            terms['velocityBackwardEulerL2'].append(area*((u-velocity_be*vx)**2+(v-velocity_be*vy)**2))
            terms['scalarTimeDiscretizationL2'].append(area*((theta_be-theta_exact)*sx*sy)**2)
            terms['velocityTimeDiscretizationL2'].append(area*(velocity_be-velocity_exact)**2*(vx*vx+vy*vy))
    scalar.require(len(areas)==info['cells']==flow['cells'],'reference cell count mismatch')
    scalar.same(math.fsum(areas),1,'reference requires full unit square',absolute=1e-10,relative=0)
    result={k:math.sqrt(math.fsum(v)/math.fsum(areas)) for k,v in terms.items()}
    result.update(cells=len(areas),h=1/math.sqrt(len(areas)),dt=dt,steps=steps,time=dt*steps,
                  nu=nu,diffusivity=diffusivity,speed=speed,velocityRelaxation=relaxation,
                  relaxationMetadata='explicit' if 'flowVelocityRelaxation' in info else 'legacy fixed 0.6',
                  carrierCheckpointByteIdentical=True,
                  definition='Area-weighted RMS at real cell centroids; velocity is vector error in m/s, scalar uses input units.',
                  scope='Backward-Euler reference is continuous in space; residual error includes spatial/carrier/iteration effects. Continuous-reference error remains the physical error diagnostic.')
    return result


def iteration_counts(prefix,flow_prefix):
    info=read(prefix,'.json');dt=info['timeStep'];steps=info['steps']
    with Path(str(prefix)+'.thermal-history.csv').open() as f:thermal=list(csv.DictReader(f))
    with Path(str(prefix)+'.history.csv').open() as f:history=list(csv.DictReader(f))
    flow=transient.read_history(Path(str(flow_prefix)+'.time-history.csv'))
    scalar.require(len(thermal)==len(flow)==steps,'accepted step coverage mismatch')
    flow_total=0;linear_total=0;corrections=0
    groups={i:[] for i in range(1,steps+1)}
    for row in history:
        step=int(row['step']);scalar.require(step in groups,'scalar history step outside trajectory')
        groups[step].append(row)
    for i,(row,carrier) in enumerate(zip(thermal,flow),1):
        scalar.require(int(row['step'])==i and int(row['accepted'])==1,'thermal step ordering/acceptance mismatch')
        scalar.same(row['time'],i*dt,'thermal step clock mismatch',absolute=1e-14,relative=1e-13)
        scalar.require(int(row['flowIterations'])==carrier['innerIterations'],'separate carrier iteration count mismatch')
        scalar.require(int(row['scalarIterations'])==len(groups[i])>0,'scalar correction coverage mismatch')
        flow_total+=carrier['innerIterations'];corrections+=len(groups[i])
        for j,r in enumerate(groups[i],1):
            scalar.require(int(r['iteration'])==j,'scalar correction ordering mismatch')
            scalar.same(r['time'],i*dt,'scalar correction clock mismatch',absolute=1e-14,relative=1e-13)
            count=int(r['linearIterations']);scalar.require(count>=0,'negative scalar linear iterations')
            linear_total+=count
    return {'flowNonlinearIterations':flow_total,'scalarCorrections':corrections,
            'scalarLinearIterations':linear_total,'steps':steps,
            'scope':'Raw accepted histories checked against separately computed carrier step counts; intermediate fields are not independently archived.'}


def audit(prefix, flow_prefix, output):
    result={'valid':False,'prefix':str(prefix),'flowPrefix':str(flow_prefix),'issues':[],
            'verifierSha256':{str(Path(p).resolve()):native.sha256_file(Path(p)) for p in
                              (__file__,native.__file__,scalar.__file__,transient.__file__)}}
    try:
        info=read(prefix,'.json')
        result['scalar']=scalar.verify(prefix)
        result['flow']=transient.verify(Path(info['mesh']),flow_prefix,output.with_name('flow-audit.json'))
        scalar.require(result['flow']['valid'],'independent carrier momentum/time audit failed')
        result['reference']=reference_errors(prefix,flow_prefix)
        result['iterations']=iteration_counts(prefix,flow_prefix)
        result['valid']=True
    except (ValueError,OSError,KeyError,OverflowError,TypeError) as exc:result['issues'].append(str(exc))
    native.write_json(output,result)
    return result


def series_checks(cases):
    scalar.require(len(cases)>=3,'spatial study requires at least three grids')
    scalar.require(all(c.get('valid') is True for c in cases),'invalid physical audit in series')
    rows=[dict(c['reference']) for c in cases]
    for row in rows[1:]:
        for key in ('dt','steps','time','nu','diffusivity','speed','velocityRelaxation'):
            scalar.require(row[key]==rows[0][key],'mixed spatial study '+key)
    decreasing={key:True for key in ('scalarBackwardEulerL2','velocityBackwardEulerL2')}
    for coarse,fine in zip(rows,rows[1:]):
        scalar.require(0<fine['h']<coarse['h'],'spatial grids must strictly refine')
        for key in decreasing:
            a,b=coarse[key],fine[key]
            scalar.require(math.isfinite(a) and math.isfinite(b) and a>=0 and b>=0,'invalid error')
            fine[key+'Order']=math.log(a/b)/math.log(coarse['h']/fine['h']) if a>0 and b>0 else None
            decreasing[key] &= b<a
    return {'valid':all(decreasing.values()),'decreases':decreasing,'series':rows,
            'scope':'Study consistency and decreasing discrete-time-reference errors only; no universal engineering error threshold.'}


def generate(args):
    root=args.output.resolve();scalar.require(not root.exists(),'output must be fresh')
    root.mkdir(parents=True)
    binaries={k:getattr(args,k).resolve(strict=True) for k in ('mesh_cli','flow_cli','transport_cli')}
    report={'valid':False,'runs':[],'cases':[],'issues':[],'meshes':[],
            'executables':{str(p):native.sha256_file(p) for p in binaries.values()},
            'runnerSha256':native.sha256_file(Path(__file__).resolve()),
            'platform':platform.platform(),'timeoutSeconds':args.timeout,
            'velocityRelaxation':args.velocity_relaxation,
            'scope':'Serial two-equation unforced thermal vortex scale/accuracy study; fixed properties and step count. Not long-time, external CFD, GUI or turbulent qualification.'}
    def save():native.write_json(root/'study.json',report)
    def execute(cmd,label):
        scalar.require(shutil.disk_usage(root).free>=1.5*1024**3,'resource stop: less than 1.5 GiB free')
        wrapped=list(map(str,cmd));wrapper=['/usr/bin/time','-l' if platform.system()=='Darwin' else '-v']
        if Path('/usr/bin/time').is_file():wrapped=wrapper+wrapped
        start=time.monotonic();record=run_process_group(wrapped,root/'logs'/label,args.timeout)
        record['elapsedSeconds']=time.monotonic()-start
        stderr=Path(record['stderr']).read_text()
        pattern=r'(\d+)\s+maximum resident set size' if platform.system()=='Darwin' else r'Maximum resident set size \(kbytes\):\s*(\d+)'
        m=re.search(pattern,stderr);record['peakRssBytes']=int(m[1])*(1 if platform.system()=='Darwin' else 1024) if m else None
        report['runs'].append(record);save()
        scalar.require(record['returncode']==0 and not record['timedOut'],label+' failed; recorded output retained')
    try:
        for across in args.cells_across:
            level=math.ceil(math.log2(across+2));padding=((1<<level)/across-1)/2
            directory=root/f'n{across}';directory.mkdir()
            if args.reuse_mesh_study:
                mesh=(args.reuse_mesh_study/f'n{across}'/'mesh'/'square.solver.cm2d').resolve(strict=True)
            else:
                mesh_prefix=directory/'mesh'/'square';mesh_prefix.parent.mkdir()
                request=native.Request(f'n{across}','manufactured',level,padding,.1,1.)
                execute(native.mesh_command(binaries['mesh_cli'],mesh_prefix,request,root),f'mesh-{across}')
                mesh=Path(str(mesh_prefix)+'.solver.cm2d')
            report['meshes'].append({'path':str(mesh),'sha256':native.sha256_file(mesh),'cellsAcross':across})
            save();prefix=directory/'thermal';flow=directory/'flow'
            execute([binaries['transport_cli'],'--mesh',mesh,'--output',prefix,'--verification','thermal-vortex',
                     '--dt',str(args.dt),'--steps',str(args.steps),'--flow-nu','.1','--diffusivity','.02',
                     '--flow-convection','limited-linear','--convection','limited-linear','--pressure-preconditioner','aggregation',
                     '--flow-velocity-relaxation',str(args.velocity_relaxation)],f'thermal-{across}')
            execute([binaries['flow_cli'],'--mesh',mesh,'--output',flow,'--case','taylor-green',
                     '--time-step',str(args.dt),'--steps',str(args.steps),'--nu','.1','--speed','1','--tolerance','1e-8',
                     '--convection','limited-linear','--pressure-preconditioner','aggregation','--profile',
                     '--velocity-relaxation',str(args.velocity_relaxation)],f'flow-{across}')
            execute([sys.executable,Path(__file__).resolve(),'--audit-prefix',prefix,'--flow-prefix',flow,'--output',directory/'audit.json'],f'audit-{across}')
            record=read(directory/'audit','.json');report['cases'].append(record);save()
            scalar.require(record['valid'],'independent thermal/flow audit failed')
            scalar.require(record['reference']['cells']==across**2,'mesh did not produce requested square grid')
        report['spatial']=series_checks(report['cases']);report['valid']=report['spatial']['valid']
        if not report['valid']:report['issues'].append('discrete-time-reference errors did not decrease')
    except (ValueError,OSError,KeyError) as exc:report['issues'].append(str(exc))
    save();return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--audit-prefix',type=Path)
    p.add_argument('--flow-prefix',type=Path);p.add_argument('--cells-across',type=int,nargs='+',default=[100,200,320])
    p.add_argument('--dt',type=float,default=.005);p.add_argument('--steps',type=int,default=2)
    p.add_argument('--timeout',type=float,default=180.)
    p.add_argument('--velocity-relaxation',type=float,default=.6)
    p.add_argument('--reuse-mesh-study',type=Path,help='Reuse already generated nN/mesh/square.solver.cm2d files; never copy or modify them')
    for name in ('mesh','flow','transport'):p.add_argument('--'+name+'-cli',type=Path,default=Path('build/cartmesh2d_'+('cli' if name=='mesh' else name+'_cli')))
    a=p.parse_args()
    if not math.isfinite(a.timeout) or not 0<a.timeout<=300:p.error('timeout must be in (0,300] seconds')
    if not math.isfinite(a.velocity_relaxation) or not 0<a.velocity_relaxation<=1:p.error('velocity relaxation must be in (0,1]')
    if a.audit_prefix:
        if not a.flow_prefix:p.error('--flow-prefix required for audit')
        result=audit(a.audit_prefix,a.flow_prefix,a.output)
    else:
        if len(a.cells_across)<3 or any(n<4 or n>1024 for n in a.cells_across) or any(b<=a for a,b in zip(a.cells_across,a.cells_across[1:])):p.error('at least three increasing cells-across values in [4,1024] required')
        if a.steps<1 or not math.isfinite(a.dt) or a.dt<=0:p.error('positive dt/steps required')
        result=generate(a)
    print(json.dumps({'valid':result['valid'],'issues':result.get('issues',[])}),flush=True)
    return 0 if result['valid'] else 1

if __name__=='__main__':raise SystemExit(main())
