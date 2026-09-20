#!/usr/bin/env python3
"""Short-window temporal sensitivity from one independently audited wake state.

This cannot remove error accumulated before the common starting state and is
not a spatial-convergence or long-time shedding-frequency qualification.
"""
import argparse,csv,json,math,shutil,subprocess,tempfile,time
from pathlib import Path
import verify_flow_trajectory as trajectory
n=trajectory.native


def study(args):
    root=args.output.resolve()
    if root.exists():raise ValueError('Choose a new output directory; preserve earlier failures')
    checked=json.loads(args.verified.read_text());source=Path(checked['source'])
    if checked.get('valid') is not True or n.sha256_file(source)!=checked['sourceSha256']:
        raise ValueError('The source trajectory is not bound to its verified record')
    original=json.loads(source.read_text());mesh=Path(original['mesh']);binary=Path(original['binary'])
    if n.sha256_file(mesh)!=checked['meshSha256'] or n.sha256_file(binary)!=checked['binarySha256']:
        raise ValueError('Mesh/executable changed after trajectory audit')
    candidates=[r for r in original['runs'] if r.get('valid') is True and r.get('returnCode')==0 and r.get('time')==args.start_time]
    if len(candidates)!=1 or not any(s.get('time')==args.start_time for s in checked['segments']):
        raise ValueError('Starting time must be a fully audited segment endpoint')
    selected=candidates[0];command=selected['command'];prefix=Path(command[command.index('--output')+1])
    if prefix.parent.parent!=source.resolve().parent:raise ValueError('Source segment path differs')
    root.mkdir(parents=True);report=dict(format='cartmesh2d-wake-time-refinement-v1',valid=False,
        scope=__doc__.strip(),sourceVerified=str(args.verified),sourceVerifiedSha256=n.sha256_file(args.verified),
        sourceTrajectorySha256=checked['sourceSha256'],mesh=str(mesh),meshSha256=checked['meshSha256'],
        binary=str(binary),binarySha256=checked['binarySha256'],driverSha256=n.sha256_file(Path(__file__)),
        startTime=args.start_time,duration=args.duration,steps=args.steps,
        budget=dict(perRunSeconds=args.timeout,totalSeconds=args.total_timeout,retainedBytes=60*1024**2,minimumFreeBytes=5*1024**3),runs=[],issues=[])
    def save():n.write_json(root/'study.json',report)
    save();began=time.monotonic()
    try:
        with tempfile.TemporaryDirectory(prefix='cm2d-wake-origin-') as temporary:
            folder=trajectory.materialize(source.resolve().parent,prefix.parent.name,Path(temporary))
            # Independently bind the common starting state to its final field,
            # rather than trusting just a recorded hash.
            origin=trajectory.transient.verify(mesh,folder/'flow',root/'origin-audit.json')
            trajectory.audit_final_checkpoint(mesh,folder/'flow')
            if origin['time']!=args.start_time or n.sha256_file(folder/'flow.checkpoint')!=selected['checkpointSha256']:
                raise ValueError('Common initial state differs from audited source')
            restart=root/'initial.checkpoint';shutil.copyfile(folder/'flow.checkpoint',restart)
            report['initialCheckpointSha256']=n.sha256_file(restart);report['controls']=origin['controls'];save()
        for steps in args.steps:
            remaining=args.total_timeout-(time.monotonic()-began)
            if remaining<=0:raise ValueError('Total wall-clock budget exhausted')
            if shutil.disk_usage(root).free<5*1024**3:raise ValueError('Less than 5 GiB free disk')
            if sum(p.stat().st_size for p in root.rglob('*') if p.is_file())>50*1024**2:
                raise ValueError('Insufficient remaining space in 60 MiB study budget')
            folder=root/f'n{steps}';folder.mkdir();prefix=folder/'flow';dt=args.duration/steps
            # Inherit every physical/discretization control, replacing only
            # the original time controller and run-owned paths. Fixed equal
            # steps make the three histories directly comparable.
            flags={};items=iter(command[1:])
            for flag in items:flags[flag]=True if flag=='--profile' else next(items)
            for flag in ('--time-step','--end-time','--max-courant','--min-time-step','--max-time-steps','--max-step-retries','--profile','--output','--restart'):
                flags.pop(flag,None)
            actual=[str(binary),*(part for key,value in flags.items() for part in (key,str(value))),
                '--output',str(prefix),'--restart',str(restart),'--time-step',str(dt),'--steps',str(steps),'--profile']
            run=dict(steps=steps,dt=dt,command=actual,valid=False);start=time.monotonic()
            with (folder/'stdout.log').open('w') as out,(folder/'stderr.log').open('w') as err:
                try:run['returnCode']=subprocess.run(actual,stdout=out,stderr=err,timeout=min(args.timeout,remaining)).returncode
                except subprocess.TimeoutExpired:run['returnCode']=None;run['reason']='Predeclared resource budget exhausted; partial state retained'
            run['elapsedSeconds']=time.monotonic()-start;report['runs'].append(run);save()
            if run['returnCode']!=0:raise ValueError(f'n{steps} did not complete; original outputs retained')
            audit=trajectory.transient.verify(mesh,prefix,folder/'audit.json');trajectory.audit_final_checkpoint(mesh,prefix)
            if audit['controls']!=report['controls'] or not math.isclose(audit['time'],args.start_time+args.duration,rel_tol=0,abs_tol=1e-9):
                raise ValueError('Trial changed physics or end time')
            if len(audit['history'])!=steps or any(r['dt']!=dt for r in audit['history']):raise ValueError('Trial step history differs')
            run.update(valid=True,acceptedTime=audit['time'],forceX=audit['history'][-1]['forceX'],forceY=audit['history'][-1]['forceY'],
                auditSha256=n.sha256_file(folder/'audit.json'),checkpointSha256=n.sha256_file(folder/'flow.checkpoint'),
                artifactSha256={p.name:n.sha256_file(p) for p in folder.iterdir() if p.is_file()})
            save();print(json.dumps({k:run[k] for k in ('steps','dt','valid','elapsedSeconds','forceX','forceY')}),flush=True)
        geometry=n.read_cm2d(mesh);measured=n.measure(geometry,1e-11,1e-9)
        fields=[n.read_cells(root/f'n{steps}/flow.cells.csv',geometry,measured,'external') for steps in args.steps]
        weights=measured.areas;differences=[]
        for coarse,fine in zip(fields,fields[1:]):
            differences.append(dict(velocityL2=n.weighted_l2([math.hypot(a['u']-b['u'],a['v']-b['v']) for a,b in zip(coarse,fine)],weights),
                pressureL2=n.weighted_l2([a['p']-b['p'] for a,b in zip(coarse,fine)],weights),
                maxVelocityDifference=max(math.hypot(a['u']-b['u'],a['v']-b['v']) for a,b in zip(coarse,fine))))
        report['successiveDifferences']=differences
        report['observedDifferenceOrders']={key:math.log(differences[0][key]/differences[1][key])/math.log(args.steps[1]/args.steps[0])
            if differences[0][key]>0 and differences[1][key]>0 else None for key in differences[0]}
        # These are observations, not an automatically passed accuracy target.
        report['valid']=True;report['accuracyQualification']='not-qualified';save()
    except Exception as error:
        report['issues'].append(str(error));save();raise
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verified',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--start-time',type=float,default=250);p.add_argument('--duration',type=float,default=10)
    p.add_argument('--steps',type=int,nargs=3,default=[200,400,800]);p.add_argument('--timeout',type=float,default=900)
    p.add_argument('--total-timeout',type=float,default=2100);a=p.parse_args()
    if not(math.isfinite(a.start_time) and a.start_time>=0 and math.isfinite(a.duration) and a.duration>0 and
        0<a.timeout<=900 and 0<a.total_timeout<=2100 and 1<=a.steps[0]<a.steps[1]<a.steps[2]<=10000 and
        a.steps[1]**2==a.steps[0]*a.steps[2]):p.error('Finite positive duration and three increasing equal-ratio step counts are required; maximum budgets 900 s/run and 2100 s total')
    r=study(a);print(json.dumps(dict(valid=r['valid'],accuracyQualification=r['accuracyQualification'],observedDifferenceOrders=r['observedDifferenceOrders'])))
