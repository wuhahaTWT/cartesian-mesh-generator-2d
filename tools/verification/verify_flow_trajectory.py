#!/usr/bin/env python3
"""Re-audit retained external-flow segments and byte-exact restart continuity.

Intermediate segment directories may be losslessly archived by the bounded
runner. Each archive and member is verified before a safe temporary extraction.
Only segment-final fields are reconstructed; intermediate monitors are records.
"""
import argparse,csv,hashlib,json,math,shlex,shutil,tarfile,tempfile
from pathlib import Path,PurePosixPath
import verify_transient_flow as transient

native=transient.native


def fail(message):raise native.VerificationError('trajectory: '+message)


def audit_final_checkpoint(mesh_path,prefix):
    """Bind the full retained v2 state to independent mesh geometry and final CSV."""
    mesh=native.read_cm2d(mesh_path);m=native.measure(mesh,1e-11,1e-9);g=native.face_geometry(mesh,m)
    summary=json.loads(Path(str(prefix)+'.json').read_text())
    cells=native.read_cells(Path(str(prefix)+'.cells.csv'),mesh,m,'external')
    faces,_=native.read_faces(Path(str(prefix)+'.faces.csv'),mesh)
    stream=iter(shlex.split(line) for line in Path(str(prefix)+'.checkpoint').read_text().splitlines() if line.strip())
    def row(label,n):
        fields=next(stream,[])
        if len(fields)!=n or fields[0]!=label:fail('malformed checkpoint '+label)
        return fields[1:]
    def equal(values,expected):
        if any(not native.close(native.finite(a,'checkpoint geometry'),b,1e-12,1e-10) for a,b in zip(values,expected)):
            fail('checkpoint geometry differs from CM2D')
    if row('CARTMESH2D_FLOW_CHECKPOINT',2)!=['2'] or row('DISCRETIZATION',2)!=['Euler-RC-v2']:
        fail('trajectory requires constant-viscosity preset checkpoint v2')
    config=row('CONFIG',8)
    if config!=['external',str(summary['nu']),str(summary['speed']),summary['convection'],summary['viscousStress'],'0',summary['outletBackflow']]:
        # Decimal formatting of nu/speed is not a different physical value.
        if (config[0],config[3],config[4],config[6])!=('external',summary['convection'],summary['viscousStress'],summary['outletBackflow']) or [native.finite(config[j],'configuration') for j in (1,2,5)]!=[summary['nu'],summary['speed'],0]:
            fail('checkpoint physical controls differ')
    if row('CELLS',2)!=[str(len(mesh.cells))]:fail('checkpoint cell count differs')
    for c in mesh.cells:
        fields=row('CELL',6+len(c.edges))
        if fields[0]!=str(c.id) or fields[4:]!=[str(len(c.edges)),*map(str,c.edges)]:fail('checkpoint cell topology differs')
        equal(fields[1:4],[*m.centroids[c.id],m.areas[c.id]])
    if row('FACES',2)!=[str(len(mesh.edges))]:fail('checkpoint face count differs')
    for e,f in zip(mesh.edges,g):
        fields=row('FACE',13)
        if fields[:4]!=[str(e.id),str(e.owner),str(e.neighbour) if e.neighbour>=0 else '-',str(e.patch)]:fail('checkpoint face topology differs')
        equal(fields[4:],[*f.centre,*f.area_vector,*f.correction,f.transmissibility,f.neighbour_weight])
    if native.finite(row('TIME',2)[0],'checkpoint time')!=summary['acceptedTime']:fail('checkpoint time differs from accepted field')
    for key,expected in [('U',[c['u'] for c in cells]),('V',[c['v'] for c in cells]),('P',[c['p'] for c in cells]),('FLUX',native.face_fluxes(faces))]:
        fields=row(key,len(expected)+2)
        if fields[0]!=str(len(expected)) or [native.finite(v,'checkpoint '+key) for v in fields[1:]]!=expected:
            fail('checkpoint '+key+' differs from final CSV')
    row('END',1)
    if next(stream,None) is not None:fail('trailing checkpoint data')


def materialize(root,name,destination):
    source=root/name
    if source.is_dir():return source
    archive=root/(name+'.tar.gz');manifest_path=root/(name+'.manifest.json')
    manifest=json.loads(manifest_path.read_text())
    if manifest.get('verified') is not True or native.sha256_file(archive)!=manifest.get('sha256'):
        fail('archive digest/retention manifest mismatch')
    files=manifest.get('files')
    if not isinstance(files,dict) or not files:fail('empty archive inventory')
    found=set()
    with tarfile.open(archive) as tar:
        for entry in tar:
            path=PurePosixPath(entry.name)
            if path.is_absolute() or '..' in path.parts or not path.parts or path.parts[0]!=name:
                fail('archive path escapes its segment')
            if entry.isdir():continue
            if not entry.isfile() or entry.name not in files or entry.name in found or entry.size>64*1024**2:
                fail('unsupported/duplicate/oversized archive member')
            data=tar.extractfile(entry).read()
            if hashlib.sha256(data).hexdigest()!=files[entry.name]:fail('archived source file digest mismatch')
            target=destination.joinpath(*path.parts);target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
            found.add(entry.name)
    if found!=set(files):fail('archive inventory does not cover all retained source files')
    return destination/name


def verify(trajectory,initial,output):
    report=json.loads(trajectory.read_text());root=trajectory.resolve().parent
    if report.get('format')!='cartmesh2d-bounded-trajectory-v1' or report.get('complete') is not True:
        fail('a completed bounded trajectory is required')
    if not isinstance(report.get('runs'),list) or not report['runs'] or not 0<=native.finite(report['initialTime'],'initial time')<native.finite(report['targetTime'],'target time'):
        fail('invalid run list or time interval')
    mesh=Path(report['mesh']);binary=Path(report['binary'])
    if native.sha256_file(mesh)!=report['meshSha256'] or native.sha256_file(binary)!=report['binarySha256']:
        fail('mesh or retained executable has changed')
    rows=[];segments=[];previous=report.get('initialCheckpointSha256');time=report['initialTime'];controls=None
    if previous is None and initial is None:fail('a recorded restart origin or an initial run is required')
    if report.get('parentTrajectory'):
        parent_path=Path(report['parentTrajectory'])
        if native.sha256_file(parent_path)!=report.get('parentTrajectorySha256'):fail('parent trajectory changed')
        parent=json.loads(parent_path.read_text())
        if (parent.get('format')!='cartmesh2d-bounded-trajectory-v1' or parent.get('complete') is not True or parent.get('targetTime')!=time or
            parent.get('meshSha256')!=report['meshSha256'] or parent.get('binarySha256')!=report['binarySha256'] or
            not parent.get('runs') or parent['runs'][-1].get('checkpointSha256')!=previous):
            fail('parent final checkpoint/time/mesh/executable does not connect to this trajectory')
    with tempfile.TemporaryDirectory(prefix='cm2d-trajectory-audit-') as temporary:
        tmp=Path(temporary)
        if initial:
            first=transient.verify(mesh,initial,tmp/'initial-audit.json')
            audit_final_checkpoint(mesh,initial)
            if first['case']!='external' or first['time']!=time or first['history'][0]['time']!=first['history'][0]['dt']:
                fail('initial run does not connect zero time to the first restart')
            controls=first['controls'];initial_sha=native.sha256_file(Path(str(initial)+'.checkpoint'))
            if previous is not None and previous!=initial_sha:fail('provided initial run is not the recorded restart origin')
            previous=initial_sha
            for r in first['history']:rows.append(dict(segment='initial',**r))
            segments.append(dict(label='initial',time=time,finalStateAudited=True,initialCondition=first['initialVortex']))
        for index,run in enumerate(report['runs']):
            target=native.finite(run['targetTime'],'target time')
            if not run.get('valid') or run.get('returnCode')!=0 or not target>time:
                fail('failed or unordered segment')
            command=run['command']
            if not isinstance(command,list) or not command or command[0]!=str(binary):fail('executable command identity differs')
            args=command[1:];flags={};i=0
            while i<len(args):
                key=args[i]
                if key in flags:fail('duplicate command option')
                if key=='--profile':flags[key]=True;i+=1
                else:
                    if i+1>=len(args):fail('truncated command')
                    flags[key]=args[i+1];i+=2
            if flags.get('--mesh')!=str(mesh) or flags.get('--case')!='external' or any(k.startswith('--initial-vortex') for k in flags):
                fail('segment changes mesh/case or reapplies initial condition')
            prefix_path=Path(flags['--output']);name=prefix_path.parent.name
            if prefix_path.parent.parent!=root or prefix_path.name!='flow' or Path(flags['--restart'])!=prefix_path.parent/'input.checkpoint':
                fail('recorded artifact paths are outside the segment')
            directory=materialize(root,name,tmp/f'segment-{index}')
            execution=json.loads((directory/'execution.json').read_text())
            if execution!=run:fail('execution record differs from trajectory manifest')
            source=directory/'input.checkpoint';final=directory/'flow.checkpoint'
            if native.sha256_file(source)!=run['restartSha256'] or (previous is not None and run['restartSha256']!=previous):
                fail('restart bytes differ from preceding accepted state')
            if native.sha256_file(final)!=run['checkpointSha256'] or native.sha256_file(directory/'audit.json')!=run['auditSha256']:
                fail('accepted checkpoint or original audit changed')
            audited=transient.verify(mesh,directory/'flow',tmp/f'reaudit-{index}.json')
            audit_final_checkpoint(mesh,directory/'flow')
            if audited['case']!='external' or audited['time']!=target or audited['initialVortex']:
                fail('wrong case/target or repeated initial condition')
            summary=json.loads((directory/'flow.json').read_text())
            expected={'--nu':summary['nu'],'--speed':summary['speed'],'--tolerance':summary['tolerance'],
                '--time-step':summary['maximumTimeStep'],'--end-time':target,'--max-courant':summary['targetCourant'],
                '--min-time-step':summary['minimumTimeStep'],'--max-time-steps':summary['maximumAcceptedSteps']}
            for key,value in expected.items():
                if native.finite(flags.get(key),'command '+key)!=value:fail('recorded command controls differ from actual output')
            for key,field in (('--convection','convection'),('--pressure-preconditioner','pressurePreconditioner')):
                if flags.get(key)!=summary[field]:fail('recorded numerical scheme differs')
            if summary.get('startTime')!=time or len(audited['history'])!=run['acceptedSteps']:
                fail('segment time origin/step count differs')
            if controls is not None and audited['controls']!=controls:fail('physical or discretization controls changed between segments')
            controls=audited['controls'];time=target;previous=run['checkpointSha256']
            for r in audited['history']:rows.append(dict(segment=name,**r))
            segments.append(dict(label=name,startTime=summary['startTime'],time=target,steps=len(audited['history']),
                elapsedSeconds=run['elapsedSeconds'],restartSha256=run['restartSha256'],checkpointSha256=previous,
                independentContinuity=audited['independentContinuity']['nativeDefinitionContinuity'],
                independentMomentumResidual=audited['momentumAudit']['cellResidual']['maxNormalized']))
            if directory!=root/name:shutil.rmtree(directory.parent)
        if time!=report['targetTime']:fail('final target not reached')
    for a,b in zip(rows,rows[1:]):
        if not math.isclose(a['time']+b['dt'],b['time'],rel_tol=1e-12,abs_tol=1e-14):fail('gap/overlap in physical time')
    output.parent.mkdir(parents=True,exist_ok=True)
    csv_path=Path(str(output)+'.time-history.csv')
    with csv_path.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=['globalStep',*rows[0]]);writer.writeheader()
        for i,row in enumerate(rows,1):writer.writerow(dict(globalStep=i,**row))
    result=dict(valid=True,scope='Each segment-final field independently audited; restart bytes identical. Intermediate monitor rows are recorded, not all reconstructed. No mesh/time independence or shedding qualification.',
        source=str(trajectory),sourceSha256=native.sha256_file(trajectory),meshSha256=report['meshSha256'],binarySha256=report['binarySha256'],
        startTime=rows[0]['time']-rows[0]['dt'],endTime=time,acceptedSteps=len(rows),controls=controls,segments=segments,
        historyFile=str(csv_path),historySha256=native.sha256_file(csv_path),finalForceX=rows[-1]['forceX'],finalForceY=rows[-1]['forceY'],
        initialRunPrefix=str(initial) if initial else None,
        parentTrajectory=report.get('parentTrajectory'),parentTrajectorySha256=report.get('parentTrajectorySha256'))
    Path(str(output)+'.json').write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--trajectory',type=Path,required=True)
    p.add_argument('--initial',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    try:r=verify(a.trajectory,a.initial,a.output)
    except (ValueError,OSError,KeyError,TypeError,tarfile.TarError) as error:
        failure=dict(valid=False,issues=[str(error)])
        a.output.parent.mkdir(parents=True,exist_ok=True);Path(str(a.output)+'.json').write_text(json.dumps(failure,indent=2)+'\n')
        print(json.dumps(failure));raise SystemExit(1)
    print(json.dumps({k:r[k] for k in ('valid','startTime','endTime','acceptedSteps')}))
