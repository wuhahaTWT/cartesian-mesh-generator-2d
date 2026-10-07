#!/usr/bin/env python3
"""Native black-box state transactions, events and controller restart.
Reads published native fields; no independent equation/topology audit backend.
"""
import argparse,csv,hashlib,json,math,os,signal,subprocess,sys,time,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'/'thermal'))
from workflow import run,configure

def vortex_initial_condition(a,root,mesh,flow,bc):
    common=[a.cli,'--mesh',mesh,'--evolve-flow','custom','--flow-boundary',flow,'--boundary',bc,
        '--initial','300','--flow-speed','.2','--flow-nu','.1','--diffusivity','.1',
        '--dt','.01','--flow-convection','limited-linear','--convection','limited-linear',
        '--flux-correction','bounded-spatial']
    seed=['--initial-vortex-x','.5','--initial-vortex-y','.5','--initial-vortex-radius','.2','--initial-vortex-speed','.05']
    def invoke(label,steps=4,extra=(),code=0):
        out=root/label;r=run([*common,'--output',out,'--steps',str(steps),*extra],root/(label+'.log'))
        assert r['code']==code,(label,r,(root/(label+'.log')).read_text()[-2000:])
        return out
    full=invoke('seed-full',extra=seed);first=invoke('seed-first',steps=2,extra=seed)
    initial=full.with_suffix('.initial.checkpoint');joint=first.with_suffix('.thermal.checkpoint');snapshot=joint.read_bytes()
    assert '\nTIME 0\n' in initial.read_text()
    assert initial.read_bytes()==first.with_suffix('.initial.checkpoint').read_bytes()
    summary=json.loads(full.with_suffix('.json').read_text())
    assert summary['initialVortex']=={'definition':'compact-cubic-v1','centre':[.5,.5],'radius':.2,'peakSpeed':.05,'checkpointSuffix':'.initial.checkpoint'}
    resumed=invoke('seed-resumed',steps=2,extra=['--restart',joint])
    for suffix in ['.thermal.checkpoint','.carrier.checkpoint','.cells.csv','.faces.csv']:
        assert full.with_suffix(suffix).read_bytes()==resumed.with_suffix(suffix).read_bytes(),suffix
    assert not resumed.with_suffix('.initial.checkpoint').exists()
    assert 'initialVortex' not in json.loads(resumed.with_suffix('.json').read_text())
    zero=invoke('seed-zero')
    assert full.with_suffix('.carrier.checkpoint').read_bytes()!=zero.with_suffix('.carrier.checkpoint').read_bytes()
    def progress(label):
        data=[json.loads(line) for line in (root/(label+'.log')).read_text().splitlines() if line.startswith('{')]
        return [{k:v for k,v in r.items() if k!='step'} for r in data if r.get('type')=='thermal-time-step']
    observed=progress('seed-full');assert len(observed)==4
    assert observed==progress('seed-first')+progress('seed-resumed')
    assert all(math.isfinite(r[key]) for r in observed for key in ['forceX','forceY'])
    failed=invoke('seed-first-failed',extra=[*seed,'--flow-max-iterations','1'],code=2)
    assert failed.with_suffix('.initial.checkpoint').exists() and not failed.with_suffix('.thermal.checkpoint').exists()
    for label,args in [('seed-repeat',[*seed,'--restart',joint]),('seed-partial',seed[:2]),
                       ('seed-duplicate',seed+seed[:2]),('seed-radius',seed[:5]+['-1']+seed[6:]),
                       ('seed-outside',['--initial-vortex-x','-1',*seed[2:]]),
                       ('seed-as-restart',['--restart',initial])]:
        out=invoke(label,extra=args,code=1);assert not out.with_suffix('.thermal.checkpoint').exists()
    assert joint.read_bytes()==snapshot

def main(a):
    base=Path(a.output).resolve();base.mkdir(parents=True,exist_ok=True)
    root=Path(tempfile.mkdtemp(prefix='run-',dir=base))
    geom=root/'square.xy';geom.write_text('0 0\n1 0\n1 1\n0 1\n')
    assert run([a.mesh_cli,geom,root/'mesh','4',str(1/14),'.1','interior',root/'openfoam','4','0'],root/'mesh.log')['code']==0
    mesh=root/'mesh.solver.cm2d';flow,bc,event,_=configure(mesh,'cavity',root)
    vortex_initial_condition(a,root,mesh,flow,bc)
    common=[a.cli,'--mesh',mesh,'--evolve-flow','custom','--flow-boundary',flow,'--boundary',bc,
        '--thermal-events',event,'--initial','300','--flow-speed','.2','--flow-nu','.1','--diffusivity','.1',
        '--dt','.1','--end-time','6','--min-dt','.000001','--flow-convection','limited-linear','--convection','limited-linear']
    def invoke(label,extra=(),code=0):
        out=root/label;r=run([*common,'--output',out,*extra],root/(label+'.log'))
        assert r['code']==code,(label,r,(root/(label+'.log')).read_text()[-2000:])
        return out
    continuous=invoke('continuous')
    check=invoke('check-only',['--restart',continuous.with_suffix('.thermal.checkpoint'),'--check-restart','on'])
    checked=json.loads((root/'check-only.log').read_text())
    assert checked['status']=='valid' and checked['time']==6 and not check.with_suffix('.json').exists()
    bad=root/'truncated.thermal.checkpoint';bad.write_text(continuous.with_suffix('.thermal.checkpoint').read_text()[:-12])
    invoke('check-truncated',['--restart',bad,'--check-restart','on'],1)
    # Frozen transport must accept the custom-boundary checkpoint written by
    # this same native coupled run, including its exact carrier flux and clock.
    carrier=continuous.with_suffix('.carrier.checkpoint')
    frozen=root/'frozen-custom'
    cmd=[a.cli,'--mesh',mesh,'--flow-checkpoint',carrier,'--boundary',bc,'--initial','300',
         '--diffusivity','.1','--dt','.01','--steps','2','--output',frozen]
    def frozen_run(label,text=None,code=0):
        args=list(cmd);out=root/label;args[args.index('--output')+1]=out
        if text is not None:
            file=root/(label+'.input.checkpoint');file.write_text(text)
            args[args.index('--flow-checkpoint')+1]=file
        result=run(args,root/(label+'.log'))
        assert result['code']==code,(label,result,(root/(label+'.log')).read_text()[-1000:])
        return out
    frozen_run('frozen-custom')
    source=carrier.read_text();assert source.startswith('CARTMESH2D_FLOW_CHECKPOINT 4')
    frozen_run('frozen-crlf',source.replace('\n','\r\n'))
    frozen_summary=json.loads(frozen.with_suffix('.json').read_text())
    assert frozen_summary['carrierTime']==6 and not frozen_summary['evolvingFlow']
    field_rows=lambda prefix:list(csv.DictReader(prefix.with_suffix('.faces.csv').open()))
    assert [r['volumeFlux'] for r in field_rows(frozen)]==[r['volumeFlux'] for r in field_rows(continuous)]
    frozen_run('frozen-bad-kind',source.replace(' moving-wall ', ' invalid-kind ',1),1)
    frozen_run('frozen-truncated',source[:source.rfind('END')],1)
    lines=source.splitlines();index=next(i for i,line in enumerate(lines) if line.startswith('CELL '))
    fields=lines[index].split();fields[2]=str(float(fields[2])+.001);lines[index]=' '.join(fields)
    frozen_run('frozen-wrong-mesh','\n'.join(lines)+'\n',1)
    published=continuous.with_suffix('.thermal.checkpoint').read_bytes()
    published_summary=continuous.with_suffix('.json').read_bytes()
    invoke('continuous',code=1)
    assert continuous.with_suffix('.thermal.checkpoint').read_bytes()==published
    assert continuous.with_suffix('.json').read_bytes()==published_summary
    part=invoke('part',['--max-time-steps','17'],1)
    joint=part.with_suffix('.thermal.checkpoint');assert joint.exists()
    snapshot=joint.read_bytes();split=invoke('split',['--restart',joint])
    assert snapshot==joint.read_bytes()
    assert continuous.with_suffix('.thermal.checkpoint').read_bytes()==split.with_suffix('.thermal.checkpoint').read_bytes()
    assert continuous.with_suffix('.cells.csv').read_bytes()==split.with_suffix('.cells.csv').read_bytes()
    rows=lambda p:list(csv.DictReader(p.open()))
    h=rows(continuous.with_suffix('.thermal-history.csv'));attempts=rows(continuous.with_suffix('.attempt-history.csv'))
    assert all(r['accepted']=='1' for r in h)
    for t in (1.37,2.43,3.23):assert t in [float(r['time']) for r in h]
    assert float(h[-1]['time'])==6
    # Actual exported native budget integrated over accepted steps, not rebuilt equations.
    heat=rows(continuous.with_suffix('.heat-history.csv'))
    integrated=sum(float(r['dt'])*(float(r['sourceIntegral'])-float(r['boundaryFlux'])) for r in heat)
    actual=float(heat[-1]['heatContent'])-300
    defect=sum(float(r['dt'])*float(r['globalBalance']) for r in heat)
    assert abs(actual-integrated-defect)<1e-9,(actual,integrated,defect)
    assert abs(defect)<1e-5
    # Changed full event law cannot silently resume.
    original=event.read_text();event.write_text(original.replace('1.37,lid,value,300,300','1.37,lid,value,300.1,300'))
    invoke('bad-event',['--restart',joint],1);event.write_text(original)
    assert snapshot==joint.read_bytes()
    exhausted=invoke('exhausted',['--time-error','on','--time-rtol','.000001','--max-step-retries','0'],2)
    assert not exhausted.with_suffix('.thermal.checkpoint').exists()
    assert len(rows(exhausted.with_suffix('.thermal-history.csv')))==0
    assert rows(exhausted.with_suffix('.attempt-history.csv'))[-1]['reason']=='time-error'
    failure=json.loads(exhausted.with_suffix('.json').read_text())
    assert failure['controller']['reason']=='time-error' and failure['controller']['errorRatio']>1
    assert failure['controller']['attempts']==1 and failure['controller']['maximumRetries']==0
    minimum=invoke('minimum',['--time-error','on','--min-dt','.1'],2)
    failure=json.loads(minimum.with_suffix('.json').read_text())
    assert failure['acceptedTime']==0 and not minimum.with_suffix('.thermal.checkpoint').exists()
    assert failure['controller']['timeStep']==failure['controller']['minimumTimeStep']==.1
    assert failure['controller']['reason']=='time-error' and failure['controller']['errorRatio']>1
    scalar=invoke('scalar-fail',['--max-corrections','1','--max-step-retries','0'],2)
    assert not scalar.with_suffix('.thermal.checkpoint').exists()
    invoke('nan',['--temperature-scale','nan','--restart',joint],1)
    # Changed maximum step is supported and reaches exactly the same physical target.
    changed=invoke('changed-step',['--restart',joint,'--dt','.05'])
    meta=json.loads(changed.with_suffix('.json').read_text());assert meta['acceptedTime']==6
    # Cooperative cancellation while solving: saved accepted transaction survives.
    live=root/'cancelled';log=(root/'cancelled.log').open('w')
    proc=subprocess.Popen([str(x) for x in [*common,'--output',live,'--end-time','1000']],stdout=log,stderr=subprocess.STDOUT)
    checkpoint=live.with_suffix('.thermal.checkpoint');deadline=time.monotonic()+10
    while time.monotonic()<deadline and not checkpoint.exists() and proc.poll() is None:time.sleep(.01)
    assert checkpoint.exists() and proc.poll() is None
    proc.send_signal(signal.SIGTERM);proc.wait(timeout=10);log.close();assert proc.returncode!=0
    saved=checkpoint.read_bytes();resumed=invoke('after-cancel',['--restart',checkpoint]);assert saved==checkpoint.read_bytes()
    # State comparison uses matching accepted timestamps; interruption has no hidden history.
    assert resumed.with_suffix('.thermal.checkpoint').read_bytes()==continuous.with_suffix('.thermal.checkpoint').read_bytes()
    replaced_path_checked=False
    if os.name=='posix':
        # Reproduce replacement of a live history file: writes to the old
        # unlinked inode otherwise succeed while the published path stays stale.
        replaced=root/'replaced-output'
        command=[str(x) for x in [*common,'--output',replaced,'--end-time','1000']]
        (root/'replaced-output.command.json').write_text(json.dumps(command,indent=2))
        proc=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        first=proc.stdout.readline();assert 'thermal-time-step' in first,first
        proc.send_signal(signal.SIGSTOP)
        os.waitpid(proc.pid,os.WUNTRACED)
        victim=replaced.with_suffix('.thermal-history.csv')
        old_stat=victim.stat();original_bytes=victim.read_bytes()
        copy=victim.with_suffix('.replacement');copy.write_bytes(original_bytes);copy.replace(victim)
        assert victim.stat().st_ino!=old_stat.st_ino
        proc.send_signal(signal.SIGCONT)
        rest=proc.communicate(timeout=10)[0]
        (root/'replaced-output.log').write_text(first+rest)
        assert proc.returncode==1 and 'output path replaced or truncated' in rest,rest[-2000:]
        assert victim.read_bytes()==original_bytes
        saved=replaced.with_suffix('.thermal.checkpoint').read_bytes()
        recovery=invoke('after-replacement',['--restart',replaced.with_suffix('.thermal.checkpoint')])
        assert replaced.with_suffix('.thermal.checkpoint').read_bytes()==saved
        assert recovery.with_suffix('.thermal.checkpoint').read_bytes()==continuous.with_suffix('.thermal.checkpoint').read_bytes()
        replaced_path_checked=True
    # Each bounded algorithm is explicit in v4 and cannot silently restart
    # with another operator. Complete event law and accepted controller survive.
    bounded_methods=[]
    for mode in ['bounded','bounded-spatial']:
        controls=['--flux-correction',mode,'--time-error','on']
        whole=invoke(mode+'-continuous',controls)
        partial=invoke(mode+'-partial',[*controls,'--max-time-steps','2'],1)
        saved=partial.with_suffix('.thermal.checkpoint').read_bytes()
        assert saved.startswith(b'CARTMESH2D_THERMAL_CHECKPOINT 4\n')
        config=next(line.split() for line in saved.decode().splitlines() if line.startswith('THERMAL_CONFIG '))
        assert config[-1]==mode
        resumed=invoke(mode+'-resumed',[*controls,'--restart',partial.with_suffix('.thermal.checkpoint')])
        assert partial.with_suffix('.thermal.checkpoint').read_bytes()==saved
        assert whole.with_suffix('.thermal.checkpoint').read_bytes()==resumed.with_suffix('.thermal.checkpoint').read_bytes()
        assert whole.with_suffix('.cells.csv').read_bytes()==resumed.with_suffix('.cells.csv').read_bytes()
        assert json.loads(whole.with_suffix('.json').read_text())['fluxCorrection']==mode
        for other in ['unrestricted','bounded','bounded-spatial']:
            if other==mode:continue
            rejected=invoke(mode+'-as-'+other,['--flux-correction',other,'--restart',partial.with_suffix('.thermal.checkpoint')],1)
            assert not rejected.with_suffix('.thermal.checkpoint').exists()
        assert partial.with_suffix('.thermal.checkpoint').read_bytes()==saved
        bounded_methods.append(mode)
    report={'readOnlyNativeCheckpointValidation':True,'frozenCustomCarrierExact':True,'frozenCrLfSupported':True,'frozenMalformedRejected':True,'boundedMethodsRestartedExactly':bounded_methods,'passed':True,'outputDirectory':str(root),'existingAcceptedOutputNotOverwritten':True,'continuousSteps':len(h),'splitSteps':17,'eventTimes':[1.37,2.43,3.23],
            'heatGain':actual,'integratedHeatGain':integrated,'integratedBudgetDefect':defect,
            'continuousSplitCheckpointIdentical':True,'cancelResumeCheckpointIdentical':True,
            'changedStepTime':meta['acceptedTime'],'replacedLiveOutputFailsClosedAndResumesIdentically':replaced_path_checked,
            'continuousCheckpointSha256':hashlib.sha256(continuous.with_suffix('.thermal.checkpoint').read_bytes()).hexdigest()}
    (root/'evidence.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cli',required=True);p.add_argument('--mesh-cli',required=True);p.add_argument('--output',required=True);main(p.parse_args())
