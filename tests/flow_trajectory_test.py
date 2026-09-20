#!/usr/bin/env python3
"""Re-audit actual restarted outputs, archived segments and forged snapshots."""
import argparse,hashlib,json,shutil,subprocess,sys,tarfile,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_transient_flow as transient
import verify_flow_trajectory as trajectory
n=transient.native


def main(cli,mesh_cli):
    with tempfile.TemporaryDirectory(prefix='cm2d-trajectory-test-') as temporary:
        root=Path(temporary).resolve();mesh=root/'mesh.solver.cm2d'
        generated=subprocess.run(list(map(str,[mesh_cli,ROOT/'examples/acceptance/circle.xy',root/'mesh',8,.25,.15,
            'exterior',root/'foam',0,0,'--size-field','--far-field-spans',10,'--wall-relative-size',.3281250000003282,
            '--background-relative-size',1.3125000000013127,'--cells-per-level',22,'--max-safe-wall-level',11])),
            capture_output=True,text=True,timeout=30)
        assert generated.returncode==0,generated.stdout+generated.stderr
        initial=root/'initial';current=root/'initial.checkpoint'
        def command(prefix,target,restart=None):
            args=[str(cli),'--mesh',str(mesh),'--output',str(prefix),'--case','external','--nu','.1','--speed','1',
                '--convection','face-limited-linear','--pressure-preconditioner','aggregation','--tolerance','1e-8',
                '--max-iterations','200','--time-step','.005','--end-time',str(target),'--max-courant','1','--min-time-step','.00001','--max-time-steps','1000']
            if restart:args+=['--restart',str(restart)]
            r=subprocess.run(args,capture_output=True,text=True,timeout=20)
            assert r.returncode==0,r.stdout+r.stderr
            return args
        command(initial,.005)
        manifest=dict(format='cartmesh2d-bounded-trajectory-v1',complete=True,initialTime=.005,targetTime=.015,
            mesh=str(mesh),meshSha256=n.sha256_file(mesh),binary=str(cli),binarySha256=n.sha256_file(cli),
            initialCheckpointSha256=n.sha256_file(current),runs=[])
        for label,target in [('first',.01),('second',.015)]:
            directory=root/label;directory.mkdir();source=directory/'input.checkpoint';shutil.copyfile(current,source)
            prefix=directory/'flow';args=command(prefix,target,source);a=transient.verify(mesh,prefix,directory/'audit.json')
            current=directory/'flow.checkpoint'
            q=dict(valid=True,returnCode=0,targetTime=target,command=args,elapsedSeconds=0.,acceptedSteps=len(a['history']),
                restartSha256=n.sha256_file(source),checkpointSha256=n.sha256_file(current),auditSha256=n.sha256_file(directory/'audit.json'))
            (directory/'execution.json').write_text(json.dumps(q));manifest['runs'].append(q)
        path=root/'trajectory.json';path.write_text(json.dumps(manifest));output=root/'result'
        def verify():return trajectory.verify(path,initial,output)
        assert verify()['acceptedSteps']==3
        stopped=json.loads(json.dumps(manifest));stopped.update(complete=False,targetTime=.02,stopReason='Retained bounded-stop fixture')
        stopped['runs'].append(dict(valid=False,returnCode=None,targetTime=.02,reason='Unverified trailing attempt'))
        stopped_path=root/'stopped.json';stopped_path.write_text(json.dumps(stopped))
        try:trajectory.verify(stopped_path,initial,root/'stopped-result')
        except ValueError as e:assert 'completed' in str(e),str(e)
        else:raise AssertionError('stopped trajectory mislabeled as completed')
        prefix=trajectory.verify(stopped_path,initial,root/'stopped-result',completed_prefix=True)
        assert prefix['valid'] and not prefix['complete'] and prefix['endTime']==.015 and prefix['requestedEndTime']==.02
        assert prefix['unauditedTail'][0]['targetTime']==.02 and prefix['acceptedSteps']==3
        del stopped['stopReason'];stopped['complete']=None;stopped_path.write_text(json.dumps(stopped))
        try:trajectory.verify(stopped_path,initial,root/'stopped-result',completed_prefix=True)
        except ValueError as e:assert 'explicitly stopped' in str(e),str(e)
        else:raise AssertionError('actively changing run accepted as a stopped prefix')
        # A separately retained study must connect to the actual parent endpoint,
        # not merely carry a self-consistent hash for an unrelated parent file.
        parent=json.loads(json.dumps(manifest));parent['targetTime']=.01;parent['runs']=parent['runs'][:1]
        parent_path=root/'parent.json';parent_path.write_text(json.dumps(parent))
        child=json.loads(json.dumps(manifest));child['initialTime']=.01;child['runs']=child['runs'][1:]
        child['initialCheckpointSha256']=parent['runs'][-1]['checkpointSha256']
        child['parentTrajectory']=str(parent_path);child['parentTrajectorySha256']=n.sha256_file(parent_path)
        child_path=root/'child.json';child_path.write_text(json.dumps(child))
        assert trajectory.verify(child_path,None,root/'child-result')['acceptedSteps']==1
        parent['targetTime']=.009;parent_path.write_text(json.dumps(parent))
        child['parentTrajectorySha256']=n.sha256_file(parent_path);child_path.write_text(json.dumps(child))
        try:trajectory.verify(child_path,None,root/'child-result')
        except ValueError as e:assert 'parent final checkpoint/time/mesh/executable' in str(e),str(e)
        else:raise AssertionError('unrelated parent endpoint accepted')
        actual_origin=manifest['initialCheckpointSha256'];manifest['initialCheckpointSha256']='0'*64;path.write_text(json.dumps(manifest))
        try:verify()
        except ValueError as e:assert 'recorded restart origin' in str(e),str(e)
        else:raise AssertionError('unrelated initial state accepted')
        manifest['initialCheckpointSha256']=actual_origin;path.write_text(json.dumps(manifest))
        first=root/'first';archive=root/'first.tar.gz';files={str(p.relative_to(root)):n.sha256_file(p) for p in first.iterdir()}
        with tarfile.open(archive,'w:gz') as tar:tar.add(first,arcname='first')
        (root/'first.manifest.json').write_text(json.dumps(dict(verified=True,sha256=n.sha256_file(archive),files=files)))
        shutil.rmtree(first);assert verify()['valid']
        # A checkpoint changed together with its claimed hash still differs from the independently read field.
        checkpoint=root/'second/flow.checkpoint';original=checkpoint.read_text()
        lines=original.splitlines();position=next(i for i,line in enumerate(lines) if line.startswith('U '));words=lines[position].split();words[2]=str(float(words[2])+.1);lines[position]=' '.join(words);checkpoint.write_text('\n'.join(lines)+'\n')
        manifest['runs'][1]['checkpointSha256']=n.sha256_file(checkpoint)
        (root/'second/execution.json').write_text(json.dumps(manifest['runs'][1]));path.write_text(json.dumps(manifest))
        try:verify()
        except ValueError as e:assert 'U differs from final CSV' in str(e),str(e)
        else:raise AssertionError('forged final checkpoint accepted')
        checkpoint.write_text(original);manifest['runs'][1]['checkpointSha256']=n.sha256_file(checkpoint)
        (root/'second/execution.json').write_text(json.dumps(manifest['runs'][1]));path.write_text(json.dumps(manifest))
        assert verify()['valid']
        source=root/'second/input.checkpoint';source.write_text(source.read_text().replace('TIME 0.01','TIME 0.02'))
        try:verify()
        except ValueError as e:assert 'restart bytes' in str(e),str(e)
        else:raise AssertionError('broken restart chain accepted')
        # Archive bytes cannot be silently substituted even when no raw directory remains.
        archive.write_bytes(archive.read_bytes()+b'changed')
        try:trajectory.materialize(root,'first',root/'extract')
        except ValueError as e:assert 'archive digest' in str(e),str(e)
        else:raise AssertionError('changed archive accepted')
        print('Actual raw/archived trajectory, checkpoint-field binding and tamper rejection passed')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cli',type=Path,required=True);p.add_argument('--mesh-cli',type=Path,required=True)
    a=p.parse_args();main(a.cli.resolve(),a.mesh_cli.resolve())
