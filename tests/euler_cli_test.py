#!/usr/bin/env python3
import argparse,csv,json,math,re,subprocess,sys,tempfile,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
import verify_euler as euler
parser=argparse.ArgumentParser();parser.add_argument('--cli',required=True);args=parser.parse_args()
with tempfile.TemporaryDirectory(prefix='cm2d-euler-') as directory:
    root=Path(directory)
    def run(mesh,name,extra=(),expected=0):
        prefix=root/name
        command=[args.cli,'--mesh',str(mesh),'--output',str(prefix),'--end-time','.2','--gas-r','1',*map(str,extra)]
        result=subprocess.run(command,capture_output=True,text=True,timeout=30)
        assert result.returncode==expected,(command,result.stdout,result.stderr)
        return prefix
    previous=None
    for n in [50,100,200]:
        mesh=root/f'n{n}.solver.cm2d';euler.rectangle(mesh,n,4)
        prefix=run(mesh,f'sod{n}')
        audit=euler.audit(mesh,prefix);assert audit['valid']
        rows=list(csv.DictReader(Path(str(prefix)+'.cells.csv').open()))
        errors={key:sum(float(r['area'])*abs(float(r[key])-euler.sod(float(r['x']),.2)[k]) for r in rows)/.1 for key,k in [('rho',0),('u',1),('p',2)]}
        print(n,errors,flush=True)
        if previous:assert all(previous[k]/errors[k]>1.1 for k in errors),(previous,errors)
        previous=errors
    assert errors['rho']<.025 and errors['p']<.02 and errors['u']<.05,errors
    # Full versus step-budget interruption/restart follows identical accepted states.
    mesh=root/'n50.solver.cm2d';first=run(mesh,'budget',['--max-steps','20'],2)
    assert euler.audit(mesh,first)['valid']
    summary=json.loads(Path(str(first)+'.json').read_text());assert not summary['targetReached'] and summary['acceptedSteps']==20
    resumed=run(mesh,'resumed',['--restart',str(first)+'.checkpoint'])
    assert euler.audit(mesh,resumed)['valid']
    assert Path(str(resumed)+'.checkpoint').read_bytes()==(root/'sod50.checkpoint').read_bytes()
    # Serialized endpoint order is arbitrary: outward normals come from the owner polygon.
    reverse=root/'reverse.solver.cm2d';lines=mesh.read_text().splitlines();start=next(i for i,s in enumerate(lines) if s.startswith('EDGES '));size=int(lines[start].split()[1])
    for index in range(start+1,start+1+size):
        tokens=lines[index].split();tokens[1],tokens[2]=tokens[2],tokens[1];lines[index]=' '.join(tokens)
    reverse.write_text('\n'.join(lines)+'\n');flipped=run(reverse,'flipped')
    assert euler.audit(reverse,flipped)['valid']
    assert Path(str(flipped)+'.checkpoint').read_bytes()==(root/'sod50.checkpoint').read_bytes()
    for extra in [['--gamma','1'],['--density','-1'],['--cfl','.5'],['--u','1'],['--case','vortex'],['--restart',str(first)+'.checkpoint','--gamma','1.3']]:
        run(mesh,'bad',extra,1)
    for speed in [-5,.3,5]:
        p=run(mesh,'uniform'+str(speed),['--case','uniform','--u',speed,'--v','.17'])
        assert euler.audit(mesh,p)['valid']
        rows=list(csv.DictReader(Path(str(p)+'.cells.csv').open()))
        assert max(max(abs(float(r['u'])-speed),abs(float(r['v'])-.17),abs(float(r['p'])-1),abs(float(r['rho'])-1)) for r in rows)<1e-12
    boundary=root/'input-custom.boundaries'
    subprocess.run([args.cli,'--mesh',str(mesh),'--case','uniform','--u','.2','--export-boundaries',str(boundary)],check=True,capture_output=True)
    custom=run(mesh,'custom',['--case','custom','--u','.2','--boundary',boundary])
    assert euler.audit(mesh,custom)['valid']
    original=boundary.read_text();rows=original.splitlines();tokens=rows[2].split();tokens[2]=str(float(tokens[2])+.01);rows[2]=' '.join(tokens);boundary.write_text('\n'.join(rows)+'\n')
    run(mesh,'bad-custom',['--case','custom','--boundary',boundary],1)
    boundary.write_text(original)
    cancelled=root/'cancelled'
    process=subprocess.Popen([args.cli,'--mesh',str(mesh),'--output',str(cancelled),'--case','uniform','--gas-r','1','--end-time','100'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    deadline=time.monotonic()+5;accepted=False
    try:
        while process.poll() is None and time.monotonic()<deadline:
            checkpoint=Path(str(cancelled)+'.checkpoint')
            if checkpoint.exists():
                match=re.search(r'^STATE (\S+) (\d+)$',checkpoint.read_text(),re.M)
                if match and float(match[1])>0:accepted=True;break
            time.sleep(.005)
        assert accepted,'no accepted state before cancellation'
        process.terminate();code=process.wait(timeout=5)
        if sys.platform=='win32':
            # TerminateProcess cannot invoke the C++ SIGTERM handler. Windows
            # preserves the last atomic checkpoint, not a fresh final summary.
            assert code!=0
        else:
            assert code==2
            assert json.loads(Path(str(cancelled)+'.json').read_text())['status']=='cancelled'
            assert euler.audit(mesh,cancelled)['valid']
        # Exercise the accepted-checkpoint recovery on every platform, including
        # the hard-stop path actually used by the Windows desktop.
        match=re.search(r'^STATE (\S+) (\d+)$',checkpoint.read_text(),re.M)
        assert match and float(match[1])>0
        restart_time=float(match[1])
        resumed_cancel=run(mesh,'cancel-resumed',['--case','uniform',
            '--restart',checkpoint,'--end-time',str(restart_time+.01)])
        assert euler.audit(mesh,resumed_cancel)['valid']
        state=json.loads(Path(str(resumed_cancel)+'.json').read_text())
        assert state['targetReached']
        rows=list(csv.DictReader(Path(str(resumed_cancel)+'.cells.csv').open()))
        assert max(max(abs(float(r['rho'])-1),abs(float(r['p'])-1),
            abs(float(r['u'])),abs(float(r['v']))) for r in rows)<1e-12
    finally:
        if process.poll() is None:process.kill();process.wait()
    mesh=root/'periodic.solver.cm2d';euler.rectangle(mesh,32,32,20,20)
    prefix=run(mesh,'vortex',['--case','vortex','--end-time','.1'])
    assert euler.audit(mesh,prefix)['valid']
    # The independent verifier must detect an altered physical field and flux.
    for suffix,field in [('.cells.csv','p'),('.faces.csv','energy')]:
        path=Path(str(prefix)+suffix);original=path.read_text();rows=list(csv.DictReader(original.splitlines()));rows[10][field]=str(float(rows[10][field])+.1)
        with path.open('w') as out:
            writer=csv.DictWriter(out,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
        rejected=False
        try:euler.audit(mesh,prefix)
        except ValueError:rejected=True
        assert rejected,'independent audit accepted tampering'
        path.write_text(original)
    print('Euler exact Sod refinement, independent conservation, positivity, restart and mutation checks passed')
