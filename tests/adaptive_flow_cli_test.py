#!/usr/bin/env python3
"""Exercise rollback, bounded retries, exact resumed trajectory and forged records."""
import argparse,csv,json,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_rotating_annulus as annulus
import verify_transient_flow as audit

def main(cli):
    with tempfile.TemporaryDirectory(prefix='cm2d-adaptive-') as directory:
        root=Path(directory);mesh=root/'ring.solver.cm2d';bc=root/'input.boundaries'
        annulus.write_polar_mesh(mesh,8,64)
        def command(args,success=True):
            r=subprocess.run(list(map(str,args)),capture_output=True,text=True,timeout=30)
            assert (r.returncode==0)==success,r.stdout+r.stderr
            return r
        command([cli,'--mesh',mesh,'--case','annulus','--speed',.5,'--export-boundaries',bc])
        def solve(name,extra=(),success=True):
            prefix=root/name
            r=command([cli,'--mesh',mesh,'--case','custom','--boundary',bc,'--nu',.1,'--speed',.5,
                '--convection','face-limited-linear','--pressure-preconditioner','aggregation','--tolerance',1e-8,
                '--max-iterations',80,'--time-step',.04,'--end-time',.08,'--max-courant',.15,'--output',prefix,*extra],success)
            return prefix,r
        def check(prefix):return audit.verify(mesh,prefix,root/(prefix.name+'-audit.json'))
        def checkpoint(prefix):return Path(str(prefix)+'.checkpoint').read_bytes()
        whole,_=solve('whole');assert check(whole)['valid']
        q=json.loads(Path(str(whole)+'.json').read_text())
        assert q['rejectedSteps']>=1 and q['completedSteps']>2 and q['time']==.08
        attempts=list(csv.DictReader(open(str(whole)+'.attempt-history.csv')))
        assert attempts[0]['reason']=='courant'
        first,stop=solve('first',['--max-time-steps',2],False)
        assert 'accepted-step budget' in stop.stderr
        assert json.loads(Path(str(first)+'.json').read_text())['status']=='failed'
        resumed,_=solve('resumed',['--restart',str(first)+'.checkpoint'])
        assert check(resumed)['valid'] and checkpoint(whole)==checkpoint(resumed)
        failed,stop=solve('failed',['--max-step-retries',0],False)
        assert 'retry budget' in stop.stderr
        failed_again,_=solve('failed-again',['--max-step-retries',0,'--restart',str(failed)+'.checkpoint'],False)
        assert checkpoint(failed)==checkpoint(failed_again)
        restored,_=solve('restored',['--restart',str(failed)+'.checkpoint'])
        assert checkpoint(restored)==checkpoint(whole)
        # Inner convergence failure, distinct from a converged but CFL-rejected trial.
        inner,_=solve('inner',['--max-iterations',45,'--end-time',.04,'--max-courant',1])
        assert check(inner)['valid']
        assert 'nonconverged' in Path(str(inner)+'.attempt-history.csv').read_text()
        # Minimum step exhaustion preserves the initial accepted state as well.
        floor,_=solve('floor',['--min-time-step',.04],False)
        assert checkpoint(floor)==checkpoint(failed)
        for args in (['--steps',2],['--min-time-step',.05],['--max-courant',0],['--max-step-retries',31],
                     ['--end-time',.08,'--restart',str(whole)+'.checkpoint']):
            solve('bad',args,False)
        # The independent reader must reject a deleted rejection and a fictitious accepted time.
        path=Path(str(whole)+'.attempt-history.csv');original=path.read_text()
        for bad in ('\n'.join(original.splitlines()[0:1]+original.splitlines()[2:])+'\n',
                    original.replace(',courant,',',accepted,',1),
                    original.replace('2,1,0,0.02','2,1,0.01,0.02',1)):
            path.write_text(bad)
            try:check(whole)
            except (ValueError,RuntimeError):pass
            else:raise AssertionError('forged adaptive attempt history accepted')
        path.write_text(original)
        assert check(whole)['valid']
        print(json.dumps(dict(valid=True,continuousRestartIdentical=True,courantRetries=q['rejectedSteps'],
            acceptedSteps=q['completedSteps'],innerRetry=True,minimumStepRollback=True)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cli',type=Path,required=True);main(p.parse_args().cli.resolve())
