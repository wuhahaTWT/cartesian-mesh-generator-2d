#!/usr/bin/env python3
"""Initial vortex geometry, discrete conservation, restart and tamper rejection."""
import argparse,json,subprocess,sys,tempfile
from pathlib import Path
from transient_flow_cli_test import rectangle
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_transient_flow as audit
import verify_rotating_annulus as annulus


def main(cli):
    with tempfile.TemporaryDirectory(prefix='cm2d-vortex-') as directory:
        root=Path(directory);mesh=root/'box.solver.cm2d';rectangle(mesh)
        seed=['--initial-vortex-x',.5,'--initial-vortex-y',.5,'--initial-vortex-radius',.3,'--initial-vortex-speed',.2]
        def solve(name,extra=(),seeded=True,success=True,case='cavity',grid=mesh):
            prefix=root/name
            command=[cli,'--mesh',grid,'--output',prefix,'--case',case,'--nu',.1,'--speed',1,
                '--time-step',.002,'--steps',2,'--tolerance',1e-9,'--max-iterations',200,
                '--convection','face-limited-linear','--pressure-preconditioner','aggregation',
                *(seed if seeded else []),*extra]
            r=subprocess.run(list(map(str,command)),text=True,capture_output=True,timeout=30)
            assert (r.returncode==0)==success,r.stderr
            return prefix,r
        def check(p,grid=mesh):return audit.verify(grid,p,Path(str(p)+'.audit.json'))
        def checkpoint(p):return Path(str(p)+'.checkpoint').read_bytes()
        whole,_=solve('whole');report=check(whole)
        assert report['initialVortex']['valid'] and report['initialVortex']['maximumCellFluxImbalance']<1e-16
        assert .19<report['initialVortex']['sampledPeakSpeed']<=.2
        first,_=solve('first',['--steps',1]);assert check(first)['valid']
        resumed,_=solve('resumed',['--steps',1,'--restart',str(first)+'.checkpoint'],seeded=False)
        assert check(resumed)['valid'] and checkpoint(whole)==checkpoint(resumed)
        _,bad=solve('double-seed',['--restart',str(first)+'.checkpoint'],success=False)
        assert 'fresh physical transient' in bad.stderr
        for extra in (seed[:2],seed+seed[:2],seed[:-1]+[float('nan')],
                      ['--initial-vortex-x',.5,'--initial-vortex-y',.5,'--initial-vortex-radius',.6,'--initial-vortex-speed',.2],
                      ['--initial-vortex-x',2,'--initial-vortex-y',.5,'--initial-vortex-radius',.1,'--initial-vortex-speed',.2]):
            solve('invalid',extra,seeded=False,success=False)
        solve('verification',case='taylor-green',success=False)
        _,tiny=solve('unresolved',['--initial-vortex-x',.5,'--initial-vortex-y',.5,
            '--initial-vortex-radius',1e-8,'--initial-vortex-speed',.2],seeded=False,success=False)
        assert 'unresolved' in tiny.stderr
        # Clockwise sign changes only the initial U/V/FLUX, not pressure or metadata.
        opposite,_=solve('opposite',seed[:-1]+[-.2],seeded=False)
        assert check(opposite)['valid']
        def vectors(p):
            return {line.split()[0]:list(map(float,line.split()[2:])) for line in Path(str(p)+'.initial.checkpoint').read_text().splitlines() if line.split()[0] in ('U','V','P','FLUX')}
        a,b=vectors(whole),vectors(opposite)
        assert all(x==-y for k in ('U','V','FLUX') for x,y in zip(a[k],b[k])) and a['P']==b['P']
        # Closed multiply connected fluid: an interior seed works, solid-hole centre fails.
        ring=root/'ring.solver.cm2d';annulus.write_polar_mesh(ring,8,64);bc=root/'ring.boundaries'
        r=subprocess.run(list(map(str,[cli,'--mesh',ring,'--case','annulus','--speed',.5,'--export-boundaries',bc])),capture_output=True,text=True,timeout=10)
        assert r.returncode==0,r.stderr
        ringseed=['--boundary',bc,'--initial-vortex-x',.75,'--initial-vortex-y',0,'--initial-vortex-radius',.1,'--initial-vortex-speed',.05]
        p,_=solve('ring',ringseed,False,case='custom',grid=ring);assert check(p,ring)['valid']
        for centre in (0,.52,1.2):
            args=ringseed.copy();args[3]=centre
            p,r=solve('ring-invalid',args,False,False,case='custom',grid=ring)
            assert 'support disk' in r.stderr,r.stderr
        # Tampering with seed metadata, checkpoint geometry or actual state must be rejected.
        summary_path=Path(str(whole)+'.json');original=summary_path.read_text();value=json.loads(original)
        value['initialVortex']['peakSpeed']*=2;summary_path.write_text(json.dumps(value))
        try:check(whole)
        except (ValueError,RuntimeError):pass
        else:raise AssertionError('forged seed metadata accepted')
        summary_path.write_text(original)
        path=Path(str(whole)+'.initial.checkpoint');original=path.read_text()
        for bad in (original.replace('TIME 0','TIME 1'),original.replace('CELL 0 ','CELL 1 ',1),
                    original.replace('FLUX 420 0','FLUX 420 0.1',1),original+'garbage\n'):
            assert bad!=original
            path.write_text(bad)
            try:check(whole)
            except (ValueError,RuntimeError):pass
            else:raise AssertionError('forged initial checkpoint accepted')
        path.write_text(original);assert check(whole)['valid']
        print(json.dumps(dict(valid=True,restartByteIdentical=True,initialAudit=report['initialVortex'],solidAndBoundaryRejection=True)))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--cli',type=Path,required=True)
    main(parser.parse_args().cli.resolve())
