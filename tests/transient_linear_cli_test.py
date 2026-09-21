#!/usr/bin/env python3
"""Every accepted physical step uses strict certification; numerical policy may change on restart."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from transient_flow_cli_test import rectangle
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
import verify_transient_flow as audit


def main(cli):
    with tempfile.TemporaryDirectory(prefix='cm2d-transient-linear-') as directory:
        root=Path(directory);mesh=root/'taylor.solver.cm2d';rectangle(mesh,16)
        def solve(name,policy,steps=3,restart=None,iterations=300,success=True):
            prefix=root/name
            command=[str(cli),'--mesh',str(mesh),'--output',str(prefix),'--case','taylor-green',
                '--nu','.01','--speed','1','--tolerance','1e-10','--max-iterations',str(iterations),
                '--time-step','.02','--steps',str(steps),'--convection','face-limited-linear',
                '--pressure-preconditioner','aggregation','--linear-policy',policy]
            if restart:command+=['--restart',str(restart)+'.checkpoint']
            result=subprocess.run(command,capture_output=True,text=True,timeout=30)
            assert (result.returncode==0)==success,result.stderr
            summary=json.loads(prefix.with_suffix('.json').read_text())
            if success:
                assert summary['converged'] and summary['strictLinearFinal']
                assert summary['strictAcceptedSteps']==summary['completedSteps']==steps
                assert summary['adaptiveLinear']==(policy=='adaptive')
                assert (summary['adaptiveLinearAttemptIterations']>0)==(policy=='adaptive')
                assert audit.verify(mesh,prefix,root/(name+'-audit.json'))['valid']
            return prefix,summary
        strict,_=solve('strict','strict');adaptive,_=solve('adaptive','adaptive')
        first,_=solve('first','adaptive',1)
        resumed,_=solve('resumed','adaptive',2,first)
        assert adaptive.with_suffix('.checkpoint').read_bytes()==resumed.with_suffix('.checkpoint').read_bytes()
        switched,_=solve('switched','strict',2,first)
        rows=lambda p:list(csv.DictReader(p.with_suffix('.cells.csv').open()))
        for other in (adaptive,switched):
            for field in ('u','v','p'):
                difference=max(abs(float(a[field])-float(b[field])) for a,b in zip(rows(strict),rows(other)))
                assert difference<1e-7,(field,difference)
        failed,q=solve('failed','adaptive',1,first,iterations=1,success=False)
        assert q['completedSteps']==q['strictAcceptedSteps']==0 and not q['converged']
        assert failed.with_suffix('.checkpoint').read_bytes()==first.with_suffix('.checkpoint').read_bytes()
        restored,_=solve('restored','adaptive',2,failed)
        assert restored.with_suffix('.checkpoint').read_bytes()==adaptive.with_suffix('.checkpoint').read_bytes()
    print('Transient adaptive linear: independent time equations, strict accepted steps, policy switch and failure/restart pass')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cli',type=Path,required=True);main(p.parse_args().cli.resolve())
