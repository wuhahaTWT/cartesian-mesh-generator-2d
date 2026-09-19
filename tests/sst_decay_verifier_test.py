#!/usr/bin/env python3
"""Real native output readback and deliberate corruption of SST decay evidence."""
import argparse
import csv
import sys
import tempfile
import subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'/'verification'))
import verify_sst_decay as verifier


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--mesh-cli',type=Path,required=True)
    p.add_argument('--probe',type=Path,required=True)
    a=p.parse_args()
    with tempfile.TemporaryDirectory() as root:
        root=Path(root)
        # Generate a real final solver partition, not a hand-written mesh CSV.
        request=verifier.native.Request('sst-test','manufactured',3,1/6,.1,1.)
        command=verifier.native.mesh_command(a.mesh_cli,root/'mesh',request,root)
        run=subprocess.run(command,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30)
        if run.returncode: raise RuntimeError(run.stdout)
        mesh=root/'mesh.solver.cm2d';prefix=root/'decay'
        subprocess.run([str(a.probe),str(mesh),str(prefix),'.2','4'],check=True,capture_output=True,timeout=30)
        verifier.audit(mesh,prefix,.2,4)
        def tamper(suffix,column,value):
            path=Path(str(prefix)+suffix); original=path.read_text()
            with path.open() as f:
                rows=list(csv.DictReader(f))
            rows[0][column]=value
            with path.open('w') as f:
                writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
            try:
                try:verifier.audit(mesh,prefix,.2,4)
                except (ValueError,KeyError):return
                raise AssertionError(f'corrupt {column} accepted')
            finally:path.write_text(original)
        for suffix,column,value in [('.cells.csv','k','-1'),('.cells.csv','omega','nan'),
                                    ('.cells.csv','previousK','.5'),('.cells.csv','area','100'),
                                    ('.faces.csv','omegaFlux','.01'),('.history.csv','time','.7'),
                                    ('.history.csv','k','.3'),('.history.csv','maxNonlinearRateResidual','1')]:
            tamper(suffix,column,value)
        for dt,steps in [(0.,4),(.2,0),(.2,3)]:
            try:verifier.audit(mesh,prefix,dt,steps)
            except ValueError:continue
            raise AssertionError('invalid time controls accepted')
    print('SST homogeneous-decay independent readback/tamper tests passed')


if __name__=='__main__':main()
