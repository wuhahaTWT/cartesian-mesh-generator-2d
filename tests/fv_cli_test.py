#!/usr/bin/env python3
"""Black-box native scalar CLI: outputs, honest nonconvergence, bad inputs."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile

p=argparse.ArgumentParser();p.add_argument('--cli',required=True);args=p.parse_args()
cli=str(Path(args.cli).resolve())
with tempfile.TemporaryDirectory(prefix='cartmesh-fv-') as name:
    root=Path(name);mesh=root/'one.solver.cm2d'
    mesh.write_text('''CM2D 1
VERTICES 4
0 0 0
1 1 0
2 1 1
3 0 1
EDGES 4
0 0 1 0 -1 2
1 1 2 0 -1 2
2 2 3 0 -1 2
3 3 0 0 -1 2
CELLS 1
0 0 0 1 4 0 1 2 3 4 0 1 2 3
AUDIT 0 0 0 0 0 0 0
END
''')
    def run(label,extra,code=0,input_mesh=mesh):
        prefix=root/label
        result=subprocess.run([cli,'--mesh',str(input_mesh),'--output',str(prefix),*extra],capture_output=True,text=True,timeout=20)
        assert result.returncode==code,(label,result.stdout,result.stderr)
        return prefix
    def load(prefix):
        return json.loads(prefix.with_suffix('.json').read_text(),parse_constant=lambda v: (_ for _ in ()).throw(ValueError(v)))
    normal=run('linear',['--problem','linear']);data=load(normal)
    assert data['converged'] and data['linfError']<1e-8
    rows=list(csv.DictReader(normal.with_suffix('.cells.csv').open()))
    assert len(rows)==1 and abs(float(rows[0]['value'])-2.5)<1e-8
    text=normal.with_suffix('.vtk').read_text()
    assert text.count('CELL_DATA ')==1 and 'SCALARS value double 1' in text
    failed=run('limited',['--problem','linear','--max-corrections','1'],2)
    assert load(failed)['converged'] is False
    for label,extra in [('overflow',['--outer-value','1e154']),('nan',['--outer-value','nan']),('negative',['--diffusivity','-1']),('unknown',['--problem','navier-stokes'])]:
        failed=run(label,extra,1);assert not failed.with_suffix('.json').exists()
    intermediate=root/'one.cm2d';intermediate.write_text(mesh.read_text())
    run('intermediate',[],1,intermediate)
    invalid=root/'invalid.solver.cm2d';invalid.write_text(mesh.read_text().replace('0 0 1 0 -1 2','0 0 1 99 -1 2'))
    run('invalid',[],1,invalid)
print('Native FV CLI outputs and explicit failures verified.')
