#!/usr/bin/env python3
"""Adaptive linear solves must end at the original strict nonlinear gates."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
from verify_pressure_openings import rectangle,pressure_boundaries
import verify_native_flow as native

parser=argparse.ArgumentParser();parser.add_argument('--cli',required=True);args=parser.parse_args()
cli=str(Path(args.cli).resolve())
with tempfile.TemporaryDirectory(prefix='cartmesh-adaptive-linear-') as name:
 root=Path(name);mesh=root/'channel.solver.cm2d';rectangle(mesh,8)
 p=subprocess.run([cli,'--mesh',str(mesh),'--case','duct','--export-boundaries',str(root/'template')],capture_output=True,text=True,timeout=20)
 assert p.returncode==0,p.stderr
 pressure_boundaries(root/'template',root/'conditions')
 outputs={}
 for policy in ['strict','adaptive','engineering','relaxation08']:
  prefix=root/policy
  cmd=[cli,'--mesh',str(mesh),'--case','custom','--boundary',str(root/'conditions'),'--output',str(prefix),
   '--nu','.1','--speed','1','--tolerance','1e-9','--max-iterations','3000','--pressure-preconditioner','aggregation',
   '--steady-acceleration','anderson','--linear-policy','adaptive' if policy=='adaptive' else 'strict',
   '--convergence','engineering' if policy=='engineering' else 'strict']
  if policy=='relaxation08':cmd+=['--velocity-relaxation','.8']
  p=subprocess.run(cmd,capture_output=True,text=True,timeout=30);assert p.returncode==0,p.stderr
  summary=json.loads(Path(str(prefix)+'.json').read_text())
  assert summary['converged'] and summary['strictLinearFinal']
  assert summary['convergenceMode']==('engineering' if policy=='engineering' else 'strict')
  assert summary['momentumResidual']<1e-9 and summary['velocityChange']<1e-9 and summary['pressureChange']<1e-9
  assert summary['continuity']<1e-8
  assert (summary['adaptiveLinearSteps']>0)==(policy=='adaptive')
  audit=native.verify_case(mesh,prefix,'custom',.1,1,native.argument_parser().parse_args(['--max-iterations','3000']))
  assert audit['valid'],audit['issues']
  outputs[policy]=list(csv.DictReader(Path(str(prefix)+'.cells.csv').open()))
 for policy in ['adaptive','engineering','relaxation08']:
  for field in ['u','v','p']:
   difference=max(abs(float(a[field])-float(b[field])) for a,b in zip(outputs['strict'],outputs[policy]))
   assert difference<1e-6,(policy,field,difference)
 # Unsupported temporal combinations are refused before any accepted state.
 p=subprocess.run([cli,'--mesh',str(mesh),'--case','channel','--output',str(root/'invalid'),
   '--linear-policy','adaptive','--time-step','.01','--steps','1'],capture_output=True,text=True,timeout=20)
 assert p.returncode!=0 and 'steady laminar' in p.stderr
print('Adaptive and strict solves pass original gates, independent equations and field comparison; transient misuse refused')
