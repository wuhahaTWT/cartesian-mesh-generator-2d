#!/usr/bin/env python3
"""Bounded serial laminar benchmark. Timing alone never certifies CFD accuracy."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import shutil
import time
from verify_thermal_scale import run_process_group


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cli',type=Path,default=Path('build/cartmesh2d_flow_cli'))
    p.add_argument('--mesh',type=Path,required=True)
    p.add_argument('--boundary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--iterations',type=int,default=15000)
    p.add_argument('--timeout',type=float,default=180)
    p.add_argument('--linear-policy',choices=['strict','adaptive'],default='adaptive')
    p.add_argument('--tolerance',type=float,default=1e-8)
    p.add_argument('--velocity-relaxation',type=float,default=.6)
    p.add_argument('--no-resource-counters',action='store_true',help='Run without platform time wrapper; peak RSS is unmeasured')
    args=p.parse_args()
    if not 0<args.velocity_relaxation<=1:p.error('velocity relaxation must be in (0,1]')
    if not 0<args.timeout<=600 or not 0<args.iterations<=100000 or not 1e-12<=args.tolerance<=1e-6:
        p.error('invalid bounded budget or strict tolerance')
    cli=args.cli.resolve(strict=True);mesh=args.mesh.resolve(strict=True);bc=args.boundary.resolve(strict=True)
    prefix=args.output.resolve();prefix.parent.mkdir(parents=True,exist_ok=True)
    if prefix.suffix:p.error('output prefix must have no suffix (log helper uses suffix replacement)')
    if any(prefix.parent.glob(prefix.name+'.*')):p.error('output prefix already exists; preserve prior evidence')
    if shutil.disk_usage(prefix.parent).free<5*1024**3:p.error('resource stop: less than 5 GiB free')
    cmd=[str(cli),'--mesh',str(mesh),'--output',str(prefix),'--case','custom','--boundary',str(bc),
         '--nu','.1','--speed','1','--tolerance',str(args.tolerance),'--max-iterations',str(args.iterations),
         '--convection','face-limited-linear','--pressure-preconditioner','aggregation',
         '--steady-acceleration','anderson','--linear-policy',args.linear_policy,'--convergence','strict','--velocity-relaxation',str(args.velocity_relaxation),'--profile']
    wrapper=[] if args.no_resource_counters else ['/usr/bin/time','-l' if platform.system()=='Darwin' else '-v']
    input_hashes={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [cli,mesh,bc,Path(__file__).resolve()]}
    start=time.monotonic();record=run_process_group(wrapper+cmd,prefix,args.timeout);record['wallSeconds']=time.monotonic()-start
    record['inputHashes']=input_hashes
    record['timeoutSeconds']=args.timeout
    record['resourceCountersMeasured']=not args.no_resource_counters
    record['scope']='Same strict stopping gates; custom incompressible laminar flow, nu=.1 and Uref=1. Independent equation/accuracy audit required separately.'
    stderr=Path(record['stderr']).read_text()
    pattern=r'(\d+)\s+maximum resident set size' if platform.system()=='Darwin' else r'Maximum resident set size \(kbytes\):\s*(\d+)'
    match=re.search(pattern,stderr);record['maximumRssBytes']=int(match[1])*(1 if platform.system()=='Darwin' else 1024) if match else None
    for suffix,key in [('.json','summary'),('.performance.json','performance')]:
        path=Path(str(prefix)+suffix)
        if path.exists():record[key]=json.loads(path.read_text())
    summary=record.get('summary',{})
    record['converged']=record['returncode']==0 and summary.get('converged') is True and summary.get('strictLinearFinal') is True
    record['accuracyQualified']=False
    Path(str(prefix)+'.measurement.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:record.get(k) for k in ['returncode','timedOut','wallSeconds','maximumRssBytes','converged','performance']}))
    return 0 if record['converged'] else 2


if __name__=='__main__':raise SystemExit(main())
