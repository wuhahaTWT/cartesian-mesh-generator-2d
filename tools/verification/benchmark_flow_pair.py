#!/usr/bin/env python3
"""Serial identical-input baseline/candidate flow timing and artifact comparison.

Uses the OS time tool for peak RSS on macOS/Linux when available. Timing and
byte equality do not replace independent physical/mesh validation. Output must
be new; every completed, failed or timed-out invocation is retained.
"""
import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import re
import signal
import subprocess
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def write(path, value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def differences(left, right):
    result = dict.fromkeys(('u','v','p'),0.)
    with Path(str(left)+'.cells.csv').open() as x,Path(str(right)+'.cells.csv').open() as y:
        for a,b in itertools.zip_longest(csv.DictReader(x),csv.DictReader(y)):
            if a is None or b is None or a['cell'] != b['cell'] or a['area'] != b['area']:
                raise ValueError('cell coverage/area/order differs')
            for key in result:
                av,bv=float(a[key]),float(b[key])
                if not math.isfinite(av) or not math.isfinite(bv):raise ValueError('nonfinite cell value')
                result[key]=max(result[key],abs(av-bv))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('baseline','candidate','mesh','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--case',choices=('external','channel','cavity','taylor-green'),required=True)
    p.add_argument('--nu',type=float,required=True);p.add_argument('--speed',type=float,default=1.)
    p.add_argument('--convection',choices=('upwind','limited-linear'),default='limited-linear')
    p.add_argument('--tolerance',type=float,default=1e-9);p.add_argument('--max-iterations',type=int,default=1500)
    p.add_argument('--dt',type=float);p.add_argument('--steps',type=int)
    p.add_argument('--repeats',type=int,default=2);p.add_argument('--timeout',type=float,default=180.)
    a=p.parse_args()
    if any(not math.isfinite(v) or v<=0 for v in (a.nu,a.speed,a.tolerance,a.timeout)) or a.repeats<1 or a.max_iterations<1:p.error('invalid controls or budgets')
    if (a.dt is None)!=(a.steps is None) or (a.dt is not None and (not math.isfinite(a.dt) or a.dt<=0 or a.steps<1)):p.error('positive dt and steps must be supplied together')
    binaries={key:getattr(a,key).resolve(strict=True) for key in ('baseline','candidate')}
    mesh=a.mesh.resolve(strict=True);a.output.mkdir(parents=True,exist_ok=False)
    report={'scope':'serial timing and output comparison only; independent physical audit separate',
            'platform':platform.platform(),'machine':platform.machine(),'mesh':str(mesh),'meshSha256':sha(mesh),
            'binarySha256':{k:sha(v) for k,v in binaries.items()},'runs':[],'pairs':[]}
    suffixes=['.cells.csv','.faces.csv','.fields.json','.vtk','.residuals.csv']
    if a.dt is not None:suffixes+=['.checkpoint','.time-history.csv']
    for repeat in range(a.repeats):
        pair={}
        # Reverse the second pair to reduce systematic order bias.
        for label in (('baseline','candidate') if repeat%2==0 else ('candidate','baseline')):
            prefix=(a.output/f'{repeat}-{label}'/'result').resolve();prefix.parent.mkdir()
            cmd=[str(binaries[label]),'--mesh',str(mesh),'--output',str(prefix),'--case',a.case,'--nu',str(a.nu),'--speed',str(a.speed),'--convection',a.convection,'--viscous-stress','symmetric','--tolerance',str(a.tolerance),'--max-iterations',str(a.max_iterations),'--profile']
            if a.dt is not None:cmd+=['--time-step',str(a.dt),'--steps',str(a.steps)]
            wrapper=[]
            if Path('/usr/bin/time').is_file():
                if platform.system()=='Darwin':wrapper=['/usr/bin/time','-l']
                elif platform.system()=='Linux':wrapper=['/usr/bin/time','-v']
            start=time.monotonic();timed_out=False
            with Path(str(prefix)+'.stdout.log').open('w') as out,Path(str(prefix)+'.stderr.log').open('w') as err:
                child=subprocess.Popen(wrapper+cmd,stdout=out,stderr=err,start_new_session=os.name=='posix')
                try:code=child.wait(timeout=a.timeout)
                except subprocess.TimeoutExpired:
                    timed_out=True
                    if os.name=='posix':
                        try:os.killpg(child.pid,signal.SIGKILL)
                        except ProcessLookupError:pass  # Child exited at the timeout boundary.
                    else:child.kill()
                    child.wait();code=None
            record={'label':label,'repeat':repeat,'command':wrapper+cmd,'prefix':str(prefix),'returnCode':code,'timedOut':timed_out,'timeoutSeconds':a.timeout,'elapsedSeconds':time.monotonic()-start}
            stderr=Path(str(prefix)+'.stderr.log').read_text()
            match=re.search(r'(\d+)\s+maximum resident set size',stderr) if platform.system()=='Darwin' else re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',stderr)
            record['peakRssBytes']=int(match.group(1))*(1 if platform.system()=='Darwin' else 1024) if match else None
            report['runs'].append(record);write(a.output/'comparison.json',report)
            if code!=0:print(json.dumps(record),flush=True);return 1
            record['summary']=json.loads(Path(str(prefix)+'.json').read_text())
            record['performance']=json.loads(Path(str(prefix)+'.performance.json').read_text())
            solve_seconds=record['performance'].get('solveSeconds')
            if not isinstance(solve_seconds,(int,float)) or isinstance(solve_seconds,bool) or not math.isfinite(solve_seconds) or solve_seconds<=0:
                raise ValueError('performance solveSeconds must be finite and positive')
            if record['summary'].get('converged') is not True:raise ValueError('zero exit but unconverged summary')
            record['sha256']={suffix:sha(str(prefix)+suffix) for suffix in suffixes};pair[label]=record
            write(a.output/'comparison.json',report);print(json.dumps({k:record[k] for k in ('label','repeat','returnCode','elapsedSeconds','peakRssBytes')}),flush=True)
        report['pairs'].append({'repeat':repeat,'maxAbsoluteFieldDifference':differences(pair['baseline']['prefix'],pair['candidate']['prefix']),
            'byteIdentical':{s:pair['baseline']['sha256'][s]==pair['candidate']['sha256'][s] for s in suffixes},
            'elapsedSpeedup':pair['baseline']['elapsedSeconds']/pair['candidate']['elapsedSeconds'],
            'solveSpeedup':pair['baseline']['performance']['solveSeconds']/pair['candidate']['performance']['solveSeconds']})
        write(a.output/'comparison.json',report)
    return 0


if __name__=='__main__':raise SystemExit(main())
