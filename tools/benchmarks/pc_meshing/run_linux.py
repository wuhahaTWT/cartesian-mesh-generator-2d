"""Measure an external solver process inside its Linux container."""
import argparse,json,os,resource,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
p.add_argument('--cwd',type=Path,required=True);p.add_argument('--budget',type=int,required=True)
p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
command=a.command[1:] if a.command and a.command[0]=='--' else a.command
start=time.monotonic()
try:
    r=subprocess.run(command,cwd=a.cwd,timeout=a.budget);code=r.returncode;status='finished'
except subprocess.TimeoutExpired:code=124;status='timeout'
usage=resource.getrusage(resource.RUSAGE_CHILDREN)
data=dict(status=status,exitCode=code,wallSeconds=time.monotonic()-start,
    cpuSeconds=usage.ru_utime+usage.ru_stime,peakChildRssMiB=usage.ru_maxrss/1024,
    scope='Linux child-process rusage; peak is largest child RSS, not simultaneous tree sum; no energy measurement')
a.output.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data),flush=True)
raise SystemExit(code)
