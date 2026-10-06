#!/usr/bin/env python3
"""Continue retained native cylinder fields; no independent PDE or remapping."""
import pathlib,json,subprocess,time
import argparse
p=argparse.ArgumentParser();p.add_argument('--output',type=pathlib.Path,default=pathlib.Path('outputs/corrected-curved-reference'));p.add_argument('--mid-step',type=float,default=4e-9);p.add_argument('--mid-restart',type=pathlib.Path);p.add_argument('--end',type=float,default=4e-7);p.add_argument('--budget',type=float,default=900);p.add_argument('--cases',nargs='+',choices=['mid','fine','fine-half'],default=['mid','fine','fine-half']);a=p.parse_args()
r=a.output;r.mkdir(parents=True,exist_ok=False);records=[]
for case,old,dt in [('mid','N192-W48/linear-resume',4e-9),('fine','N192-W96/linear',4e-9),('fine-half','N192-W96/linear-half',2e-9)]:
 if case not in a.cases:continue
 if case=='mid':dt=a.mid_step
 old=pathlib.Path('outputs/curved-cylinder')/old;cmd=json.loads(old.with_suffix('.command.json').read_text());
 if '--restart' in cmd:i=cmd.index('--restart');cmd=cmd[:i]
 for opt,value in [('--output',str(r/case)),('--end-time',str(a.end)),('--max-step',str(dt)),('--max-seconds',str(a.budget))]:cmd[cmd.index(opt)+1]=value
 cmd+=['--restart',str(a.mid_restart if case=='mid' and a.mid_restart else old.with_suffix('.checkpoint'))];(r/(case+'.command.json')).write_text(json.dumps(cmd,indent=2));t=time.monotonic()
 with (r/(case+'.log')).open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 records.append({'case':case,'command':cmd,'seconds':time.monotonic()-t,'returncode':p.returncode,'source':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()});(r/'campaign.json').write_text(json.dumps(records,indent=2));print(json.dumps(records[-1]),flush=True)
