"""Actual Fluent GPU launch probe, no fallback or licence changes."""
import json,time
from pathlib import Path
import ansys.fluent.core as f
out=Path('outputs/pc-mesh-study/gpu/probe2d');out.mkdir(parents=True,exist_ok=True)
r={'requestedDimension':2,'requestedGpu':[0]};s=None;t=time.monotonic()
try:
 s=f.launch_fluent(fluent_path=r'E:\ANSYS Inc\v261\fluent\ntbin\win64\fluent.exe',dimension=2,precision='double',gpu=[0],processor_count=1,ui_mode='no_gui',start_timeout=70,cwd=str(out.resolve()),start_transcript=False,start_watchdog=True,env={'REMOTING_SERVER_ADDRESS':'127.0.0.1'})
 s.transcript.start(file_name=str((out/'transcript.log').resolve()))
 r['version']=str(s.get_fluent_version());r['launchPassed']=True
 r['solverVersion']=str(s.scheme.eval('(cx-version)'))
except Exception as e:r['error']=repr(e)
finally:
 if s:
  try:s.exit(timeout=10,wait=10);r['closed']=True
  except Exception as e:r['exitError']=repr(e)
 r['seconds']=time.monotonic()-t;(out/'result.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(r,flush=True)
