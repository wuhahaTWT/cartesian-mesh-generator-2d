"""Shared-device telemetry; not attributed solely to this CFD process."""
from pathlib import Path
import subprocess,time
out=Path('outputs/pc-mesh-study/gpu-v2');out.mkdir(parents=True,exist_ok=True)
with (out/'device-telemetry.csv').open('w',encoding='utf-8') as f:
 f.write('timestamp, utilization.gpu [%], memory.used [MiB], power.draw [W], temperature.gpu, pstate\n')
 deadline=time.monotonic()+1800
 while time.monotonic()<deadline and not (out/'telemetry.stop').exists():
  r=subprocess.run(['nvidia-smi','--query-gpu=timestamp,utilization.gpu,memory.used,power.draw,temperature.gpu,pstate','--format=csv,noheader'],capture_output=True,text=True)
  if r.returncode==0:f.write(r.stdout);f.flush()
  time.sleep(2)
