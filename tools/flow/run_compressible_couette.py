#!/usr/bin/env python3
"""Prepare polygon inputs and invoke native Couette reference; no PDE evaluation."""
import argparse, json, math, pathlib, subprocess, time
p=argparse.ArgumentParser()
p.add_argument('--output',required=True); p.add_argument('--segments',type=int,default=64)
p.add_argument('--level',type=int,default=5); p.add_argument('--phase',type=float,default=.375)
p.add_argument('--initial',choices=['exact','rest','thermal-exact','thermal-rest'],default='exact')
p.add_argument('--end',type=float,default=100); p.add_argument('--dt',type=float,default=.05)
p.add_argument('--budget',type=float,default=240); p.add_argument('--checkpoint')
a=p.parse_args(); out=pathlib.Path(a.output); out.mkdir(parents=True,exist_ok=False)
if a.segments<8: p.error('segments must be >=8')
xy=out/'annulus.xy'
with xy.open('w') as f:
    for radius in (2,1):
        for j in range(a.segments):
            theta=2*math.pi*(j+a.phase)/a.segments
            f.write(f'{radius*math.cos(theta):.17g} {radius*math.sin(theta):.17g}\n')
        f.write('\n')
commands=[]
def run(cmd,name):
    start=time.monotonic()
    with (out/(name+'.log')).open('w') as log: result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
    commands.append(dict(command=cmd,seconds=time.monotonic()-start,returncode=result.returncode))
    (out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    if result.returncode: raise SystemExit(result.returncode)
run(['build/cartmesh2d_cli',str(xy),str(out/'mesh'),str(a.level),'.125','.1','interior',str(out/'foam'),str(a.level-1),'0'],'mesh')
cmd=['build/cartmesh2d_compressible_couette_benchmark',str(out/'mesh.solver.cm2d'),str(out/a.initial),a.initial,str(a.end),str(a.dt),str(a.budget)]
if a.checkpoint: cmd.append(a.checkpoint)
run(cmd,'solve')
