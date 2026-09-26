#!/usr/bin/env python3
"""Wall flux precision, cold-start steady Couette, temporal refinement and audit.

Analytic seed files are test inputs, never production accepted checkpoints.
Spatial order remains two; scalar point recovery is not a high-order cell-
average reconstruction. The harmonic operator test covers non-polynomial heat.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
import verify_euler as e
p=argparse.ArgumentParser();p.add_argument('--cli',required=True);p.add_argument('--output',type=Path);args=p.parse_args()
temp=None if args.output else tempfile.TemporaryDirectory(prefix='cm2d-wall-accuracy-')
root=args.output or Path(temp.name);root.mkdir(parents=True,exist_ok=True)
report={'couette':[],'temporal':[],'coldRefinement':[]}

def execute(mesh,name,extra=(),expected=0):
    out=root/name;command=[args.cli,'--mesh',str(mesh),'--output',str(out),'--gas-r','1','--viscosity','.1','--conductivity','.1','--flux','hllc','--order','2','--wall-gradient','quadratic','--end-time','.1',*map(str,extra)]
    start=time.perf_counter();result=subprocess.run(command,capture_output=True,text=True,timeout=150)
    assert result.returncode==expected,(command,result.stdout[-1000:],result.stderr)
    return out,time.perf_counter()-start

def rows(prefix,suffix):
    with Path(str(prefix)+suffix).open() as f:return list(csv.DictReader(f))
def summary(prefix):return json.loads(Path(str(prefix)+'.json').read_text())
def state(prefix):
    return [[float(row[key]) for key in ('rho','rhoU','rhoV','rhoE')] for row in rows(prefix,'.cells.csv')]
def distance(a,b):return math.sqrt(math.fsum((x-y)**2 for row,col in zip(a,b) for x,y in zip(row,col))/len(a))

for n in (8,16,32):
    mesh=root/f'couette-{n}.solver.cm2d';e.rectangle(mesh,2,n,width=1,height=1)
    template,_=execute(mesh,f'template-{n}',['--case','thermal-wave','--wall-gradient','linear','--end-time',1e-6])
    lines=Path(str(template)+'.boundaries').read_text().splitlines();lines[0]='CM2D_EULER_BOUNDARY 3'
    for i in range(2,len(lines)-1):
        r=shlex.split(lines[i])[:15]
        if float(r[5])!=0:
            r[6]='no-slip-wall';r[11]='-';r[12]='top' if float(r[3])>.5 else 'bottom';r[13]='temperature';r[14]='1';r+=['.5' if float(r[3])>.5 else '0','0']
        else:r+=['0','0']
        r[12]=json.dumps(r[12]);lines[i]=' '.join(r)
    bc=root/f'couette-{n}.boundary';bc.write_text('\n'.join(lines)+'\n')
    initial,_=execute(mesh,f'initial-{n}',['--case','custom','--boundary',bc,'--end-time',1e-6])
    data=Path(str(initial)+'.checkpoint').read_text().splitlines();index=next(i for i,s in enumerate(data) if s.startswith('STATE '));data[index]='STATE 0 0'
    g=e.geometry.read_cm2d(mesh);measure=e.geometry.measure(g,1e-11,1e-10)
    for i,(x,y) in enumerate(measure.centroids):data[index+1+i]=' '.join(format(v,'.17g') for v in e.conservative((1/(1+.125*y*(1-y)),.5*y,0,1),1.4))
    seed=root/f'analytic-input-{n}.checkpoint';seed.write_text('\n'.join(data)+'\n')
    for scheme in ('linear','quadratic'):
        out,seconds=execute(mesh,f'{scheme}-{n}',['--case','custom','--boundary',bc,'--restart',seed,'--wall-gradient',scheme])
        audit=e.audit(mesh,out);s=summary(out);r=rows(out,'.cells.csv')
        err=sum(abs(float(c['temperature'])-1-.125*float(c['y'])*(1-float(c['y']))) for c in r)/len(r)
        flux_error=abs(s['boundaryHeat']-.025)/.025
        report['couette'].append(dict(n=n,scheme=scheme,temperatureL1=err,wallHeatRelativeError=flux_error,wallWorkRelativeError=abs(s['boundaryViscousWork']+.025)/.025,seconds=seconds,audit=audit))
        if scheme=='quadratic':assert err<4096*math.ulp(1.) and flux_error<4096*math.ulp(1.),report['couette'][-1]
    if n in (8,16,32):
        cold_data=data.copy()
        for i in range(len(measure.areas)):cold_data[index+1+i]='1 0 0 2.5000000000000004'
        cold_seed=root/f'cold-input-{n}.checkpoint';cold_seed.write_text('\n'.join(cold_data)+'\n')
        cold,seconds=execute(mesh,f'cold-steady-{n}',['--case','custom','--boundary',bc,'--restart',cold_seed,'--end-time',60,'--max-steps',1000000])
        audit=e.audit(mesh,cold);s=summary(cold);r=rows(cold,'.cells.csv');h=rows(cold,'.history.csv')
        error=max(abs(float(c['temperature'])-1-.125*float(c['y'])*(1-float(c['y']))) for c in r)
        speed_error=max(abs(float(c['u'])-.5*float(c['y'])) for c in r)
        # Continuum peak temperature rise is 0.03125 K. The 1e-6 relative goal
        # separates genuine steady approach from the prior short seeded test.
        assert error/.03125<1e-6 and speed_error/.5<1e-6,(error,speed_error)
        qdot=abs(s['boundaryHeat']+s['boundaryViscousWork'])/.025
        assert qdot<1e-6,qdot
        # Additionally evaluate dU/dt, rather than inferring steady state solely
        # from target time or exact global cancellation of two boundary totals.
        change=max(abs(float(c[key])-float(c[old]))/s['lastStep'] for c in r for key,old in [('rho','previousRho'),('rhoU','previousRhoU'),('rhoV','previousRhoV'),('rhoE','previousRhoE')])
        assert change<1e-7,change
        energy0=1/.4;incoming=-math.fsum(float(c['dt'])*(float(c['boundaryHeat'])+float(c['boundaryViscousWork'])) for c in h);gain=float(h[-1]['totalEnergy'])-energy0
        budget=abs(gain-incoming)/(energy0+float(h[-1]['totalEnergy'])+abs(incoming));assert budget<8192*math.ulp(1.),budget
        a=.125;b=1+a/4;exact_pressure=math.sqrt(a*b)/(2*math.atanh(.5*math.sqrt(a/b)))
        pressure_error=max(abs(float(c['p'])-exact_pressure)/exact_pressure for c in r)
        item=dict(continuumPressure=exact_pressure,pressureRelativeLinf=pressure_error,n=n,seconds=seconds,steps=s['steps'],temperatureLinf=error,velocityLinf=speed_error,wallPowerImbalanceRelative=qdot,maximumConservativeTimeDerivative=change,cumulativeEnergyBalanceRelative=budget,audit=audit)
        report['coldRefinement'].append(item)
        if n==16:report['coldStart']=item
    if n==8:
        cold_data=data.copy()
        for i in range(len(measure.areas)):cold_data[index+1+i]='1 0 0 2.5000000000000004'
        cold_seed=root/'time-input.checkpoint';cold_seed.write_text('\n'.join(cold_data)+'\n')
        temporal=[]
        for dt in (.001,.0005,.00025,.00003125,.000015625):
            out,_=execute(mesh,f'time-{dt}',['--case','custom','--boundary',bc,'--restart',cold_seed,'--end-time','.02','--max-step',dt])
            e.audit(mesh,out);history=rows(out,'.history.csv');assert all(1e-14<=float(c['dt'])<=dt for c in history);assert summary(out)['time']==.02;temporal.append((dt,state(out),out))
        ref=temporal[-1][1];ref_delta=distance(ref,temporal[-2][1])
        for dt,field,_ in temporal[:3]:report['temporal'].append(dict(dt=dt,conservativeRmsError=distance(field,ref)))
        errors=[v['conservativeRmsError'] for v in report['temporal']]
        order=math.log2(errors[-2]/errors[-1]);assert order>1.8 and ref_delta<.1*errors[-1],(order,errors,ref_delta)
        report['temporalOrder']=order;report['temporalReferenceDifference']=ref_delta
        # Restart at an accepted step, with the same end-time truncation. The
        # wall reconstruction is numerical control, not a changed physical BC.
        full=temporal[0][2];split,_=execute(mesh,'split',['--case','custom','--boundary',bc,'--restart',cold_seed,'--end-time','.02','--max-step','.001','--max-steps','8'],2)
        resumed,_=execute(mesh,'resumed',['--case','custom','--boundary',bc,'--restart',str(split)+'.checkpoint','--end-time','.02','--max-step','.001'])
        assert Path(str(full)+'.checkpoint').read_bytes()==Path(str(resumed)+'.checkpoint').read_bytes()
        switched,_=execute(mesh,'switched',['--case','custom','--boundary',bc,'--restart',str(split)+'.checkpoint','--end-time','.02','--max-step','.001','--wall-gradient','linear'])
        e.audit(mesh,switched);report['restartByteIdentical']=True;report['numericalSwitchResumable']=True
        # Metadata alone must not be enough to certify a quadratic run.
        file=Path(str(full)+'.json');original=file.read_text()
        for key,value in [('wallGradient','linear'),('quadraticHeatWalls',0),('quadraticViscousWalls',0)]:
            changed=json.loads(original);changed[key]=value;file.write_text(json.dumps(changed))
            try:e.audit(mesh,full)
            except ValueError:pass
            else:raise AssertionError('modified wall diagnostics accepted')
            finally:file.write_text(original)
        for suffix,key in [('.faces.csv','heatFlux'),('.faces.csv','viscousMomentumX'),('.faces.csv','viscousWork'),('.cells.csv','heatRate'),('.cells.csv','viscousRate')]:
            file=Path(str(full)+suffix);original=file.read_text();r=rows(full,suffix);r[0][key]=str(float(r[0][key])+1)
            with file.open('w') as f:w=csv.DictWriter(f,fieldnames=r[0].keys());w.writeheader();w.writerows(r)
            try:e.audit(mesh,full)
            except ValueError:pass
            else:raise AssertionError('modified wall flux/rate accepted')
            finally:file.write_text(original)
        report['mutationsRejected']=True
    print('wall precision n=',n,flush=True)
pressure_order=math.log2(report['coldRefinement'][-2]['pressureRelativeLinf']/report['coldRefinement'][-1]['pressureRelativeLinf'])
assert pressure_order>1.8 and report['coldRefinement'][-1]['pressureRelativeLinf']<.0001,report['coldRefinement']
report['pressureOrder']=pressure_order
(root/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print('Wall precision verification passed',flush=True)
