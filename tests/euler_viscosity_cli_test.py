#!/usr/bin/env python3
"""Compressible laminar transport: continuum modes, Couette wall power, restart.

Refinement checks are small-signal/analytic consistency tests. They do not confer
engineering accuracy on arbitrary Cut-cell boundary layers.
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
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
import verify_euler as e
from verify_heat_conduction import linear_euler_fourier_mode
p=argparse.ArgumentParser();p.add_argument('--cli',required=True);p.add_argument('--output',type=Path);args=p.parse_args()
temp=None if args.output else tempfile.TemporaryDirectory(prefix='cm2d-viscous-')
root=args.output or Path(temp.name);root.mkdir(parents=True,exist_ok=True)

def run(mesh,name,extra=(),expected=0):
    prefix=root/name
    command=[args.cli,'--mesh',str(mesh),'--output',str(prefix),'--case','sealed','--wall-model','no-slip',
        '--viscosity','.1','--conductivity','.1','--gas-r','1','--flux','hllc','--order','2','--end-time','.02',*map(str,extra)]
    proc=subprocess.run(command,capture_output=True,text=True,timeout=90)
    assert proc.returncode==expected,(command,proc.stdout,proc.stderr)
    return prefix

def rows(prefix,suffix):
    with Path(str(prefix)+suffix).open() as f:return list(csv.DictReader(f))

def mutate(mesh,prefix,suffix,field,index=-1):
    file=Path(str(prefix)+suffix);original=file.read_text();r=rows(prefix,suffix);r[index][field]=str(float(r[index][field])+1)
    try:
        with file.open('w') as f:w=csv.DictWriter(f,fieldnames=r[0].keys());w.writeheader();w.writerows(r)
        try:e.audit(mesh,prefix)
        except ValueError:pass
        else:raise AssertionError('mutation accepted: '+field)
    finally:file.write_text(original)

mesh=root/'box.solver.cm2d';e.rectangle(mesh,16,4)
closed=run(mesh,'closed',['--u',.2,'--v',-.1]);report={'closed':e.audit(mesh,closed),'refinement':{},'couette':[]}
h=rows(closed,'.history.csv');energy0=.1*(1/.4+.025)
assert abs(float(h[-1]['totalEnergy'])-energy0)/(2*energy0)<1024*math.ulp(1.)
assert all(float(r['boundaryHeat'])==float(r['boundaryViscousWork'])==float(r['boundaryEnergy'])==0 for r in h)
for suffix,key in [('.faces.csv','viscousMomentumX'),('.faces.csv','viscousMomentumY'),('.faces.csv','viscousWork'),('.cells.csv','viscousRate'),('.history.csv','viscousCourant'),('.history.csv','boundaryViscousWork')]:mutate(mesh,closed,suffix,key)
mutate(mesh,closed,'.history.csv','boundaryViscousWork',0)
# Dimensional similarity: x=L*x*, u=V*u*, rho=D*rho*, p=D V^2 p*,
# mu=D V L mu*, k=D R V L k*, t=L/V t*. Same gas gamma.
scaled_mesh=root/'scaled.solver.cm2d';e.rectangle(scaled_mesh,16,4,width=.02,height=.002)
scaled=run(scaled_mesh,'scaled',['--gas-r',287.05,'--density',2,'--pressure',20000,'--u',20,'--v',-10,
    '--viscosity',.4,'--conductivity',.4*287.05,'--end-time',.000004])
report['dimensionalAudit']=e.audit(scaled_mesh,scaled);scaling_error=0
for a,b in zip(rows(closed,'.cells.csv'),rows(scaled,'.cells.csv')):
    for key,scale in [('rho',2),('u',100),('v',100),('p',20000),('temperature',10000/287.05)]:
        x,y=float(a[key]),float(b[key])/scale;scaling_error=max(scaling_error,abs(x-y)/max(1,abs(x),abs(y)))
assert scaling_error<4096*math.ulp(1.),scaling_error
report['dimensionalSimilarityError']=scaling_error
# Every selected convective/time method must include the same constitutive law.
report['methods']={}
for flux,order in [('rusanov',1),('rusanov',2),('hllc',1),('hllc',2)]:
    result=run(mesh,f'{flux}-{order}',['--u',.2,'--flux',flux,'--order',order]);report['methods'][f'{flux}-{order}']=e.audit(mesh,result)
# Stop with the accepted-step budget; target-time truncation must not alter dt.
first=run(mesh,'split',['--u',.2,'--v',-.1,'--max-steps',len(h)//2+1],2)
resumed=run(mesh,'resumed',['--u',.2,'--v',-.1,'--restart',str(first)+'.checkpoint'])
assert Path(str(resumed)+'.checkpoint').read_bytes()==Path(str(closed)+'.checkpoint').read_bytes()
for changed in [('--viscosity',.2),('--wall-model','slip'),('--conductivity',.2)]:run(mesh,'bad-restart',['--u',.2,'--v',-.1,'--restart',str(first)+'.checkpoint',*changed],1)
for bad in [('--viscosity',-1),('--viscosity',0),('--wall-model','rough')]:run(mesh,'invalid',bad,1)
# Legacy zero-mu Fourier (v2) and inviscid (v1) remain resumable and byte-identical
# with explicit zero viscosity vs an omitted option; the old executable comparison
# is retained separately at delivery, rather than bundled in the normal tests.
for k,version in [(0,1),(.1,2)]:
    r=run(mesh,f'legacy{version}',['--wall-model','slip','--viscosity',0,'--conductivity',k,'--end-time','.001'])
    assert Path(str(r)+'.checkpoint').read_text().startswith(f'CM2D_EULER_CHECKPOINT {version}\n')
    run(mesh,f'legacy{version}-resume',['--wall-model','slip','--viscosity',0,'--conductivity',k,'--restart',str(r)+'.checkpoint'])
# Linearized continuum NS: shear uses nu, longitudinal velocity uses 4 nu / 3.
# Normalize errors by the imposed small amplitude; no native coefficient enters
# either exact solution. Three refinements should approach second order (>1.7).
for problem in ('shear-wave','thermal-wave'):
    refinement=[]
    for n in (16,32,64):
        m=root/f'{problem}-{n}.solver.cm2d';e.rectangle(m,n,4)
        result=run(m,f'{problem}-{n}',['--case',problem,'--wall-model','slip','--end-time',.1])
        audit=e.audit(m,result);r=rows(result,'.cells.csv');errors={}
        if problem=='shear-wave':
            amplitude=1e-4*math.sqrt(1.4);decay=math.exp(-.1*(2*math.pi)**2*.1)
            errors['v']=sum(abs(float(c['v'])-amplitude*decay*math.cos(2*math.pi*float(c['x']))) for c in r)/len(r)/amplitude
        else:
            A,B,C=linear_euler_fourier_mode(.1,.1,viscosity=.1)
            for field,base,amplitude,phase in [('rho',1,A,math.cos),('u',0,B,math.sin),('temperature',1,C,math.cos)]:
                errors[field]=sum(abs(float(c[field])-base-amplitude*phase(2*math.pi*float(c['x']))) for c in r)/len(r)/1e-5
        refinement.append(dict(n=n,errors=errors,audit=audit))
    orders={key:math.log(refinement[-2]['errors'][key]/refinement[-1]['errors'][key],2) for key in refinement[-1]['errors']}
    print(problem,orders,flush=True);assert min(orders.values())>1.7,refinement
    report['refinement'][problem]={'levels':refinement,'orders':orders}
# Moving tangential walls with periodic x. Exact steady Couette solution:
# u=Uy/H, T=Tw+mu U^2/(2k)*(y/H)*(1-y/H), p=const, rho=p/(RT).
# Start from this analytic field to test discretization consistency. This is
# not a test of convergence from rest; target time is not steady convergence.
for n in (8,16,32):
    m=root/f'couette-{n}.solver.cm2d';e.rectangle(m,4,n,width=.2,height=1.)
    exported=run(m,f'couette-template-{n}',['--case','thermal-wave','--wall-model','slip','--end-time','.00001'])
    lines=Path(str(exported)+'.boundaries').read_text().splitlines();lines[0]='CM2D_EULER_BOUNDARY 3'
    for index in range(2,len(lines)-1):
        r=shlex.split(lines[index]);r=r[:15]
        if abs(float(r[5]))>0: # horizontal wall, normal y
            r[6]='no-slip-wall';r[11]='-';r[12]='top' if float(r[3])>.5 else 'bottom';r[13]='temperature';r[14]='1'
            r+=['.5' if float(r[3])>.5 else '0','0']
        else:r+=['0','0']
        r[12]=json.dumps(r[12]);lines[index]=' '.join(r)
    boundary=root/f'couette-{n}.boundary';boundary.write_text('\n'.join(lines)+'\n')
    # Custom walls are fully specified by the file, not CLI preset switches.
    prefix=root/f'couette-initial-{n}'
    cmd=[args.cli,'--mesh',str(m),'--output',str(prefix),'--case','custom','--boundary',str(boundary),'--gas-r','1','--viscosity','.1','--conductivity','.1','--flux','hllc','--order','2','--end-time','.00001']
    proc=subprocess.run(cmd,capture_output=True,text=True,timeout=20);assert proc.returncode==0,proc.stderr
    data=Path(str(prefix)+'.checkpoint').read_text().splitlines();index=next(i for i,s in enumerate(data) if s.startswith('STATE '));data[index]='STATE 0 0'
    g=e.geometry.read_cm2d(m);measure=e.geometry.measure(g,1e-11,1e-10)
    initial_energy=0
    for i,(x,y) in enumerate(measure.centroids):
        t=1+.5*.5/2*y*(1-y);q=e.conservative((1/t,.5*y,0,1),1.4);data[index+1+i]=' '.join(format(v,'.17g') for v in q);initial_energy+=measure.areas[i]*q[3]
    start=root/f'couette-{n}.checkpoint';start.write_text('\n'.join(data)+'\n')
    out=root/f'couette-{n}';cmd[cmd.index('--output')+1]=str(out);cmd[cmd.index('--end-time')+1]='.1';cmd+=['--restart',str(start)]
    proc=subprocess.run(cmd,capture_output=True,text=True,timeout=40);assert proc.returncode==0,proc.stderr
    audit=e.audit(m,out);r=rows(out,'.cells.csv');hist=rows(out,'.history.csv')
    error=sum(abs(float(c['temperature'])-(1+.125*float(c['y'])*(1-float(c['y'])))) for c in r)/len(r)
    gain=float(hist[-1]['totalEnergy'])-initial_energy
    incoming=-math.fsum(float(c['dt'])*(float(c['boundaryHeat'])+float(c['boundaryViscousWork'])) for c in hist)
    budget=abs(gain-incoming)/(initial_energy+float(hist[-1]['totalEnergy'])+abs(incoming))
    assert budget<1024*math.ulp(1.) and audit['boundaryViscousWork']<0
    report['couette'].append(dict(n=n,temperatureL1=error,energyBalanceRelative=budget,audit=audit))
    if n==8:
        altered=boundary.read_text().replace(' .5 0',' .6 0');boundary.write_text(altered)
        bad=subprocess.run(cmd,capture_output=True,text=True,timeout=10);assert bad.returncode==1 and 'incompatible' in bad.stderr
        boundary.write_text('\n'.join(lines)+'\n')
couette_order=math.log(report['couette'][-2]['temperatureL1']/report['couette'][-1]['temperatureL1'],2)
print('Couette consistency order',couette_order,flush=True);assert couette_order>1.7,report['couette'];report['couetteOrder']=couette_order
report.update(restartByteIdentical=True,legacyVersionsAccepted=[1,2],mutationsRejected=True)
(root/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
print('Compressible viscous transport checks passed',flush=True)
