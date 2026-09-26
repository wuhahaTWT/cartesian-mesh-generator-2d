#!/usr/bin/env python3
"""Total-energy Fourier transport: continuum coupling, wall budgets, restart, audit.

The refinement gate measures normalized small-signal error, not an engineering
accuracy certificate. Conservation/identity checks use roundoff-scaled budgets.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
import verify_euler as euler
from verify_heat_conduction import linear_euler_fourier_mode

parser=argparse.ArgumentParser()
parser.add_argument('--cli',required=True)
parser.add_argument('--output',type=Path,help='retain reproducible evidence instead of a temporary directory')
args=parser.parse_args()
temporary=None if args.output else tempfile.TemporaryDirectory(prefix='cm2d-conduction-')
root=args.output or Path(temporary.name)
root.mkdir(parents=True,exist_ok=True)


def run(mesh,name,extra=(),expected=0):
    prefix=root/name
    command=[args.cli,'--mesh',str(mesh),'--output',str(prefix),'--case','sealed',
             '--gas-r','1','--conductivity','.5','--flux','hllc','--order','2','--end-time','.02',*map(str,extra)]
    result=subprocess.run(command,capture_output=True,text=True,timeout=60)
    assert result.returncode==expected,(command,result.stdout,result.stderr)
    return prefix


def rows(prefix,suffix):
    with Path(str(prefix)+suffix).open() as stream:return list(csv.DictReader(stream))


def summary(prefix):return json.loads(Path(str(prefix)+'.json').read_text())


def rejects_mutation(mesh,prefix,suffix,field,index=-1):
    file=Path(str(prefix)+suffix);original=file.read_text()
    data=list(csv.DictReader(original.splitlines()))
    data[index][field]=str(float(data[index][field])+1)
    try:
        with file.open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=data[0].keys());writer.writeheader();writer.writerows(data)
        try:euler.audit(mesh,prefix)
        except ValueError:pass
        else:raise AssertionError('independent audit accepted modified '+field)
    finally:file.write_text(original)


# Constant field is an exact heat-operator nullspace, including insulated walls.
mesh=root/'box.solver.cm2d';euler.rectangle(mesh,16,4)
insulated=run(mesh,'insulated');report={'insulated':euler.audit(mesh,insulated),'walls':{},'refinement':[]}
assert all(float(r['heatFlux'])==0 for r in rows(insulated,'.faces.csv'))
assert all(float(r['rho'])==1 and float(r['p'])==1 for r in rows(insulated,'.cells.csv'))
for kind,value in [('temperature',2),('flux',-1)]:
    prefix=run(mesh,kind,['--wall-thermal',kind,'--wall-value',value])
    report['walls'][kind]=euler.audit(mesh,prefix)
    history=rows(prefix,'.history.csv');s=summary(prefix)
    # Sealed enclosure: integrated total energy + all outward wall heat = E(0).
    injected=-math.fsum(float(r['dt'])*float(r['boundaryHeat']) for r in history)
    energy0=.1/(1.4-1);energy=float(history[-1]['totalEnergy'])
    relative=abs(energy-energy0-injected)/(abs(energy)+energy0+abs(injected))
    assert relative<1024*math.ulp(1.),relative
    assert injected>0 and min(float(r['temperature']) for r in rows(prefix,'.cells.csv'))>1
    assert max(abs(float(r['u']))+abs(float(r['v'])) for r in rows(prefix,'.cells.csv'))>0
    if kind=='flux':assert abs(injected-2.2*.02)<256*math.ulp(1.),injected
    report['walls'][kind].update(integratedHeatInput=injected,totalEnergyGain=energy-energy0,historyEnergyBalanceRelative=relative)
    # Last stage heat data, stability bounds and history must all be checked.
    for suffix,field in [('.faces.csv','heatFlux'),('.faces.csv','convectiveEnergy'),
                         ('.cells.csv','heatRate'),('.history.csv','combinedCourant'),('.history.csv','boundaryHeat')]:
        rejects_mutation(mesh,prefix,suffix,field)
    rejects_mutation(mesh,prefix,'.history.csv','boundaryHeat',0)
    rejects_mutation(mesh,prefix,'.history.csv','totalEnergy',1)

# Both convective flux choices and time integrators transport the same wall heat.
report['methods']={}
for scheme,order in [('rusanov',1),('rusanov',2),('hllc',1)]:
    p=run(mesh,f'flux-{scheme}{order}',['--wall-thermal','flux','--wall-value',-1,'--flux',scheme,'--order',order])
    a=euler.audit(mesh,p);h=rows(p,'.history.csv')
    assert abs(float(h[-1]['totalEnergy'])-(.1/(1.4-1)+2.2*.02))<1024*math.ulp(1.)
    report['methods'][f'{scheme}{order}']=a

# Dimensional similarity scales length, time, density, EOS and conductivity
# together. This catches confusing k with diffusivity and using cp instead of cv.
length_scale=.02;rho_scale=1.225;pressure_scale=101325.;gas_r=287.05
speed_scale=math.sqrt(pressure_scale/rho_scale)
temperature_scale=pressure_scale/(rho_scale*gas_r)
time_scale=length_scale/speed_scale
si_mesh=root/'si.solver.cm2d';euler.rectangle(si_mesh,16,4,length_scale,.1*length_scale)
si=run(si_mesh,'si-temperature',['--density',rho_scale,'--pressure',pressure_scale,'--gas-r',gas_r,
    '--end-time',.02*time_scale,'--conductivity',.5*rho_scale*gas_r*speed_scale*length_scale,
    '--wall-thermal','temperature','--wall-value',2*temperature_scale])
report['dimensionalSimilarity']=euler.audit(si_mesh,si)
error=0
for actual,base in zip(rows(si,'.cells.csv'),rows(root/'temperature','.cells.csv')):
    for key,scale in [('rho',rho_scale),('p',pressure_scale),('u',speed_scale),('v',speed_scale),('temperature',temperature_scale)]:
        error=max(error,abs(float(actual[key])/scale-float(base[key]))/max(1,abs(float(base[key]))))
assert error<8192*math.ulp(1.),error
report['dimensionalSimilarity']['maximumNormalizedFieldDifference']=error
print('Euler-Fourier dimensional similarity',error,flush=True)

# Restart binds conductivity, wall type and wall value; a split trajectory is
# byte-identical to an uninterrupted run. Numerical controls are still mutable.
first=run(mesh,'budget',['--wall-thermal','temperature','--wall-value',2,'--max-steps',5],2)
assert euler.audit(mesh,first)['valid']
resumed=run(mesh,'resumed',['--wall-thermal','temperature','--wall-value',2,'--restart',str(first)+'.checkpoint'])
assert euler.audit(mesh,resumed)['valid']
assert Path(str(resumed)+'.checkpoint').read_bytes()==(root/'temperature.checkpoint').read_bytes()
for changed in [('--conductivity','.6'),('--wall-value','2.1'),('--wall-thermal','flux'),('--conductivity','0')]:
    run(mesh,'bad-restart',['--wall-thermal','temperature','--wall-value',2,'--restart',str(first)+'.checkpoint',*changed],1)
for invalid in [('--conductivity',-1),('--wall-thermal','temperature','--wall-value',0),
                ('--conductivity',0,'--wall-thermal','flux','--wall-value',-1)]:
    run(mesh,'bad-input',invalid,1)
legacy=run(mesh,'legacy',['--conductivity',0,'--max-steps',1],2)
assert Path(str(legacy)+'.checkpoint').read_text().startswith('CM2D_EULER_CHECKPOINT 1\n')
run(mesh,'legacy-resumed',['--conductivity',0,'--restart',str(legacy)+'.checkpoint'])

# A boundary file may mix independent thermal conditions on physical patches.
# A steady linear temperature gradient is covered in the skew-grid native test;
# here exercise actual file parsing, shared flux audit and custom restart binding.
boundary=root/'input-mixed.boundaries'
command=[args.cli,'--mesh',str(mesh),'--case','sealed','--gas-r','1','--conductivity','.5','--export-boundaries',str(boundary)]
subprocess.run(command,check=True,capture_output=True)
lines=boundary.read_text().splitlines()
for i,line in enumerate(lines[2:-1],2):
    tokens=line.split()
    if float(tokens[4])<0:tokens[-2:]=['temperature','2']
    elif float(tokens[4])>0:tokens[-2:]=['flux','.2']
    lines[i]=' '.join(tokens)
boundary.write_text('\n'.join(lines)+'\n')
custom=run(mesh,'mixed',['--case','custom','--boundary',boundary])
report['mixed']=euler.audit(mesh,custom)
original=Path(str(custom)+'.boundaries').read_text()
try:
    Path(str(custom)+'.boundaries').write_text(original.replace('temperature 2','temperature 3'))
    try:euler.audit(mesh,custom)
    except ValueError:pass
    else:raise AssertionError('independent audit accepted wrong thermal boundary')
finally:Path(str(custom)+'.boundaries').write_text(original)

# Continuum isobaric disturbance generates acoustic motion while heat diffuses.
# 1e-5 amplitude makes nonlinear corrections O(1e-5) after amplitude normalization,
# below these spatial errors. The independent 3x3 continuum system uses cv, not cp.
end=.08;amplitude=1e-5;mode=linear_euler_fourier_mode(end,.5)
assert abs(mode[1])>amplitude/20,'reference did not generate acoustic motion'
previous=None
for nx in [16,32,64]:
    mesh=root/f'wave{nx}.solver.cm2d';euler.rectangle(mesh,nx,4,1,.25)
    prefix=run(mesh,f'wave{nx}',['--case','thermal-wave','--end-time',end])
    audit=euler.audit(mesh,prefix);data=rows(prefix,'.cells.csv')
    errors={name:0. for name in ['rho','u','temperature']}
    for r in data:
        x=float(r['x']);area=float(r['area'])
        exact={'rho':1+mode[0]*math.cos(2*math.pi*x),'u':mode[1]*math.sin(2*math.pi*x),
               'temperature':1+mode[2]*math.cos(2*math.pi*x)}
        for name in errors:errors[name]+=area*abs(float(r[name])-exact[name])/(.25*amplitude)
    orders={key:math.log2(previous[key]/errors[key]) for key in errors} if previous else {}
    # Convergence in all three coupled quantities is the feature gate; min 1.5
    # separates second-order behavior with limiting/finite-grid effects from a
    # first-order transport regression. 16/32/64 makes this a bounded local test.
    if nx==64:assert min(orders.values())>1.5,(errors,orders)
    previous=errors
    print('Euler-Fourier mode',nx,errors,orders,flush=True)
    history=rows(prefix,'.history.csv')
    energy=[float(r['totalEnergy']) for r in history]
    assert max(energy)-min(energy)<1024*math.ulp(1.)*max(energy),'periodic energy drift'
    report['refinement'].append(dict(nx=nx,normalizedL1=errors,observedOrder=orders,audit=audit))
report.update(restartByteIdentical=True,legacyCheckpointAccepted=True,thermalMutationsRejected=True,
              continuumMode=dict(time=end,amplitude=amplitude,conductivity=.5,amplitudes=mode),
              scope='constant Fourier total-energy transport; bounded continuum small-signal and discrete checks, not general engineering qualification')
(root/'conduction-validation.json').write_text(json.dumps(report,indent=2)+'\n')
print('Euler-Fourier wall budgets, independent fluxes, coupled refinement and physical restart checks passed',flush=True)
if temporary:temporary.cleanup()
