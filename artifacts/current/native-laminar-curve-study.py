"""Native CLI circular-Couette joint geometry/mesh refinement; field postprocessing only."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

ROOT = Path('outputs/cloud-laminar')

def metrics(prefix):
    summary = json.loads(Path(str(prefix)+'.json').read_text())
    if not summary.get('converged'):
        return {'converged': False, 'summary': summary}
    rows = list(csv.DictReader(open(str(prefix)+'.cells.csv')))
    area = sum(float(r['area']) for r in rows)
    samples = []
    for r in rows:
        x,y,a,u,v,p = (float(r[k]) for k in ['x','y','area','u','v','p'])
        r2=x*x+y*y; q=(1/r2-1)/3
        pe=(r2/18-math.log(r2)/9-1/(18*r2))*summary['momentumInertia']
        samples.append((a, math.hypot(u+y*q,v-x*q), p-pe, math.sqrt(r2)))
    offset=sum(a*p for a,u,p,r in samples)/area
    interior=[s for s in samples if .6 <= s[3] <= .9]
    rotor=next(x for x in summary['namedWallLoads'] if x['name']=='rotor')
    # Laplacian-only traction excludes the transpose gradient; its torque is different.
    exact_torque=-4*math.pi*.1/3 if summary['viscousStress']=='symmetric' else -2*math.pi*.1*5/12
    return dict(converged=True, cells=len(rows), evaluations=summary['iterations'],
        momentumInertia=summary['momentumInertia'],viscousStress=summary['viscousStress'],
        fluidArea_m2=area, circleAreaDifference_m2=area-math.pi*.75,
        velocityRms_m_s=math.sqrt(sum(a*u*u for a,u,p,r in samples)/area),
        velocityMax_m_s=max(u for a,u,p,r in samples),
        pressureRms_m2_s2=math.sqrt(sum(a*(p-offset)**2 for a,u,p,r in samples)/area),
        pressureMax_m2_s2=max(abs(p-offset) for a,u,p,r in samples),
        pressureGaugeOffset_m2_s2=offset,
        fixedInteriorPressureRms_m2_s2=math.sqrt(sum(a*(p-offset)**2 for a,u,p,r in interior)/sum(a for a,u,p,r in interior)),
        rotorTorque=rotor['torque'],exactRotorTorque=exact_torque,
        rotorTorqueRelativeError=abs(rotor['torque']-exact_torque)/abs(exact_torque),
        continuity=summary['continuity'],momentumResidual=summary['momentumResidual'],
        cellFieldSha256=hashlib.sha256(Path(str(prefix)+'.cells.csv').read_bytes()).hexdigest())

def run(command, log, timeout):
    started=time.monotonic()
    with open(log,'w') as out:
        try:
            r=subprocess.run(command,stdout=out,stderr=subprocess.STDOUT,timeout=timeout)
            code=r.returncode;status='completed'
        except subprocess.TimeoutExpired:
            code=None;status='timed-out-no-solver-qualification'
    return dict(command=command,returncode=code,status=status,seconds=time.monotonic()-started,log=str(log))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--levels',nargs='+',type=int,default=[4,5,6,7]);parser.add_argument('--schemes',nargs='+',default=['upwind','face-limited-linear']);args=parser.parse_args()
    ROOT.mkdir(parents=True,exist_ok=True);records=[]
    for level in args.levels:
        n=2**(level+2);prefix=ROOT/f'annulus-joint-{level}';source=Path(str(prefix)+'.xy')
        with source.open('w') as f:
            for radius in [1.,.5]:
                for i in range(n):f.write(f'{radius*math.cos(2*math.pi*i/n):.17g} {radius*math.sin(2*math.pi*i/n):.17g}\n')
                f.write('\n')
        row=dict(level=level,segmentsPerRing=n,inputSha256=hashlib.sha256(source.read_bytes()).hexdigest())
        row['mesh']=run(['build/cartmesh2d_cli',str(source),str(prefix),str(level),'.15','.1','interior',str(prefix)+'-foam',str(level)],str(prefix)+'-mesh.log',300)
        row['solves']=[]
        if row['mesh']['returncode']==0:
            boundary=Path(str(prefix)+'.boundaries')
            row['boundary']=run(['build/cartmesh2d_flow_cli','--mesh',str(prefix)+'.solver.cm2d','--case','annulus','--speed','.5','--export-boundaries',str(boundary)],str(prefix)+'-boundary.log',30)
            if row['boundary']['returncode']!=0:
                records.append(row)
                (ROOT/'curve-joint-refinement.json').write_text(json.dumps(records,indent=2)+'\n')
                continue
            for scheme in args.schemes:
                out=Path(str(prefix)+'-'+scheme)
                result=run(['build/cartmesh2d_flow_cli','--mesh',str(prefix)+'.solver.cm2d','--case','custom','--boundary',str(boundary),'--nu','.1','--speed','.5','--convection',scheme,'--tolerance','1e-8','--max-iterations','1800','--output',str(out)],str(out)+'.log',240)
                result['scheme']=scheme
                if result['returncode']==0:result['metrics']=metrics(out)
                row['solves'].append(result)
                print(json.dumps(dict(level=level,**result)),flush=True)
        records.append(row)
        (ROOT/'curve-joint-refinement.json').write_text(json.dumps(records,indent=2)+'\n')
