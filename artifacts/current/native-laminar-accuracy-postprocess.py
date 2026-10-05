#!/usr/bin/env python3
"""Compare saved native fields to analytic/published references; no PDE audit."""
import bisect
import csv
import hashlib
import json
import math
from pathlib import Path
ROOT = Path('outputs/cloud-laminar')

def read_fields(prefix):
    with Path(str(prefix)+'.cells.csv').open() as stream:
        return [{k: float(v) for k,v in row.items()} for row in csv.DictReader(stream)]

def compare_fields(first,second):
    a,b=read_fields(first),read_fields(second)
    assert len(a)==len(b)
    for x,y in zip(a,b):
        assert all(x[k]==y[k] for k in ['x','y','area'])
    area=sum(x['area'] for x in a)
    shift=sum(x['area']*(y['p']-x['p']) for x,y in zip(a,b))/area
    return {'velocityMax_m_s':max(math.hypot(y['u']-x['u'],y['v']-x['v']) for x,y in zip(a,b)),
            'pressureMaxGaugeAligned_m2_s2':max(abs(y['p']-x['p']-shift) for x,y in zip(a,b)),
            'pressureGaugeShift_m2_s2':shift}

def cavity_comparisons(rows):
    reference_path=Path('tools/verification/references/cavity-marchi-2009-re100.json')
    reference=json.loads(reference_path.read_text());out=[]
    for row in rows:
        r=row['result']
        if r['case']!='cavity': continue
        n=r['n'];fields=read_fields(row['command'][-1]);axis=[0]+[(i+.5)/n for i in range(n)]+[1]
        def sample(x,y,key):
            i=min(n,bisect.bisect_right(axis,x)-1);j=min(n,bisect.bisect_right(axis,y)-1)
            a=(x-axis[i])/(axis[i+1]-axis[i]);b=(y-axis[j])/(axis[j+1]-axis[j])
            def value(ix,iy):
                if ix in [0,n+1] or iy==0: return 0
                if iy==n+1: return 1 if key=='u' else 0
                return fields[(iy-1)*n+ix-1][key]
            return sum(value(i+dx,j+dy)*wx*wy for dx,wx in [(0,1-a),(1,a)] for dy,wy in [(0,1-b),(1,b)])
        for key in ['u','v']:
            samples=[]
            for point in reference[key]:
                q=point['coordinate'];pred=sample(.5,q,key) if key=='u' else sample(q,.5,key)
                samples.append({'coordinate':q,'computed':pred,'reference':point['value'],'error':pred-point['value']})
            out.append({'n':n,'scheme':r['scheme'],'component':key,'rms_m_s':math.sqrt(sum(x['error']**2 for x in samples)/len(samples)),
                        'max_m_s':max(abs(x['error']) for x in samples),'samples':samples})
    return {'source':str(reference_path),'sourceSha256':hashlib.sha256(reference_path.read_bytes()).hexdigest(),
            'interpolation':'bilinear cell-centre interpolation plus prescribed wall values; no extrapolation',
            'normalization':'lid speed = 1 m/s; dimensionless error equals numeric m/s error','comparisons':out}

def annulus_comparison(prefix):
    rows=read_fields(prefix);samples=[];area=0;offset=0;v2=0;vmax=0
    ri,ro,speed,nu=.5,1.,.5,.1
    B=speed*ri*ro*ro/(ro*ro-ri*ri);A=-B/(ro*ro)
    for row in rows:
        x,y,a,u,v,p=[row[k] for k in ['x','y','area','u','v','p']];r=math.hypot(x,y)
        ut=A*r+B/r;ue=-ut*y/r;ve=ut*x/r;pe=.5*A*A*r*r+2*A*B*math.log(r)-.5*B*B/(r*r)
        er=math.hypot(u-ue,v-ve);samples.append((a,p,pe,row));area+=a;offset+=a*(p-pe);v2+=a*er*er;vmax=max(vmax,er)
    offset/=area;p2=sum(a*(p-pe-offset)**2 for a,p,pe,_ in samples)
    ranked=sorted(samples,key=lambda item:abs(item[1]-item[2]-offset),reverse=True)
    summary=json.loads(Path(str(prefix)+'.json').read_text());loads=summary['namedWallLoads'];torque=next(x['torque'] for x in loads if x['name']=='rotor')
    exact=-4*math.pi*nu*B
    return {'prefix':str(prefix),'cells':len(rows),'converged':summary['converged'],'evaluations':summary['iterations'],
            'reference':'circular Couette, ri=.5 m, ro=1 m, inner speed=.5 m/s, nu=.1 m2/s',
            'velocityFormula':'u_theta=A*r+B/r; A=-1/3 s^-1, B=1/3 m2/s',
            'pressureFormula':'p/rho=A^2*r^2/2+2*A*B*log(r)-B^2/(2*r^2)',
            'qualification':'polygon walls compared with circular reference; both geometry and discretization error included',
            'fluidArea_m2':area,'velocityRms_m_s':math.sqrt(v2/area),'velocityMax_m_s':vmax,
            'velocityRms_over_innerSpeed':math.sqrt(v2/area)/speed,'pressureRms_m2_s2':math.sqrt(p2/area),
            'pressureMax_m2_s2':abs(ranked[0][1]-ranked[0][2]-offset),'pressureRms_over_speedSquared':math.sqrt(p2/area)/speed**2,
            'pressureGaugeOffset_m2_s2':offset,'rotorTorque':torque,'exactRotorTorque':exact,'rotorTorqueRelativeError':abs(torque-exact)/abs(exact),
            'worstPressureCells':[dict(item[3],pressureError_m2_s2=item[1]-item[2]-offset) for item in ranked[:8]]}

rows=json.loads((ROOT/'accuracy-runs.json').read_text())
cavity=cavity_comparisons(rows)
(ROOT/'cavity-reference-comparison.json').write_text(json.dumps(cavity,indent=2)+'\n')
annuli=[annulus_comparison(ROOT/name) for name in ['annulus-default','annulus-tighter','annulus-l5-flow','annulus-l6-flow'] if (ROOT/(name+'.cells.csv')).exists()]
controls=[]
for row in json.loads((ROOT/'accuracy-controls.json').read_text()):
    r=row['result'];baseline=ROOT/f"{r['case']}-{r['scheme']}-{r['n']}"
    controls.append({'case':r['case'],'scheme':r['scheme'],'tolerance':r['tolerance'],'converged':r['converged'],
                     'evaluations':r['evaluations'],'seconds':r['seconds'],**compare_fields(baseline,row['command'][-1])})
costs=[]
if (ROOT/'cost-runs.json').exists():
    for row in json.loads((ROOT/'cost-runs.json').read_text()):
        r=row['result'];baseline=ROOT/f"{r['case']}-{r['scheme']}-{r['n']}"
        default=json.loads(Path(str(baseline)+'.json').read_text())
        costs.append({'case':r['case'],'scheme':r['scheme'],'simple':r,'default':default,
                      'fieldComparison':compare_fields(baseline,row['command'][-1]),
                      'matchedConvergedSolutions':bool(r['converged'] and default['converged']),
                      'timingQualification':'complete native solve and complete subprocess recorded; runs overlapped other cloud tasks; no universal speedup claim'})
report={'status':'accuracy-investigation-in-progress-not-mature','baseCommit':'f2e0469','units':{'velocity':'m/s','pressure':'p/rho in m2/s2','continuity':'native dimensionless residual'},
        'thresholds':'No new acceptance threshold is introduced. Errors are measured against references and mesh trends; original solver gates unchanged.',
        'nativeGridRuns':rows,'cavity':cavity,'annulus':annuli,'iterationControlComparisons':controls,'completeCostComparisons':costs,
        'annulusControlComparison':compare_fields(ROOT/'annulus-default',ROOT/'annulus-tighter'),
        'limitations':['Curved-wall local pressure errors remain unresolved.','Single circular-Couette comparison does not by itself establish mesh-independent curved-wall accuracy.','Headless Linux App is separate from physical macOS GUI and packaged-platform qualifications.']}
quality_path=ROOT/'annulus-l6.failed.solver-quality.json'
if quality_path.exists():
    quality=json.loads(quality_path.read_text())
    report['coarseCurveEntrance']={'level6':{'status':'rejected-by-original-Solver-quality','issueCount':quality['issue_count'],'policy':quality['policy'],'metrics':quality['metrics'],'failureMesh':'outputs/cloud-laminar/annulus-l6.failed.solver.cm2d','log':'outputs/cloud-laminar/mesh-l6-solver.log'},'level5':json.loads((ROOT/'mesh-l5-cancelled.json').read_text())}
report['sourcesSha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path('artifacts/current/native-laminar-accuracy-probe.cpp'),Path(__file__),Path('tools/verification/references/cavity-marchi-2009-re100.json')]}
Path('artifacts/current/native-laminar-cloud-accuracy.json').write_text(json.dumps(report,indent=2)+'\n')
