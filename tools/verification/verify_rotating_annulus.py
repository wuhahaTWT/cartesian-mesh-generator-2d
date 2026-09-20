#!/usr/bin/env python3
"""Concentric circular Couette verification on independently constructed polygons.

These ring quadrilaterals isolate the flow solver; they are NOT Cut-cell output.
Ri=.5, Ro=1, inner Omega=1, outer stationary, nu=.1, no body force.
The polygon walls converge to circles as radial and angular spacing refine.
Reference: Richard Fitzpatrick, UT Austin, Taylor-Couette Flow, eq.10.35-40.
"""
import argparse
import json
import shutil
import time
import verify_native_flow as native

REFERENCE='https://farside.ph.utexas.edu/teaching/336L/Fluidhtml/node137.html'
import math
from pathlib import Path

def write_polar_mesh(path,nr,nt):
 vertices=[((.5+.5*j/nr)*math.cos(2*math.pi*i/nt),(.5+.5*j/nr)*math.sin(2*math.pi*i/nt)) for j in range(nr+1) for i in range(nt)]
 edges=[];cells=[];pairs={}
 for j in range(nr):
  for i in range(nt):
   ids=[j*nt+i,(j+1)*nt+i,(j+1)*nt+(i+1)%nt,j*nt+(i+1)%nt];incident=[]
   for a,b in zip(ids,ids[1:]+ids[:1]):
    key=tuple(sorted((a,b)))
    if key not in pairs:
     pairs[key]=len(edges);edges.append([len(edges),a,b,len(cells),-1,1])
    else:edges[pairs[key]][4:]=[len(cells),0]
    incident.append(pairs[key])
   area=.5*math.fsum(vertices[a][0]*vertices[b][1]-vertices[b][0]*vertices[a][1] for a,b in zip(ids,ids[1:]+ids[:1]))
   cells.append([len(cells),0,0,area,4,*ids,4,*incident])
 lines=['CM2D 1',f'VERTICES {len(vertices)}']+[f'{i} {x:.17g} {y:.17g}' for i,(x,y) in enumerate(vertices)]
 lines += [f'EDGES {len(edges)}',*(' '.join(map(str,e)) for e in edges),f'CELLS {len(cells)}',*(' '.join(map(str,c)) for c in cells),'AUDIT 0 0 0 0 0 0 0','END']
 path.write_text('\n'.join(lines)+'\n')


def metrics(mesh_path,prefix):
    mesh=native.read_cm2d(mesh_path);measured=native.measure(mesh,1e-11,1e-9)
    cells=native.read_cells(Path(str(prefix)+'.cells.csv'),mesh,measured,'custom')
    summary=json.loads(Path(str(prefix)+'.json').read_text())
    a,b=-1/3,1/3
    def pressure(r): return .5*a*a*r*r+2*a*b*math.log(r)-.5*b*b/(r*r)
    gauge=pressure(math.hypot(cells[0]['x'],cells[0]['y']))
    weights=[c['area'] for c in cells];velocity=[];radial=[];p_error=[]
    for c in cells:
        x,y=c['x'],c['y'];r=math.hypot(x,y);ut=a*r+b/r
        velocity.append(math.hypot(c['u']+ut*y/r,c['v']-ut*x/r)/.5)
        radial.append((c['u']*x+c['v']*y)/r/.5)
        p_error.append((c['p']-pressure(r)+gauge)/.25)
    loads={p['name']:p for p in summary['namedWallLoads']}
    expected=4*math.pi*.1/3
    return dict(cells=len(cells),velocityL2=native.weighted_l2(velocity,weights),
        radialVelocityL2=native.weighted_l2(radial,weights),pressureL2=native.weighted_l2(p_error,weights),
        rotorTorque=loads['rotor']['torque'],housingTorque=loads['housing']['torque'],
        rotorTorqueError=abs(loads['rotor']['torque']+expected)/expected,
        housingTorqueError=abs(loads['housing']['torque']-expected)/expected,
        torqueImbalance=abs(loads['rotor']['torque']+loads['housing']['torque'])/expected,
        iterations=summary['iterations'])


def study(args):
    root=args.output.resolve()
    if root.exists():raise ValueError('Choose a fresh output directory; existing failures are retained.')
    root.mkdir(parents=True)
    cli=args.cli.resolve(strict=True)
    report=dict(valid=False,reference=REFERENCE,meshType='concentric quadrilaterals, not Cut-cell',
        definition=dict(innerRadius=.5,outerRadius=1,innerOmega=1,outerOmega=0,nu=.1,
            velocityTheta='-r/3+1/(3*r)',pressure='r*r/18-2*log(r)/9-1/(18*r*r)+C',
            pressureGauge='cell 0',torque=4*math.pi*.1/3),
        predeclaredLimits=dict(fineVelocityL2=.005,finePressureL2=.005,fineTorqueRelative=.01,torqueImbalance=.005,
            tighterVelocityMaxChange=1e-6),controls=dict(tolerance=args.tolerance,tightTolerance=args.tolerance/100,
            maxIterations=6000,timeout=args.timeout),sourceSha256={str(p):native.sha256_file(p) for p in [cli,Path(__file__)]},
        runs=[],cases=[],issues=[])
    def save():native.write_json(root/'summary.json',report)
    def run(cmd,label):
        if shutil.disk_usage(root).free<5*1024**3:raise ValueError('Under 5 GiB free disk')
        start=time.monotonic();r=native.run(list(map(str,cmd)),root/'logs'/label,args.timeout)
        r['seconds']=time.monotonic()-start;report['runs'].append(r);save()
        print(label,r['status'],round(r['seconds'],2),flush=True)
        if r['status']!='passed':raise ValueError(label+' failed; original logs retained')
    def solve(mesh,bc,prefix,label,tol):
        run([cli,'--mesh',mesh,'--case','custom','--boundary',bc,'--output',prefix,'--nu',.1,'--speed',.5,
            '--convection','face-limited-linear','--pressure-preconditioner','aggregation',
            '--max-iterations',6000,'--tolerance',tol],label)
        audit=native.verify_case(mesh,prefix,'custom',.1,.5,native.argument_parser().parse_args(['--max-iterations','6000']))
        native.write_json(Path(str(prefix)+'.audit.json'),audit)
        if not audit['valid']:raise ValueError(label+': '+str(audit['issues']))
        return audit
    try:
        for nr in args.radial_cells:
            directory=root/f'n{nr}';directory.mkdir();mesh=directory/'ring.solver.cm2d'
            write_polar_mesh(mesh,nr,8*nr);bc=directory/'input.boundaries';prefix=directory/'flow'
            run([cli,'--mesh',mesh,'--case','annulus','--speed',.5,'--export-boundaries',bc],f'n{nr}-boundary')
            audit=solve(mesh,bc,prefix,f'n{nr}-flow',args.tolerance)
            item=dict(radialCells=nr,facets=8*nr,mesh=str(mesh),prefix=str(prefix),metrics=metrics(mesh,prefix),
                independentValid=True,artifactSha256=audit['artifactSha256'],meshSha256=audit['meshSha256'])
            report['cases'].append(item);save();print(item['metrics'],flush=True)
        keys=['velocityL2','pressureL2','rotorTorqueError','housingTorqueError']
        report['observedOrders']={k:[math.log(c['metrics'][k]/f['metrics'][k])/math.log(f['radialCells']/c['radialCells'])
            for c,f in zip(report['cases'],report['cases'][1:])] for k in keys}
        fine=report['cases'][-1]['metrics']
        for k,limit in [('velocityL2',.005),('pressureL2',.005),('rotorTorqueError',.01),('housingTorqueError',.01),('torqueImbalance',.005)]:
            if fine[k]>limit:report['issues'].append(k+' misses predeclared fine-grid target')
        for k in keys:
            if any(o<=0 for o in report['observedOrders'][k]):report['issues'].append(k+' does not decrease')
        tight=directory/'tight';solve(mesh,bc,tight,'fine-tight',args.tolerance/100)
        m=native.read_cm2d(mesh);measured=native.measure(m,1e-11,1e-9)
        base=native.read_cells(Path(str(prefix)+'.cells.csv'),m,measured,'custom')
        rows=native.read_cells(Path(str(tight)+'.cells.csv'),m,measured,'custom')
        sensitivity={k:max(abs(a[k]-b[k]) for a,b in zip(base,rows)) for k in ('u','v','p')}
        report['iterationSensitivity']=dict(maxAbs=sensitivity,metrics=metrics(mesh,tight))
        if max(sensitivity[k] for k in ('u','v'))>1e-6:report['issues'].append('iterative velocity uncertainty exceeds target')
        report['valid']=not report['issues'];save()
    except Exception as e:
        report['issues'].append(str(e));save();raise
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--radial-cells',type=int,nargs='+',default=[8,16,32]);p.add_argument('--timeout',type=int,default=180)
    p.add_argument('--cli',type=Path,default=native.REPO/'build/cartmesh2d_flow_cli');p.add_argument('--tolerance',type=float,default=1e-9);a=p.parse_args()
    if len(a.radial_cells)<3 or sorted(set(a.radial_cells))!=a.radial_cells or min(a.radial_cells)<4 or max(a.radial_cells)>64 or not 0<a.timeout<=600 or not 0<a.tolerance<=1e-8:p.error('Need >=3 increasing radial counts in 4..64 and timeout in 1..600')
    r=study(a);print(json.dumps(dict(valid=r['valid'],issues=r['issues'])));raise SystemExit(0 if r['valid'] else 1)
