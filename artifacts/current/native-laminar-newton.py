#!/usr/bin/env python3
"""Native analytic Newton with actual-residual backtracking, research only.

Python orchestrates native equations and linear algebra; it does not implement
another PDE. U=L=1, pressure U^2. The merit is half the squared L2 norm of all
free uncondensed equations in their existing dimensionless polynomial basis.
Original equation/divergence AND state-change convergence gates remain.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root=Path(__file__).resolve().parents[2]
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mesh');p.add_argument('n',type=int);p.add_argument('problem')
    p.add_argument('viscosity',type=float);p.add_argument('equation',choices=['ns','stokes'])
    p.add_argument('boundary',choices=['closed','open','pressure','traction','pseudo-traction','normal-stress'])
    p.add_argument('output',type=Path)
    p.add_argument('--transport',type=Path,default=root/'build/native-laminar-newton')
    p.add_argument('--block',type=Path,default=root/'build/native-laminar-block-precondition')
    p.add_argument('--fields',type=Path,default=root/'build/native-laminar-state-compare')
    p.add_argument('--backend',choices=['ilu0','dense-reference'],default='ilu0')
    p.add_argument('--preconditioner',choices=['picard','newton'],default='picard')
    p.add_argument('--pressure-preconditioner',choices=['viscous','local-oseen'],default='viscous')
    p.add_argument('--transport-speed',type=float,default=0.)
    p.add_argument('--initial-state',default='zero');p.add_argument('--iterations',type=int,default=40)
    p.add_argument('--backtracks',type=int,default=15);p.add_argument('--armijo',type=float,default=1e-4)
    p.add_argument('--no-line-search',action='store_true')
    p.add_argument('--tolerance',type=float,default=1e-9);p.add_argument('--linear-tolerance',type=float,default=1e-13)
    p.add_argument('--restarts',type=int,default=50);p.add_argument('--quadrature',type=int,choices=range(4,13),default=6)
    a=p.parse_args()
    if any(not math.isfinite(v) or v<=0 for v in [a.viscosity,a.tolerance,a.linear_tolerance,a.armijo]) or a.armijo>=.5:
        p.error('positive finite viscosity/tolerances and 0<Armijo<.5 required')
    if min(a.iterations,a.backtracks,a.restarts)<1 or a.linear_tolerance>.01:p.error('invalid iteration budget')
    if not math.isfinite(a.transport_speed) or a.transport_speed<0 or ((a.pressure_preconditioner=='local-oseen')!=(a.transport_speed>0)):
        p.error('local-oseen requires positive transport speed; viscous mode requires zero')
    if a.pressure_preconditioner=='local-oseen' and (a.boundary=='closed' or a.backend!='ilu0'):
        p.error('local-oseen requires an open-boundary ILU0 solve')
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    controls=[a.mesh,str(a.n),a.problem,str(a.viscosity),a.equation,a.boundary,str(a.quadrature)]
    start=time.monotonic()
    record={'method':'native analytic Newton','controls':controls,'completed':False,'pid':os.getpid(),'backend':a.backend,'velocityPreconditioner':a.preconditioner,'pressurePreconditioner':a.pressure_preconditioner,'transportSpeed':a.transport_speed,'lineSearch':not a.no_line_search,'armijo':a.armijo,'maximumBacktracks':a.backtracks,'maximumIterations':a.iterations,'maximumRestarts':a.restarts,'initialState':a.initial_state,'iterationTolerance':a.tolerance,'linearTolerance':a.linear_tolerance,'iterations':[],'commands':[],'lastAcceptedIterate':None,'finalState':None,'physicalCheckpoint':False}

    def save():
        record['totalSeconds']=time.monotonic()-start
        temporary=a.output/'result.json.tmp';temporary.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');temporary.replace(a.output/'result.json')

    def run(command,label,allow_failure=False):
        cmd=list(map(str,command));t=time.monotonic();r=subprocess.run(cmd,text=True,capture_output=True)
        (a.output/(label+'.stdout')).write_text(r.stdout);(a.output/(label+'.stderr')).write_text(r.stderr)
        record['commands'].append({'label':label,'argv':cmd,'exitCode':r.returncode,'seconds':time.monotonic()-t});save()
        if r.returncode and not allow_failure:raise RuntimeError(label+' failed; no candidate acceptance')
        return r

    def native(command,label):
        value=json.loads(run(command,label).stdout)
        (a.output/(label+'.json')).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
        return value

    def metrics(check):
        values=[check[k] for k in ('internalResidual','freeFaceResidual','maxDivergence','uncondensedResidualL2')]
        if not all(math.isfinite(v) and v>=0 for v in values):raise RuntimeError('nonfinite native equation/merit metric')
        return max(values[:3]),values[3]

    try:
        binaries={'transport':a.transport,'fields':a.fields}
        if a.backend!='dense-reference':binaries['block']=a.block
        record['nativeBinaries']={k:{'path':str(v),'sha256':sha(v)} for k,v in binaries.items()}
        names=['native-laminar-newton.cpp','native-laminar-newton.py','native-laminar-p1-oseen.cpp','native-laminar-p1-transport.cpp','native-laminar-open-boundary.cpp','native-laminar-p1-stress.cpp','native-laminar-hybrid-stokes-p1.cpp','native-laminar-block-precondition.cpp','native-laminar-state-compare.cpp']
        record['sourceSha256']={n:sha(root/'artifacts/current'/n) for n in names}
        if Path(a.mesh).is_file():record['meshSha256']=sha(a.mesh)
        if a.initial_state!='zero':record['initialStateSha256']=sha(a.initial_state)
        geometry=a.output/'geometry';record['geometry']=native([a.fields,'geometry',a.mesh,a.n,geometry],'geometry')
        current=a.output/'boundary-consistent-seed.state'
        record['seed']=native([a.transport,'seed',a.mesh,a.n,a.problem,a.boundary,a.quadrature,a.initial_state,current],'seed')
        check=native([a.transport,'check',*controls,current,a.output/'initial'],'initial.check')
        record['initialCheck']=check;_,oldnorm=metrics(check);save()
        for iteration in range(a.iterations):
            label=f'iteration-{iteration:03d}';prefix=a.output/label
            step={'iteration':iteration,'previous':str(current),'previousSha256':sha(current),'previousMerit':.5*oldnorm*oldnorm,'trials':[]};record['iterations'].append(step);save()
            native_args=[*controls,current,prefix]
            assembly=native([a.transport,'assemble',*native_args],label+'.assemble');step['assembly']=assembly
            if a.backend=='dense-reference':
                import numpy as np
                size=assembly['unknowns']
                if size>2175:raise RuntimeError('dense reference limited to 2175 unknowns')
                z=np.fromfile(str(prefix)+'.entries',dtype=[('i','<i8'),('j','<i8'),('v','<f8')]);matrix=np.zeros((size,size));np.add.at(matrix,(z['i'],z['j']),z['v'])
                rhs=np.fromfile(str(prefix)+'.rhs',dtype='<f8');x=np.linalg.solve(matrix,rhs)
                for _ in range(2):x+=np.linalg.solve(matrix,rhs-matrix@x)
                residual=rhs-matrix@x;norm=float(np.linalg.norm(rhs));rel=float(np.linalg.norm(residual))/(norm if norm else 1.)
                ok=bool(np.isfinite(x).all() and math.isfinite(rel) and rel<=a.linear_tolerance)
                x.astype('<f8').tofile(str(prefix)+('.solution' if ok else '.candidate'))
                step['linear']={'converged':ok,'relative_residual':rel,'method':'dense-reference'}
                del z,matrix,rhs,x,residual
                code=0 if ok else 2
            else:
                preconditioner_args=['--velocity-preconditioner',str(prefix)+'.picard'] if a.preconditioner=='picard' else []
                if preconditioner_args:step['preconditionerEntriesSha256']=sha(str(prefix)+'.picard.entries')
                pressure_mode='outlet-oseen' if a.pressure_preconditioner=='local-oseen' else 'gauge' if a.boundary=='closed' else 'outlet'
                r=run([a.block,prefix,str(geometry)+'.cells.csv','ilu0',pressure_mode,prefix,a.linear_tolerance,a.restarts,a.viscosity,a.transport_speed,*preconditioner_args],label+'.linear',True)
                code=r.returncode
                if Path(str(prefix)+'.json').exists():step['linear']=json.loads(Path(str(prefix)+'.json').read_text())
            step['linearExitCode']=code;save()
            if code or not step['linear']['converged']:
                record['failure']='linear solve failed or exhausted its budget; no candidate recovery';save();return 2
            recovery=native([a.transport,'recover',*native_args],label+'.recover');step['recovery']=recovery
            if not recovery['linearizedNewton'] or recovery['physicalDiagnostics']:raise RuntimeError('linearized candidate was mislabeled')
            candidate=Path(str(prefix)+'.recover.state');accepted=None
            for bt in range(1 if a.no_line_search else a.backtracks):
                alpha=2.**(-bt);trial_label=f'{label}-trial-{bt:02d}';trial_prefix=a.output/trial_label
                if bt==0:
                    trial_state=candidate;change=recovery['stateCoefficientChange']
                else:
                    trial_state=Path(str(trial_prefix)+'.state')
                    blend=native([a.transport,'blend',a.mesh,a.n,current,candidate,alpha,trial_state],trial_label+'.blend');change=blend['stateCoefficientChange']
                check=native([a.transport,'check',*controls,trial_state,trial_prefix],trial_label+'.check')
                equation,norm=metrics(check)
                if not math.isfinite(change) or change<0:raise RuntimeError('invalid state change')
                converged=equation<=a.tolerance and (a.equation=='stokes' or change<=a.tolerance)
                merit=.5*norm*norm;bound=(.5-a.armijo*alpha)*oldnorm*oldnorm
                sufficient=merit<=bound
                accept=sufficient or converged or a.no_line_search
                trial={'alpha':alpha,'state':str(trial_state),'stateSha256':sha(trial_state),'stateCoefficientChange':change,'check':check,'equationResidualMax':equation,'merit':merit,'armijoBound':bound,'sufficientDecrease':sufficient,'strictConvergence':converged,'accepted':accept}
                step['trials'].append(trial);save()
                if accept:
                    accepted=trial;current=trial_state;oldnorm=norm;record['lastAcceptedIterate']=str(current);step['acceptedState']=str(current);save();break
            if accepted is None:
                record['failure']='line-search budget exhausted; all trial fields retained, last accepted iterate unchanged';save();return 2
            print(json.dumps({'iteration':iteration,'alpha':accepted['alpha'],'residual':accepted['equationResidualMax'],'change':accepted['stateCoefficientChange'],'trials':len(step['trials'])}),flush=True)
            if accepted['strictConvergence']:
                record['completed']=True;record['finalState']=str(current);record['finalStateSha256']=sha(current);save();return 0
        record['failure']='nonlinear iteration budget exhausted; no final acceptance';save();return 2
    except KeyboardInterrupt:
        record['failure']='interrupted; no final acceptance';save();return 130
    except Exception as exc:
        record['failure']=str(exc);save();raise


if __name__=='__main__':
    raise SystemExit(main())
