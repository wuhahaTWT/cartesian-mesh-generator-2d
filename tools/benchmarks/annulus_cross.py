#!/usr/bin/env python3
"""Matched original-polygon annulus: real native and external FVM runs.

Only orchestration, interchange and analytic comparison; no duplicate FV PDE
assembly. Circle Couette errors are diagnostics for the polygon/facet model.
"""
from pathlib import Path
import argparse
import csv
import json
import math
import re
import shutil
import sys
from openfoam import prepare, dictionary
from cutcell_cross import run, dump

ROOT=Path(__file__).resolve().parents[2]
STUDY=ROOT/'outputs/annulus-cross'
sys.path.insert(0,str(ROOT/'tools/flow'))
import native_mesh as nm
import openfoam_data as foam
SIMPLE_IMAGE='opencfd/openfoam-run@sha256:4229997e74defb81548222d511b8e3b95b98305e5df41b8e88b031813fe47eeb'
COUPLED_IMAGE='pmariappan/foam-extend-4.0@sha256:63ab10ff5d5886d0428c76e56a50181d4e078a024ce4482c712e5b6e63893143'

def cases():
    return {
      'original':(ROOT/'build/laminar-fixtures/annulus-l7.solver.cm2d',ROOT/'outputs/laminar-foundation/foam-annulus-l7'),
      'geos96':(STUDY/'geos96.solver.cm2d',STUDY/'geos96.foam'),
      'geos128':(STUDY/'geos-canonical128.solver.cm2d',STUDY/'geos-canonical128.foam'),
      'hybrid':(STUDY/'preserved-hybrid.solver.cm2d',STUDY/'preserved-hybrid.foam'),
      'cfmesh':(STUDY/'cfmesh-native.solver.cm2d',STUDY/'cfmesh-precise'),
    }

def prepare_case(name,solver,destination):
    mesh,source=cases()[name]
    prepare('annulus',source,mesh,destination,'upwind',4000)
    if solver=='simple':return
    dictionary(destination,'constant','transportProperties','transportModel Newtonian;\nnu nu [0 2 -1 0 0 0 0] .1;')
    dictionary(destination,'constant','RASProperties','RASModel laminar; turbulence off; printCoeffs on;')
    dictionary(destination,'system','fvSchemes','''ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) Gauss upwind; div((nuEff*dev(T(grad(U))))) Gauss linear; div(U) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }''')
    dictionary(destination,'system','fvSolution','''solvers {
 Up { solver BiCGStab; preconditioner Cholesky; tolerance 1e-11; relTol 0; minIter 1; maxIter 2000; }
 p { solver CG; preconditioner DIC; tolerance 1e-11; relTol 0; }
 U { solver BiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
}
blockSolver { convergence 1e-8; pRefCell 0; pRefValue 0; }
fieldBounds { U 1e100; p -1e100 1e100; }
relaxationFactors { U 0.999; p 1; }''')
    dictionary(destination,'system','controlDict','''application pUCoupledFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 2000; deltaT 1;
writeControl timeStep; writeInterval 25; purgeWrite 1; writeFormat ascii; writePrecision 17;
writeCompression off; timeFormat general; timePrecision 10; runTimeModifiable false;''')

def external(name,solver,directory,budget):
    image=SIMPLE_IMAGE if solver=='simple' else COUPLED_IMAGE
    env='/usr/lib/openfoam/openfoam2606/etc/bashrc' if solver=='simple' else '/root/foam/foam-extend-4.0/etc/bashrc'
    command=['docker','run','--rm','--user','0','--entrypoint','/bin/bash']
    if solver=='coupled':command+=['--platform','linux/amd64']
    command+=['-v',str(ROOT)+':/work','-w','/work/'+str(directory.relative_to(ROOT)),image,'-lc',
      f'source {env} >/dev/null 2>&1; checkMesh > checkMesh.log 2>&1; '
      'checkMesh -allGeometry -allTopology > checkMesh-extended.log 2>&1; '
      f'timeout {budget} '+('simpleFoam' if solver=='simple' else 'pUCoupledFoam')+' > solve.log 2>&1']
    return run(command,directory/'container.log',budget+30)

def exact(x,y):
    r=math.hypot(x,y);t=(1/r-r)/3
    return -y/r*t,x/r*t,(r*r-1/(r*r))/18-2/9*math.log(r)

def errors(rows,boundary_model='facet-wall'):
    area=math.fsum(r['area'] for r in rows)
    refs=[exact(r['x'],r['y']) for r in rows]
    gauge=math.fsum(r['area']*(r['p']-p[2]) for r,p in zip(rows,refs))/area
    velocity=math.sqrt(math.fsum(r['area']*((r['u']-p[0])**2+(r['v']-p[1])**2) for r,p in zip(rows,refs))/area)/.5
    pressure=math.sqrt(math.fsum(r['area']*(r['p']-p[2]-gauge)**2 for r,p in zip(rows,refs))/area)/.25
    radial=math.sqrt(math.fsum(r['area']*((r['x']*r['u']+r['y']*r['v'])/math.hypot(r['x'],r['y']))**2 for r in rows)/area)/.5
    pressure_span=-exact(.5,0)[2]
    point_pressure=[abs(r['p']-p[2]-gauge) for r,p in zip(rows,refs)]
    worst=max(range(len(rows)),key=lambda i:point_pressure[i])
    return dict(velocityL2OverUref=velocity,pressureL2OverUrefSquared=pressure,radialVelocityL2OverUref=radial,
                pressureL2OverAnalyticPressureSpan=pressure*.25/pressure_span,
                pressureMaxOverAnalyticPressureSpan=point_pressure[worst]/pressure_span,
                pressureWorstCell=int(rows[worst].get('cell',worst)),pressureWorstArea=rows[worst]['area'],
                velocityMaxOverUref=max(math.hypot(r['u']-p[0],r['v']-p[1]) for r,p in zip(rows,refs))/.5,
                removedPressureGauge=gauge,reference=('circular Couette restricted to polygon domain; prescribed analytic face means'
                    if boundary_model=='analytic-trace' else 'true circular Couette; polygon geometry and prescribed facet trace differ'),
                fieldDefinition='cell-centre value; volume-weighted errors; no spatial qualification inferred')

def collect_external(name,solver,directory):
    mesh,_=cases()[name];m=nm.read_cm2d(mesh);g=nm.measure(m)
    log=(directory/'solve.log').read_text()
    declared=re.findall(r'SIMPLE solution converged in (\d+) iterations',log) if solver=='simple' else re.findall(r'latestTime = ([0-9.eE+-]+)',log)
    times=re.findall(r'^Time = ([0-9.eE+-]+)\s*$',log,re.M)
    # OpenCFD announces "Floating point exception trapping enabled" at startup.
    # That announcement is not a solver exception.
    finite_log='FOAM FATAL' not in log and not re.search(r'^Floating point exception(?:\s|$)',log,re.M)
    bounding=bool(re.search(r'\.\s+Bounding\.',log))
    final=directory/declared[-1] if declared else None
    reported_convergence=bool(declared) if solver=='simple' else 'reached convergence criterion:' in log
    converged=bool(reported_convergence and declared and times and float(declared[-1])==float(times[-1]) and finite_log and not bounding
                   and (final/'U').exists() and (final/'p').exists())
    run_file=directory/'container.run.json'
    if run_file.exists():converged=converged and json.loads(run_file.read_text())['exitCode']==0
    boundary_model='analytic-trace' if (directory/'analytic-trace.json').exists() else 'facet-wall'
    record=dict(case=name,solver=solver,cells=len(m.cells),converged=converged,boundaryModel=boundary_model,
                status='converged' if converged else 'not-converged',
                iterations=float(declared[-1]) if declared else float(times[-1]) if times else None,
                boundingOccurred=bounding,meshSha256=nm.sha256_file(mesh),
                standardCheckMesh='Mesh OK' in (directory/'checkMesh.log').read_text() if (directory/'checkMesh.log').exists() else None,
                extendedCheckMesh=('Mesh OK' in (directory/'checkMesh-extended.log').read_text()
                    if (directory/'checkMesh-extended.log').exists() else None),
                boundaryFaceCounts={p['name']:p['nFaces'] for p in foam.read_boundary(directory/'constant/polyMesh/boundary')},
                normalization='external solver-specific; not equivalent to native strict gates')
    if converged:
        u=foam.vector_field(final/'U');p=foam.scalar_field(final/'p')
        if len(u)==1:u*=len(m.cells)
        if len(p)==1:p*=len(m.cells)
        assert len(u)==len(p)==len(m.cells)
        rows=[dict(cell=i,x=c[0],y=c[1],area=a,u=v[0],v=v[1],p=pr) for i,(c,a,v,pr) in enumerate(zip(g.centroids,g.areas,u,p))]
        assert all(math.isfinite(v) for r in rows for v in r.values())
        record['circleErrors']=errors(rows,boundary_model)
        with (directory/'comparison.cells.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=['cell','x','y','area','u','v','p']);writer.writeheader();writer.writerows(rows)
    continuity=re.findall(r'time step continuity errors : (.*)',log)
    record['lastContinuity']=continuity[-1] if continuity else None
    dump(directory/'result.json',record);return record

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase',choices=['prepare','run','collect'],required=True)
    parser.add_argument('--solver',choices=['simple','coupled'],required=True)
    parser.add_argument('--names',nargs='+',choices=list(cases()),required=True)
    parser.add_argument('--budget',type=int,default=180)
    parser.add_argument('--directory',type=Path,help='Explicit case directory for one case, e.g. an existing preserved setup')
    args=parser.parse_args()
    if args.directory and len(args.names)!=1:parser.error('explicit directory requires one case')
    for name in args.names:
        directory=(args.directory.resolve() if args.directory else STUDY/(name+'-external-'+args.solver))
        if args.phase=='prepare':prepare_case(name,args.solver,directory)
        elif args.phase=='run':print(name,external(name,args.solver,directory,args.budget),flush=True)
        else:print(json.dumps(collect_external(name,args.solver,directory)),flush=True)
