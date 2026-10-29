#!/usr/bin/env python3
"""Controlled native/AMReX geometry x native/OpenFOAM FVM experiment.

Orchestration, mesh/field comparison and exact-solution error only; no second
implementation of the finite-volume equations. See DEVELOPMENT_CN.md.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/flow'))
import native_mesh as nm
import openfoam_data as foam
from openfoam import dictionary

OFFSETS={'mild':.2375,'small':.246875,'tiny':.24996875}
IMAGE='opencfd/openfoam-run@sha256:4229997e74defb81548222d511b8e3b95b98305e5df41b8e88b031813fe47eeb'

def dump(path,data):
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')

def run(command,log,budget=60):
    start=time.monotonic()
    with log.open('w') as stream:
        try:
            result=subprocess.run(list(map(str,command)),cwd=ROOT,stdout=stream,
                                  stderr=subprocess.STDOUT,timeout=budget)
            code=result.returncode
        except subprocess.TimeoutExpired:
            code='timeout'
    record=dict(command=list(map(str,command)),exitCode=code,wallSeconds=time.monotonic()-start,budgetSeconds=budget)
    dump(log.with_suffix('.run.json'),record)
    return record

def cases(root):
    for severity,b in OFFSETS.items():
        for source in ('native','amrex'):
            for merged in (False,True):
                name=f'{source}-{severity}-'+('merged' if merged else 'raw')
                yield name,root/name,b,source,merged

def generate(root,adapter):
    root.mkdir(parents=True,exist_ok=False)
    for name,p,b,source,merged in cases(root):
        r=run([adapter,source,16,b,p,int(merged)],p.with_suffix('.generate.log'))
        if r['exitCode'] not in (0,3):raise RuntimeError(f'geometry failed: {name}')
    comparisons=[]
    for severity in OFFSETS:
        for mode in ('raw','merged'):
            left=root/f'native-{severity}-{mode}.solver.cm2d'
            right=root/f'amrex-{severity}-{mode}.solver.cm2d'
            a,b=nm.read_cm2d(left),nm.read_cm2d(right)
            ma,mb=nm.measure(a),nm.measure(b)
            # Match polygons by sorted rounded coordinates for comparison only.
            # No rounding enters generation, connectivity, boundary data or solve.
            def records(mesh,measure):
                return {tuple(sorted((round(mesh.vertices[v][0],12),round(mesh.vertices[v][1],12))
                                     for v in c.vertices)):(measure.areas[c.id],measure.centroids[c.id])
                        for c in mesh.cells}
            ra,rb=records(a,ma),records(b,mb)
            common=ra.keys()&rb.keys()
            comparisons.append(dict(severity=severity,mode=mode,cellsNative=len(a.cells),cellsAMReX=len(b.cells),
                matchedPolygons=len(common),allPolygonsMatched=ra.keys()==rb.keys(),
                maxMatchedAreaDifference=max(abs(ra[k][0]-rb[k][0]) for k in common),
                maxMatchedCentreDifference=max(math.dist(ra[k][1],rb[k][1]) for k in common)))
    dump(root/'geometry-comparison.json',comparisons)

def native(root):
    for name,p,b,source,merged in cases(root):
        metrics=json.loads(p.with_suffix('.metrics.json').read_text())
        if metrics['qualityIssues']:continue
        mesh=nm.read_cm2d(Path(str(p)+'.solver.cm2d'));measure=nm.measure(mesh)
        initial=Path(str(p)+'.rest.csv')
        if initial.exists():raise ValueError('existing initial iterate')
        with initial.open('w') as f:
            writer=csv.writer(f);writer.writerow(['cell','x','y','u','v','p'])
            for i,(x,y) in enumerate(measure.centroids):writer.writerow([i,x,y,0,0,0])
        for coupling in ('simple','coupled'):
            output=Path(str(p)+'-rest-'+coupling)
            if output.with_suffix('.run.json').exists():raise ValueError('existing native run')
            command=[ROOT/'build/cartmesh2d_flow_cli','--mesh',str(p)+'.solver.cm2d',
                     '--case','custom','--boundary',str(p)+'.boundaries','--nu','.1','--speed','1',
                     '--convection','upwind','--viscous-stress','symmetric','--coupling',coupling,
                     '--linear-policy','strict','--max-iterations','1500','--tolerance','1e-8',
                     '--initial-guess',initial,'--profile','--output',output]
            run(command,output.with_suffix('.log'))

def prepare(root):
    for name,p,b,source,merged in cases(root):
        source_mesh=Path(str(p)+'.foam/constant/polyMesh')
        target=Path(str(p)+'-openfoam');target.mkdir(exist_ok=False)
        mesh_dir=target/'constant/polyMesh';mesh_dir.mkdir(parents=True)
        points=foam.read_points(source_mesh/'points');faces=foam.read_faces(source_mesh/'faces')
        owners=foam.read_labels(source_mesh/'owner');neighbours=foam.read_labels(source_mesh/'neighbour')
        groups={}
        for fid in range(len(neighbours),len(faces)):
            face=faces[fid];area=foam.face_area(face,points)
            x,y=[sum(points[v][axis] for v in face)/len(face) for axis in (0,1)]
            if abs(area[2])>0:kind='empty';u=(0,0,0)
            elif abs(y-.5*x-b)<1e-12:kind='wall';u=(0,0,0)
            else:
                shear=(y-.5*x-b)/1.25;u=(shear,.5*shear,0)
                kind='inlet' if u[0]*area[0]+u[1]*area[1]<0 else 'velocityOutlet'
            groups.setdefault(kind,[]).append((fid,u))
        order=list(range(len(neighbours)))+[fid for group in groups.values() for fid,u in group]
        if sorted(order)!=list(range(len(faces))):raise ValueError('nonbijective face permutation')
        for f in ('points','neighbour'):shutil.copyfile(source_mesh/f,mesh_dir/f)
        (mesh_dir/'faces').write_text(foam.header('faceList','constant/polyMesh','faces')+str(len(faces))+'\n(\n'+
            ''.join(str(len(faces[i]))+'('+' '.join(map(str,faces[i]))+')\n' for i in order)+')\n')
        (mesh_dir/'owner').write_text(foam.header('labelList','constant/polyMesh','owner')+str(len(owners))+'\n(\n'+
            ''.join(str(owners[i])+'\n' for i in order)+')\n')
        boundary=[];ub=[];pb=[];start=len(neighbours)
        for kind,group in groups.items():
            patch='empty' if kind=='empty' else 'wall' if kind=='wall' else 'patch'
            boundary.append(f'{kind} {{ type {patch}; nFaces {len(group)}; startFace {start}; }}\n');start+=len(group)
            if kind=='empty':uentry=pentry='type empty;'
            else:
                uentry='type fixedValue; value nonuniform List<vector>\n'+str(len(group))+'\n(\n'+''.join(
                    '('+' '.join(f'{v:.17g}' for v in u)+')\n' for fid,u in group)+');'
                pentry='type zeroGradient;'
            ub.append(f'{kind} {{ {uentry} }}\n');pb.append(f'{kind} {{ {pentry} }}\n')
        (mesh_dir/'boundary').write_text(foam.header('polyBoundaryMesh','constant/polyMesh','boundary')+
            str(len(groups))+'\n(\n'+''.join(boundary)+')\n')
        (target/'0').mkdir()
        for field,cls,dim,initial,bc in [('U','volVectorField','0 1 -1 0 0 0 0','(0 0 0)',ub),
                                       ('p','volScalarField','0 2 -2 0 0 0 0','0',pb)]:
            (target/'0'/field).write_text(foam.header(cls,'0',field)+f'dimensions [{dim}];\ninternalField uniform {initial};\nboundaryField {{\n'+''.join(bc)+'}\n')
        dictionary(target,'constant','transportProperties','transportModel Newtonian;\nnu [0 2 -1 0 0 0 0] .1;')
        dictionary(target,'constant','turbulenceProperties','simulationType laminar;')
        dictionary(target,'system','fvSchemes','''ddtSchemes { default steadyState; }
gradSchemes { default leastSquares; }
divSchemes { default none; div(phi,U) bounded Gauss upwind; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
wallDist { method meshWave; }''')
        dictionary(target,'system','fvSolution','''solvers
{
 p { solver PCG; preconditioner DIC; tolerance 1e-11; relTol 0; maxIter 3000; }
 U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0; }
}
SIMPLE { nNonOrthogonalCorrectors 3; pRefCell 0; pRefValue 0; residualControl { p 1e-8; U 1e-8; } }
relaxationFactors { fields { p 0.25; } equations { U 0.6; } }''')
        dictionary(target,'system','controlDict','''application simpleFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 1500; deltaT 1;
writeControl timeStep; writeInterval 1500; purgeWrite 1; writeFormat ascii; writePrecision 17;
writeCompression off; timeFormat general; timePrecision 10; runTimeModifiable false;''')
        dump(target/'conversion.json',dict(mesh=str(p)+'.solver.cm2d',sha256=nm.sha256_file(Path(str(p)+'.solver.cm2d')),
            cellNumberingUnchanged=True,pointsUnchanged=True,internalFacesUnchanged=True,boundaryPermutation=order[len(neighbours):]))

def run_foam(root):
    for name,p,b,source,merged in cases(root):
        target=Path(str(p)+'-openfoam');relative=target.relative_to(ROOT)
        command=['docker','run','--rm','--user','0','--entrypoint','/bin/bash',
                 '-v',str(ROOT)+':/work','-w','/work/'+str(relative),IMAGE,'-lc',
                 'source /usr/lib/openfoam/openfoam2606/etc/bashrc; '
                 'checkMesh > checkMesh.log 2>&1; '
                 'checkMesh -allGeometry -allTopology > checkMesh-extended.log 2>&1; '
                 'timeout 45 simpleFoam > solve.log 2>&1']
        run(command,target/'container.log',60)

def error_metrics(rows,b):
    area=sum(r['area'] for r in rows)
    gauge=sum(r['area']*r['p'] for r in rows)/area
    return dict(velocityL2=math.sqrt(sum(r['area']*((r['u']-(r['y']-.5*r['x']-b)/1.25)**2+
        (r['v']-.5*(r['y']-.5*r['x']-b)/1.25)**2) for r in rows)/area),
        pressureL2=math.sqrt(sum(r['area']*(r['p']-gauge)**2 for r in rows)/area))

def collect(root):
    records=[]
    for name,p,b,source,merged in cases(root):
        metrics=json.loads(p.with_suffix('.metrics.json').read_text())
        record=dict(name=name,geometry=metrics,native={},openfoam={})
        for coupling in ('simple','coupled'):
            output=Path(str(p)+'-rest-'+coupling);summary=output.with_suffix('.json')
            if metrics['qualityIssues']:record['native'][coupling]=dict(status='rejected-by-original-quality-gate');continue
            runfile=output.with_suffix('.run.json')
            result=dict(run=json.loads(runfile.read_text())) if runfile.exists() else {}
            if summary.exists():
                data=json.loads(summary.read_text())
                result.update({k:data[k] for k in ('status','iterations','converged','momentumResidual','continuity','velocityChange','pressureChange','strictLinearFinal','globalRelativeImbalance')})
                with output.with_suffix('.cells.csv').open() as f:rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
                result['exactErrors']=error_metrics(rows,b)
            else:result.update(status='failed',lastLog=output.with_suffix('.log').read_text().splitlines()[-3:])
            record['native'][coupling]=result
        target=Path(str(p)+'-openfoam');log=(target/'solve.log').read_text()
        declared=re.findall(r'SIMPLE solution converged in (\d+) iterations',log)
        times=re.findall(r'^Time = ([0-9.eE+-]+)\s*$',log,re.M)
        container_run=json.loads((target/'container.run.json').read_text())
        completed=bool(declared and times and int(declared[-1])==float(times[-1])
                       and container_run['exitCode']==0 and 'FOAM FATAL' not in log
                       and (target/declared[-1]/'U').exists() and (target/declared[-1]/'p').exists())
        result=dict(status='converged' if completed else 'not-converged',
            declaredIterations=int(declared[-1]) if declared else None,lastIteration=float(times[-1]) if times else None,
            standardCheckMesh='Mesh OK' in (target/'checkMesh.log').read_text(),
            extendedCheckMesh='Mesh OK' in (target/'checkMesh-extended.log').read_text(),
            run=container_run)
        final=target/(declared[-1] if declared else '1500')
        if (final/'U').exists():
            mesh=nm.read_cm2d(Path(str(p)+'.solver.cm2d'));measure=nm.measure(mesh)
            u=foam.vector_field(final/'U');pressure=foam.scalar_field(final/'p')
            if len(u)==1:u*=len(mesh.cells)
            if len(pressure)==1:pressure*=len(mesh.cells)
            if len(u)!=len(mesh.cells) or len(pressure)!=len(mesh.cells):raise ValueError('OpenFOAM field count mismatch')
            rows=[dict(x=c[0],y=c[1],area=a,u=v[0],v=v[1],p=pr) for c,a,v,pr in zip(measure.centroids,measure.areas,u,pressure)]
            if all(math.isfinite(v) for row in rows for v in row.values()):result['exactErrors']=error_metrics(rows,b)
            with (target/'comparison.cells.csv').open('w') as f:
                writer=csv.DictWriter(f,fieldnames=['x','y','area','u','v','p']);writer.writeheader();writer.writerows(rows)
        # OpenFOAM residuals/continuity use its own normalization, retained verbatim.
        residuals=re.findall(r'^.*Solving for (Ux|Uy|p), Initial residual = ([^,]+), Final residual = ([^,]+), No Iterations (\d+)',log,re.M)
        result['lastLinearResiduals']=[dict(field=f,initial=float(i),final=float(v),iterations=int(n)) for f,i,v,n in residuals[-6:]]
        continuity=re.findall(r'time step continuity errors : (.*)',log)
        result['lastContinuityReport']=continuity[-1] if continuity else None
        record['openfoam']=result;records.append(record)
    dump(root/'results.json',records)
    for r in records:
        print(r['name'],r['geometry']['rawMinFraction'],r['geometry']['qualityIssues'],
              [(k,v['status'],v.get('iterations')) for k,v in r['native'].items()],r['openfoam']['status'],r['openfoam']['declaredIterations'])

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase',required=True,choices=['generate','native','prepare-foam','foam','collect'])
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--adapter',type=Path)
    args=parser.parse_args();args.root=args.root.resolve()
    if args.phase=='generate':
        if args.adapter is None:parser.error('generate requires --adapter')
        generate(args.root,args.adapter.resolve())
    elif args.phase=='native':native(args.root)
    elif args.phase=='prepare-foam':prepare(args.root)
    elif args.phase=='foam':run_foam(args.root)
    else:collect(args.root)
