#!/usr/bin/env python3
"""Prepare and measure a same-cell OpenFOAM comparison for laminar.py.

Only boundary face ordering/patch names change. Points, internal faces,
owner/neighbour incidence and cell numbering are preserved.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import re
import shutil
import sys
import laminar
sys.path.insert(0, str(laminar.ROOT/'tools/flow'))
import openfoam_data as foam


def dictionary(root, location, name, body):
    path=root/location/name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(foam.header('dictionary',location,name)+body+'\n')


def prepare(case, source, mesh_file, root, scheme, iterations):
    source=source/'constant/polyMesh'
    target=root/'constant/polyMesh'
    target.mkdir(parents=True, exist_ok=True)
    points=foam.read_points(source/'points')
    faces=foam.read_faces(source/'faces')
    owners=foam.read_labels(source/'owner')
    neighbours=foam.read_labels(source/'neighbour')
    groups={}
    for i in range(len(neighbours),len(faces)):
        face=faces[i];area=foam.face_area(face,points)
        x,y=[sum(points[j][axis] for j in face)/len(face) for axis in (0,1)]
        if abs(area[2])>0:
            kind,name,u,v,p='empty','frontAndBack',0.,0.,0.
        else:
            kind,name,u,v,p=laminar.boundary_values(case,x,y,*area[:2])
        group=groups.setdefault(name,dict(kind=kind,faces=[],velocities=[]))
        if group['kind']!=kind:raise ValueError('mixed boundary kind in '+name)
        group['faces'].append(i);group['velocities'].append((u,v,0))
    order=list(range(len(neighbours)))+[i for g in groups.values() for i in g['faces']]
    if sorted(order)!=list(range(len(faces))):raise ValueError('face permutation is not bijective')
    shutil.copyfile(source/'points',target/'points')
    shutil.copyfile(source/'neighbour',target/'neighbour')
    (target/'faces').write_text(foam.header('faceList','constant/polyMesh','faces')+
        f'{len(faces)}\n(\n'+''.join(str(len(faces[i]))+'('+' '.join(map(str,faces[i]))+')\n' for i in order)+')\n')
    (target/'owner').write_text(foam.header('labelList','constant/polyMesh','owner')+
        f'{len(owners)}\n(\n'+''.join(str(owners[i])+'\n' for i in order)+')\n')
    boundary=[];u_boundary=[];p_boundary=[];start=len(neighbours)
    for name,g in groups.items():
        kind=g['kind'];count=len(g['faces'])
        patch_type='empty' if kind=='empty' else 'wall' if 'wall' in kind else 'patch'
        boundary.append(f'{name}\n{{ type {patch_type}; nFaces {count}; startFace {start}; }}\n')
        start+=count
        if kind=='empty':ub=pb='type empty;'
        elif kind=='pressure-outlet':ub='type zeroGradient;';pb='type fixedValue; value uniform 0;'
        else:
            ub='type fixedValue; value nonuniform List<vector>\n'+str(count)+'\n(\n'+''.join('('+' '.join(f'{v:.17g}' for v in row)+')\n' for row in g['velocities'])+');'
            pb='type zeroGradient;'
        u_boundary.append(f'{name}\n{{ {ub} }}\n');p_boundary.append(f'{name}\n{{ {pb} }}\n')
    (target/'boundary').write_text(foam.header('polyBoundaryMesh','constant/polyMesh','boundary')+f'{len(groups)}\n(\n'+''.join(boundary)+')\n')
    (root/'0').mkdir(exist_ok=True)
    for name,cls,dimensions,initial,bc in [('U','volVectorField','0 1 -1 0 0 0 0','(0 0 0)',u_boundary),('p','volScalarField','0 2 -2 0 0 0 0','0',p_boundary)]:
        (root/'0'/name).write_text(foam.header(cls,'0',name)+f'dimensions [{dimensions}];\ninternalField uniform {initial};\nboundaryField\n{{\n'+''.join(bc)+'}\n')
    dictionary(root,'constant','transportProperties',f'transportModel Newtonian;\nnu [0 2 -1 0 0 0 0] {laminar.controls(case)["nu"]};')
    dictionary(root,'constant','turbulenceProperties','simulationType laminar;')
    convection='upwind' if scheme=='upwind' else 'linearUpwind grad(U)'
    dictionary(root,'system','fvSchemes',f'''ddtSchemes {{ default steadyState; }}
gradSchemes {{ default leastSquares; }}
divSchemes {{ default none; div(phi,U) bounded Gauss {convection}; div((nuEff*dev2(T(grad(U))))) Gauss linear; }}
laplacianSchemes {{ default Gauss linear corrected; }}
interpolationSchemes {{ default linear; }}
snGradSchemes {{ default corrected; }}
wallDist {{ method meshWave; }}''')
    dictionary(root,'system','fvSolution','''solvers
{
 p { solver PCG; preconditioner DIC; tolerance 1e-11; relTol 0; maxIter 3000; }
 U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0; }
}
SIMPLE { nNonOrthogonalCorrectors 3; pRefCell 0; pRefValue 0; residualControl { p 1e-8; U 1e-8; } }
relaxationFactors { fields { p 0.25; } equations { U 0.6; } }''')
    forces=''
    if case in ('dfg20','annulus'):
        patches='cylinder' if case=='dfg20' else 'rotor housing'
        forces=f'''functions {{ wallLoads {{ type forces; libs (forces); patches ({patches}); rho rhoInf; rhoInf 1; CofR (0 0 0); writeControl timeStep; writeInterval 1; log false; }} }}'''
    dictionary(root,'system','controlDict',f'''application simpleFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {iterations}; deltaT 1;
writeControl timeStep; writeInterval {iterations}; purgeWrite 1; writeFormat ascii; writePrecision 17; writeCompression off;
timeFormat general; timePrecision 10; runTimeModifiable false;
{forces}''')
    laminar.write_json(root/'comparison.json',dict(case=case,scheme=scheme,mesh=str(mesh_file.resolve()),
        meshSha256=laminar.meshio.sha256_file(mesh_file),source=str(source.resolve()),
        pointsUnchanged=True,internalFacesUnchanged=True,cellNumberingUnchanged=True,
        boundaryPermutation=order[len(neighbours):],convection=convection,
        gradient='leastSquares',note='OpenFOAM linearUpwind and native limited-linear are distinct discretizations.'))


def collect(root):
    config=json.loads((root/'comparison.json').read_text())
    times=sorted((p for p in root.iterdir() if p.is_dir() and re.fullmatch(r'[0-9]+(?:\.[0-9]+)?',p.name)),key=lambda p:float(p.name))
    final=times[-1];mesh=laminar.meshio.read_cm2d(Path(config['mesh']));measured=laminar.meshio.measure(mesh)
    u=foam.vector_field(final/'U');p=foam.scalar_field(final/'p')
    if len(u)==1:u*=len(mesh.cells)
    if len(p)==1:p*=len(mesh.cells)
    if len(u)!=len(mesh.cells) or len(p)!=len(mesh.cells):raise ValueError('field/cell mismatch')
    prefix=root/'comparison'
    with Path(str(prefix)+'.cells.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['cell','x','y','area','u','v','p'])
        for i,(centre,area) in enumerate(zip(measured.centroids,measured.areas)):
            writer.writerow([i,*centre,area,*u[i][:2],p[i]])
    log=(root/'solve.log').read_text()
    converged='SIMPLE solution converged' in log
    summary=dict(namedWallLoads=[])
    if config['case']=='dfg20':
        force_files=list((root/'postProcessing/wallLoads').glob('*/force.dat'))
        if not force_files:raise ValueError('missing OpenFOAM cylinder forces')
        last=[line for line in force_files[-1].read_text().splitlines() if line and not line.startswith('#')][-1]
        numbers=[float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',last)]
        # force.dat: time, total vector, pressure vector, viscous vector.
        depth=max(v[2] for v in foam.read_points(root/'constant/polyMesh/points'))-min(v[2] for v in foam.read_points(root/'constant/polyMesh/points'))
        summary['namedWallLoads']=[dict(name='cylinder',forceX=numbers[1]/depth,forceY=numbers[2]/depth)]
    result=dict(case=config['case'],scheme=config['scheme'],iterations=float(final.name),converged=converged,
        metrics=laminar.evaluate(config['case'],prefix,summary),meshSha256=config['meshSha256'])
    laminar.write_json(root/'result.json',result)
    print({**result,'metrics':{k:v for k,v in result['metrics'].items() if isinstance(v,(float,int))}})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--collect',action='store_true')
    parser.add_argument('--source',type=Path)
    parser.add_argument('--mesh',type=Path)
    parser.add_argument('--case',choices=laminar.CASES)
    parser.add_argument('--scheme',choices=['upwind','limited-linear'],default='limited-linear')
    parser.add_argument('--iterations',type=int,default=16000)
    args=parser.parse_args()
    if args.collect:collect(args.output)
    else:
        if args.source is None or args.mesh is None or args.case is None:parser.error('prepare requires source, mesh and case')
        prepare(args.case,args.source,args.mesh,args.output,args.scheme,args.iterations)
