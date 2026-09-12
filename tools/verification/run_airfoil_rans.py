#!/usr/bin/env python3
"""Run a bounded SST trial on an existing six-patch 2D airfoil mesh.

No remeshing or quality threshold changes. Close outer boundaries make this a
confined-domain numerical trial, not a free-air aerodynamic validation.
"""
import argparse
import hashlib
import json
import math
import re
import uuid
from pathlib import Path
import shutil

from openfoam_nozzle_flow import header, read_points, read_boundary
from run_nozzle_flow import run_stage

IMAGE = 'opencfd/openfoam-run:2606'

def save(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n')

def mesh_hashes(directory):
    return {n: hashlib.sha256((directory/n).read_bytes()).hexdigest()
            for n in ('points','faces','owner','neighbour','boundary')}

def configure(source, output):
    output.mkdir(parents=True, exist_ok=False)
    case=output/'case'
    (case/'constant').mkdir(parents=True)
    shutil.copytree(source/'constant/polyMesh',case/'constant/polyMesh')
    for folder in ('0','system'): (case/folder).mkdir()
    patches=read_boundary(case/'constant/polyMesh/boundary')
    assert {p['name'] for p in patches} == {'wall_0','left','right','top','bottom','frontAndBack'}
    assert next(p for p in patches if p['name']=='wall_0')['type']=='wall'
    assert next(p for p in patches if p['name']=='frontAndBack')['type']=='empty'
    points=read_points(case/'constant/polyMesh/points')
    bounds=[(min(p[i] for p in points),max(p[i] for p in points)) for i in range(3)]
    speed=15.; nu=1.5e-5; chord=1.; intensity=.01; viscosity_ratio=10.
    k=1.5*(speed*intensity)**2; omega=k/(viscosity_ratio*nu)
    length=math.sqrt(k)/(.09**.25*omega)
    def dictionary(name, body):
        location,obj=name.split('/')
        (case/name).write_text(header('dictionary',location,obj)+body+'\n')
    dictionary('constant/transportProperties',f'transportModel Newtonian;\nnu [0 2 -1 0 0 0 0] {nu};')
    dictionary('constant/turbulenceProperties','simulationType RAS;\nRAS { RASModel kOmegaSST; turbulence on; printCoeffs on; }')
    def field(name,dimensions,internal,conditions,vector=False):
        cls='volVectorField' if vector else 'volScalarField'
        text=header(cls,'0',name)+f'dimensions {dimensions};\ninternalField uniform {internal};\nboundaryField\n{{\n'
        for patch,body in conditions.items(): text+=f' {patch} {{ {body} }}\n'
        (case/'0'/name).write_text(text+' frontAndBack { type empty; }\n}\n')
    field('U','[0 1 -1 0 0 0 0]','(15 0 0)',{
        'left':'type fixedValue; value uniform (15 0 0);',
        'right':'type pressureInletOutletVelocity; value uniform (15 0 0);',
        'top':'type slip;', 'bottom':'type slip;', 'wall_0':'type noSlip;'},True)
    field('p','[0 2 -2 0 0 0 0]',0,{'left':'type zeroGradient;',
        'right':'type fixedValue; value uniform 0;', 'top':'type zeroGradient;',
        'bottom':'type zeroGradient;', 'wall_0':'type zeroGradient;'})
    for name,dim,val,wall in [('k','[0 2 -2 0 0 0 0]',k,'kqRWallFunction'),
                             ('omega','[0 0 -1 0 0 0 0]',omega,'omegaWallFunction')]:
        field(name,dim,val,{'left':f'type fixedValue; value uniform {val};',
          'right':f'type inletOutlet; inletValue uniform {val}; value uniform {val};',
          'top':'type zeroGradient;', 'bottom':'type zeroGradient;',
          'wall_0':f'type {wall}; value uniform {val};'})
    field('nut','[0 2 -1 0 0 0 0]',0,{'left':'type calculated; value uniform 0;',
        'right':'type calculated; value uniform 0;', 'top':'type calculated; value uniform 0;',
        'bottom':'type calculated; value uniform 0;',
        'wall_0':'type nutUSpaldingWallFunction; value uniform 0;'})
    dictionary('system/fvSolution','''solvers
{
 p { solver GAMG; tolerance 1e-9; relTol 0.01; smoother GaussSeidel; }
 "(U|k|omega)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-9; relTol 0.05; }
}
SIMPLE { nNonOrthogonalCorrectors 2; consistent yes; residualControl { p 1e-6; U 1e-6; k 1e-6; omega 1e-6; } }
relaxationFactors { fields { p 0.25; } equations { U 0.5; k 0.5; omega 0.5; } }''')
    config={'source':str(source.resolve()),'image':IMAGE,'cells':next(p['nFaces']//2 for p in patches if p['name']=='frontAndBack'),
        'meshSha256':mesh_hashes(case/'constant/polyMesh'),'bounds':bounds,
        'velocity':speed,'referenceChord':chord,'nu':nu,'Re':speed*chord/nu,
        'inletIntensity':intensity,'inletLengthScale':length,'inletNutOverNuEstimate':viscosity_ratio,'k':k,'omega':omega,
        'angleOfAttackDegrees':0,'model':'kOmegaSST','wallTreatment':'nutUSpaldingWallFunction + omegaWallFunction + kqRWallFunction',
        'outerBoundary':'left uniform inlet; right p=0 outlet; top/bottom impermeable slip',
        'scope':'Existing close-domain mesh; steady fully turbulent RANS trial, not free-air validation or grid independence.'}
    save(output/'configuration.json',config)
    (case/'airfoil.foam').touch()
    return case,config

def phase_settings(case, phase, end):
    upwind=phase=='startup'
    (case/'system/fvSchemes').write_text(header('dictionary','system','fvSchemes')+f'''
ddtSchemes {{ default steadyState; }}
gradSchemes {{ default cellLimited Gauss linear 1; }}
divSchemes {{ default none; div(phi,U) bounded Gauss {'upwind' if upwind else 'linearUpwind grad(U)'}; div(phi,k) bounded Gauss upwind; div(phi,omega) bounded Gauss upwind; div((nuEff*dev2(T(grad(U))))) Gauss linear; }}
laplacianSchemes {{ default Gauss linear limited 0.5; }}
interpolationSchemes {{ default linear; }}
snGradSchemes {{ default limited 0.5; }}
wallDist {{ method meshWave; }}
''')
    funcs='''
 yPlus { type yPlus; libs ("libfieldFunctionObjects.so"); executeControl timeStep; executeInterval 25; writeControl timeStep; writeInterval 25; }
 forces { type forceCoeffs; libs ("libforces.so"); patches (wall_0); rho rhoInf; rhoInf 1.225; CofR (0.25 0 0); liftDir (0 1 0); dragDir (1 0 0); pitchAxis (0 0 1); magUInf 15; lRef 1; Aref 0.02; writeControl timeStep; writeInterval 10; }
 extrema { type fieldMinMax; libs ("libfieldFunctionObjects.so"); fields (U p k omega nut); mode magnitude; writeControl timeStep; writeInterval 25; }
'''
    for patch in ('left','right','top','bottom'):
        funcs+=f' flux_{patch} {{ type surfaceFieldValue; libs ("libfieldFunctionObjects.so"); regionType patch; name {patch}; operation sum; fields (phi); writeFields false; writeControl timeStep; writeInterval 10; }}\n'
    (case/'system/controlDict').write_text(header('dictionary','system','controlDict')+f'''
application simpleFoam;
startFrom {'startTime' if upwind else 'latestTime'};
startTime 0; stopAt endTime; endTime {end}; deltaT 1;
writeControl timeStep; writeInterval 100; purgeWrite 2;
writeFormat ascii; writePrecision 12; writeCompression off;
timeFormat general; timePrecision 8; runTimeModifiable false;
functions {{ {funcs} }}
''')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--resume',action='store_true',help='Continue this configured case without changing its mesh or solver tolerances')
    parser.add_argument('--end-iteration',type=int,default=1500)
    args=parser.parse_args()
    if args.end_iteration<=200: parser.error('--end-iteration must exceed 200')
    if args.resume:
        if args.prepare_only: parser.error('--resume and --prepare-only cannot be combined')
        case=args.output/'case'
        config=json.loads((args.output/'configuration.json').read_text())
        results=json.loads((args.output/'stages.json').read_text())
        assert results['check-standard']['mesh_ok']
        assert config['meshSha256']==mesh_hashes(args.source/'constant/polyMesh')==mesh_hashes(case/'constant/polyMesh')
        iterations=[float(p.name) for p in case.iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name)]
        if not iterations or args.end_iteration<=max(iterations): parser.error('Resume end must exceed the latest saved iteration')
    else:
        case,config=configure(args.source,args.output)
        phase_settings(case,'startup',200)
        if args.prepare_only: return
        results={}
    save(args.output/'completion.json',{'solverProcessesCompleted':False,'status':'running','requestedEndIteration':args.end_iteration})
    def run(label,program,extra=(),timeout=120):
        name='cartmesh-airfoil-'+label+'-'+uuid.uuid4().hex[:8]
        cmd=['docker','run','--rm','--name',name,'--network','none','--cpus','2',
             '-v',f'{case.resolve()}:/home/openfoam/workingDir/case',IMAGE,program,
             '-case','/home/openfoam/workingDir/case',*extra]
        stage=run_stage(args.output/'stages'/label,cmd,timeout,name)
        if program=='checkMesh':
            log=(args.output/'stages'/label/'stdout.log').read_text()
            stage['mesh_ok']='Mesh OK.' in log and 'Failed ' not in log
            failed=re.search(r'Failed (\d+) mesh checks',log)
            stage['failed_mesh_checks']=int(failed.group(1)) if failed else 0
        results[label]=stage
        save(args.output/'stages.json',results)
        return stage
    if not args.resume:
        std=run('check-standard','checkMesh')
        log=(args.output/'stages/check-standard/stdout.log').read_text()
        if std['status']!='passed' or 'Mesh OK.' not in log:
            raise RuntimeError('Standard checkMesh failed; solver not started')
        run('check-expanded','checkMesh',('-allGeometry','-allTopology'))
    # Only momentum convection is higher order; k and omega retain upwind.
    phases=[(f'continuation-{args.end_iteration}',args.end_iteration)] if args.resume else [('startup',200),('second-order',args.end_iteration)]
    for name,end in phases:
        phase_settings(case,name,end)
        stage=run(name,'simpleFoam',timeout=900)
        shutil.copyfile(case/'system/fvSchemes',args.output/'stages'/name/'fvSchemes')
        shutil.copyfile(case/'system/controlDict',args.output/'stages'/name/'controlDict')
        if stage['status']!='passed': raise RuntimeError(f'{name} failed; inspect preserved log')
    assert config['meshSha256']==mesh_hashes(args.source/'constant/polyMesh')==mesh_hashes(case/'constant/polyMesh')
    for label,program,extra in [
        ('write-centres','postProcess',('-func','writeCellCentres','-latestTime')),
        ('final-yplus','simpleFoam',('-postProcess','-func','yPlus','-latestTime'))]:
        stage=run(f'{label}-{args.end_iteration}',program,extra)
        if stage['status']!='passed': raise RuntimeError(f'{label} post-processing failed')
    save(args.output/'completion.json',{'solverProcessesCompleted':True,'meshUnchanged':True,
        'convergence':'Must evaluate residuals, force history, continuity, yPlus and fields; exit status is not convergence.'})

if __name__=='__main__': main()
