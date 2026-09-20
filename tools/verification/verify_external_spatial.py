#!/usr/bin/env python3
"""Re20 external-flow spatial sensitivity on one fixed polygonal circle.

Three resolutions refine the wall and background together, keeping the domain,
32-segment geometry, physical controls and boundary conditions fixed. Differences
in integrated traction are observations, not errors against a continuum circle
or a qualification of spatial accuracy. The finest grid is also solved with a
100-times tighter stopping tolerance to expose remaining iteration sensitivity.
"""
import argparse
import json
import math
import shutil
import time
from pathlib import Path
import verify_native_flow as native


def study(args):
    root=args.output.resolve()
    if root.exists():raise ValueError('Use a new output directory; preserve earlier failures')
    mesh_cli=args.mesh_cli.resolve(strict=True);flow_cli=args.flow_cli.resolve(strict=True)
    geometry=args.geometry.resolve(strict=True);root.mkdir(parents=True)
    report=dict(format='cartmesh2d-external-spatial-sensitivity-v1',valid=False,
        accuracyQualification='not-qualified',scope=__doc__.strip(),
        sourceSha256={str(p):native.sha256_file(p) for p in (mesh_cli,flow_cli,geometry,Path(__file__))},
        controls=dict(nu=.1,speed=1,referenceLength=2,farFieldSpans=10,cellsPerLevel=3,smallAlpha=.1,
            convection='face-limited-linear',pressurePreconditioner='aggregation',tolerance=1e-8,
            tightTolerance=1e-10,maximumIterations=4000),
        budget=dict(perCommandSeconds=180,totalSeconds=1000,retainedBytes=60*1024**2,minimumFreeBytes=5*1024**3),
        runs=[],cases=[],issues=[])
    began=time.monotonic()
    retained_geometry=root/'input.xy';shutil.copyfile(geometry,retained_geometry)
    def save():native.write_json(root/'study.json',report)
    def run(command,label):
        remaining=1000-(time.monotonic()-began)
        if remaining<=0:raise ValueError('Total 1000 second study budget exhausted')
        if shutil.disk_usage(root).free<5*1024**3:raise ValueError('Under 5 GiB free disk')
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file())>50*1024**2:
            raise ValueError('Insufficient space remaining in the 60 MiB study budget')
        start=time.monotonic();result=native.run([str(x) for x in command],root/'logs'/label,min(180,remaining))
        result['elapsedSeconds']=time.monotonic()-start;report['runs'].append(result);save()
        print(label,result['status'],round(result['elapsedSeconds'],3),flush=True)
        if result['status']!='passed':raise ValueError(label+' failed; original inputs and outputs retained')
    def solve(mesh,prefix,tolerance,label):
        run([flow_cli,'--mesh',mesh,'--output',prefix,'--case','external','--nu',.1,'--speed',1,
            '--convection','face-limited-linear','--pressure-preconditioner','aggregation',
            '--tolerance',tolerance,'--max-iterations',4000],label)
        audit=native.verify_case(mesh,prefix,'external',.1,1,native.argument_parser().parse_args(['--max-iterations','4000']))
        native.write_json(Path(str(prefix)+'.audit.json'),audit)
        if not audit['valid']:raise ValueError(label+' failed independent audit: '+str(audit['issues']))
        summary=json.loads(Path(str(prefix)+'.json').read_text())
        if summary['tolerance']!=tolerance or summary['convection']!='face-limited-linear':raise ValueError('Numerical controls changed')
        return summary,audit
    save()
    try:
        common=None
        for index,wall in enumerate((.125,.0625,.03125)):
            folder=root/f'level{index}';folder.mkdir();prefix=folder/'mesh';mesh=folder/'mesh.solver.cm2d'
            run([mesh_cli,retained_geometry,prefix,8,.25,.1,'exterior',folder/'foam',0,0,
                '--size-field','--reference-length',2,'--wall-relative-size',wall,
                '--background-relative-size',wall*8,'--far-field-spans',10,'--cells-per-level',3],f'level{index}-mesh')
            parsed=native.read_cm2d(mesh);measured=native.measure(parsed,1e-11,1e-9)
            reference=native.circular_obstacle_reference(parsed)
            if not reference['circleReferenceApplicable'] or reference['segments']!=32 or abs(reference['radius']-1)>1e-8:
                raise ValueError('Re20 controls require the same unit-radius 32-segment circle')
            domain=[min(v[0] for v in parsed.vertices),min(v[1] for v in parsed.vertices),
                max(v[0] for v in parsed.vertices),max(v[1] for v in parsed.vertices)]
            identity=(domain,reference['segments'])
            if common is None:common=identity
            if identity!=common:raise ValueError('Domain or polygon segment count changed')
            quality=json.loads((folder/'foam/solver_quality.json').read_text())
            if quality['valid'] is not True:raise ValueError('Solver-quality gate did not pass')
            summary,audit=solve(mesh,folder/'flow',1e-8,f'level{index}-flow')
            report['cases'].append(dict(label=folder.name,cells=len(parsed.cells),wallRelativeSize=wall,
                actualCharacteristicH=measured.characteristic_h,area=math.fsum(measured.areas),domain=domain,
                geometryReference=reference,solverQuality=quality,
                resolution=json.loads(Path(str(prefix)+'.resolution.json').read_text()),
                summary=summary,meshSha256=native.sha256_file(mesh),auditSha256=native.sha256_file(folder/'flow.audit.json')))
            save()
        fine=root/'level2';mesh=fine/'mesh.solver.cm2d';tight,tight_audit=solve(mesh,fine/'tight',1e-10,'level2-tight')
        parsed=native.read_cm2d(mesh);measured=native.measure(parsed,1e-11,1e-9)
        a=native.read_cells(fine/'flow.cells.csv',parsed,measured,'external');b=native.read_cells(fine/'tight.cells.csv',parsed,measured,'external')
        report['iterationSensitivity']=dict(maxVelocityChange=max(math.hypot(x['u']-y['u'],x['v']-y['v']) for x,y in zip(a,b)),
            maxPressureChange=max(abs(x['p']-y['p']) for x,y in zip(a,b)),
            forceXChange=abs(tight['forceX']-report['cases'][-1]['summary']['forceX']),
            forceYChange=abs(tight['forceY']-report['cases'][-1]['summary']['forceY']),
            tightSummary=tight,tightAuditSha256=native.sha256_file(fine/'tight.audit.json'))
        report['successiveTractionDifferences']=[dict(forceX=abs(a['summary']['forceX']-b['summary']['forceX']),
            forceY=abs(a['summary']['forceY']-b['summary']['forceY'])) for a,b in zip(report['cases'],report['cases'][1:])]
        report['artifactSha256']={str(p.relative_to(root)):native.sha256_file(p) for p in root.rglob('*') if p.is_file() and p.name!='study.json'}
        report['valid']=True;save()
    except Exception as error:
        report['issues'].append(str(error));save();raise
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--geometry',type=Path,default=Path('examples/acceptance/circle.xy'))
    parser.add_argument('--mesh-cli',type=Path,default=Path('build/cartmesh2d_cli'))
    parser.add_argument('--flow-cli',type=Path,default=Path('outputs/native-flow/binaries/flow-0.4.20'))
    result=study(parser.parse_args());print(json.dumps(dict(valid=result['valid'],cases=[c['cells'] for c in result['cases']],accuracyQualification=result['accuracyQualification'])))
