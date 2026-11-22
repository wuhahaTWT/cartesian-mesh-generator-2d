#!/usr/bin/env python3
"""Reproducible physical benchmarks for the native incompressible solver.

No acceptance decision is inferred from a process exit code alone. Every run
retains its mesh, fields, residual history, controls and elapsed time. Results
from a nonconverged run remain diagnostic and are excluded from observed orders.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import shutil
from pathlib import Path
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "flow"))
import native_mesh as meshio

CASES = ("cavity100", "cavity400", "cavity1000", "dfg20", "annulus", "poiseuille")
# Ghia, Ghia & Shin (1982), Tables I and II, DOI 10.1016/0021-9991(82)90058-4.
# Coordinates are the printed table coordinates; no silent point substitution.
GHIA_Y = [1, .9766, .9688, .9609, .9531, .8516, .7344, .6172, .5, .4531,
          .2813, .1719, .1016, .0703, .0625, .0547, 0]
GHIA_X = [1, .9688, .9609, .9531, .9453, .9063, .8594, .8047, .5, .2344,
          .2266, .1563, .0938, .0781, .0703, .0625, 0]
GHIA_U = {
    100: [1,.84123,.78871,.73722,.68717,.23151,.00332,-.13641,-.20581,-.21090,
          -.15662,-.10150,-.06434,-.04775,-.04192,-.03717,0],
    400: [1,.75837,.68439,.61756,.55892,.29093,.16256,.02135,-.11477,-.17119,
          -.32726,-.24299,-.14612,-.10338,-.09266,-.08186,0],
    1000: [1,.65928,.57492,.51117,.46604,.33304,.18719,.05702,-.06080,-.10648,
           -.27805,-.38289,-.29730,-.22220,-.20196,-.18109,0],
}
GHIA_V = {
    100: [0,-.05906,-.07391,-.08864,-.10313,-.16914,-.22445,-.24533,.05454,
          .17527,.17507,.16077,.12317,.10890,.10091,.09233,0],
    400: [0,-.12146,-.15663,-.19254,-.22847,-.23827,-.44993,-.38598,.05188,
          .30174,.30203,.28124,.22965,.20920,.19713,.18360,0],
    1000: [0,-.21388,-.27669,-.33714,-.39188,-.51550,-.42665,-.31966,.02526,
           .32235,.33075,.37095,.32627,.30353,.29012,.27485,0],
}
DFG = dict(Cd=5.57953523384, Cl=.010618948146, pressureDrop=.11752016697)
DFG_SOURCE = "https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html"


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def run(command, log, timeout):
    start = time.perf_counter()
    with log.open("w") as stream:
        try:
            p = subprocess.run([str(x) for x in command], stdout=stream,
                               stderr=subprocess.STDOUT, timeout=timeout)
            status = dict(returnCode=p.returncode, timedOut=False)
        except subprocess.TimeoutExpired:
            status = dict(returnCode=None, timedOut=True)
    return dict(command=list(map(str, command)), seconds=time.perf_counter()-start, **status)


def circle(n, radius, centre=(0., 0.)):
    result = []
    for i in range(n):
        # Analytic cardinal vertices avoid trigonometric roundoff at axes.
        if (4*i) % n == 0:
            x, y = ((1,0),(0,1),(-1,0),(0,-1))[(4*i)//n]
        else:
            angle = 2*math.pi*i/n
            x, y = math.cos(angle), math.sin(angle)
        result.append((centre[0]+radius*x, centre[1]+radius*y))
    return result


def geometry(case, level, segments=None):
    segments = segments if segments is not None else 4*2**level
    if case.startswith("cavity"):
        return [[(0,0),(1,0),(1,1),(0,1)]]
    if case == "poiseuille":
        return [[(0,0),(2,0),(2,1),(0,1)]]
    if case == "dfg20":
        return [[(0,0),(2.2,0),(2.2,.41),(0,.41)],
                circle(segments, .05, (.2,.2))]
    return [circle(segments, 1), circle(segments, .5)]


def controls(case):
    if case.startswith("cavity"):
        return dict(nu=1/int(case[6:]), speed=1.)
    if case == "dfg20":
        return dict(nu=.001, speed=.3)
    return dict(nu=.1, speed=.5 if case == "annulus" else 1.)


def boundary_values(case, x, y, sx, sy):
    if case == "annulus":
        if math.hypot(x,y) < .75:
            length=math.hypot(sx,sy)
            return "smooth-moving-wall", "rotor", .5*sy/length, -.5*sx/length, 0.
        return "wall", "housing", 0., 0., 0.
    if case.startswith("cavity"):
        return ("moving-wall", "lid", 1., 0., 0.) if abs(y-1)<1e-10 else ("wall", "walls", 0.,0.,0.)
    width, height = (2.2,.41) if case == "dfg20" else (2.,1.)
    speed=controls(case)["speed"]
    if abs(x) < 1e-10:
        return "velocity-inlet", "inlet", 4*speed*y*(height-y)/height**2, 0., 0.
    if abs(x-width) < 1e-10:
        return "pressure-outlet", "outlet", 0., 0., 0.
    name="cylinder" if case == "dfg20" and math.hypot(x-.2,y-.2)<.051 else "walls"
    return "wall", name, 0.,0.,0.


def prepare(case, level, root, mesh_cli, timeout, layout='native', segments=None):
    directory=root / "meshes" / f"{case}-l{level}"
    directory.mkdir(parents=True, exist_ok=True)
    manifest=directory / "mesh.json"
    if manifest.exists():
        record=json.loads(manifest.read_text())
        if record.get('layout','native') != layout or record.get('geometrySegments') != segments:
            raise ValueError('mesh layout or geometry resolution changed; use a fresh output directory')
        return record
    xy=directory / "boundary.xy"
    xy.write_text("\n\n".join("\n".join(f"{x:.17g} {y:.17g}" for x,y in loop)
                              for loop in geometry(case,level,segments))+"\n")
    prefix=directory / "mesh"
    # One background cell outside each side gives nested interior cell sizes.
    padding=1/(2**level-2)
    command=[mesh_cli,xy,prefix,level,padding,.1,"interior",directory/"openfoam",level,0]
    if layout=='square':
        command += ['--size-field','--far-field-spans',padding,
                    '--wall-cells-per-span',2**level-2,'--cells-per-level',0,'--far-level',level]
    generation=run(command,directory/"generation.log",timeout)
    record=dict(case=case,level=level,layout=layout,geometrySegments=segments,generation=generation,path=str(prefix)+".solver.cm2d",
                boundary=str(directory/"flow.boundaries"),geometry=str(xy),openfoam=str(directory/"openfoam"))
    record['geometrySha256']=meshio.sha256_file(xy)
    if case in ('annulus','dfg20'):
        facets=segments if segments is not None else 4*2**level
        radii=[1.,.5] if case=='annulus' else [.05]
        record['circleApproximation']=[dict(radius=r,facets=facets,
            maximumSagitta=r*(1-math.cos(math.pi/facets)),
            circleAreaMinusPolygonArea=r*r*(math.pi-.5*facets*math.sin(2*math.pi/facets))) for r in radii]
    if generation["returnCode"] == 0:
        mesh=meshio.read_cm2d(Path(record["path"]))
        measured=meshio.measure(mesh)
        faces=meshio.face_geometry(mesh,measured)
        boundaries=[(edge,face) for edge,face in zip(mesh.edges,faces) if edge.neighbour<0]
        with Path(record["boundary"]).open("w") as stream:
            stream.write(f"CARTMESH2D_FLOW_BOUNDARIES 1\nCOUNTS {len(mesh.cells)} {len(mesh.edges)} {len(boundaries)}\n")
            for edge,face in boundaries:
                x,y=face.centre;sx,sy=face.area_vector
                kind,name,u,v,p=boundary_values(case,x,y,sx,sy)
                stream.write(f'BOUNDARY {edge.id} {edge.owner} {x:.17g} {y:.17g} {sx:.17g} {sy:.17g} {kind} "{name}" {u:.17g} {v:.17g} {p:.17g}\n')
            stream.write("END\n")
        record.update(cells=len(mesh.cells),h=measured.characteristic_h,fluidArea=measured.total_area,
                      meshSha256=meshio.sha256_file(Path(record["path"])))
    write_json(manifest,record)
    return record


def rms(values):
    return math.sqrt(math.fsum(x*x for x in values)/len(values))


def evaluate(case, prefix, summary):
    rows=meshio.read_cells(str(prefix)+".cells.csv")
    result={}
    if case.startswith("cavity"):
        reynolds=int(case[6:]);samples=[]
        for field,coordinates,reference in [("u",GHIA_Y,GHIA_U[reynolds]),("v",GHIA_X,GHIA_V[reynolds])]:
            for coordinate,expected in zip(coordinates,reference):
                x,y=(.5,coordinate) if field=="u" else (coordinate,.5)
                walls=[dict(x=.5,y=0,u=0),dict(x=.5,y=1,u=1)] if field=="u" else [dict(x=0,y=.5,v=0),dict(x=1,y=.5,v=0)]
                actual=meshio.affine_sample(rows,x,y,field,boundary=walls)
                samples.append(dict(field=field,coordinate=coordinate,actual=actual,reference=expected,error=actual-expected))
        result.update(ghiaRmse=rms([r['error'] for r in samples]),ghiaMax=max(abs(r['error']) for r in samples),samples=samples)
        if reynolds==400:
            result['referenceNotes']=[dict(field='v',coordinate=.9063,reference=-.23827,
                note='Suspicious printed Ghia Table II entry retained in every metric; no substitution or exclusion.')]
        if reynolds==100:
            ref=json.loads((ROOT/'tools/verification/references/cavity-marchi-2009-re100.json').read_text())
            errors=[]
            for field in ('u','v'):
                for r in ref[field]:
                    x,y=(.5,r['coordinate']) if field=='u' else (r['coordinate'],.5)
                    errors.append(meshio.affine_sample(rows,x,y,field)-r['value'])
            result['marchiRmse']=rms(errors)
    elif case=="dfg20":
        load=next(x for x in summary['namedWallLoads'] if x['name']=='cylinder')
        result.update(Cd=load['forceX']/.002,Cl=load['forceY']/.002)
        result['pressureDrop']=meshio.affine_sample(rows,.15,.2,'p',count=16)-meshio.affine_sample(rows,.25,.2,'p',count=16)
        for k in DFG:
            result[k+'Error']=result[k]-DFG[k]
            result[k+'RelativeError']=abs(result[k]-DFG[k])/abs(DFG[k])
        result['reference']=dict(**DFG,source=DFG_SOURCE,coefficientMeanVelocity=.2)
    else:
        errors_u=[];errors_p=[];weights=[]
        if case=='annulus':
            A=-1/3;B=1/3
            pressure=lambda r: .5*A*A*r*r+2*A*B*math.log(r)-.5*B*B/(r*r)
            first=rows[0];gauge=pressure(math.hypot(first['x'],first['y']))
        for r in rows:
            x,y=r['x'],r['y'];weight=r.get('area',1.)
            if case=='annulus':
                radius=math.hypot(x,y);tangent=A*radius+B/radius
                u,v=-y*tangent/radius,x*tangent/radius
                p=pressure(radius)-gauge
            else:
                u,v=4*y*(1-y),0.;p=.8*(2-x)
            errors_u.append((r['u']-u)**2+(r['v']-v)**2)
            errors_p.append((r['p']-p)**2);weights.append(weight)
        area=math.fsum(weights)
        result['velocityL2']=math.sqrt(math.fsum(w*e for w,e in zip(weights,errors_u))/area)/controls(case)['speed']
        result['pressureL2']=math.sqrt(math.fsum(w*e for w,e in zip(weights,errors_p))/area)/controls(case)['speed']**2
        if case=='annulus':
            exact=4*math.pi*.1/3
            for load in summary['namedWallLoads']:
                if load['name'] in ('rotor','housing'):
                    result[load['name']+'Torque']=load['torque']
                    result[load['name']+'TorqueError']=abs(abs(load['torque'])-exact)/exact
    return result


def solve(case, level, scheme, args):
    directory=args.output / f"{case}-l{level}-{scheme}"
    directory.mkdir(parents=True,exist_ok=True)
    record_file=directory/'result.json'
    if record_file.exists() and not args.rerun:
        return json.loads(record_file.read_text())
    mesh=json.loads((args.output/'meshes'/f'{case}-l{level}'/'mesh.json').read_text())
    record=dict(case=case,level=level,scheme=scheme,mesh=mesh,status='mesh-failed',converged=False)
    if mesh['generation']['returnCode'] == 0:
        prefix=directory/'flow';physical=controls(case)
        command=[args.flow_cli,'--mesh',mesh['path'],'--case','custom','--boundary',mesh['boundary'],
                 '--output',prefix,'--nu',physical['nu'],'--speed',physical['speed'],
                 '--convection',scheme,'--max-iterations',args.iterations,'--tolerance',args.tolerance,
                 '--velocity-relaxation',args.relaxation,'--pressure-corrections',args.corrections,
                 '--pressure-preconditioner',args.preconditioner,'--profile']
        if args.adaptive:command+=['--linear-policy','adaptive']
        if args.anderson:command+=['--steady-acceleration','anderson']
        if args.extra:command+=args.extra
        record['run']=run(command,directory/'solve.log',args.timeout)
        summary=Path(str(prefix)+'.json')
        if summary.exists():
            data=json.loads(summary.read_text());record['summary']=data
            record['converged']=data.get('converged',False) and record['run']['returnCode']==0
            record['status']='converged' if record['converged'] else 'not-converged'
            record['metrics']=evaluate(case,prefix,data)
        else:record['status']='solve-failed'
        record['binarySha256']=args.binary_hash
    write_json(record_file,record)
    metrics={k:v for k,v in record.get('metrics',{}).items() if isinstance(v,(float,int))}
    print(f"{case} l{level} {scheme}: {record['status']} {metrics}",flush=True)
    return record


def report(records, output):
    rows=[]
    for r in records:
        row={k:r[k] for k in ['case','level','scheme','status','converged']}
        row.update(cells=r['mesh'].get('cells'),h=r['mesh'].get('h'),seconds=r.get('run',{}).get('seconds'),
                   iterations=r.get('summary',{}).get('iterations'))
        row.update({k:v for k,v in r.get('metrics',{}).items() if isinstance(v,(float,int))})
        rows.append(row)
    for row in rows:
        prior=next((r for r in rows if r['case']==row['case'] and r['scheme']==row['scheme'] and r['level']==row['level']-1),None)
        if prior and prior['converged'] and row['converged']:
            for key in ['ghiaRmse','marchiRmse','velocityL2','pressureL2','CdRelativeError','rotorTorqueError']:
                if row.get(key,0)>0 and prior.get(key,0)>0:
                    row[key+'Order']=math.log(prior[key]/row[key])/math.log(prior['h']/row['h'])
    write_json(output/'summary.json',rows)
    columns=list(dict.fromkeys(key for row in rows for key in row))
    with (output/'summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,lineterminator='\n');writer.writeheader();writer.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(13,8),constrained_layout=True)
    for case,ax in zip(CASES,axes.flat):
        for scheme in sorted({r['scheme'] for r in rows}):
            values=sorted([r for r in rows if r['case']==case and r['scheme']==scheme and r['converged']],key=lambda r:r['h'])
            key='ghiaRmse' if case.startswith('cavity') else 'CdRelativeError' if case=='dfg20' else 'pressureL2'
            if values:ax.loglog([r['h'] for r in values],[r[key] for r in values],'o-',label=scheme)
        ax.set_title(case);ax.set_xlabel('sqrt(area / cells)');ax.set_ylabel('error');ax.grid(True,which='both',alpha=.25)
        if ax.lines:ax.legend(fontsize=8)
    fig.savefig(output/'convergence.png',dpi=180);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cases',nargs='+',choices=CASES,default=CASES)
    parser.add_argument('--levels',nargs='+',type=int,default=[5,6,7])
    parser.add_argument('--schemes',nargs='+',default=['upwind','limited-linear'])
    parser.add_argument('--mesh-cli',type=Path,default=ROOT/'build/cartmesh2d_cli')
    parser.add_argument('--mesh-layout',choices=['native','square'],default='native')
    parser.add_argument('--geometry-segments',type=int,help='fixed polygon facet count for each circle across all mesh levels; default is 4*2**level')
    parser.add_argument('--preconditioner',choices=['jacobi','ic0','aggregation','cholesky'],
                        default='cholesky' if platform.system()=='Darwin' else 'ic0')
    parser.add_argument('--flow-cli',type=Path,default=ROOT/'build/cartmesh2d_flow_cli')
    parser.add_argument('--iterations',type=int,default=16000)
    parser.add_argument('--tolerance',type=float,default=1e-8)
    parser.add_argument('--relaxation',type=float,default=.6)
    parser.add_argument('--corrections',type=int,default=4)
    parser.add_argument('--timeout',type=float,default=1800)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--adaptive',action='store_true')
    parser.add_argument('--anderson',action='store_true')
    parser.add_argument('--rerun',action='store_true')
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--extra',nargs=argparse.REMAINDER,default=[])
    args=parser.parse_args()
    if args.geometry_segments is not None and (args.geometry_segments<8 or args.geometry_segments%4):
        parser.error('--geometry-segments must be a multiple of four, at least eight')
    args.output=args.output.resolve();args.flow_cli=args.flow_cli.resolve();args.mesh_cli=args.mesh_cli.resolve()
    args.output.mkdir(parents=True,exist_ok=True)
    args.binary_hash=meshio.sha256_file(args.flow_cli)
    # A queued run must not pick up a later rebuild of build/ midway through
    # the study. Keep one executable per content hash in the study directory.
    binary=args.output/'binaries'/('flow-'+args.binary_hash+args.flow_cli.suffix)
    binary.parent.mkdir(exist_ok=True)
    if not binary.exists():shutil.copy2(args.flow_cli,binary)
    if meshio.sha256_file(binary)!=args.binary_hash:raise ValueError('benchmark executable snapshot hash mismatch')
    args.flow_cli=binary
    for case in args.cases:
        for level in args.levels:prepare(case,level,args.output,args.mesh_cli,args.timeout,args.mesh_layout,args.geometry_segments)
    if args.prepare_only:return
    jobs=[(case,level,scheme) for case in args.cases for level in args.levels for scheme in args.schemes]
    records=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for record in pool.map(lambda item:solve(*item,args),jobs):
            records.append(record);write_json(args.output/'results.json',records)
    report(records,args.output)


if __name__=='__main__':main()
