"""Same-cell 2D Fluent interchange, setup and real solver comparisons.

Mesh format follows Fluent's documented face-based ASCII format. No physical
wall is moved; one wall zone per original facet supports its exact constant
tangential velocity. Field data and settings are retained, including failures.
"""
from pathlib import Path
import argparse, hashlib, json, math, re, sys, time
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tools/flow'))
import native_mesh as nm

def export_mesh(meshfile,destination):
    m=nm.read_cm2d(meshfile);g=nm.measure(m);faces=nm.face_geometry(m,g)
    groups={3:[]};velocities={};names={3:'interior'}
    # Deduplicate source wall direction at 1e-10 radians for zone assignment
    # only; coordinates, short edges and per-face connectivity are unchanged.
    wall_keys={}
    for e,f in zip(m.edges,faces):
        if e.neighbour>=0:zone=3
        else:
            x,y=f.centre;inner=math.hypot(x,y)<.75
            sx,sy=f.area_vector;length=math.hypot(sx,sy)
            u,v=(.5*sy/length,-.5*sx/length) if inner else (0.,0.)
            # Inner outward normal points into solid; (sy,-sx) gives CCW.
            if inner:
                key=(round(u,10),round(v,10))
                if key not in wall_keys:wall_keys[key]=4+len(wall_keys)
                zone=wall_keys[key];velocities[zone]=(u,v);names[zone]=f'rotor_{zone}'
            else:zone=1000;velocities[zone]=(0.,0.);names[zone]='housing'
        # In Fluent's 2D convention the directed segment keeps c0 on its
        # left. Preserve the native owner-outward (dy,-dx) orientation.
        # The manual's 3D face-loop right-hand wording must not be applied
        # as a reversal of this 2D segment convention.
        a,b=e.v0,e.v1;pa,pb=m.vertices[a],m.vertices[b]
        if (pb[1]-pa[1])*f.area_vector[0]-(pb[0]-pa[0])*f.area_vector[1]<0:a,b=b,a
        groups.setdefault(zone,[]).append((a+1,b+1,e.owner+1,e.neighbour+1))
    H=lambda i:format(i,'x')
    with destination.open('x') as out:
        out.write('(0 "CartMesh2D same-cell 2D reference")\n(2 2)\n')
        out.write(f'(10 (0 1 {H(len(m.vertices))} 0))\n(10 (1 1 {H(len(m.vertices))} 1 2)(\n')
        out.writelines(f'{x:.17g} {y:.17g}\n' for x,y in m.vertices);out.write('))\n')
        out.write(f'(12 (0 1 {H(len(m.cells))} 0))\n(12 (2 1 {H(len(m.cells))} 1 0)(\n')
        out.write(' '.join('1' if len(c.vertices)==3 else '3' if len(c.vertices)==4 else '7' for c in m.cells)+'\n))\n')
        out.write(f'(13 (0 1 {H(len(m.edges))} 0))\n');start=1
        for zone,rows in sorted(groups.items()):
            out.write(f'(13 ({H(zone)} {H(start)} {H(start+len(rows)-1)} {2 if zone==3 else 3} 2)(\n')
            out.writelines(' '.join(map(H,row))+'\n' for row in rows);out.write('))\n');start+=len(rows)
        out.write('(45 (2 fluid fluid)())\n')
        for zone in sorted(groups):out.write(f'(45 ({zone} {"interior" if zone==3 else "wall"} {names[zone]})())\n')
    data=dict(cells=len(m.cells),faces=len(m.edges),vertices=len(m.vertices),area=g.total_area,
        meshSha256=hashlib.sha256(meshfile.read_bytes()).hexdigest(),wallVelocities={names[z]:v for z,v in velocities.items()},
        boundaryModel='Original facet-wall motion, rho=1, nu=0.1',
        source='https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/tgd_user_format_gridsect.html')
    destination.with_suffix('.interchange.json').write_text(json.dumps(data,indent=2)+'\n')
    return data

def execute(args):
    import ansys.fluent.core as pyfluent
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);record=dict(status='starting',threads=args.cores,gpuUsed=False,
        dimension=2,precision='double',rho=1,nu=.1,tolerance=args.tolerance,mesh=str(args.mesh),
        coupling=args.coupling,
        convergenceScope='Fluent scaled residuals; different from native strict combined stopping gate')
    started=time.monotonic();session=None
    def save():
        record['totalSeconds']=time.monotonic()-started
        (out/'result.json').write_text(json.dumps(record,indent=2,default=str)+'\n')
    try:
        data=export_mesh(args.mesh,out/'input.msh');record['interchange']=data;save()
        session=pyfluent.launch_fluent(fluent_path=r'E:\ANSYS Inc\v261\fluent\ntbin\win64\fluent.exe',mode='solver',
            dimension=2,precision='double',processor_count=args.cores,ui_mode='no_gui',start_timeout=120,
            cwd=str(out),cleanup_on_exit=True,start_watchdog=True,start_container=False,
            start_transcript=False,env={'REMOTING_SERVER_ADDRESS':'127.0.0.1'})
        session.transcript.start(file_name=str(out/'fluent-transcript.log'))
        record['version']=str(session.get_fluent_version());record['launchSeconds']=time.monotonic()-started;save()
        s=session.settings;s.file.read_mesh(file_name=str(out/'input.msh'));session.tui.mesh.check()
        record['meshRead']=True
        transcript=(out/'fluent-transcript.log').read_text(encoding='utf-8',errors='replace')
        minimums=re.findall(r'minimum volume\s*\(m3\):\s*([-+0-9.eE]+)',transcript)
        record['meshCheckPassed']=bool(minimums and math.isfinite(float(minimums[-1])) and float(minimums[-1])>0
            and 'Checking mesh' in transcript and 'Done.' in transcript and 'Mesh check failed' not in transcript
            and 'left-handed faces detected' not in transcript)
        if not record['meshCheckPassed']:raise RuntimeError('Fluent mesh check failed; no repair or solving is allowed')
        record['initialState']={'walls':s.setup.boundary_conditions.wall.get_object_names(),
            'material':s.setup.materials.fluid.get_state(),'methods':s.solution.methods.get_state(),
            'wallExample':s.setup.boundary_conditions.wall['housing'].get_state(),
            'exportSurfaces':s.file.export.ascii.surface_name_list.get_attr('allowed-values')}
        save()
        if args.probe:
            record['status']='mesh-probe-passed';save();return
        s.setup.models.viscous.model='laminar'
        material=s.setup.materials.fluid['air'];material.density.value=1.;material.viscosity.value=.1
        for name,velocity in data['wallVelocities'].items():
            w=s.setup.boundary_conditions.wall[name].momentum
            if name=='housing':continue
            w.wall_motion='Moving Wall';w.velocity_spec='Components';w.relative=False
            w.velocity_components[0].value=velocity[0];w.velocity_components[1].value=velocity[1]
        s.solution.methods.p_v_coupling.flow_scheme='SIMPLE' if args.coupling=='simple' else 'Coupled'
        s.solution.methods.discretization_scheme={'pressure':'standard','mom':'first-order-upwind'}
        if args.coupling=='simple':s.solution.controls.under_relaxation={'pressure':.25,'mom':.6}
        for name in ('continuity','x-velocity','y-velocity'):
            s.solution.monitor.residual.equations[name].absolute_criteria=args.tolerance
        s.solution.initialization.initialization_type='standard'
        s.solution.initialization.standard_initialize()
        record['configuredState']={'walls':s.setup.boundary_conditions.wall.get_state(),
            'material':material.get_state(),'methods':s.solution.methods.get_state(),
            'controls':s.solution.controls.get_state(),'residuals':s.solution.monitor.residual.get_state()};save()
        t=time.monotonic();s.solution.run_calculation.iterate(iter_count=args.iterations)
        record['solveSeconds']=time.monotonic()-t
        text=(out/'fluent-transcript.log').read_text(encoding='utf-8',errors='replace')
        history=[]
        for row in re.findall(r'^\s*(\d+)\s+([0-9.eE+\-]+)\s+([0-9.eE+\-]+)\s+([0-9.eE+\-]+)',text,re.M):
            iteration,c,u,v=row
            history.append(dict(iteration=int(iteration),continuity=float(c),xVelocity=float(u),yVelocity=float(v)))
        (out/'residuals.json').write_text(json.dumps(history,indent=2)+'\n')
        record.update(iterations=history[-1]['iteration'] if history else None,
            residualConverged=bool(history and all(history[-1][k]<=args.tolerance for k in ['continuity','xVelocity','yVelocity'])),
            finalResiduals=history[-1] if history else None,cells=data['cells']);save()
        s.file.write_case_data(file_name=str(out/'accepted.cas.h5'))
        record['caseDataWritten']=True;save()
        s.file.export.ascii(file_name=str(out/'cells.csv'),surface_name_list=[],delimiter='comma',
            quantities=['x-velocity','y-velocity','pressure'],location='cell-center')
        record.update(status='finished',caseDataWritten=True)
    except Exception as error:
        record.update(status='failed',error=repr(error));save();raise
    finally:
        if session is not None:
            try:session.exit(timeout=15,wait=15);record['sessionClosed']=True
            except Exception as error:record['closeError']=repr(error)
        save()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mesh',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--cores',type=int,default=2)
    p.add_argument('--iterations',type=int,default=5000);p.add_argument('--tolerance',type=float,default=1e-8)
    p.add_argument('--coupling',choices=['simple','coupled'],default='simple')
    p.add_argument('--probe',action='store_true');execute(p.parse_args())
