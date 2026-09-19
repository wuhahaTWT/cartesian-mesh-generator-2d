#!/usr/bin/env python3
"""Independent geometry, constitutive face flux and scalar balance readback.

No external CFD run or engineering certification is implied. Geometry utilities
are shared with the independent flow reader, never with the native solver.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import time as clock
from pathlib import Path
import verify_native_flow as native


def require(ok, message):
    if not ok:
        raise ValueError(message)


def finite(value):
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError('invalid scalar numeric evidence') from exc
    require(math.isfinite(result), 'nonfinite scalar evidence')
    return result


def same(actual, expected, label, absolute=1e-10, relative=1e-8):
    require(math.isclose(finite(actual), finite(expected), abs_tol=absolute, rel_tol=relative), label)


def verify(prefix: Path) -> dict:
    meta_path = Path(str(prefix)+'.json')
    info = json.loads(meta_path.read_text())
    require(info['format'] == 'cartmesh2d-scalar-transport-v1', 'unsupported scalar summary')
    require(info['converged'] is True and info['status'] == 'converged', 'scalar did not converge')
    mesh_path = Path(info['mesh'])
    mesh = native.read_cm2d(mesh_path)
    measured = native.measure(mesh, 1e-10, 1e-9)
    require(not measured.issues, str(measured.issues))
    geo = native.face_geometry(mesh, measured)
    def rows(suffix):
        with Path(str(prefix)+suffix).open() as stream:
            return list(csv.DictReader(stream))
    cells, faces = rows('.cells.csv'), rows('.faces.csv')
    require(len(cells) == len(mesh.cells) == info['cells'], 'cell count mismatch')
    require(len(faces) == len(mesh.edges) == info['faces'], 'face count mismatch')
    k, dt, time = map(finite, (info['diffusivity'], info['timeStep'], info['time']))
    require(k > 0 and dt >= 0 and time >= 0, 'invalid scalar controls')
    require(info['convection'] in ('upwind', 'limited-linear'), 'unknown convection')
    for key in ('relativeTolerance', 'absoluteTolerance', 'cellTolerance'):
        require(0 < finite(info[key]) <= (1e-9 if key != 'absoluteTolerance' else 1e-12), 'verification does not permit relaxed stops')
    mode = info['verification']
    require(mode in ('', 'sine', 'variable-sine', 'decay', 'thermal-vortex'), 'unknown verification mode')
    u = [finite(c['value']) for c in cells]
    q = [finite(f['volumeFlux']) for f in faces]
    bc = [None] * len(faces)
    inputs = [mesh_path, meta_path, Path(__file__), Path(native.__file__)]
    model = info.get('diffusivityModel', 'uniform')
    require(model in ('uniform', 'face-values'), 'unknown diffusivity model')
    face_k = [k] * len(faces)
    if model == 'face-values':
        require(not info.get('evolvingFlow') and not info.get('restart'), 'variable diffusivity cannot use joint restart/evolving flow')
        if mode == 'variable-sine':
            require(not info.get('diffusivityFile'), 'manufactured diffusivity cannot use a file')
            require(all(x > -1 for x,y in mesh.vertices), 'nonpositive manufactured domain diffusivity')
            face_k = [k*(1+g.centre[0]) for g in geo]
        else:
            require(not mode and bool(info.get('diffusivityFile')), 'face diffusivity requires explicit input')
            field_path = Path(info['diffusivityFile']); inputs.append(field_path)
            face_k = [None] * len(faces)
            with field_path.open() as stream:
                lines = stream.read().splitlines()
                require(all(lines), 'empty diffusivity CSV row')
                reader = csv.DictReader(lines)
                require(reader.fieldnames == ['face','diffusivity'], 'invalid diffusivity CSV header')
                for row in reader:
                    require(set(row) == {'face','diffusivity'}, 'invalid diffusivity row')
                    idx = finite(row['face'])
                    require(idx.is_integer() and 0 <= idx < len(faces), 'invalid diffusivity face ID')
                    idx = int(idx)
                    require(face_k[idx] is None, 'duplicate diffusivity face')
                    face_k[idx] = finite(row['diffusivity'])
            require(all(d is not None for d in face_k), 'missing diffusivity face')
        require(all(math.isfinite(d) and d > 0 for d in face_k), 'nonpositive face diffusivity')
        for row,d in zip(faces,face_k):
            same(row['diffusivity'],d,'face diffusivity mismatch',absolute=0.,relative=1e-13)
    else:
        require(mode != 'variable-sine' and not info.get('diffusivityFile'), 'missing face diffusivity model')
        require(all('diffusivity' not in row for row in faces), 'unexpected variable coefficients in uniform output')
    if mode in ('sine','variable-sine'):
        require(dt == 0 and time == 0, 'sine verification must be steady')
    if mode:
        pi = math.pi
        def exact(x, y):
            return math.sin(pi*x)*math.sin(pi*y)*math.exp(-2*pi*pi*k*time if mode not in ('sine','variable-sine') else 0)
        speed = finite(info['verificationSpeed']) if mode in ('sine','variable-sine') else 0
        for e, g in zip(mesh.edges, geo):
            if mode != 'thermal-vortex':
                same(q[e.id], speed*g.area_vector[0], 'manufactured carrier mismatch')
            if e.neighbour < 0:
                bc[e.id] = ('value', exact(*g.centre), None)
    else:
        bc_path, carrier_path = Path(info['boundaryFile']), Path(info['carrierCheckpoint'])
        inputs += [bc_path, carrier_path]
        with bc_path.open() as stream:
            for row in csv.DictReader(stream):
                idx = int(row['face'])
                require(0 <= idx < len(faces) and mesh.edges[idx].neighbour < 0 and bc[idx] is None, 'invalid duplicate boundary row')
                require(row['type'] in ('value', 'flux'), 'unknown boundary type')
                bc[idx] = (row['type'], finite(row['value']), finite(row['inflowValue']) if row['inflowValue'] else None)
        flux_lines = [line.split() for line in carrier_path.read_text().splitlines() if line.startswith('FLUX ')]
        require(len(flux_lines) == 1 and int(flux_lines[0][1]) == len(faces), 'checkpoint flux count mismatch')
        carrier_q = list(map(finite, flux_lines[0][2:]))
        require(carrier_q == q, 'scalar carrier differs from imported checkpoint')
    if info.get('evolvingFlow'):
        require(dt > 0 and time > 0, 'invalid evolving flow clock')
        same(info['acceptedTime'],time,'joint accepted clock mismatch',absolute=1e-14,relative=1e-14)
        same(info['carrierTime'],time,'lagged carrier time',absolute=1e-14,relative=1e-14)
        carrier_path = Path(info['carrierCheckpoint'])
        joint_path = Path(str(prefix)+'.thermal.checkpoint')
        inputs += [carrier_path,joint_path,Path(str(prefix)+'.thermal-history.csv')]
        carrier_text=carrier_path.read_text();joint=joint_path.read_text()
        parts=joint.split('\nFLOW\n',1)
        require(len(parts)==2 and parts[1]==carrier_text,'joint carrier differs from diagnostic checkpoint')
        scalar_lines=[line.split() for line in joint.splitlines() if line.startswith('SCALAR ')]
        require(len(scalar_lines)==1 and int(scalar_lines[0][1])==len(u), 'joint scalar dimensions mismatch')
        require(list(map(finite,scalar_lines[0][2:]))==u,'joint scalar differs from output')
        flux_lines=[line.split() for line in carrier_text.splitlines() if line.startswith('FLUX ')]
        require(len(flux_lines)==1 and int(flux_lines[0][1])==len(q), 'joint carrier flux count mismatch')
        require(list(map(finite,flux_lines[0][2:]))==q,'scalar used a different carrier flux')
        times=[line.split() for line in carrier_text.splitlines() if line.startswith('TIME ')]
        require(len(times)==1,'joint checkpoint clock missing')
        same(times[0][1],time,'joint checkpoint time mismatch',absolute=1e-14,relative=1e-14)
        with Path(str(prefix)+'.thermal-history.csv').open() as stream:
            thermal_history=list(csv.DictReader(stream))
        require(bool(thermal_history) and all(int(row['accepted'])==1 for row in thermal_history),'unaccepted thermal history')
        same(thermal_history[-1]['time'],time,'thermal history clock mismatch',absolute=1e-14,relative=1e-14)
        for a,b in zip(thermal_history,thermal_history[1:]):
            same(float(b['time'])-float(a['time']),dt,'thermal history step mismatch',absolute=1e-14,relative=1e-12)
    for e in mesh.edges:
        if e.neighbour < 0:
            require(bc[e.id] is not None, 'missing boundary condition')
            if q[e.id] < 0:
                require(bc[e.id][0] == 'value' or bc[e.id][2] is not None, 'missing inflow scalar')
    # Independently solve the two-coordinate least-squares normal equations.
    gradients = []
    for cell, centre in zip(mesh.cells, measured.centroids):
        xx = xy = yy = bx = by = 0.
        for fid in cell.edges:
            e, g = mesh.edges[fid], geo[fid]
            j = e.neighbour if e.owner == cell.id else e.owner
            if j >= 0 or bc[fid][0] == 'value':
                other = measured.centroids[j] if j >= 0 else g.centre
                dx, dy = other[0]-centre[0], other[1]-centre[1]
                length = math.hypot(dx, dy)
                delta = ((u[j] if j >= 0 else bc[fid][1])-u[cell.id])/length
                dx, dy = dx/length, dy/length
            else:
                sx, sy = g.area_vector
                length = math.hypot(sx, sy)
                dx, dy, delta = sx/length, sy/length, -bc[fid][1]/face_k[fid]
            xx += dx*dx; xy += dx*dy; yy += dy*dy
            bx += dx*delta; by += dy*delta
        det = xx*yy-xy*xy
        require(det > 0, 'rank deficient independent stencil')
        gradients.append(((yy*bx-xy*by)/det, (xx*by-xy*bx)/det))
    fixed = [b is not None and (b[0] == 'value' or q[i] < 0) for i,b in enumerate(bc)]
    boundary = [(b[1] if b[0] == 'value' else (b[2] or 0)) if b else 0 for b in bc]
    limiter = native.face_limiter(mesh, measured, u, gradients, boundary, fixed) if info['convection'] == 'limited-linear' else [0.] * len(u)
    residual, base, diagonal, continuity, flow_scale, temporal, sources = ([0.] * len(u) for _ in range(7))
    errors = []
    for i, (row, centre, area) in enumerate(zip(cells, measured.centroids, measured.areas)):
        require(int(row['cell']) == i, 'cell ordering mismatch')
        same(row['x'], centre[0], 'cell centre x mismatch'); same(row['y'], centre[1], 'cell centre y mismatch')
        same(row['area'], area, 'cell area mismatch')
        if mode in ('sine','variable-sine'):
            x,y = centre
            d = k*(1+x) if mode == 'variable-sine' else k
            grad_d = k if mode == 'variable-sine' else 0.
            src = 2*d*math.pi**2*exact(x,y)+(speed-grad_d)*math.pi*math.cos(math.pi*x)*math.sin(math.pi*y)
        else:
            src = 0. if mode in ('decay','thermal-vortex') else finite(info['constantSource'])
        sources[i] = src*area
        same(row['sourceIntegral'], sources[i], 'source integral mismatch')
        temporal[i] = area*(u[i]-finite(row['previous']))/dt if dt else 0.
        same(row['temporalIntegral'], temporal[i], 'time term mismatch')
        base[i] = sources[i]+(area*finite(row['previous'])/dt if dt else 0.)
        diagonal[i] = area/dt if dt else 0.
        residual[i] = temporal[i]-sources[i]
        if mode:
            same(row['exact'], exact(*centre), 'manufactured exact field mismatch')
            errors.append(u[i]-exact(*centre))
    boundary_fluxes, constitutive_errors = [], []
    for fid, (row, e, g) in enumerate(zip(faces, mesh.edges, geo)):
        i,j = e.owner,e.neighbour
        require((int(row['face']), int(row['owner']), int(row['neighbour'])) == (fid,i,j), 'face incidence mismatch')
        flow = q[fid]; d = face_k[fid]*g.transmissibility
        continuity[i] += flow; flow_scale[i] += abs(flow)
        if j >= 0:
            continuity[j] -= flow; flow_scale[j] += abs(flow)
            diagonal[i] += d+max(flow,0); diagonal[j] += d+max(-flow,0)
            gf = tuple((1-g.neighbour_weight)*gradients[i][axis]+g.neighbour_weight*gradients[j][axis] for axis in (0,1))
        else:
            gf = gradients[i]
            if bc[fid][0] == 'value':
                diagonal[i] += d; base[i] += d*bc[fid][1]
            else:
                base[i] -= bc[fid][1]*math.hypot(*g.area_vector)
            if flow < 0:
                base[i] -= flow*boundary[fid]
            else:
                diagonal[i] += flow
        if j < 0 and bc[fid][0] == 'flux':
            diffusion = bc[fid][1]*math.hypot(*g.area_vector)
        else:
            diffusion = d*(u[i]-(u[j] if j >= 0 else bc[fid][1]))-face_k[fid]*sum(a*b for a,b in zip(gf,g.correction))
        up = j if j >= 0 and flow < 0 else i
        advected = u[up]+limiter[up]*sum(gradients[up][axis]*(g.centre[axis]-measured.centroids[up][axis]) for axis in (0,1))
        if j < 0 and flow < 0:
            advected = boundary[fid]
        advection = flow*advected
        same(row['advectiveFlux'], advection, 'advective constitutive flux mismatch')
        same(row['diffusiveFlux'], diffusion, 'diffusive constitutive flux mismatch')
        constitutive_errors += [abs(finite(row['advectiveFlux'])-advection), abs(finite(row['diffusiveFlux'])-diffusion)]
        total = advection+diffusion
        residual[i] += total
        if j >= 0:
            residual[j] -= total
        else:
            boundary_fluxes.append(total)
    for defect,scale in zip(continuity, flow_scale):
        require(abs(defect) <= 1e-12+1e-8*scale, 'carrier violates cell continuity')
    norm = math.sqrt(math.fsum(r*r for r in residual))
    stop = info['absoluteTolerance']+info['relativeTolerance']*math.sqrt(math.fsum(b*b for b in base))
    scaled = max(abs(r)/d for r,d in zip(residual,diagonal))
    # Roundoff allowance for the independent geometry reconstruction only.
    require(norm <= stop+1e-11 and scaled <= info['cellTolerance']+1e-11, 'independent full equation is not converged')
    boundary_total = math.fsum(boundary_fluxes)
    global_balance = math.fsum(temporal)+boundary_total-math.fsum(sources)
    for key,value in dict(boundaryFlux=boundary_total, temporalIntegral=math.fsum(temporal), sourceIntegral=math.fsum(sources), globalBalance=global_balance).items():
        same(info[key], value, key+' report mismatch', absolute=1e-9)
    report = dict(valid=True, cells=len(u), faces=len(faces), residualNorm=norm, residualStop=stop,
        maxDiagonalScaledImbalance=scaled, maxCarrierImbalance=max(map(abs,continuity)),
        globalBalance=global_balance, boundaryAbsoluteFlux=math.fsum(map(abs,boundary_fluxes)),
        minValue=min(u), maxValue=max(u), maxConstitutiveDifference=max(constitutive_errors),
        l2Error=native.weighted_l2(errors, measured.areas) if mode else None,
        diffusivityModel=model, minFaceDiffusivity=min(face_k), maxFaceDiffusivity=max(face_k),
        scope='independent actual geometry, constitutive flux, local/global scalar and time balance; not external CFD or accuracy certification')
    inputs += [Path(str(prefix)+s) for s in ('.cells.csv','.faces.csv','.history.csv')]
    report['sha256'] = {str(p):native.sha256_file(p) for p in inputs}
    return report

def generate(root: Path, mesh_cli: Path, transport_cli: Path, variable=False, shear=0.) -> dict:
    """Fresh native meshes; fixed analytic problem and common-time dt sequence."""
    require(not root.exists(), 'generation output must be a new directory')
    require(math.isfinite(shear) and 0 <= shear <= .5 and (variable or shear == 0),
            'shear requires variable study and must lie in [0,.5]')
    root.mkdir(parents=True)
    summary = dict(valid=False, scope='prescribed variable face diffusion spatial refinement; no turbulence' if variable else 'scalar spatial/time refinement only; no coupled thermal flow certification',
                   meshRuns=[], runs=[], spatial={}, temporal=[], executables={}, shear=shear)
    for exe in (mesh_cli, transport_cli):
        summary['executables'][str(exe)] = native.sha256_file(exe)
    def execute(command, label):
        start = clock.monotonic()
        record = native.run(command, root/'logs'/label, 90.)
        record['wallSeconds'] = clock.monotonic()-start
        record['command'] = command
        require(record.get('returncode') == 0 and not record.get('timedOut'), f'{label} native run failed: {record}')
        return record
    try:
        meshes = []
        for level in (4,5,6):
            label = f'l{level}'
            prefix = root/'meshes'/label/'mesh'
            prefix.parent.mkdir(parents=True)
            request = native.Request(label,'manufactured',level,1/((1<<level)-2),.1,1.)
            mesh_cmd = native.mesh_command(mesh_cli,prefix,request,root)
            if shear:
                Path(mesh_cmd[1]).write_text(f'0 0\n1 0\n{1+shear:.17g} 1\n{shear:.17g} 1\n')
            summary['meshRuns'].append(execute(mesh_cmd, 'mesh-'+label))
            mesh = Path(str(prefix)+'.solver.cm2d'); meshes.append(mesh)
            for scheme in ('upwind','limited-linear'):
                out = root/scheme/label
                command = [str(transport_cli),'--mesh',str(mesh),'--output',str(out),'--verification','variable-sine' if variable else 'sine',
                           '--speed','.35','--diffusivity','.08','--convection',scheme]
                summary['runs'].append(execute(command,scheme+'-'+label))
                audit = verify(out); audit['prefix']=str(out)
                summary['spatial'].setdefault(scheme,[]).append(audit)
        for steps in (() if variable else (5,10,20)):
            out=root/'decay'/str(steps)
            command=[str(transport_cli),'--mesh',str(meshes[1]),'--output',str(out),'--verification','decay',
                     '--diffusivity','.2','--dt',str(.5/steps),'--steps',str(steps)]
            summary['runs'].append(execute(command,'decay-'+str(steps)))
            audit=verify(out);audit.update(prefix=str(out),time=.5,dt=.5/steps)
            summary['temporal'].append(audit)
        for name,series in list(summary['spatial'].items())+[('time',summary['temporal'])]:
            for coarse,fine in zip(series,series[1:]):
                require(0<fine['l2Error']<coarse['l2Error'], name+' error did not decrease')
                ratio=coarse['dt']/fine['dt'] if name=='time' else math.sqrt(fine['cells']/coarse['cells'])
                fine['observedOrder']=math.log(coarse['l2Error']/fine['l2Error'])/math.log(ratio)
        summary['valid']=True
    except (ValueError, OSError) as exc:
        summary['error']=str(exc)
    (root/'audit.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--prefix',type=Path)
    group.add_argument('--generate',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--mesh-cli',type=Path,default=Path('build/cartmesh2d_cli'))
    parser.add_argument('--transport-cli',type=Path,default=Path('build/cartmesh2d_transport_cli'))
    parser.add_argument('--variable-diffusivity',action='store_true',help='generate steady variable-sine spatial series instead of constant-D spatial/time series')
    parser.add_argument('--shear',type=float,default=0.,help='variable study domain: x -> x + shear*y, range [0,.5]; genuine Cartesian cut cells')
    args = parser.parse_args()
    if args.variable_diffusivity and not args.generate:parser.error('--variable-diffusivity requires --generate')
    if args.shear and not (args.variable_diffusivity and args.generate):parser.error('--shear requires --variable-diffusivity --generate')
    if args.generate:
        result=generate(args.generate.resolve(),args.mesh_cli.resolve(),args.transport_cli.resolve(),args.variable_diffusivity,args.shear)
    else:
        if not args.output: parser.error('--prefix requires --output')
        try:
            result = verify(args.prefix)
        except (ValueError,KeyError,OSError,ZeroDivisionError) as exc:
            result = dict(valid=False,error=str(exc))
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('sha256','spatial','temporal','meshRuns','runs')},indent=2))
    raise SystemExit(0 if result['valid'] else 1)
