#!/usr/bin/env python3
"""Independent ideal-gas Euler polygon/flux/time audit and exact Sod reference.

Only the mesh reader and polygon measurements are reused from the incompressible
verifier. No solver code, native residual or native pass flag is an oracle.
Reference: https://www.clawpack.org/riemann_book/html/Euler.html
"""
import argparse
import csv
import json
import math
from pathlib import Path
import shlex
import verify_native_flow as geometry


def require(ok, message):
    if not ok:
        raise ValueError(message)


def primitive(q, gamma):
    rho, mx, my, energy = q
    require(all(math.isfinite(x) for x in q) and rho > 0, 'invalid conservative state')
    u, v = mx / rho, my / rho
    pressure = (gamma - 1) * math.fsum([energy, -.5 * mx * u, -.5 * my * v])
    require(pressure > 0 and math.isfinite(pressure), 'nonpositive pressure')
    return rho, u, v, pressure


def conservative(q, gamma):
    rho, u, v, p = q
    return rho, rho * u, rho * v, p / (gamma - 1) + .5 * rho * (u*u + v*v)


def acoustic(q, gamma):
    return math.sqrt(gamma * q[3] / q[0])


def ghost(q, boundary, normal, gamma):
    rho, u, v, p = q
    nx, ny = normal
    vn = u * nx + v * ny
    if boundary['kind'] == 'slip-wall':
        return rho, u - 2*vn*nx, v - 2*vn*ny, p
    if boundary['kind'] == 'transmissive':
        return q
    require(boundary['kind'] == 'farfield', 'unsupported physical boundary')
    far = boundary['reference']
    a = acoustic(q, gamma)
    if vn >= a:
        return q
    if vn <= -a:
        return far
    vfar = far[1]*nx + far[2]*ny
    rout, rin = vn + 2*a/(gamma-1), vfar - 2*acoustic(far, gamma)/(gamma-1)
    vb, ab = (rout+rin)/2, (gamma-1)*(rout-rin)/4
    require(ab > 0, 'farfield nonpositive sound speed')
    upstream = far if vb < 0 else q
    entropy = upstream[3] / upstream[0]**gamma
    rhob = (ab**2 / (gamma*entropy))**(1/(gamma-1))
    old = upstream[1]*nx + upstream[2]*ny
    return rhob, upstream[1]+(vb-old)*nx, upstream[2]+(vb-old)*ny, entropy*rhob**gamma


def flux(q, primitive_q, normal):
    rho, u, v, p = primitive_q
    nx, ny = normal
    un = u*nx + v*ny
    return rho*un, q[1]*un+p*nx, q[2]*un+p*ny, (q[3]+p)*un


def numerical_flux(left, right, normal, gamma, scheme, restoration=1):
    """Independent HLLC star-state formula (Toro eq. 71), or Rusanov."""
    l, r = primitive(left, gamma), primitive(right, gamma)
    nx, ny = normal
    vl, vr = l[1]*nx+l[2]*ny, r[1]*nx+r[2]*ny
    al, ar = acoustic(l, gamma), acoustic(r, gamma)
    fl, fr = flux(left, l, normal), flux(right, r, normal)
    speed = max(abs(vl)+al, abs(vr)+ar)
    fallback = False
    if scheme == 'hllc':
        wl, wr = math.sqrt(l[0]), math.sqrt(r[0])
        ux, uy = ((wl*l[k]+wr*r[k])/(wl+wr) for k in (1, 2))
        h = (wl*(left[3]+l[3])/l[0]+wr*(right[3]+r[3])/r[0])/(wl+wr)
        a2 = (gamma-1)*(h-(ux*ux+uy*uy)/2)
        usable = math.isfinite(a2) and a2 > 0
        sound = math.sqrt(a2) if usable else 0
        sl, sr = min(vl-al, vr-ar, ux*nx+uy*ny-sound), max(vl+al, vr+ar, ux*nx+uy*ny+sound)
        if usable:
            speed = max(speed, abs(sl), abs(sr))
            if sl >= 0:
                return fl, speed, False
            if sr <= 0:
                return fr, speed, False
        try:
            dl, dr = l[0]*(sl-vl), r[0]*(sr-vr)
            sm = (r[3]-l[3]+dl*vl-dr*vr)/(dl-dr)
            usable = usable and math.isfinite(sm) and sl < sm < sr
            stars = []
            for q, u, sk, vk in [(l, left, sl, vl), (r, right, sr, vr)]:
                density = q[0]*(sk-vk)/(sk-sm)
                star = (density, density*(q[1]+(sm-vk)*nx), density*(q[2]+(sm-vk)*ny),
                        density*(u[3]/q[0]+(sm-vk)*(sm+q[3]/(q[0]*(sk-vk)))))
                primitive(star, gamma)
                usable = usable and q[3]+q[0]*(sk-vk)*(sm-vk) > 0
                stars.append(star)
            if usable:
                f, u, star, wave = (fl, left, stars[0], sl) if sm >= 0 else (fr, right, stars[1], sr)
                hllc = [f[k]+wave*(star[k]-u[k]) for k in range(4)]
                hlle = [(sr*fl[k]-sl*fr[k]+sl*sr*(right[k]-left[k]))/(sr-sl) for k in range(4)]
                return tuple(restoration*a+(1-restoration)*b for a,b in zip(hllc,hlle)), speed, False
        except (ValueError, ZeroDivisionError, OverflowError):
            pass
        fallback = True
    return tuple((fl[k]+fr[k]-speed*(right[k]-left[k]))/2 for k in range(4)), speed, fallback


def spatial_reference(mesh, measured, bc, states, gamma, scheme, order):
    """Rebuild face states from polygons, with no exported slopes/stages as oracle."""
    primitives = [primitive(q, gamma) for q in states]
    face_geometry = []
    for edge in mesh.edges:
        a, b = owner_segment(mesh, edge)
        sx, sy = b[1]-a[1], a[0]-b[0]
        length = math.hypot(sx, sy)
        face_geometry.append(((a[0]+b[0])/2, (a[1]+b[1])/2, sx/length, sy/length, length))
    contact_weight = [1.]*len(states)
    if scheme == 'hllc':
        for edge in mesh.edges:
            x, y, nx, ny, length = face_geometry[edge.id]
            j = edge.neighbour
            if j < 0 and bc[edge.id]['kind'] == 'periodic': j = mesh.edges[bc[edge.id]['partner']].owner
            pl = primitives[edge.owner][3]
            pr = primitives[j][3] if j >= 0 else ghost(primitives[edge.owner],bc[edge.id],(nx,ny),gamma)[3]
            weight = (min(pl,pr)/max(pl,pr))**3
            contact_weight[edge.owner] = min(contact_weight[edge.owner],weight)
            if j >= 0: contact_weight[j] = min(contact_weight[j],weight)
    gradients = [[[0., 0.] for _ in range(4)] for _ in states]
    reconstruction_fallback = 0
    if order == 2:
        for i, cell in enumerate(mesh.cells):
            centre = measured.centroids[i]
            ex, ey = face_geometry[cell.edges[0]][2:4]
            def local(q): return q[0], ex*q[1]+ey*q[2], -ey*q[1]+ex*q[2], q[3]
            q = local(primitives[i])
            stencil = []
            for face_id in cell.edges:
                edge = mesh.edges[face_id]
                x, y, nx, ny, _ = face_geometry[face_id]
                if edge.neighbour >= 0:
                    j = edge.neighbour if edge.owner == i else edge.owner
                    point, neighbour = measured.centroids[j], primitives[j]
                    dx, dy = point[0]-centre[0], point[1]-centre[1]
                elif bc[face_id]['kind'] == 'periodic':
                    partner = bc[face_id]['partner']; j = mesh.edges[partner].owner
                    px, py = face_geometry[partner][:2]; point = measured.centroids[j]
                    dx, dy = point[0]+x-px-centre[0], point[1]+y-py-centre[1]
                    neighbour = primitives[j]
                else:
                    distance = (x-centre[0])*nx+(y-centre[1])*ny
                    dx, dy = 2*distance*nx, 2*distance*ny
                    neighbour = ghost(primitives[i], bc[face_id], (nx, ny), gamma)
                stencil.append((dx, dy, 1/(dx*dx+dy*dy), local(neighbour)))
            xx = math.fsum(w*dx*dx for dx, dy, w, v in stencil)
            xy = math.fsum(w*dx*dy for dx, dy, w, v in stencil)
            yy = math.fsum(w*dy*dy for dx, dy, w, v in stencil)
            det = xx*yy-xy*xy
            if det <= 64*math.ulp(1.)*(xx+yy)**2:
                reconstruction_fallback += 1
                continue
            for k in range(4):
                bx = math.fsum(w*dx*(v[k]-q[k]) for dx, dy, w, v in stencil)
                by = math.fsum(w*dy*(v[k]-q[k]) for dx, dy, w, v in stencil)
                gx, gy = (yy*bx-xy*by)/det, (xx*by-xy*bx)/det
                low, high = min([q[k]]+[v[k] for dx, dy, w, v in stencil]), max([q[k]]+[v[k] for dx, dy, w, v in stencil])
                theta = 1.
                for face_id in cell.edges:
                    x, y = face_geometry[face_id][:2]
                    delta = gx*(x-centre[0])+gy*(y-centre[1])
                    if delta > 0: theta = min(theta, (high-q[k])/delta)
                    elif delta < 0: theta = min(theta, (low-q[k])/delta)
                theta = min(1., max(0., theta)); gradients[i][k] = [gx*theta, gy*theta]
            if any(q[k]+gradients[i][k][0]*(face_geometry[f][0]-centre[0])+gradients[i][k][1]*(face_geometry[f][1]-centre[1]) <= 0 for f in cell.edges for k in (0, 3)):
                gradients[i] = [[0., 0.] for _ in range(4)]; reconstruction_fallback += 1
            gn, gt = gradients[i][1], gradients[i][2]
            gradients[i][1] = [ex*gn[k]-ey*gt[k] for k in (0, 1)]
            gradients[i][2] = [ey*gn[k]+ex*gt[k] for k in (0, 1)]
    def face_state(i, face_id):
        if order == 1: return states[i]
        x, y = face_geometry[face_id][:2]; cx, cy = measured.centroids[i]
        q = [primitives[i][k]+gradients[i][k][0]*(x-cx)+gradients[i][k][1]*(y-cy) for k in range(4)]
        return conservative(q, gamma)
    fluxes, speeds, fallback, scales = [], [], [], []
    for edge in mesh.edges:
        face_id = edge.id; x, y, nx, ny, length = face_geometry[face_id]
        left = face_state(edge.owner, face_id)
        if edge.neighbour >= 0: right = face_state(edge.neighbour, face_id)
        elif bc[face_id]['kind'] == 'periodic':
            partner = bc[face_id]['partner']; right = face_state(mesh.edges[partner].owner, partner)
        else: right = conservative(ghost(primitive(left, gamma), bc[face_id], (nx, ny), gamma), gamma)
        neighbour = edge.neighbour
        if neighbour < 0 and bc[face_id]['kind'] == 'periodic': neighbour = mesh.edges[bc[face_id]['partner']].owner
        restoration = min(contact_weight[edge.owner],contact_weight[neighbour]) if neighbour >= 0 else contact_weight[edge.owner]
        f, speed, used_fallback = numerical_flux(left, right, (nx, ny), gamma, scheme, restoration)
        l, r = primitive(left, gamma), primitive(right, gamma)
        if bc.get(face_id, {}).get('kind') == 'slip-wall':
            # Independent scalar mirror-Riemann pressure formula; both HLLC and
            # HLLE give p + rho*un*(un+S), with symmetric bounding speed S.
            un = l[1]*nx+l[2]*ny
            pressure = l[3]+l[0]*un*(un+speed)
            f = (0., pressure*nx, pressure*ny, 0.)
        # Componentwise dimensional characteristic flux, not the possibly
        # vanishing result of cancellation. This scales covariantly under a
        # change of density/velocity units and includes the actual face length.
        velocity = max(math.hypot(*q[1:3])+acoustic(q,gamma) for q in (l,r))
        density = max(l[0],r[0]); enthalpy = max(left[3]+l[3],right[3]+r[3])
        scales.append([length*density*velocity, length*density*velocity**2,
                       length*density*velocity**2, length*enthalpy*velocity])
        fluxes.append([length*v for v in f]); speeds.append(speed); fallback.append(int(used_fallback))
    return fluxes, speeds, fallback, reconstruction_fallback, min(contact_weight), scales


def stage_update(mesh, measured, states, fluxes, dt):
    residual = [[0.]*4 for _ in states]
    for edge, f in zip(mesh.edges, fluxes):
        for k in range(4):
            residual[edge.owner][k] += f[k]
            if edge.neighbour >= 0: residual[edge.neighbour][k] -= f[k]
    return [[q[k]-dt*r[k]/area for k in range(4)] for q, r, area in zip(states, residual, measured.areas)]


def close(actual, expected, message, relative=2e-10, absolute=2e-13):
    require(math.isfinite(actual) and math.isfinite(expected) and
            abs(actual-expected) <= absolute + relative*max(abs(actual), abs(expected)), message)


def owner_segment(mesh, edge):
    cell=mesh.cells[edge.owner];local=cell.edges.index(edge.id)
    return mesh.vertices[cell.vertices[local]],mesh.vertices[cell.vertices[(local+1)%len(cell.vertices)]]


def read_boundaries(path, mesh):
    lines=path.read_text().splitlines()
    require(lines[0]=='CM2D_EULER_BOUNDARY 1' and lines[-1]=='END', 'boundary format')
    result={}
    require(lines[1]==f'MESH {len(mesh.cells)} {len(mesh.edges)}', 'boundary mesh binding')
    for line in lines[2:-1]:
        row=shlex.split(line);require(len(row)==13, 'boundary row size')
        face=int(row[0]);require(0<=face<len(mesh.edges) and mesh.edges[face].neighbour<0 and face not in result, 'boundary coverage')
        edge=mesh.edges[face];a,b=owner_segment(mesh,edge)
        require(int(row[1])==edge.owner,'boundary owner binding')
        for actual,expected in zip(map(float,row[2:6]),[(a[0]+b[0])/2,(a[1]+b[1])/2,b[1]-a[1],a[0]-b[0]]):close(actual,expected,'boundary geometry binding')
        result[face]=dict(kind=row[6],reference=tuple(map(float,row[7:11])),partner=-1 if row[11]=='-' else int(row[11]),name=row[12])
    require(set(result)=={e.id for e in mesh.edges if e.neighbour<0}, 'incomplete boundary coverage')
    for face,b in result.items():
        if b['kind']=='periodic':
            require(b['partner'] in result and result[b['partner']]['partner']==face and result[b['partner']]['name']==b['name'], 'nonreciprocal periodic boundary')
        else:require(b['partner']==-1,'nonperiodic partner')
    return result


def audit(mesh_path, prefix):
    mesh=geometry.read_cm2d(Path(mesh_path));measured=geometry.measure(mesh,1e-11,1e-10)
    require(not measured.issues, str(measured.issues))
    prefix=str(prefix);summary=json.loads(Path(prefix+'.json').read_text())
    scheme, order = summary.get('fluxScheme', 'rusanov'), summary.get('order', 1)
    require(scheme in ('rusanov', 'hllc') and order in (1, 2), 'unsupported numerical controls')
    method = ('first-order ' if order == 1 else 'limited-linear ')+('HLLC-HLLE' if scheme == 'hllc' else 'Rusanov')+(' / forward Euler' if order == 1 else ' / SSPRK2')
    require(summary['method'] == method, 'unsupported method')
    gamma,gas_r=summary['gamma'],summary['gasConstant'];require(gamma>1 and gas_r>0,'gas properties')
    dt=summary['lastStep'];require(dt>0,'no accepted step to audit')
    bc=read_boundaries(Path(prefix+'.boundaries'),mesh)
    cells=list(csv.DictReader(Path(prefix+'.cells.csv').open()));faces=list(csv.DictReader(Path(prefix+'.faces.csv').open()))
    require(len(cells)==len(mesh.cells) and len(faces)==len(mesh.edges),'field size')
    current=[];old=[]
    for i,(row,c,a) in enumerate(zip(cells,measured.centroids,measured.areas)):
        require(int(row['cell'])==i,'cell id')
        for field,value in zip(['x','y','area'],[*c,a]):close(float(row[field]),value,'cell geometry '+field)
        q=tuple(float(row[k]) for k in ['rho','rhoU','rhoV','rhoE']);previous=tuple(float(row[k]) for k in ['previousRho','previousRhoU','previousRhoV','previousRhoE'])
        pq=primitive(q,gamma);primitive(previous,gamma)
        for field,value in zip(['rho','u','v','p'],pq):close(float(row[field]),value,'ideal-gas reconstruction '+field)
        close(float(row['temperature']),pq[3]/(pq[0]*gas_r),'temperature EOS')
        close(float(row['mach']),math.hypot(pq[1],pq[2])/acoustic(pq,gamma),'Mach EOS')
        current.append(q);old.append(previous)
    residual=[[0.]*4 for _ in cells];absolute=[[0.]*4 for _ in cells];spectral=[0.]*len(cells);boundary=[[] for _ in range(4)]
    expected_flux, expected_speed, expected_fallback, reconstruction_fallback, minimum_restoration, flux_scales = spatial_reference(mesh, measured, bc, old, gamma, scheme, order)
    if order == 2:
        stage = stage_update(mesh, measured, old, expected_flux, dt)
        for q in stage: primitive(q, gamma)
        f2, s2, b2, r2, restoration2, scales2 = spatial_reference(mesh, measured, bc, stage, gamma, scheme, order)
        for q in stage_update(mesh, measured, stage, f2, dt): primitive(q, gamma)
        expected_flux = [[(a+b)/2 for a, b in zip(left, right)] for left, right in zip(expected_flux, f2)]
        expected_speed = [max(a,b) for a,b in zip(expected_speed,s2)]
        flux_scales = [[max(a,b) for a,b in zip(s,t)] for s,t in zip(flux_scales,scales2)]
        expected_fallback = [a | (b << 1) for a,b in zip(expected_fallback,b2)]
        reconstruction_fallback += r2
        minimum_restoration = min(minimum_restoration,restoration2)
    close(summary.get('lastMinimumContactRestoration',1), minimum_restoration, 'pressure sensor weight')
    fallback_count = sum(mask.bit_count() for edge, mask in zip(mesh.edges, expected_fallback)
                         if bc.get(edge.id, {}).get('partner', -1) < 0 or edge.id < bc[edge.id]['partner'])
    require(summary.get('lastHllcFallbackEvaluations', 0) == fallback_count, 'HLLC fallback count')
    require(summary.get('lastReconstructionFallbackCells', 0) == reconstruction_fallback, 'reconstruction fallback count')
    max_flux_error=0
    for edge,row in zip(mesh.edges,faces):
        require(int(row['face'])==edge.id and int(row['owner'])==edge.owner and int(row['neighbour'])==edge.neighbour,'face incidence')
        a,b=owner_segment(mesh,edge);sx,sy=b[1]-a[1],a[0]-b[0]
        centre=((a[0]+b[0])/2,(a[1]+b[1])/2);length=math.hypot(sx,sy);normal=(sx/length,sy/length)
        for key,value in zip(['x','y','sx','sy'],[*centre,sx,sy]):close(float(row[key]),value,'face geometry '+key)
        boundary_condition=bc.get(edge.id)
        require(row['kind']==('internal' if boundary_condition is None else boundary_condition['kind']), 'face physical kind')
        neighbour=edge.neighbour
        if boundary_condition and boundary_condition['kind']=='periodic':neighbour=mesh.edges[boundary_condition['partner']].owner
        require(int(row['partner'])==(-1 if boundary_condition is None else boundary_condition['partner']),'periodic face export')
        speed = expected_speed[edge.id]
        close(float(row['waveSpeed']),speed,'face acoustic speed')
        require(int(row.get('hllcFallbackStages', 0)) == expected_fallback[edge.id], 'HLLC fallback stage mask')
        expected=expected_flux[edge.id]
        values=[float(row[k]) for k in ['mass','momentumX','momentumY','energy']]
        for k in range(4):
            error = abs(values[k]-expected[k])/flux_scales[edge.id][k]
            # A floating-point equation check (512 binary64 eps), not a physical
            # accuracy threshold. The same dimensionless gate applies in SI and
            # nondimensional runs, including nearly zero wall/transverse fluxes.
            require(math.isfinite(error) and error <= 512*math.ulp(1.),
                    f'independent numerical flux/stage equation: face {edge.id}, component {k}, normalized error {error}')
            max_flux_error=max(max_flux_error,error)
            residual[edge.owner][k]+=values[k];absolute[edge.owner][k]+=abs(values[k])
            if edge.neighbour>=0:residual[edge.neighbour][k]-=values[k];absolute[edge.neighbour][k]+=abs(values[k])
            else:boundary[k].append(values[k])
        spectral[edge.owner]+=speed*length
        if edge.neighbour>=0:spectral[edge.neighbour]+=speed*length
        if boundary_condition and boundary_condition['kind']=='slip-wall':
            require(values[0] == 0 and values[3] == 0, 'mirror wall mass and energy must be exactly zero')
            tangential = abs(-normal[1]*values[1]+normal[0]*values[2])/flux_scales[edge.id][1]
            require(tangential <= 512*math.ulp(1.), 'inviscid wall tangential traction')
        if boundary_condition and boundary_condition['kind']=='periodic':
            partner=faces[boundary_condition['partner']]
            for k,key in enumerate(['mass','momentumX','momentumY','energy']):require(values[k]==-float(partner[key]),'periodic flux pair differs')
    max_local=0
    for i,area in enumerate(measured.areas):
        for k in range(4):
            error=area*(current[i][k]-old[i][k])+dt*residual[i][k]
            scale=area*(abs(current[i][k])+abs(old[i][k]))+dt*absolute[i][k]
            max_local=max(max_local,abs(error)/(scale or 1))
    require(max_local<1e-12,'cell conservative update failed')
    cfl=max(dt*s/a for s,a in zip(spectral,measured.areas));require(cfl<=summary['cflLimit']*(1+1e-12),'acoustic CFL limit')
    global_errors=[]
    for k in range(4):
        terms=[a*(q[k]-p[k]) for a,q,p in zip(measured.areas,current,old)]+[dt*f for f in boundary[k]]
        scale=math.fsum(a*(abs(q[k])+abs(p[k])) for a,q,p in zip(measured.areas,current,old))+dt*math.fsum(abs(f) for f in boundary[k])
        global_errors.append(abs(math.fsum(terms))/(scale or 1))
    require(max(global_errors)<1e-12,'global conservation failed')
    history=list(csv.DictReader(Path(prefix+'.history.csv').open()));require(len(history)==summary['acceptedSteps'] and history,'accepted history size')
    previous_time=summary['initialTime'];previous_step=summary['steps']-summary['acceptedSteps']
    for row in history:
        values={k:float(v) for k,v in row.items()};require(all(math.isfinite(x) for x in values.values()),'nonfinite history')
        require(int(row['step'])==previous_step+1 and values['time']>previous_time,'history not sequential')
        close(values['time']-previous_time,values['dt'],'history physical time',absolute=1e-14)
        require(values['minimumDensity']>0 and values['minimumPressure']>0 and values['acousticCourant']<=summary['cflLimit']*(1+1e-12),'history positivity/CFL')
        previous_time=values['time'];previous_step+=1
    for k, key in enumerate(('mass','momentumX','momentumY','totalEnergy')):
        close(float(history[-1][key]), math.fsum(a*q[k] for a,q in zip(measured.areas,current)), 'history final integral '+key)
    close(float(history[-1]['minimumDensity']), min(q[0] for q in current), 'history final minimum density')
    close(float(history[-1]['minimumPressure']), min(primitive(q,gamma)[3] for q in current), 'history final minimum pressure')
    close(float(history[-1]['acousticCourant']), cfl, 'history final acoustic CFL')
    for key in ('hllcFallbackEvaluations', 'reconstructionFallbackCells'):
        require(summary.get(key, 0) == sum(int(row.get(key, 0)) for row in history), 'history fallback total')
    close(summary.get('minimumContactRestoration',1), min(float(row.get('minimumContactRestoration',1)) for row in history),'history pressure sensor')
    close(previous_time,summary['time'],'summary physical time')
    if summary['targetReached']:require(summary['time']==summary['requestedEndTime'] and summary['status']=='target_reached','false target-time completion')
    return dict(valid=True,minimumContactRestoration=minimum_restoration,fluxScheme=scheme,order=order,hllcFallbackEvaluations=fallback_count,reconstructionFallbackCells=reconstruction_fallback,cells=len(cells),faces=len(faces),time=summary['time'],acousticCourant=cfl,
        maximumCellBalanceRelative=max_local,globalBalanceRelative=global_errors,maximumFluxRelative=max_flux_error,
        fluxErrorNormalization="face length times characteristic rho*V, rho*V^2, (rhoE+p)*V; V=max(|velocity|+sound); 512 binary64 eps",
        minimumDensity=min(q[0] for q in current),minimumPressure=min(primitive(q,gamma)[3] for q in current),
        note='last accepted step four-equation/flux/EOS audit plus every-row time/CFL/positivity checks; not independent replay of all steps')


def sod(x,time,gamma=1.4,split=.5):
    """Exact self-similar Sod Riemann state, independently bracketed star pressure."""
    left=(1.,0.,1.);right=(.125,0.,.1)
    def wave(p,state):
        rho,u,pk=state;a=math.sqrt(gamma*pk/rho)
        if p>pk:return (p-pk)*math.sqrt((2/((gamma+1)*rho))/(p+(gamma-1)/(gamma+1)*pk))
        return 2*a/(gamma-1)*((p/pk)**((gamma-1)/(2*gamma))-1)
    low,high=.1,1.
    for _ in range(100):
        middle=(low+high)/2
        if wave(middle,left)+wave(middle,right)>0:high=middle
        else:low=middle
    ps=(low+high)/2;us=(wave(ps,right)-wave(ps,left))/2
    if time==0:return (1.,0.,1.) if x<split else (.125,0.,.1)
    xi=(x-split)/time;al=math.sqrt(gamma);ar=math.sqrt(gamma*.1/.125);astar=al*ps**((gamma-1)/(2*gamma))
    if xi<=-al:return left
    if xi<=us-astar:
        u=2/(gamma+1)*(al+xi);a=2/(gamma+1)*(al-(gamma-1)*xi/2)
        return (a/al)**(2/(gamma-1)),u,(a/al)**(2*gamma/(gamma-1))
    if xi<=us:return ps**(1/gamma),us,ps
    shock=ar*math.sqrt((gamma+1)/(2*gamma)*ps/.1+(gamma-1)/(2*gamma))
    if xi<shock:
        ratio=ps/.1;k=(gamma-1)/(gamma+1)
        return .125*(ratio+k)/(k*ratio+1),us,ps
    return right


def rectangle(path,nx,ny,width=1,height=.1):
    vertices=[(width*i/nx,height*j/ny) for j in range(ny+1) for i in range(nx+1)];edges=[];cells=[];pairs={}
    for j in range(ny):
        for i in range(nx):
            a=j*(nx+1)+i;ids=[a,a+1,a+nx+2,a+nx+1];incident=[]
            for x,y in zip(ids,ids[1:]+ids[:1]):
                pair=tuple(sorted((x,y)))
                if pair not in pairs:pairs[pair]=len(edges);edges.append([len(edges),x,y,len(cells),-1,2])
                else:edges[pairs[pair]][4:]=[len(cells),0]
                incident.append(pairs[pair])
            cells.append([len(cells),0,0,width*height/nx/ny,4,*ids,4,*incident])
    path.write_text('\n'.join(['CM2D 1',f'VERTICES {len(vertices)}',*[f'{i} {x:.17g} {y:.17g}' for i,(x,y) in enumerate(vertices)],
        f'EDGES {len(edges)}',*[' '.join(map(str,e)) for e in edges],f'CELLS {len(cells)}',*[' '.join(map(str,c)) for c in cells],
        'AUDIT 0 0 0 0 0 0 0','END','']))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--mesh',type=Path,required=True);parser.add_argument('--prefix',type=Path,required=True)
    args=parser.parse_args();print(json.dumps(audit(args.mesh,args.prefix),indent=2))
