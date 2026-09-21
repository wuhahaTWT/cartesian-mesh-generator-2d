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
    require(summary['method']=='first-order Rusanov / forward Euler', 'unsupported method')
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
        ql=old[edge.owner];pl=primitive(ql,gamma)
        if neighbour>=0:qr=old[neighbour];pr=primitive(qr,gamma)
        else:pr=ghost(pl,boundary_condition,normal,gamma);qr=conservative(pr,gamma)
        speed=max(abs(pl[1]*normal[0]+pl[2]*normal[1])+acoustic(pl,gamma),abs(pr[1]*normal[0]+pr[2]*normal[1])+acoustic(pr,gamma))
        close(float(row['waveSpeed']),speed,'face acoustic speed')
        fl,fr=flux(ql,pl,normal),flux(qr,pr,normal)
        expected=[.5*length*(fl[k]+fr[k]-speed*(qr[k]-ql[k])) for k in range(4)]
        values=[float(row[k]) for k in ['mass','momentumX','momentumY','energy']]
        for k in range(4):
            close(values[k],expected[k],'Rusanov flux equation')
            max_flux_error=max(max_flux_error,abs(values[k]-expected[k])/(1+abs(expected[k])))
            residual[edge.owner][k]+=values[k];absolute[edge.owner][k]+=abs(values[k])
            if edge.neighbour>=0:residual[edge.neighbour][k]-=values[k];absolute[edge.neighbour][k]+=abs(values[k])
            else:boundary[k].append(values[k])
        spectral[edge.owner]+=speed*length
        if edge.neighbour>=0:spectral[edge.neighbour]+=speed*length
        if boundary_condition and boundary_condition['kind']=='slip-wall':
            close(values[0],0,'impermeable wall mass',absolute=1e-13);close(values[3],0,'adiabatic inviscid wall energy',absolute=1e-13)
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
    close(previous_time,summary['time'],'summary physical time')
    if summary['targetReached']:require(summary['time']==summary['requestedEndTime'] and summary['status']=='target_reached','false target-time completion')
    return dict(valid=True,cells=len(cells),faces=len(faces),time=summary['time'],acousticCourant=cfl,
        maximumCellBalanceRelative=max_local,globalBalanceRelative=global_errors,maximumFluxRelative=max_flux_error,
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
