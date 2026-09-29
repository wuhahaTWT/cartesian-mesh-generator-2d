#!/usr/bin/env python3
"""Independent read-back of the Cartesian MAC immersed analysis files.

Audits the exported staggered fields, steady equations and wall interpolation.
This does not qualify a finite-penalty wall as an exact no-slip boundary.
Uses only the Python standard library; visualization is a separate consumer.
"""
import argparse
import csv
import json
import math
from pathlib import Path


class VerificationError(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise VerificationError(message)


def table(path):
    with path.open(newline='', encoding='utf-8') as source:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(source)]
    require(all(math.isfinite(v) for row in rows for v in row.values()), f'nonfinite field: {path.name}')
    return rows


def verify(directory, require_converged=False):
    root = Path(directory)
    summary = json.loads((root/'summary.json').read_text())
    require(summary['format'] == 'cartmesh2d-immersed-flow-v1', 'incorrect format')
    require(summary['product_solver_ready'] is False, 'analysis lattice cannot claim product fluid topology')
    g, c = summary['grid'], summary['controls']
    nx, ny, hx, hy = g['nx'], g['ny'], g['dx'], g['dy']
    require(nx >= 4 and ny >= 4 and hx > 0 and hy > 0, 'invalid grid')
    ref = c['drive']*g['height']**2/(12*c['viscosity'])
    scale = c['drive']
    cells, us, vs = table(root/'cells.csv'), table(root/'u.csv'), table(root/'v.csv')
    walls = table(root/'walls.csv')
    for rows, count, xshift, yshift in ((cells,nx*ny,.5,.5),(us,nx*ny,0,.5),(vs,nx*(ny+1),.5,0)):
        require(len(rows) == count, 'field cardinality mismatch')
        for k, row in enumerate(rows):
            i, j = k % nx, k // nx
            require(row['i'] == i and row['j'] == j, 'duplicate/missing/reordered grid entry')
            require(abs(row['x']-(i+xshift)*hx) <= 1e-12*g['length'] and
                    abs(row['y']-(j+yshift)*hy) <= 1e-12*g['height'], 'field coordinate mismatch')
            require(0 <= row['solid_mask'] <= 1, 'invalid solid mask')
    u, v, p = [r['u'] for r in us], [r['v'] for r in vs], [r['pressure_fluctuation'] for r in cells]
    def U(i,j):
        if j < 0:
            return -U(i,0)
        if j == ny:
            return -U(i,ny-1)
        return u[j*nx+i%nx]
    def V(i,j):
        return v[j*nx+i%nx]
    def P(i,j):
        return p[j*nx+i%nx]
    # Reconstruct the weighted-adjoint force from exported multipliers rather
    # than trusting the Eulerian force or a solver-reported momentum residual.
    marker_path=root/'wall-markers.csv'
    markers=table(marker_path) if marker_path.exists() else []
    surface=c.get('wall_method','brinkman')=='surface-penalty'
    require(surface or not markers, 'wall multipliers in Brinkman-only mode')
    force_u,force_v=[0.]*len(u),[0.]*len(v)
    wall_law=marker_speed=0.
    for row in markers:
        require(row['segment_weight']>0 and row['weight']>0,'invalid surface weight')
        require(abs(row['weight']-row['segment_weight']*min(hx,hy)/(hx*hy))<1e-12*row['weight'],'surface quadrature weight mismatch')
        velocities=[]
        for is_u,field,forces,key in ((True,u,force_u,'multiplier_u'),(False,v,force_v,'multiplier_v')):
            x=row['x']/hx-(0 if is_u else .5)
            y=row['y']/hy-(.5 if is_u else 0)
            i,j=math.floor(x),math.floor(y);a,b=x-i,y-j
            require(0<=i<nx-1 and 0<j<ny-1,'surface stencil outside interior lattice')
            sample=0.
            for di,dj,weight in ((0,0,(1-a)*(1-b)),(1,0,a*(1-b)),(0,1,(1-a)*b),(1,1,a*b)):
                k=(j+dj)*nx+i+di
                sample+=weight*field[k]
                forces[k]-=weight*math.sqrt(row['weight'])*row[key]
            velocities.append(sample)
            wall_law=max(wall_law,abs(sample-c['wall_penalty_time']*row[key]/math.sqrt(row['weight']))/ref)
        marker_speed=max(marker_speed,math.hypot(*velocities))
    for rows,forces in ((us,force_u),(vs,force_v)):
        for row,force in zip(rows,forces):
            require(abs(row.get('wall_force',0)-force)<=1e-11*max(scale,abs(force)), 'spread wall force mismatch')
    require(wall_law<=c['continuity_tolerance'],'coupled wall law failed')
    surface_power=hx*hy*(sum(a*b for a,b in zip(force_u,u))+sum(a*b for a,b in zip(force_v,v)))
    expected_power=-hx*hy*c.get('wall_penalty_time',0)*sum(row['multiplier_u']**2+row['multiplier_v']**2 for row in markers)
    require(abs(surface_power-expected_power)<=1e-8*max(scale*ref*g['length']*g['height'],abs(expected_power)), 'surface adjoint power identity failed')
    def transported(flux, a, b):
        return flux*(a if flux >= 0 else b)
    continuity, momentum, fluxes = 0., 0., []
    error2, exact2, maxspeed = 0., 0., 0.
    classifications = [0,0,0]
    for j in range(ny):
        for i in range(nx):
            k = j*nx+i
            cls = cells[k]['classification']
            require(cls in (0,1,2), 'bad geometric classification')
            classifications[int(cls)] += 1
            div = (U(i+1,j)-U(i,j))/hx+(V(i,j+1)-V(i,j))/hy
            require(abs(div-cells[k]['divergence']) <= 1e-12*ref/g['height'], 'stored divergence mismatch')
            continuity = max(continuity, abs(div)*g['height']/ref)
            uc, vc = (U(i,j)+U(i+1,j))/2, (V(i,j)+V(i,j+1))/2
            require(abs(uc-cells[k]['u']) <= 1e-13*ref and abs(vc-cells[k]['v']) <= 1e-13*ref, 'cell velocity mismatch')
            maxspeed = max(maxspeed, math.hypot(uc,vc))
            lap = (U(i-1,j)-2*U(i,j)+U(i+1,j))/hx**2+(U(i,j-1)-2*U(i,j)+U(i,j+1))/hy**2
            advx = (transported((U(i,j)+U(i+1,j))/2,U(i,j),U(i+1,j))-
                    transported((U(i-1,j)+U(i,j))/2,U(i-1,j),U(i,j)))/hx
            advy = (transported((V(i-1,j+1)+V(i,j+1))/2,U(i,j),U(i,j+1))-
                    transported((V(i-1,j)+V(i,j))/2,U(i,j-1),U(i,j)))/hy
            residual = c['viscosity']*lap-advx-advy-(P(i,j)-P(i-1,j))/hx-us[k]['solid_mask']*U(i,j)/c['penalty_time']+c['drive']+force_u[k]
            momentum = max(momentum, abs(residual)/scale)
            y = (j+.5)*hy
            exact = c['drive']*y*(g['height']-y)/(2*c['viscosity'])
            error2 += (uc-exact)**2+vc**2
            exact2 += exact**2
    for j in range(1,ny):
        for i in range(nx):
            lap = (V(i-1,j)-2*V(i,j)+V(i+1,j))/hx**2+(V(i,j-1)-2*V(i,j)+V(i,j+1))/hy**2
            advx = (transported((U(i+1,j-1)+U(i+1,j))/2,V(i,j),V(i+1,j))-
                    transported((U(i,j-1)+U(i,j))/2,V(i-1,j),V(i,j)))/hx
            advy = (transported((V(i,j)+V(i,j+1))/2,V(i,j),V(i,j+1))-
                    transported((V(i,j-1)+V(i,j))/2,V(i,j-1),V(i,j)))/hy
            residual = c['viscosity']*lap-advx-advy-(P(i,j)-P(i,j-1))/hy-vs[j*nx+i]['solid_mask']*V(i,j)/c['penalty_time']+force_v[j*nx+i]
            momentum = max(momentum, abs(residual)/scale)
    require(all(V(i,0) == 0 and V(i,ny) == 0 for i in range(nx)), 'domain wall leakage')
    fluxes = [sum(U(i,j)*hy for j in range(ny)) for i in range(nx)]
    flux_spread = (max(fluxes)-min(fluxes))/(ref*g['height'])
    def interpolate(point, is_u):
        gx, gy = point[0]/hx-(0 if is_u else .5), point[1]/hy-(.5 if is_u else 0)
        i, j = math.floor(gx), math.floor(gy)
        a, b = gx-i, gy-j
        field = U if is_u else V
        return (1-b)*((1-a)*field(i,j)+a*field(i+1,j))+b*((1-a)*field(i,j+1)+a*field(i+1,j+1))
    wall_max = normal_max = tangent_max = 0.
    wall_flux_net = wall_flux_abs = 0.
    for row in walls:
        point = row['x'], row['y']
        wu, wv = interpolate(point,True), interpolate(point,False)
        require(abs(wu-row['u']) <= 1e-12*ref and abs(wv-row['v']) <= 1e-12*ref, 'wall interpolation mismatch')
        require(abs(math.hypot(row['nx'],row['ny'])-1) < 1e-12, 'wall normal is not unit length')
        wall_flux_net+=(wu*row['nx']+wv*row['ny'])*row['segment_weight']
        wall_flux_abs+=abs(wu*row['nx']+wv*row['ny'])*row['segment_weight']
        wall_max = max(wall_max,math.hypot(wu,wv))
        normal_max = max(normal_max,abs(wu*row['nx']+wv*row['ny']))
        tangent_max = max(tangent_max,abs(wv*row['nx']-wu*row['ny']))
    if markers and c.get('wall_quadrature')=='piecewise-gauss3':
        # Structural interpolation bound, not a wall-accuracy acceptance gate:
        # max_t sum |L_i(t)| = 7/3 for the three Gauss nodes on a quadratic piece.
        require(wall_max<=(7/3)*marker_speed+1e-12*ref, 'wall diagnostic exceeds piecewise Gauss interpolation envelope')
    audit = dict(momentum=momentum,continuity=continuity,flux_spread=flux_spread,max_speed=maxspeed,
                 wall_speed_max=wall_max,wall_normal_speed_max=normal_max,wall_tangential_speed_max=tangent_max)
    if 'wall_method' in c:
        audit.update(wall_constraint_residual=wall_law,marker_speed_max=marker_speed,
                     surface_drag_per_density=-hx*hy*sum(force_u),surface_lift_per_density=-hx*hy*sum(force_v),
                     surface_power_per_density=surface_power)
    if 'wall_normal_flux_net' in summary['metrics']:
        channel_drag=2*c['viscosity']*hx/hy*sum(U(i,0)+U(i,ny-1) for i in range(nx))
        penalty_drag=hx*hy*sum(row['solid_mask']*row['u']/c['penalty_time'] for row in us)
        drive_force=c['drive']*g['length']*g['height']
        audit.update(wall_normal_flux_net=wall_flux_net,wall_normal_flux_abs=wall_flux_abs,
                     channel_wall_drag_per_density=channel_drag,
                     force_balance=abs(drive_force-penalty_drag+hx*hy*sum(force_u)-channel_drag)/drive_force)
    for key, value in audit.items():
        require(abs(value-summary['metrics'][key]) <= 1e-10*max(1,abs(value)), f'summary mismatch: {key}')
    require(classifications == g['classification_counts'], 'classification count mismatch')
    require(abs(g['solid_area_geometry']+g['fluid_area_geometry']-g['length']*g['height']) <= 1e-12*g['length']*g['height'], 'geometric area balance')
    require(continuity <= c['continuity_tolerance'], 'accepted state is not divergence-free to requested tolerance')
    if require_converged:
        require(summary['converged'] is True and summary['stop_reason'] == 'steady-converged', 'not converged')
    if summary['converged']:
        require(momentum <= c['steady_tolerance'] and summary['metrics']['field_change'] <= c['steady_tolerance'], 'false steady convergence')
    require(summary['converged'] or summary['stop_reason'] in ('iteration-limit','candidate-failed','cancelled'), 'unknown stop status')
    history = table(root/'history.csv')
    require(len(history) == summary['steps'], 'accepted history length mismatch')
    for index, row in enumerate(history):
        require(row['step'] == index+1 and row['dt'] > 0, 'invalid step history')
        require(row['continuity'] <= c['continuity_tolerance'], 'history contains rejected state')
        require(row.get('wall_constraint_residual',0)<=c['continuity_tolerance'],'history contains rejected wall law')
        require(row['pressure_linear_residual'] <= 1.05*(1e-13+c['linear_tolerance']*row['pressure_linear_rhs_norm']), 'linear true residual failed')
    if history:
        require(abs(history[-1]['pseudo_time']-summary['pseudo_time']) < 1e-12*max(1,summary['pseudo_time']), 'accepted time mismatch')
    require(summary['state_status'] == ('last-accepted-iterate' if history else 'initial-only'), 'initial state called an accepted state')
    # Read VTK counts, coordinates and all exported fields, independent of CSV.
    tokens = (root/'field.vtk').read_text().split()
    dims = tokens.index('DIMENSIONS')
    require(list(map(int,tokens[dims+1:dims+4])) == [nx+1,ny+1,1], 'VTK dimensions mismatch')
    for name,count,spacing in (('X_COORDINATES',nx+1,hx),('Y_COORDINATES',ny+1,hy),('Z_COORDINATES',1,0)):
        k=tokens.index(name);require(int(tokens[k+1])==count,'VTK coordinate count mismatch')
        coordinates=list(map(float,tokens[k+3:k+3+count]))
        require(all(math.isfinite(x) and abs(x-i*spacing)<=1e-12*max(g['length'],g['height']) for i,x in enumerate(coordinates)),'VTK coordinates mismatch')
    k=tokens.index('CELL_DATA');require(int(tokens[k+1])==nx*ny,'VTK cell count mismatch')
    k=tokens.index('VECTORS')+3
    values=list(map(float,tokens[k:k+3*nx*ny]))
    require(all(math.isfinite(x) for x in values),'nonfinite VTK vector')
    for index,row in enumerate(cells):
        require(values[3*index:3*index+3] == [row['u'],row['v'],0.], 'VTK vector mismatch')
    for name,column in (('pressure_fluctuation','pressure_fluctuation'),('divergence','divergence'),('classification','classification'),('solid_mask','solid_mask')):
        k=tokens.index(name)+5
        values=list(map(float,tokens[k:k+nx*ny]))
        require(values == [row[column] for row in cells],f'VTK scalar mismatch: {name}')
    loops, points = [], []
    for line in (root/'boundary.xy').read_text().splitlines()+['']:
        if line.strip():
            values=tuple(map(float,line.split()))
            require(len(values)==2 and all(math.isfinite(x) for x in values),'invalid boundary coordinate')
            points.append(values)
        elif points:
            require(len(points)>=3,'invalid boundary loop')
            loops.append(points)
            points=[]
    def inside(point, loop):
        x,y=point
        parity=False
        for a,b in zip(loop,loop[1:]+loop[:1]):
            if (a[1]>y)!=(b[1]>y) and x < a[0]+(y-a[1])*(b[0]-a[0])/(b[1]-a[1]):
                parity=not parity
        return parity
    def distance(point, a, b):
        dx,dy=b[0]-a[0],b[1]-a[1]
        require(dx*dx+dy*dy>0,'duplicate boundary edge')
        t=max(0,min(1,((point[0]-a[0])*dx+(point[1]-a[1])*dy)/(dx*dx+dy*dy)))
        return math.hypot(point[0]-a[0]-t*dx,point[1]-a[1]-t*dy)
    area=0.
    segments=[]
    segment_normals=[]
    for index,loop in enumerate(loops):
        depth=sum(inside(loop[0],other) for k,other in enumerate(loops) if k!=index)
        origin=loop[0]
        twice=sum((a[0]-origin[0])*(b[1]-origin[1])-(a[1]-origin[1])*(b[0]-origin[0]) for a,b in zip(loop,loop[1:]+loop[:1]))
        area+=(-1 if depth%2 else 1)*abs(twice)/2
        sign=(1 if twice>0 else -1)*(-1 if depth%2 else 1)
        for a,b in zip(loop,loop[1:]+loop[:1]):
            length=math.hypot(b[0]-a[0],b[1]-a[1])
            segments.append((a,b))
            segment_normals.append((sign*(b[1]-a[1])/length,sign*(a[0]-b[0])/length))
    require(abs(area-g['solid_area_geometry'])<=1e-11*max(area,hx*hy),'independent polygon area mismatch')
    require(bool(walls)==bool(segments),'wall samples missing or unexpected')
    expected_markers=[]
    if surface:
        for a,b in segments:
            length=math.hypot(b[0]-a[0],b[1]-a[1])
            if c.get('wall_quadrature','midpoint')=='midpoint':
                count=max(1,math.ceil(length/min(hx,hy)))
                for k in range(count):
                    fraction=(k+.5)/count
                    expected_markers.append((a[0]+fraction*(b[0]-a[0]),a[1]+fraction*(b[1]-a[1]),length/count))
                continue
            require(c['wall_quadrature']=='piecewise-gauss3','unknown wall quadrature')
            cuts=[0.,1.]
            for axis,spacing in ((0,hx/2),(1,hy/2)):
                lo,hi=sorted((a[axis],b[axis]))
                if lo==hi: continue
                for gridline in range(math.floor(lo/spacing)+1,math.ceil(hi/spacing)):
                    fraction=(gridline*spacing-a[axis])/(b[axis]-a[axis])
                    if 0<fraction<1: cuts.append(fraction)
            cuts=sorted(set(cuts))
            for left,right in zip(cuts,cuts[1:]):
                for node,weight in ((.5*(1-math.sqrt(.6)),5/18),(.5,4/9),(.5*(1+math.sqrt(.6)),5/18)):
                    fraction=left+(right-left)*node
                    expected_markers.append((a[0]+fraction*(b[0]-a[0]),a[1]+fraction*(b[1]-a[1]),length*(right-left)*weight))
    require(len(markers)==len(expected_markers),'surface quadrature count mismatch')
    for row,expected in zip(markers,expected_markers):
        require(all(abs(row[key]-value)<=1e-12*max(g['length'],g['height']) for key,value in zip(('x','y','segment_weight'),expected)), 'surface markers changed original geometry')
    for wall in walls:
        close=[normal for (a,b),normal in zip(segments,segment_normals) if distance((wall['x'],wall['y']),a,b)<=1e-11*max(g['length'],g['height'])]
        require(close,'wall sample is not on the original polygon')
        require(any(math.hypot(wall['nx']-normal[0],wall['ny']-normal[1])<=1e-11 for normal in close),'wall normal orientation mismatch')
    # Check every nonzero mask plus a distributed sample of the outside field.
    # The signed distance is reconstructed from the exported original segments.
    width=c['mask_half_width_cells']*min(hx,hy)
    for rows in (us,vs,cells):
        stride=max(1,len(rows)//100)
        for index,row in enumerate(rows):
            if row['solid_mask']==0 and index%stride:
                continue
            point=(row['x'],row['y'])
            if segments:
                d=min(distance(point,a,b) for a,b in segments)
                if sum(inside(point,loop) for loop in loops)%2:
                    d=-d
                ratio=max(-1,min(1,d/width))
                expected=.5*(1-ratio-math.sin(math.pi*ratio)/math.pi)
            else:
                expected=0.
            require(abs(expected-row['solid_mask'])<=1e-10,'signed-distance mask mismatch')
    audit.update(independent_equation_audit='pass',staggered_and_vtk_readback='pass',steps=summary['steps'],stop_reason=summary['stop_reason'])
    if not walls:
        audit['channel_relative_l2']=math.sqrt(error2/exact2)
        require(abs(audit['channel_relative_l2']-summary['metrics']['channel_relative_l2'])<1e-12,'channel error mismatch')
    return audit


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--require-converged',action='store_true')
    args=parser.parse_args()
    try:
        print(json.dumps(verify(args.directory,args.require_converged),indent=2))
    except (VerificationError,ValueError,KeyError,OSError) as exc:
        parser.exit(1,f'immersed verification failed: {exc}\n')


if __name__ == '__main__':
    main()
