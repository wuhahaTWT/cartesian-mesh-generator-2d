#!/usr/bin/env python3
"""Prepare analytic Dirichlet TRACE verification on an existing polygon mesh.

Different boundary model from the original constant-speed facet walls. The
known circular Couette field, extended to this polygon domain, solves the PDE.
Its boundary face means are integrated with Gauss-Legendre quadrature and its
normal flux uses the exact streamfunction difference. No FV equations here.
"""
from pathlib import Path
import argparse,csv,json,math,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'flow'))
import native_mesh as nm
import openfoam_data as foam
from annulus_cross import exact,dump
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--mesh',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--external-case',type=Path)
args=parser.parse_args()
if args.output.exists():parser.error('output exists')
m=nm.read_cm2d(args.mesh);g=nm.measure(m);fg=nm.face_geometry(m,g)
nodes,weights=np.polynomial.legendre.leggauss(16)
psi=lambda p:(p[0]*p[0]+p[1]*p[1])/6-math.log(math.hypot(*p))/3
records=[];flux_sum=0.;max_adjust=0.;rounded=0
for e,f in zip(m.edges,fg):
    if e.neighbour>=0:continue
    # Outward face vector (dy,-dx) determines the owner-oriented endpoints.
    sx,sy=f.area_vector;x,y=f.centre
    a=(x+.5*sy,y-.5*sx);b=(x-.5*sy,y+.5*sx)
    u=v=0.
    for t,w in zip(nodes,weights):
        p=(a[0]+(t+1)/2*(b[0]-a[0]),a[1]+(t+1)/2*(b[1]-a[1]))
        eu,ev,_=exact(*p);u+=w*eu/2;v+=w*ev/2
    q=psi(b)-psi(a)
    # Flux units m²/s. Only cancellation/role-classification roundoff is
    # resolved: 128 eps of streamfunction or Uref*edge length; not a solver
    # threshold. The exact mathematical closed-boundary flux telescopes.
    budget=128*sys.float_info.epsilon*max(abs(psi(a)),abs(psi(b)),.5*math.hypot(sx,sy))
    if abs(q)<=budget:q=0.;rounded+=1
    adjust=(q-u*sx-v*sy)/(sx*sx+sy*sy);u+=adjust*sx;v+=adjust*sy
    max_adjust=max(max_adjust,abs(adjust)*math.hypot(sx,sy));flux_sum+=u*sx+v*sy
    if q<0:kind='velocity-inlet';name='analytic-inlet'
    elif q>0:kind='velocity-outlet';name='analytic-outlet'
    else:kind='smooth-moving-wall';name='analytic-tangent'
    records.append((e.id,e.owner,x,y,sx,sy,kind,name,u,v,0.))
with args.output.open('w') as out:
    out.write(f'CARTMESH2D_FLOW_BOUNDARIES 1\nCOUNTS {len(m.cells)} {len(m.edges)} {len(records)}\n')
    for row in records:
        e,o,x,y,sx,sy,kind,name,u,v,p=row
        out.write(f'BOUNDARY {e} {o} {x:.17g} {y:.17g} {sx:.17g} {sy:.17g} {kind} "{name}" {u:.17g} {v:.17g} {p}\n')
    out.write('END\n')
dump(args.output.with_suffix('.trace.json'),dict(model='analytic trace verification; not original facet-wall motion',
    meshSha256=nm.sha256_file(args.mesh),gaussPoints=16,roundedZeroFluxFaces=rounded,netPrescribedFlux=flux_sum,
    maxQuadratureVelocityAdjustment=max_adjust,source='u_theta=(1/r-r)/3, streamfunction=r²/6-ln(r)/3'))
if args.external_case:
    base=args.external_case/'constant/polyMesh';points=foam.read_points(base/'points')
    faces=foam.read_faces(base/'faces');owners=foam.read_labels(base/'owner');patches=foam.read_boundary(base/'boundary')
    depth=max(p[2] for p in points)-min(p[2] for p in points)
    key=lambda o,x,y,sx,sy:(o,round(x,12),round(y,12),round(sx,12),round(sy,12))
    values={key(o,x,y,sx,sy):(u,v,0.) for e,o,x,y,sx,sy,kind,name,u,v,p in records}
    assert len(values)==len(records)
    fields=[];mapped=0
    for patch in patches:
        if patch['type']=='empty':entry='type empty;'
        else:
            rows=[]
            for i in range(patch['startFace'],patch['startFace']+patch['nFaces']):
                face=faces[i];area=foam.face_area(face,points)
                x,y=[sum(points[v][axis] for v in face)/len(face) for axis in (0,1)]
                rows.append(values[key(owners[i],x,y,area[0]/depth,area[1]/depth)]);mapped+=1
            entry='type fixedValue; value nonuniform List<vector>\n'+str(len(rows))+'\n(\n'+''.join(
                '('+' '.join(format(v,'.17g') for v in row)+')\n' for row in rows)+');'
        fields.append(patch['name']+'\n{ '+entry+' }\n')
    assert mapped==len(records)
    (args.external_case/'0/U').write_text(foam.header('volVectorField','0','U')+
        'dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0 0 0);\nboundaryField\n{\n'+''.join(fields)+'}\n')
    dump(args.external_case/'analytic-trace.json',dict(mappedFaces=mapped,boundarySha256=nm.sha256_file(args.output),
        model='analytic face-mean Dirichlet verification, not the original moving-wall model'))
print(args.output)
