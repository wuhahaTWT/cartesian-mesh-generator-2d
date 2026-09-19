#!/usr/bin/env python3
"""Independent wall-segment, affine-gradient and frozen SST equation audit.

Fixed probe definition: k=.02(1+.2x+.3y), omega=4(1+.1x+.2y), U=(y,0),
nu=1e-5, dt=.01, all scalar/velocity boundary values prescribed analytically.
Only y=0 boundary segments are selected for the distance diagnostic. This is
NOT a physical wall-bounded turbulence case or a converged nonlinear RANS step.
"""
import csv
import json
import math
from pathlib import Path
import verify_native_flow as native
from verify_sst_decay import require, finite, sha


def k_value(x,y): return .02*(1+.2*x+.3*y)
def w_value(x,y): return 4*(1+.1*x+.2*y)


def point_segment(p,a,b):
    x,y=b[0]-a[0],b[1]-a[1]
    t=max(0.,min(1.,((p[0]-a[0])*x+(p[1]-a[1])*y)/(x*x+y*y)))
    return math.hypot(p[0]-a[0]-t*x,p[1]-a[1]-t*y)


def coefficient(k,w,d):
    # Independently use published SST-2003m equations with exact affine
    # gradients (.004,.006),(.4,.8), S=1, and the geometric nearest wall.
    nu=1e-5;cross=2*.856/w*(.004*.4+.006*.8)
    f1=math.tanh(min(max(math.sqrt(k)/(.09*w*d),500*nu/(d*d*w)),
                       4*.856*k/(max(cross,1e-10)*d*d))**4)
    f2=math.tanh(max(2*math.sqrt(k)/(.09*w*d),500*nu/(d*d*w))**2)
    nt=.31*k/max(.31*w,f2)
    beta=f1*.075+(1-f1)*.0828;gamma=f1*(5/9)+(1-f1)*.44
    production=min(nt,10*.09*w*k)
    return dict(F1=f1,F2=f2,nuT=nt,Dk=nu+(f1*.85+1-f1)*nt,
                Dw=nu+(f1*.5+(1-f1)*.856)*nt,
                sourceK=production,sourceW=gamma*production/nt+(1-f1)*cross,
                lossK=.09*w,lossW=beta*w)


def audit(mesh_path,prefix):
    mesh=native.read_cm2d(mesh_path);m=native.measure(mesh,1e-10,1e-9)
    require(not m.issues,str(m.issues));geo=native.face_geometry(mesh,m)
    meta=json.loads(Path(str(prefix)+'.json').read_text())
    require(meta['scope']=='one-frozen-SST-transport-iteration' and meta['nu']==1e-5 and meta['dt']==.01,
            'changed probe definition')
    require(meta['kInnerConverged'] is True and meta['omegaInnerConverged'] is True,'unconverged inner equations')
    def rows(suffix):
        with Path(str(prefix)+suffix).open() as f:return list(csv.DictReader(f))
    cells,faces=rows('.cells.csv'),rows('.faces.csv')
    require(len(cells)==len(mesh.cells)==meta['cells'] and len(faces)==len(mesh.edges),'incomplete fields')
    walls=[e.id for e,g in zip(mesh.edges,geo) if e.neighbour<0 and abs(g.centre[1])<1e-12]
    require(walls,'missing diagnostic bottom wall')
    parameters=[];max_d_error=max_g_error=0.
    for i,(row,centre,area) in enumerate(zip(cells,m.centroids,m.areas)):
        require(int(row['cell'])==i,'cell order mismatch')
        for key,value in zip(('x','y','area'),(*centre,area)):
            require(math.isclose(finite(row[key]),value,abs_tol=1e-11,rel_tol=1e-9),'geometry mismatch')
        distances={fid:point_segment(centre,mesh.vertices[mesh.edges[fid].v0],mesh.vertices[mesh.edges[fid].v1]) for fid in walls}
        d=min(distances.values());d_error=abs(finite(row['distance'])-d);max_d_error=max(max_d_error,d_error)
        nearest=int(row['nearestFace'])
        require(nearest in walls and abs(distances[nearest]-d)<1e-12 and d_error<1e-12,'nearest finite wall segment mismatch')
        k,w=k_value(*centre),w_value(*centre)
        for key,value in [('oldK',k),('oldOmega',w),('gradKx',.004),('gradKy',.006),('gradWx',.4),('gradWy',.8),('strain',1.)]:
            error=abs(finite(row[key])-value);max_g_error=max(max_g_error,error)
            require(error<1e-10,'prescribed field/affine gradient or strain mismatch')
        co=coefficient(k,w,d);parameters.append(co)
        for key,value in co.items():
            require(math.isclose(finite(row[key]),value,rel_tol=1e-9,abs_tol=1e-12),'SST model coefficient mismatch: '+key)
    diagnostics={}
    for name,initial,dk,src,loss,advkey,diffkey in [('k',k_value,'Dk','sourceK','lossK','kAdvection','kDiffusion'),
                                                 ('omega',w_value,'Dw','sourceW','lossW','omegaAdvection','omegaDiffusion')]:
        values=[finite(row[name]) for row in cells]
        require(all(x>0 for x in values),'inadmissible turbulence values')
        boundary=[initial(*g.centre) for g in geo];fixed=[e.neighbour<0 for e in mesh.edges]
        gradients=native.reconstruct_gradient(mesh,m,geo,values,boundary,fixed)
        residual=[area*((v-initial(*p))/.01+c[loss]*v-c[src]) for area,v,p,c in zip(m.areas,values,m.centroids,parameters)]
        diagonal=[area*(100+c[loss]) for area,c in zip(m.areas,parameters)]
        base=[area*(100*initial(*p)+c[src]) for area,p,c in zip(m.areas,m.centroids,parameters)]
        max_flux_error=0.
        for i,(edge,g,row) in enumerate(zip(mesh.edges,geo,faces)):
            require(int(row['face'])==i and int(row['wall'])==int(i in walls),'face ordering/wall selection mismatch')
            owner,j=edge.owner,edge.neighbour;weight=g.neighbour_weight
            q=g.centre[1]*g.area_vector[0]
            require(abs(finite(row['volumeFlux'])-q)<1e-12,'carrier flux mismatch')
            D=parameters[owner][dk] if j<0 else (1-weight)*parameters[owner][dk]+weight*parameters[j][dk]
            gf=gradients[owner] if j<0 else tuple((1-weight)*a+weight*b for a,b in zip(gradients[owner],gradients[j]))
            diff=D*g.transmissibility*(values[owner]-(boundary[i] if j<0 else values[j]))-D*sum(a*b for a,b in zip(gf,g.correction))
            adv=q*(boundary[i] if j<0 and q<0 else values[j] if j>=0 and q<0 else values[owner])
            err=max(abs(finite(row[advkey])-adv),abs(finite(row[diffkey])-diff));max_flux_error=max(max_flux_error,err)
            require(err<1e-10,'shared face constitutive mismatch')
            flux=finite(row[advkey])+finite(row[diffkey]);residual[owner]+=flux
            diagonal[owner]+=D*g.transmissibility+max(q,0)
            if j>=0:residual[j]-=flux;diagonal[j]+=D*g.transmissibility+max(-q,0)
            else:base[owner]+=D*g.transmissibility*boundary[i]-min(q,0)*boundary[i]
        maximum=max(abs(r)/diag for r,diag in zip(residual,diagonal))
        require(maximum<2e-11,'frozen source/sink/time/flux cell balance failed')
        norm=math.sqrt(math.fsum(r*r for r in residual));stop=1e-13+1e-11*math.sqrt(math.fsum(b*b for b in base))
        require(norm<=2*stop,'independent frozen residual norm failed')
        diagnostics[name]=dict(maxDiagonalScaledResidual=maximum,residualNorm=norm,nativeNormTarget=stop,
                               globalBalance=math.fsum(residual),maxConstitutiveError=max_flux_error)
    require(0<meta['segmentTests']<=len(walls)*len(cells),'invalid spatial index work count')
    return dict(valid=True,scope=meta['scope'],mesh=str(mesh_path),meshSha256=sha(mesh_path),prefix=str(prefix),cells=len(cells),
                wallSegments=len(walls),segmentTests=meta['segmentTests'],bruteForceSegmentTests=len(walls)*len(cells),
                maxDistanceError=max_d_error,maxAffineGradientError=max_g_error,balance=diagnostics,
                files={s:sha(str(prefix)+s) for s in ('.cells.csv','.faces.csv','.json')})


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=Path,required=True)
    parser.add_argument('--prefix',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=audit(args.mesh.resolve(),args.prefix.resolve())
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
