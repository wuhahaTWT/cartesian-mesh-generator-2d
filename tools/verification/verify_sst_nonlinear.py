#!/usr/bin/env python3
"""Audit nonlinear SST at RETURNED fields on prescribed shear carrier.

Fixed nu=1e-5, dt=.02, old k=.001, omega=2, U=(y,0); resolved bottom
wall k=0 and omega=60nu/(.075*owner-normal-spacing^2). Outflow is zero
normal scalar gradient, other nonwall boundaries prescribe old freestream.
This is not a momentum-coupled RANS or physical boundary-layer validation.
"""
import argparse,csv,json,math
from pathlib import Path
import verify_native_flow as native
from verify_sst_decay import require,finite,sha
from verify_sst_spatial import point_segment


def coefficients(k,w,d,gk,gw):
    cross0=2*.856/w*sum(a*b for a,b in zip(gk,gw));nu=1e-5
    arg1=min(max(math.sqrt(k)/(.09*w*d),500*nu/(d*d*w)),4*.856*k/(max(cross0,1e-10)*d*d))
    f1=math.tanh(arg1**4);f2=math.tanh(max(2*math.sqrt(k)/(.09*w*d),500*nu/(d*d*w))**2)
    nt=.31*k/max(.31*w,f2);P=min(nt,10*.09*w*k)
    beta=f1*.075+(1-f1)*.0828;gamma=f1*5/9+(1-f1)*.44
    cross=(1-f1)*cross0
    return dict(F1=f1,F2=f2,nuT=nt,Dk=nu+(f1*.85+1-f1)*nt,Dw=nu+(f1*.5+(1-f1)*.856)*nt,
                sourceK=P,sourceW=gamma*P/nt+max(cross,0),lossK=.09*w,lossW=beta*w+max(-cross,0)/w), (P-.09*w*k,gamma*P/nt-beta*w*w+cross)


def audit(mesh_path,prefix):
    mesh=native.read_cm2d(mesh_path);m=native.measure(mesh,1e-10,1e-9)
    require(not m.issues,str(m.issues));geo=native.face_geometry(mesh,m)
    meta=json.loads(Path(str(prefix)+'.json').read_text())
    require(meta['scope']=='nonlinear-SST-fixed-carrier' and meta['converged'] is True and
            meta['nu']==1e-5 and meta['dt']==.02 and meta['initialK']==.001 and meta['initialOmega']==2,'changed nonlinear probe')
    def rows(suffix):
        with Path(str(prefix)+suffix).open() as f:return list(csv.DictReader(f))
    cells,faces,history=rows('.cells.csv'),rows('.faces.csv'),rows('.history.csv')
    require(len(cells)==len(mesh.cells)==meta['cells'] and len(faces)==len(mesh.edges),'incomplete fields')
    require(len(history)==meta['iterations']+1 and 0<meta['iterations']<=500,'incomplete nonlinear history')
    for i,row in enumerate(history):
        require(int(row['iteration'])==i,'history iteration mismatch')
        require(all(finite(row[key])>=0 for key in ('kNorm','omegaNorm','kCellResidual','omegaCellResidual')),'invalid history residual')
    k=[finite(r['k']) for r in cells];w=[finite(r['omega']) for r in cells]
    require(all(v>0 for v in k+w),'inadmissible final fields')
    walls=[e.id for e,g in zip(mesh.edges,geo) if e.neighbour<0 and abs(g.centre[1])<1e-12]
    require(walls,'missing wall')
    kb=[.001]*len(faces);wb=[2.]*len(faces);fixed=[];q=[]
    for i,(e,g,row) in enumerate(zip(mesh.edges,geo,faces)):
        require(int(row['face'])==i and int(row['wall'])==int(i in walls),'face order/wall mismatch')
        q.append(0. if i in walls else g.centre[1]*g.area_vector[0]);fixed.append(e.neighbour<0 and q[-1]<=0)
        if i in walls:
            centre=m.centroids[e.owner];normalDistance=sum((a-b)*s for a,b,s in zip(g.centre,centre,g.area_vector))/math.hypot(*g.area_vector)
            kb[i]=0.;wb[i]=60*1e-5/(.075*normalDistance**2)
        require(abs(finite(row['volumeFlux'])-q[i])<1e-12,'carrier mismatch')
        if e.neighbour<0:
            for key,v in [('kBoundary',kb[i] if fixed[i] else 0),('omegaBoundary',wb[i] if fixed[i] else 0)]:
                require(math.isclose(finite(row[key]),v,rel_tol=1e-11,abs_tol=1e-13),'resolved wall/inlet/outlet value mismatch')
    kg=native.reconstruct_gradient(mesh,m,geo,k,kb,fixed);wg=native.reconstruct_gradient(mesh,m,geo,w,wb,fixed)
    co=[];sources=[];max_co_error=0.;max_gradient_error=0.;max_roundoff_allowance=0.
    for i,(row,centre,area) in enumerate(zip(cells,m.centroids,m.areas)):
        require(int(row['cell'])==i,'cell order mismatch')
        for key,value in zip(('x','y','area'),(*centre,area)):
            require(math.isclose(finite(row[key]),value,rel_tol=1e-9,abs_tol=1e-11),'geometry mismatch')
        d=min(point_segment(centre,mesh.vertices[mesh.edges[f].v0],mesh.vertices[mesh.edges[f].v1]) for f in walls)
        require(abs(finite(row['distance'])-d)<1e-12,'wall distance mismatch')
        for keys,expected in [(('gradKx','gradKy'),kg[i]),(('gradWx','gradWy'),wg[i])]:
            actual=tuple(finite(row[key]) for key in keys)
            # Native long-double polygon moments and independent compensated
            # double moments can differ by an ulp. Large wall-normal omega
            # gradients amplify that into the nearly zero tangential component.
            # Budget only roundoff relative to the VECTOR, retaining the original
            # per-component relative/absolute checks and all PDE stopping gates.
            roundoff=128*math.ulp(1.)*max(1.,math.hypot(*actual),math.hypot(*expected))
            max_roundoff_allowance=max(max_roundoff_allowance,roundoff)
            for key,a,b in zip(keys,actual,expected):
                error=abs(a-b);max_gradient_error=max(max_gradient_error,error)
                require(error<=max(1e-11,1e-8*max(abs(a),abs(b)))+roundoff,
                        f'current-field gradient mismatch: cell {i} {key}: exported={a:.17g}, independent={b:.17g}')
        require(math.isclose(finite(row['strain']),1.,rel_tol=1e-8,abs_tol=1e-11),'strain mismatch')
        c,source=coefficients(k[i],w[i],d,kg[i],wg[i]);co.append(c);sources.append(source)
        for key,value in c.items():
            error=abs(finite(row[key])-value);max_co_error=max(max_co_error,error)
            require(math.isclose(finite(row[key]),value,rel_tol=1e-8,abs_tol=1e-11),'stale or incorrect closure coefficient: '+key)
    diagnostics={}
    for index,name,values,grad,bc,old,Dkey,Skey,Lkey in [(0,'k',k,kg,kb,.001,'Dk','sourceK','lossK'),(1,'omega',w,wg,wb,2.,'Dw','sourceW','lossW')]:
        residual=[area*((value-old)/.02-source[index]) for area,value,source in zip(m.areas,values,sources)]
        diag=[area*(50+c[Lkey]) for area,c in zip(m.areas,co)]
        base=[area*(50*old+c[Skey]) for area,c in zip(m.areas,co)]
        max_flux_error=0.
        for i,(e,g,row) in enumerate(zip(mesh.edges,geo,faces)):
            a,b=e.owner,e.neighbour;weight=g.neighbour_weight
            D=1e-5 if i in walls else co[a][Dkey] if b<0 else (1-weight)*co[a][Dkey]+weight*co[b][Dkey]
            gf=grad[a] if b<0 else tuple((1-weight)*x+weight*y for x,y in zip(grad[a],grad[b]))
            adv=q[i]*(bc[i] if b<0 and q[i]<0 else values[b] if b>=0 and q[i]<0 else values[a])
            diff=0 if b<0 and not fixed[i] else D*g.transmissibility*(values[a]-(bc[i] if b<0 else values[b]))-D*sum(x*y for x,y in zip(gf,g.correction))
            error=max(abs(finite(row[name+'Advection'])-adv),abs(finite(row[name+'Diffusion'])-diff));max_flux_error=max(max_flux_error,error)
            require(error<1e-10,'current-coefficient constitutive face flux mismatch')
            flux=finite(row[name+'Advection'])+finite(row[name+'Diffusion']);residual[a]+=flux
            diag[a]+=max(q[i],0)+(D*g.transmissibility if b>=0 or fixed[i] else 0)
            if b>=0:residual[b]-=flux;diag[b]+=max(-q[i],0)+D*g.transmissibility
            else:
                base[a]-=min(q[i],0)*bc[i]
                if fixed[i]:base[a]+=D*g.transmissibility*bc[i]
        maximum=max(abs(r)/di for r,di in zip(residual,diag));norm=math.sqrt(math.fsum(r*r for r in residual))
        target=1e-13+1e-10*math.sqrt(math.fsum(x*x for x in base))
        require(maximum<=2e-10 and norm<=2*target,'original nonlinear equation residual failed')
        require(math.isclose(norm,finite(history[-1][name+'Norm']),abs_tol=1e-11,rel_tol=1e-5),'history does not describe returned fields')
        require(math.isclose(maximum,finite(history[-1][name+'CellResidual']),abs_tol=1e-13,rel_tol=1e-5),'scaled history does not describe returned fields')
        diagnostics[name]=dict(maxCellResidual=maximum,residualNorm=norm,nativeNormTarget=target,globalBalance=math.fsum(residual),maxFaceFluxError=max_flux_error)
    return dict(valid=True,scope=meta['scope'],cells=len(cells),iterations=meta['iterations'],mesh=str(mesh_path),meshSha256=sha(mesh_path),prefix=str(prefix),
                maxCoefficientError=max_co_error,maxGradientError=max_gradient_error,maxGradientRoundoffAllowance=max_roundoff_allowance,
                balance=diagnostics,files={suffix:sha(str(prefix)+suffix) for suffix in ('.cells.csv','.faces.csv','.json','.history.csv')})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mesh',type=Path,required=True);p.add_argument('--prefix',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=audit(a.mesh.resolve(),a.prefix.resolve());a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
