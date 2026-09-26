"""Independent polygon-based Fourier operator; no native coefficients are read.

Least-squares reconstruction includes normal-derivative Neumann constraints.
The complete affine face Jacobian supplies a dimensional explicit-step bound.
"""
import math
from collections import defaultdict


def require(ok,message):
    if not ok:raise ValueError("heat reference: "+message)


class HeatReference:
    def __init__(self, mesh, measured, boundaries, conductivity):
        require(math.isfinite(conductivity) and conductivity>0,"invalid conductivity")
        self.mesh, self.measured, self.k = mesh, measured, conductivity
        self.faces = []
        for edge in mesh.edges:
            cell=mesh.cells[edge.owner];local=cell.edges.index(edge.id)
            a,b=(mesh.vertices[cell.vertices[i%len(cell.vertices)]] for i in (local,local+1))
            sx,sy=b[1]-a[1],a[0]-b[0];length=math.hypot(sx,sy)
            self.faces.append(dict(owner=edge.owner,neighbour=edge.neighbour,partner=-1,x=(a[0]+b[0])/2,y=(a[1]+b[1])/2,
                sx=sx,sy=sy,length=length,kind='insulated',value=0))
        for edge,f in zip(mesh.edges,self.faces):
            if edge.neighbour<0:
                b=boundaries[edge.id];f.update(kind=b.get('thermalKind','insulated'),value=b.get('thermalValue',0),partner=b['partner'])
                if f['partner']>=0:f['neighbour']=mesh.edges[f['partner']].owner
            cx,cy=measured.centroids[f['owner']]
            if f['neighbour']>=0:
                x,y=measured.centroids[f['neighbour']]
                if f['partner']>=0:
                    paired=self.faces[f['partner']];x+=f['x']-paired['x'];y+=f['y']-paired['y']
            else:x,y=f['x'],f['y']
            dx,dy=x-cx,y-cy;sd=f['sx']*dx+f['sy']*dy
            require(sd>0,"nonpositive normal distance")
            f['a']=f['length']**2/sd;f['cx'],f['cy']=f['sx']-f['a']*dx,f['sy']-f['a']*dy
            f['w']=(f['sx']*(f['x']-cx)+f['sy']*(f['y']-cy))/sd
        self.gradients=[]
        for i,cell in enumerate(mesh.cells):
            constraints=[];cx,cy=measured.centroids[i]
            for face_id in cell.edges:
                f=self.faces[face_id];edge=mesh.edges[face_id]
                j=edge.neighbour if edge.owner==i else edge.owner
                if edge.neighbour<0:j=f['neighbour']
                if j>=0:
                    x,y=measured.centroids[j]
                    if f['partner']>=0:
                        p=self.faces[f['partner']];x+=f['x']-p['x'];y+=f['y']-p['y']
                    dx,dy=x-cx,y-cy;kind='cell';value=0
                elif f['kind']=='temperature':dx,dy=f['x']-cx,f['y']-cy;kind='temperature';value=f['value']
                else:dx,dy=f['sx']/f['length'],f['sy']/f['length'];kind='derivative';value=-f['value']/conductivity
                r2=dx*dx+dy*dy;require(r2>0,"zero gradient distance")
                constraints.append((dx,dy,1/r2,kind,j,value))
            xx=math.fsum(dx*dx*w for dx,dy,w,*_ in constraints);xy=math.fsum(dx*dy*w for dx,dy,w,*_ in constraints);yy=math.fsum(dy*dy*w for dx,dy,w,*_ in constraints)
            det=xx*yy-xy*xy;require(det>64*math.ulp(1.)*(xx+yy)**2,"rank-deficient gradient")
            self.gradients.append([(w*(yy*dx-xy*dy)/det,w*(xx*dy-xy*dx)/det,kind,j,value) for dx,dy,w,kind,j,value in constraints])
        self.coefficients=[];self.constants=[];rows=[defaultdict(float) for _ in mesh.cells]
        for f in self.faces:
            c=defaultdict(float);constant=0.;i,j=f['owner'],f['neighbour']
            if j>=0 or f['kind']=='temperature':
                c[i]+=conductivity*f['a']
                if j>=0:c[j]-=conductivity*f['a']
                else:constant-=conductivity*f['a']*f['value']
                for owner,weight in [(i,1-f['w'] if j>=0 else 1)]+([(j,f['w'])] if j>=0 else []):
                    for wx,wy,kind,other,value in self.gradients[owner]:
                        coefficient=-conductivity*weight*(wx*f['cx']+wy*f['cy'])
                        if kind=='derivative':constant+=coefficient*value
                        else:
                            c[owner]-=coefficient
                            if other>=0:c[other]+=coefficient
                            else:constant+=coefficient*value
            elif f['kind']=='flux':constant=f['value']*f['length']
            self.coefficients.append(c);self.constants.append(constant)
        for face_id,(f,c) in enumerate(zip(self.faces,self.coefficients)):
            if f['partner']>=0 and face_id>f['partner']:continue
            for j,value in c.items():
                rows[f['owner']][j]+=value
                if f['neighbour']>=0:rows[f['neighbour']][j]-=value
        self.row_norm=[.5*math.fsum(abs(v) for v in row.values()) for row in rows]
        self.non_monotone=sum(any(j!=i and value>128*math.ulp(1.)*self.row_norm[i] for j,value in row.items()) for i,row in enumerate(rows))

    def evaluate(self, temperature, capacity):
        require(len(temperature)==len(self.mesh.cells)==len(capacity),"field size")
        require(all(math.isfinite(t) and t>0 for t in [*temperature,*capacity]),"nonpositive/nonfinite temperature or capacity")
        gradients=[]
        for i,samples in enumerate(self.gradients):
            delta=[value if kind=='derivative' else (temperature[j] if j>=0 else value)-temperature[i] for wx,wy,kind,j,value in samples]
            gradients.append((math.fsum(s[0]*d for s,d in zip(samples,delta)),math.fsum(s[1]*d for s,d in zip(samples,delta))))
        flux=[]
        for f in self.faces:
            i,j=f['owner'],f['neighbour'];gx,gy=gradients[i]
            if j>=0 or f['kind']=='temperature':
                if j>=0:gx,gy=((1-f['w'])*a+f['w']*b for a,b in zip((gx,gy),gradients[j]))
                value=-self.k*(f['a']*((temperature[j] if j>=0 else f['value'])-temperature[i])+gx*f['cx']+gy*f['cy'])
            else:value=f['value']*f['length'] if f['kind']=='flux' else 0.
            flux.append(value)
        # Pair once, using the same outward-owner convention as the topology.
        for face_id,f in enumerate(self.faces):
            if f['partner']>face_id:flux[f['partner']]=-flux[face_id]
        scales=[max(temperature)*math.fsum(abs(v) for v in c.values())+abs(a) for c,a in zip(self.coefficients,self.constants)]
        rates=[r/(a*c) for r,a,c in zip(self.row_norm,self.measured.areas,capacity)]
        return flux,rates,scales


def linear_euler_fourier_mode(time, conductivity, gamma=1.4, gas_r=1., rho=1., temperature=1., wavelength=1., amplitude=1e-5):
    """Continuum linearized Euler-Fourier mode, independent of spatial stencils.

    drho=A*cos(wx), u=B*sin(wx), dT=C*cos(wx). Initially isobaric
    to first order: A=-rho*amplitude, B=0, C=T*amplitude.
    Evaluate exp(time*M) by scaling/squaring of its convergent Taylor series.
    The nonlinear PDE differs by O(amplitude**2); this is a small-signal test.
    """
    w=2*math.pi/wavelength
    alpha=conductivity/(rho*gas_r/(gamma-1))
    m=[[0.,-rho*w,0.], [gas_r*temperature*w/rho,0.,gas_r*w],
       [0.,-(gamma-1)*temperature*w,-alpha*w*w]]
    norm=max(sum(abs(x*time) for x in row) for row in m)
    squarings=max(0,math.ceil(math.log2(norm/.5))) if norm else 0
    h=time/(2**squarings)
    a=[[x*h for x in row] for row in m]
    identity=lambda:[[float(i==j) for j in range(3)] for i in range(3)]
    def multiply(left,right):
        return [[math.fsum(left[i][k]*right[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    result=identity();term=identity()
    for n in range(1,64):
        term=[[x/n for x in row] for row in multiply(term,a)]
        result=[[x+y for x,y in zip(row,delta)] for row,delta in zip(result,term)]
        if max(abs(x) for row in term for x in row)<math.ulp(1.)/8:break
    else:raise ValueError('matrix exponential did not converge')
    for _ in range(squarings):result=multiply(result,result)
    initial=(-rho*amplitude,0.,temperature*amplitude)
    return tuple(math.fsum(x*y for x,y in zip(row,initial)) for row in result)
