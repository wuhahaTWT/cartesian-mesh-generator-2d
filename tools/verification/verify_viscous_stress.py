"""Independent planar Newtonian stress and work from the polygon file.

Sparse linear forms differentiate the complete velocity-to-momentum operator.
No native face gradients, coefficients, Jacobians or accepted flags are read.
"""
import math
from collections import defaultdict
from verify_wall_gradient import wall_stencils


class Linear:
    def __init__(self, terms=None): self.terms=dict(terms or {})
    def __add__(self, other):
        if not isinstance(other,Linear):other=Linear({-1:other})
        result=defaultdict(float,self.terms)
        for i,v in other.terms.items():result[i]+=v
        return Linear(result)
    __radd__=__add__
    def __rsub__(self, other):return -self+other
    def __mul__(self, scalar):return Linear({i:v*scalar for i,v in self.terms.items()})
    __rmul__=__mul__
    def __neg__(self):return self*-1
    def __sub__(self, other):return self+-other
    def evaluate(self,u):return math.fsum(v*(u[i] if i>=0 else 1.) for i,v in self.terms.items())
    def scale(self,u):return math.fsum(abs(v*(u[i] if i>=0 else 1.)) for i,v in self.terms.items())


def dot(a,b):return sum(x*y for x,y in zip(a,b))


class ViscousReference:
    def __init__(self,mesh,measured,bc,mu,wall_gradient="linear"):
        self.mesh,self.measured,self.mu=mesh,measured,mu
        self.faces=[]
        for e in mesh.edges:
            c=mesh.cells[e.owner];index=c.edges.index(e.id)
            a,b=(mesh.vertices[c.vertices[j%len(c.vertices)]] for j in (index,index+1))
            area=(b[1]-a[1],a[0]-b[0]);length=math.hypot(*area)
            self.faces.append(dict(owner=e.owner,neighbour=e.neighbour,partner=-1,
                point=tuple((x+y)/2 for x,y in zip(a,b)),area=area,normal=tuple(x/length for x in area),kind='internal',wallVelocity=(0,0)))
        for i,f in enumerate(self.faces):
            if i in bc:f.update({key:bc[i][key] for key in ('partner','kind','wallVelocity')})
            if f['partner']>=0:f['neighbour']=mesh.edges[f['partner']].owner
            p=measured.centroids[f['owner']];f['op']=tuple(x-y for x,y in zip(f['point'],p))
            f['np']=(0,0)
            if f['neighbour']>=0:
                q=measured.centroids[f['neighbour']];point=self.faces[f['partner']]['point'] if f['partner']>=0 else f['point']
                f['np']=tuple(x-y for x,y in zip(point,q))
            f['d']=tuple(x-y for x,y in zip(f['op'],f['np']));f['dn']=dot(f['d'],f['normal'])
            f['weight']=dot(f['op'],f['normal'])/f['dn'] if f['neighbour']>=0 else 0
        velocity=[[Linear({2*i+k:1}) for k in range(2)] for i in range(len(mesh.cells))]
        gradients=[]
        for i,c in enumerate(mesh.cells):
            constraints=[];p=measured.centroids[i]
            for face_id in c.edges:
                e=mesh.edges[face_id];f=self.faces[face_id]
                j=(e.neighbour if e.owner==i else e.owner) if e.neighbour>=0 else f['neighbour']
                if j>=0:
                    d=f['d'] if f['partner']>=0 else tuple(x-y for x,y in zip(measured.centroids[j],p))
                    delta=[x-y for x,y in zip(velocity[j],velocity[i])]
                elif f['kind']=='no-slip-wall':
                    d=f['op'];delta=[x-y for x,y in zip(f['wallVelocity'],velocity[i])]
                else:
                    dn=dot(f['op'],f['normal']);d=tuple(dn*x for x in f['normal'])
                    normal_speed=dot(velocity[i],f['normal'])
                    delta=[-normal_speed*x for x in f['normal']] if f['kind']=='slip-wall' else [Linear(),Linear()]
                r2=dot(d,d);constraints.append((d,1/r2,delta))
            xx=math.fsum(d[0]**2*w for d,w,_ in constraints);xy=math.fsum(d[0]*d[1]*w for d,w,_ in constraints);yy=math.fsum(d[1]**2*w for d,w,_ in constraints)
            determinant=xx*yy-xy*xy
            if determinant<=64*math.ulp(1.)*(xx+yy)**2:raise ValueError('singular viscous reference stencil')
            gradients.append([[sum((delta[k]*(w*(yy*d[0]-xy*d[1])/determinant) for d,w,delta in constraints),Linear()),
                               sum((delta[k]*(w*(xx*d[1]-xy*d[0])/determinant) for d,w,delta in constraints),Linear())] for k in range(2)])
        self.wall_stencils=wall_stencils(mesh,measured,bc,{i for i,b in bc.items() if b['kind']=='no-slip-wall'}) if wall_gradient=='quadratic' else {}
        self.momentum=[];self.face_velocity=[];rows=[[Linear(),Linear()] for _ in mesh.cells]
        for face_id,f in enumerate(self.faces):
            i,j=f['owner'],f['neighbour'];w=f['weight'];n=f['normal'];g=[[x for x in row] for row in gradients[i]]
            uf=[Linear({-1:x}) for x in f['wallVelocity']]
            if j>=0:
                g=[[(1-w)*gradients[i][k][a]+w*gradients[j][k][a] for a in range(2)] for k in range(2)]
                uf=[(1-w)*(velocity[i][k]+dot(gradients[i][k],f['op']))+w*(velocity[j][k]+dot(gradients[j][k],f['np'])) for k in range(2)]
            if j>=0 or f['kind']=='no-slip-wall':
                other=velocity[j] if j>=0 else uf
                g=[[g[k][a]+(other[k]-velocity[i][k]-dot(g[k],f['d']))*(n[a]/f['dn']) for a in range(2)] for k in range(2)]
            if face_id in self.wall_stencils:
                g=[[sum(((Linear({-1:self.faces[index]['wallVelocity'][k]}) if boundary else velocity[index][k])-uf[k])*w[a] for index,boundary,w in self.wall_stencils[face_id]) for a in range(2)] for k in range(2)]
            div=g[0][0]+g[1][1]
            stress=[[(2*g[0][0]-2/3*div)*mu,(g[0][1]+g[1][0])*mu],[(g[0][1]+g[1][0])*mu,(2*g[1][1]-2/3*div)*mu]]
            traction=[dot(row,f['area']) for row in stress]
            if j<0:
                if f['kind']=='slip-wall':
                    normal=dot(traction,n);traction=[normal*x for x in n];uf=[Linear(),Linear()]
                elif f['kind']!='no-slip-wall':traction=[Linear(),Linear()]
            momentum=[-x for x in traction];self.momentum.append(momentum);self.face_velocity.append(uf)
            if f['partner']<0 or face_id<f['partner']:
                for k in range(2):
                    rows[i][k]=rows[i][k]+momentum[k]
                    if j>=0:rows[j][k]=rows[j][k]-momentum[k]
        self.row_norm=[.5*max(math.fsum(abs(v) for j,v in row.terms.items() if j>=0) for row in pair) for pair in rows]

    def evaluate(self,velocity,density):
        u=[x for pair in velocity for x in pair];flux=[];scales=[]
        for momentum,uf in zip(self.momentum,self.face_velocity):
            m=[x.evaluate(u) for x in momentum];v=[x.evaluate(u) for x in uf]
            flux.append([*m,dot(m,v)])
            scale=[x.scale(u) for x in momentum];scales.append([*scale,sum(x*y.scale(u) for x,y in zip(scale,uf))])
        for i,f in enumerate(self.faces):
            if f['partner']>i:flux[f['partner']]=[-x for x in flux[i]];scales[f['partner']]=scales[i]
        rates=[r/(a*d) for r,a,d in zip(self.row_norm,self.measured.areas,density)]
        return flux,rates,scales
