"""Independent point-polynomial wall recovery from original polygon geometry.

Uses reorthogonalized Gram-Schmidt, not the native Householder implementation.
Only geometry and prescribed boundary types enter preparation. No native
coefficients, gradients, or reported wall counts are used as reference data.
"""
import math


def dot(a,b):return math.fsum(x*y for x,y in zip(a,b))


def inverse_gradient(points,normal):
    tangent=(-normal[1],normal[0]);hn=max(abs(dot(p[2],normal)) for p in points);ht=max(abs(dot(p[2],tangent)) for p in points)
    if len(points)<5 or hn==0 or ht==0:return None
    a=[];weights=[]
    for _,_,offset in points:
        x,y=dot(offset,normal)/hn,dot(offset,tangent)/ht;w=1/math.hypot(x,y)
        weights.append(w);a.append([w*v for v in (x,y,x*x/2,x*y,y*y/2)])
    norms=[math.sqrt(math.fsum(row[j]**2 for row in a)) for j in range(5)]
    if min(norms)<=0:return None
    columns=[[row[j]/norms[j] for row in a] for j in range(5)];permutation=list(range(5));q=[];r=[[0.]*5 for _ in range(5)]
    for k in range(5):
        pivot=max(range(k,5),key=lambda j:dot(columns[j],columns[j]));columns[k],columns[pivot]=columns[pivot],columns[k];permutation[k],permutation[pivot]=permutation[pivot],permutation[k]
        for previous in range(k):r[previous][k],r[previous][pivot]=r[previous][pivot],r[previous][k]
        size=math.sqrt(dot(columns[k],columns[k]))
        if size<=math.sqrt(math.ulp(1.)):return None
        q.append([v/size for v in columns[k]]);r[k][k]=size
        for j in range(k+1,5):
            for _ in range(2):
                projection=dot(q[k],columns[j]);r[k][j]+=projection
                columns[j]=[x-projection*y for x,y in zip(columns[j],q[k])]
    result=[]
    for i,p in enumerate(points):
        x=[0.]*5;coefficient=[0.]*5
        for j in reversed(range(5)):x[j]=(q[j][i]-math.fsum(r[j][k]*x[k] for k in range(j+1,5)))/r[j][j]
        for j in range(5):coefficient[permutation[j]]=x[j]*weights[i]/norms[permutation[j]]
        result.append((p[0],p[1],tuple(coefficient[0]*n/hn+coefficient[1]*t/ht for n,t in zip(normal,tangent))))
    return result


def wall_stencils(mesh,measured,bc,prescribed):
    centres=[];normals=[];lengths=[]
    for edge in mesh.edges:
        cell=mesh.cells[edge.owner];local=cell.edges.index(edge.id)
        a,b=(mesh.vertices[cell.vertices[i%len(cell.vertices)]] for i in (local,local+1))
        centres.append(tuple((x+y)/2 for x,y in zip(a,b)));s=(b[1]-a[1],a[0]-b[0]);lengths.append(math.hypot(*s));normals.append(tuple(x/lengths[-1] for x in s))
    stencils={}
    for target in sorted(prescribed):
        owner=mesh.edges[target].owner;centre=measured.centroids[owner];face=centres[target]
        image_cells=[(owner,tuple(x-y for x,y in zip(centre,face)))];points=[];begin,end=0,1
        scale=max([lengths[target]]+[math.dist(centres[f],centre) for f in mesh.cells[owner].edges]);eps=128*math.ulp(1.)*scale
        def add_point(index,boundary,offset):
            if math.hypot(*offset)<=eps:return
            if any(index==i and boundary==b and math.dist(offset,d)<=eps for i,b,d in points):return
            points.append((index,boundary,offset))
        for ring in range(6):
            for cell_id,offset in image_cells[begin:end]:
                cell=mesh.cells[cell_id];c=measured.centroids[cell_id];add_point(cell_id,False,offset)
                for f in cell.edges:
                    if f in prescribed:add_point(f,True,tuple(d+x-y for d,x,y in zip(offset,centres[f],c)))
            fitted=inverse_gradient(points,normals[target]) if ring>=2 else None
            if fitted is not None:break
            if ring==5:raise ValueError('quadratic reference rank failure')
            for cell_id,offset in image_cells[begin:end]:
                cell=mesh.cells[cell_id];c=measured.centroids[cell_id]
                for f in cell.edges:
                    edge=mesh.edges[f];other=edge.neighbour if edge.owner==cell_id else edge.owner
                    if edge.neighbour>=0:delta=tuple(x-y for x,y in zip(measured.centroids[other],c))
                    elif bc[f]['partner']>=0:
                        partner=bc[f]['partner'];other=mesh.edges[partner].owner
                        delta=tuple((x-y)+(z-w) for x,y,z,w in zip(centres[f],c,measured.centroids[other],centres[partner]))
                    else:continue
                    candidate=(other,tuple(x+y for x,y in zip(offset,delta)))
                    if not any(other==i and math.dist(candidate[1],d)<=eps for i,d in image_cells):image_cells.append(candidate)
                    if len(image_cells)>512:raise ValueError('excessive quadratic reference stencil')
            begin,end=end,len(image_cells)
            if begin==end:
                fitted=inverse_gradient(points,normals[target])
                if fitted is None:raise ValueError('quadratic reference rank failure')
                break
        stencils[target]=fitted
    return stencils
