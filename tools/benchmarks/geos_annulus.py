#!/usr/bin/env python3
"""Independent mature GEOS polygon-overlay reference, not a packaged mesher.

Cut the SAME original annulus (256 facets on each wall) with exact 2D polygon intersections.
No cell merging/quality repair occurs here; those shared native stages are
reported separately by annulus_mesh_cross.cpp. No FV equations are implemented.
"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
import shapely
from shapely.geometry import Polygon, box

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--geometry',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--n',type=int,default=128)
parser.add_argument('--base-mesh',type=Path)
parser.add_argument('--protect-radius',type=float)
args=parser.parse_args()
if args.output.exists():parser.error('output exists')
loops=[]
for block in args.geometry.read_text().strip().split('\n\n'):
    loops.append([tuple(map(float,line.split())) for line in block.splitlines() if line.strip()])
assert len(loops)==2
region=Polygon(loops[0],[loops[1]])
assert region.is_valid
cells=[];h=2.6/args.n;canonicalized=0;maximum_shift=0.
protected=0
if args.base_mesh:
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'flow'))
    import native_mesh as nm
    mesh=nm.read_cm2d(args.base_mesh)
    def sources():
        for c in mesh.cells:
            vertices=[mesh.vertices[v] for v in c.vertices]
            yield c.id,Polygon(vertices),vertices
else:
    def sources():
        for j in range(args.n):
          for i in range(args.n):
            xmin=-1.3+i*h;ymin=-1.3+j*h;xmax=-1.3+(i+1)*h;ymax=-1.3+(j+1)*h
            yield j*args.n+i,box(xmin,ymin,xmax,ymax),None
for source,square,original_vertices in sources():
        xmin,ymin,xmax,ymax=square.bounds
        is_protected=bool(original_vertices and args.protect_radius and
                         max(math.hypot(x,y) for x,y in original_vertices)<args.protect_radius)
        if is_protected:
            assert region.covers(square)
            cells.append((source,-1,xmin,ymin,xmax,ymax,original_vertices));protected+=1;continue
        if original_vertices and region.covers(square):
            cells.append((source,-2,xmin,ymin,xmax,ymax,original_vertices));continue
        part=region.intersection(square)
        if part.is_empty:continue
        parts=[part] if part.geom_type=='Polygon' else list(part.geoms)
        for p in parts:
            if p.geom_type!='Polygon':
                # Point/line contacts have zero fluid area and are not cells.
                assert p.area==0;continue
            assert not p.interiors and p.is_valid and p.area>0
            vertices=[]
            # Canonical Cartesian-corner identities: the legacy trig cardinal
            # points differ from zero by ~1e-16. Resolve ONLY construction
            # roundoff (64 ulps of the unit coordinate scale), not geometry
            # tolerance or a minimum usable edge. Never smooth the wall.
            budget=64*math.ulp(max(1.,abs(xmin),abs(ymin),abs(xmax),abs(ymax)))
            for x,y in list(p.exterior.coords)[:-1]:
                a,b=x,y
                for value in (xmin,xmax):
                    if abs(a-value)<=budget:a=value
                for value in (ymin,ymax):
                    if abs(b-value)<=budget:b=value
                maximum_shift=max(maximum_shift,math.hypot(a-x,b-y))
                canonicalized+=int((a,b)!=(x,y))
                if not vertices or vertices[-1]!=(a,b):vertices.append((a,b))
            if len(vertices)>1 and vertices[0]==vertices[-1]:vertices.pop()
            assert len(vertices)>=3
            cells.append((source,int(not region.covers(square)),xmin,ymin,xmax,ymax,vertices))
with args.output.open('w') as f:
    f.write(str(len(cells))+'\n')
    for source,kind,xmin,ymin,xmax,ymax,vertices in cells:
        values=[source,kind,xmin,ymin,xmax,ymax,len(vertices)]+[v for p in vertices for v in p]
        f.write(' '.join(str(v) if isinstance(v,int) else format(v,'.17g') for v in values)+'\n')
metadata=dict(engine='GEOS polygon intersection through Shapely',shapelyVersion=shapely.__version__,geosVersion=shapely.geos_version_string,
    geometrySha256=hashlib.sha256(args.geometry.read_bytes()).hexdigest(),n=args.n,cells=len(cells),area=region.area,
    canonicalizedVertices=canonicalized,maximumRoundoffDisplacement=maximum_shift,
    baseMeshSha256=hashlib.sha256(args.base_mesh.read_bytes()).hexdigest() if args.base_mesh else None,
    protectedBaseCells=protected,protectRadius=args.protect_radius,
    note='Independent 2D cutting kernel reference; native topology, small-cell and solver-quality stages remain shared. This is not AMReX or a third-party end-to-end mesh-generator qualification.')
args.output.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')
print(metadata)
