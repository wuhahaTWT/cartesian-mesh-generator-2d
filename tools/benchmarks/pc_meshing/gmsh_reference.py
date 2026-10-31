"""Mature Gmsh meshes on the frozen straight-sided annulus, with no wall snapping.

Writes Gmsh 2.2 interchange plus polygons for the native ready-import consumer.
Mapped quads are a geometry-specific reference, not a general Cut-cell repair.
"""
from pathlib import Path
import argparse, hashlib, json, math, resource, sys, time

p=argparse.ArgumentParser()
p.add_argument('--geometry',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--mode',choices=['mapped','triangles'],default='mapped')
p.add_argument('--radial',type=int,default=32)
p.add_argument('--subdivisions',type=int,default=2)
p.add_argument('--h',type=float,default=.025)
a=p.parse_args();sys.argv[:]=[sys.argv[0]]
import gmsh
import numpy as np

loops=[[tuple(map(float,s.split())) for s in block.splitlines() if s.strip()]
       for block in a.geometry.read_text(encoding='utf-8').strip().split('\n\n')]
assert len(loops)==2 and len(loops[0])==len(loops[1])==256
# Match rays by angle only to build this annular reference; no point is moved.
loops=[sorted(loop,key=lambda p:math.atan2(p[1],p[0])) for loop in loops]
loops.sort(key=lambda loop:sum(math.hypot(*p) for p in loop))
started=time.monotonic();gmsh.initialize();gmsh.option.setNumber('General.Terminal',0)
gmsh.model.add('frozen_polygon_annulus');geo=gmsh.model.geo;n=len(loops[0])
pts=[[geo.addPoint(x,y,0,a.h) for x,y in loop] for loop in loops]
inner=[geo.addLine(pts[0][i],pts[0][(i+1)%n]) for i in range(n)]
outer=[geo.addLine(pts[1][i],pts[1][(i+1)%n]) for i in range(n)]
surfaces=[]
if a.mode=='mapped':
    rays=[geo.addLine(pts[0][i],pts[1][i]) for i in range(n)]
    for i in range(n):
        loop=geo.addCurveLoop([rays[i],outer[i],-rays[(i+1)%n],-inner[i]])
        surfaces.append(geo.addPlaneSurface([loop]))
    geo.synchronize()
    for c in inner+outer:gmsh.model.mesh.setTransfiniteCurve(c,a.subdivisions+1)
    for c in rays:gmsh.model.mesh.setTransfiniteCurve(c,a.radial+1)
    for s in surfaces:
        gmsh.model.mesh.setTransfiniteSurface(s);gmsh.model.mesh.setRecombine(2,s)
else:
    surfaces=[geo.addPlaneSurface([geo.addCurveLoop(outer),geo.addCurveLoop(inner)])]
    geo.synchronize();gmsh.option.setNumber('Mesh.Algorithm',6)
    gmsh.option.setNumber('Mesh.MeshSizeMin',a.h);gmsh.option.setNumber('Mesh.MeshSizeMax',a.h)
gmsh.model.addPhysicalGroup(1,inner,1);gmsh.model.setPhysicalName(1,1,'rotor')
gmsh.model.addPhysicalGroup(1,outer,2);gmsh.model.setPhysicalName(1,2,'housing')
gmsh.model.addPhysicalGroup(2,surfaces,3);gmsh.model.setPhysicalName(2,3,'fluid')
gmsh.model.mesh.generate(2)
if a.mode=='triangles':gmsh.model.mesh.optimize('Relocate2D',niter=3)
mesh_seconds=time.monotonic()-started
tags,coords,_=gmsh.model.mesh.getNodes();coords=np.asarray(coords).reshape(-1,3)
positions={int(tag):p[:2] for tag,p in zip(tags,coords)}
types,element_tags,element_nodes=gmsh.model.mesh.getElements(2)
polygons=[]
for kind,ids,nodes in zip(types,element_tags,element_nodes):
    width=gmsh.model.mesh.getElementProperties(kind)[3]
    assert width in (3,4)
    for row in np.asarray(nodes).reshape(-1,width):polygons.append(np.array([positions[int(i)] for i in row]))
a.output.parent.mkdir(parents=True,exist_ok=True)
with a.output.with_suffix('.polygons').open('x') as f:
    f.write(str(len(polygons))+'\n')
    for i,poly in enumerate(polygons):
        lo=poly.min(axis=0);hi=poly.max(axis=0)
        values=[i,0,*lo,*hi,len(poly),*poly.ravel()]
        f.write(' '.join(format(float(v),'.17g') if not isinstance(v,int) else str(v) for v in values)+'\n')
gmsh.option.setNumber('Mesh.MshFileVersion',2.2);gmsh.write(str(a.output.with_suffix('.msh')))
record=dict(engine='Gmsh',version=gmsh.__version__,mode=a.mode,cells=len(polygons),
    radial=a.radial,subdivisionsPerOriginalFacet=a.subdivisions,h=a.h,
    meshSeconds=mesh_seconds,totalSeconds=time.monotonic()-started,
    peakRssMiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
    cpuSeconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,
    geometrySha256=hashlib.sha256(a.geometry.read_bytes()).hexdigest(),
    boundary='Original 256 straight segments per wall retained; mapped interior uses radial interpolation.',
    source='https://gmsh.info/doc/texinfo/gmsh.html')
a.output.with_suffix('.gmsh.json').write_text(json.dumps(record,indent=2)+'\n')
gmsh.finalize();print(json.dumps(record))
