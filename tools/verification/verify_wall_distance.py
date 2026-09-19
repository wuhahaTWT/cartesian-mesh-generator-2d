#!/usr/bin/env python3
"""Independently audit ALL cells in an axis-aligned rectangular wall-distance probe.

The native algorithm handles arbitrary selected finite segments. This scale
probe selects all outer rectangle edges; this verifier rejects other domains.
No SST solve, CFD accuracy or wall-law qualification is implied by its timing.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import verify_native_flow as native
from verify_sst_decay import require,finite,sha
from verify_sst_spatial import point_segment


def audit(mesh_path,prefix):
    mesh=native.read_cm2d(mesh_path);m=native.measure(mesh,1e-10,1e-9)
    require(not m.issues,str(m.issues))
    x0,y0,x1,y1=m.bounds;tol=1e-12*max(x1-x0,y1-y0)
    require(tol>0,'degenerate rectangle')
    walls={e.id for e in mesh.edges if e.neighbour<0}
    for fid in walls:
        e=mesh.edges[fid];a,b=mesh.vertices[e.v0],mesh.vertices[e.v1]
        require(any(abs(a[axis]-value)<=tol and abs(b[axis]-value)<=tol for axis,value in [(0,x0),(0,x1),(1,y0),(1,y1)]),
                'not an axis-aligned rectangular boundary')
    require(math.isclose(m.total_area,(x1-x0)*(y1-y0),rel_tol=1e-10,abs_tol=1e-12),'rectangle area mismatch')
    meta=json.loads(Path(str(prefix)+'.json').read_text())
    require(meta['cells']==len(mesh.cells) and meta['wallSegments']==len(walls),'probe geometry counts mismatch')
    require(0<meta['segmentTests']<=len(mesh.cells)*len(walls),'invalid segment work count')
    require(finite(meta['distanceSeconds'])>=0,'invalid timing')
    with Path(str(prefix)+'.csv').open() as f:rows=list(csv.DictReader(f))
    require(len(rows)==len(mesh.cells),'missing distances')
    max_error=0.
    for i,(row,p) in enumerate(zip(rows,m.centroids)):
        require(int(row['cell'])==i,'cell ordering mismatch')
        fid=int(row['nearestFace']);require(fid in walls,'nearest face is not a wall')
        expected=min(p[0]-x0,x1-p[0],p[1]-y0,y1-p[1]);actual=finite(row['distance'])
        edge=mesh.edges[fid];nearest=point_segment(p,mesh.vertices[edge.v0],mesh.vertices[edge.v1])
        max_error=max(max_error,abs(actual-expected),abs(nearest-expected))
        require(actual>0 and abs(actual-expected)<=tol and abs(nearest-expected)<=tol,'cell wall distance/nearest segment mismatch')
    return dict(valid=True,scope='all-cell rectangular distance only',mesh=str(mesh_path),meshSha256=sha(mesh_path),
                prefix=str(prefix),allCellsAudited=len(rows),maxError=max_error,**meta,
                bruteForceSegmentTests=len(mesh.cells)*len(walls),
                csvSha256=sha(str(prefix)+'.csv'),metaSha256=sha(str(prefix)+'.json'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=Path,required=True);parser.add_argument('--prefix',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=audit(args.mesh.resolve(),args.prefix.resolve())
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
