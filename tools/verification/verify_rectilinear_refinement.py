#!/usr/bin/env python3
"""Verify exact nested wall-normal midpoint refinement of complete rectangles.

Geometry-only diagnostic. This does not replace native Solver quality or a
flow-equation audit. All original x coordinates and a declared bottom band
are fixed; every other original y interval must be bisected exactly once.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import verify_native_flow as native


def require(ok, message):
    if not ok:
        raise ValueError(message)


def tensor_axes(polygons):
    require(polygons, 'empty grid')
    require(all(math.isfinite(v) for p in polygons for point in p for v in point), 'nonfinite coordinate')
    xs=sorted({p[0] for polygon in polygons for p in polygon})
    ys=sorted({p[1] for polygon in polygons for p in polygon})
    require(len(xs)>=2 and len(ys)>=2, 'degenerate grid')
    xi,yi={x:i for i,x in enumerate(xs)},{y:i for i,y in enumerate(ys)}
    cells=set()
    for polygon in polygons:
        left,right=min(p[0] for p in polygon),max(p[0] for p in polygon)
        bottom,top=min(p[1] for p in polygon),max(p[1] for p in polygon)
        require(len(polygon)==4 and set(polygon)=={(left,bottom),(right,bottom),(right,top),(left,top)},
                'non-rectangular cell')
        require(all((a[0]==b[0]) != (a[1]==b[1]) for a,b in zip(polygon,polygon[1:]+polygon[:1])),
                'cell edges must follow the rectangle boundary')
        i,j=xi[left],yi[bottom]
        require(xi[right]==i+1 and yi[top]==j+1 and (i,j) not in cells, 'non-tensor or duplicate cell')
        cells.add((i,j))
    require(len(cells)==(len(xs)-1)*(len(ys)-1), 'incomplete rectangular partition')
    return xs,ys


def check_refinement(original, refined, preserve):
    x,y=tensor_axes(original);xx,yy=tensor_axes(refined)
    require(type(preserve) is int and 1<=preserve<len(y)-1, 'preserved band must leave rows to refine')
    require(x==xx, 'streamwise coordinates changed')
    expected=y[:preserve+1]
    for a,b in zip(y[preserve:],y[preserve+1:]):
        midpoint=float((Fraction(a)+Fraction(b))/2)
        require(a<midpoint<b, 'midpoint is not representable inside parent interval')
        expected.extend((midpoint,b))
    require(yy==expected, 'normal coordinates do not match the declared midpoint refinement')
    height=y[preserve]
    def lower(polygons):
        return Counter(tuple(sorted(p)) for p in polygons if max(v[1] for v in p)<=height)
    require(lower(original)==lower(refined), 'preserved cell polygons differ')
    width=Fraction(x[-1])-Fraction(x[0])
    old_area=width*(Fraction(y[-1])-Fraction(y[0]))
    new_area=sum(width*(Fraction(b)-Fraction(a)) for a,b in zip(yy,yy[1:]))
    require(old_area==new_area, 'exact partition area differs')
    return {'geometryValid':True,'scope':'Exact rectangular midpoint refinement only; no Solver or flow qualification',
            'originalCells':len(original),'refinedCells':len(refined),'nx':len(x)-1,
            'originalNy':len(y)-1,'refinedNy':len(yy)-1,'preservedRows':preserve,
            'preservedBandTop':height,'preservedCells':sum(lower(original).values()),
            'streamwiseNodesIdentical':True,'originalNormalNodesRetained':True,
            'lowerCellPolygonsIdentical':True,'eachUpperParentSplitOnce':True,
            'areaIdenticalAsExactRational':True,'domainArea':float(old_area),
            'areaDefinition':'Exact rational equality for parsed binary64 coordinates, not source CAD/decimal geometry',
            'originalY':y,'refinedY':yy}


def load(path):
    mesh=native.read_cm2d(path)
    measurement=native.measure(mesh,1e-10,1e-9)
    require(not measurement.issues,str(measurement.issues))
    return [[mesh.vertices[v] for v in c.vertices] for c in mesh.cells]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',type=Path,required=True)
    p.add_argument('--refined',type=Path,required=True)
    p.add_argument('--preserve-first-rows',type=int,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    result=check_refinement(load(args.original),load(args.refined),args.preserve_first_rows)
    result['files']=[{'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
                     for path in (args.original,args.refined)]
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('originalY','refinedY','files')}))


if __name__=='__main__':main()
