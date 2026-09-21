#!/usr/bin/env python3
"""Bilinear coarse-to-fine initial iterate for complete Cartesian reference grids.

This is not a physical restart, a general Cut-cell mapper, or an accepted fine
solution. The fine solver must re-establish all boundary/flux/equation gates.
"""
import argparse,csv,hashlib,json,math,time
from pathlib import Path
from bisect import bisect_left
import sys
from verify_native_flow import read_cm2d


def centres(mesh):
    result=[]
    for cell in mesh.cells:
        points=[mesh.vertices[v] for v in cell.vertices]
        xs=sorted(set(p[0] for p in points));ys=sorted(set(p[1] for p in points))
        if len(points)!=4 or len(xs)!=2 or len(ys)!=2 or set(points)!={(x,y) for x in xs for y in ys}:
            raise ValueError('Reference mapper requires complete axis-aligned rectangles')
        result.append((xs[0]+(xs[1]-xs[0])/2,ys[0]+(ys[1]-ys[0])/2))
    return result


def interpolate(source_mesh,source_rows,target_mesh):
    source=centres(source_mesh);target=centres(target_mesh)
    xs=sorted(set(p[0] for p in source));ys=sorted(set(p[1] for p in source))
    if len(xs)<2 or len(ys)<2 or len(xs)*len(ys)!=len(source):raise ValueError('Coarse mesh must be a complete tensor grid')
    if len(source_rows)!=len(source):raise ValueError('Source fields do not cover mesh')
    values={}
    for i,(point,row) in enumerate(zip(source,source_rows)):
        if int(row['cell'])!=i:raise ValueError('Source field cell order differs')
        scale=max(abs(point[0]),abs(point[1]),1.)
        if any(abs(float(row[k])-point[j])>32*sys.float_info.epsilon*scale for j,k in enumerate(['x','y'])):
            raise ValueError('Source field coordinates differ from mesh')
        ix=bisect_left(xs,point[0]);iy=bisect_left(ys,point[1])
        if (ix,iy) in values:raise ValueError('Coarse grid overlaps')
        value=[float(row[k]) for k in ['u','v','p']]
        if not all(map(math.isfinite,value)):raise ValueError('Nonfinite source field')
        values[ix,iy]=value
    # Verify identical outer domain; extrapolate only the half-cell boundary rim.
    bounds=lambda mesh:[func(p[k] for p in mesh.vertices) for func in [min,max] for k in [0,1]]
    sb=bounds(source_mesh);tb=bounds(target_mesh)
    if any(abs(a-b)>32*sys.float_info.epsilon*max(1.,*map(abs,sb)) for a,b in zip(sb,tb)):
        raise ValueError('Coarse and target domains differ')
    result=[]
    for x,y in target:
        ix=max(0,min(len(xs)-2,bisect_left(xs,x)-1));iy=max(0,min(len(ys)-2,bisect_left(ys,y)-1))
        tx=(x-xs[ix])/(xs[ix+1]-xs[ix]);ty=(y-ys[iy])/(ys[iy+1]-ys[iy])
        value=[(1-ty)*((1-tx)*values[ix,iy][k]+tx*values[ix+1,iy][k])
               +ty*((1-tx)*values[ix,iy+1][k]+tx*values[ix+1,iy+1][k]) for k in range(3)]
        if not all(map(math.isfinite,value)):raise ValueError('Interpolation overflow')
        result.append(value)
    return target,result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ['source-mesh','source-prefix','target-mesh','output']:p.add_argument('--'+flag,type=Path,required=True)
    a=p.parse_args();start=time.monotonic()
    if a.output.exists():p.error('Output exists; preserve earlier evidence')
    summary_path=Path(str(a.source_prefix)+'.json');fields_path=Path(str(a.source_prefix)+'.cells.csv')
    summary=json.loads(summary_path.read_text())
    if not summary.get('converged') or not summary.get('strictLinearFinal'):p.error('Source must be converged with strict final linear certification')
    source=read_cm2d(a.source_mesh);target=read_cm2d(a.target_mesh)
    points,values=interpolate(source,list(csv.DictReader(fields_path.open())),target)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w') as f:
        writer=csv.writer(f);writer.writerow(['cell','x','y','u','v','p'])
        writer.writerows([i,*point,*value] for i,(point,value) in enumerate(zip(points,values)))
    report={'scope':__doc__,'cells':len(points),'seconds':time.monotonic()-start,'fineSolutionAccepted':False,
        'inputHashes':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [a.source_mesh,a.target_mesh,summary_path,fields_path,Path(__file__)]},
        'outputSha256':hashlib.sha256(a.output.read_bytes()).hexdigest()}
    Path(str(a.output)+'.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
