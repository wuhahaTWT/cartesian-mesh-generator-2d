#!/usr/bin/env python3
"""Independent OpenFOAM face-plane concavity diagnostic (not overall mesh PASS)."""
from __future__ import annotations
import argparse, hashlib, json, math, re
from pathlib import Path

PLANAR_COS_ANGLE = 1.0e-6

def _list_body(path: Path):
    text = re.sub(r'/\*.*?\*/|//[^\n]*', '', path.read_text(), flags=re.S)
    m = re.search(r'(?:^|\n)\s*(\d+)\s*\(', text)
    if not m: raise ValueError(f"missing Foam list: {path}")
    start=text.find('(', m.start()); depth=0
    for i in range(start,len(text)):
        if text[i]=='(': depth+=1
        elif text[i]==')':
            depth-=1
            if depth==0:
                if text[i+1:].strip(): raise ValueError(f'extra content after Foam list: {path}')
                return int(m.group(1)), text[start+1:i]
    raise ValueError(f"unbalanced Foam list: {path}")

def read_vectors(path: Path):
    declared,body=_list_body(path)
    records=re.findall(r'\([^()]*\)',body)
    if re.sub(r'\([^()]*\)', '', body).strip(): raise ValueError('garbage in point list')
    if len(records)!=declared or any(len(re.findall(r'[^\s]+',r[1:-1]))!=3 for r in records): raise ValueError('point list count/record mismatch')
    out=[tuple(map(float,re.findall(r'[^\s]+',r[1:-1]))) for r in records]
    if any(not math.isfinite(x) for p in out for x in p): raise ValueError('non-finite point')
    return out

def read_faces(path: Path):
    declared,body=_list_body(path)
    records=re.findall(r'(\d+)\s*\(([^()]*)\)',body)
    if re.sub(r'\d+\s*\([^()]*\)', '', body).strip(): raise ValueError('garbage in face list')
    if len(records)!=declared: raise ValueError('face list count mismatch')
    out=[]
    for n,s in records:
        vals=s.split()
        if any(not re.fullmatch(r'[-+]?\d+',x) for x in vals): raise ValueError('invalid face token')
        if len(vals)!=int(n) or len(set(vals))!=len(vals) or len(vals)<3: raise ValueError('invalid face record')
        out.append([int(x) for x in vals])
    return out

def read_labels(path: Path):
    declared,body=_list_body(path); vals=body.split()
    if any(not re.fullmatch(r'[-+]?\d+',x) for x in vals): raise ValueError('invalid label token')
    if len(vals)!=declared: raise ValueError('label list count mismatch')
    return [int(x) for x in vals]

def sub(a,b): return tuple(a[i]-b[i] for i in range(3))
def add(a,b): return tuple(a[i]+b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def mag(a): return math.sqrt(dot(a,a))

def face_geometry(face, points):
    """Return (centre, area vector), using local-origin Newell geometry.

    The centre is equivalent in construction to OpenFOAM face::centre's
    area-weighted central triangle decomposition (face.C:513-575); this is an
    independent implementation, not a claim of bitwise identity.
    """
    if len(face)<3 or len(set(face))!=len(face): raise ValueError('empty, short, or repeated face')
    if any(i<0 or i>=len(points) for i in face): raise ValueError('face vertex index out of range')
    p=[points[i] for i in face]; origin=p[0]; q=[sub(x,origin) for x in p]
    if any(not math.isfinite(x) for point in p for x in point): raise ValueError('non-finite point')
    cp=tuple(math.fsum(x[k] for x in q)/len(q) for k in range(3)); sum_a=0.; sum_ac=(0.,0.,0.)
    for u,v in zip(q,q[1:]+q[:1]):
        ta=mag(cross(sub(u,cp),sub(v,cp))); sum_a+=ta
        sum_ac=add(sum_ac,tuple(ta*x for x in add(add(u,v),cp)))
    centre=add(origin,tuple(x/(3*sum_a) for x in sum_ac)) if sum_a>1e-300 else add(origin,cp)
    area=(0.,0.,0.)
    for u,v in zip(q,q[1:]+q[:1]): area=add(area,cross(u,v))
    if mag(area)<=1e-300: raise ValueError('degenerate face')
    return centre, area

def classify(cells, points, faces, owner):
    if len(owner)!=len(faces): raise ValueError('owner/face count mismatch')
    if any(x<0 or x>=len(cells) for x in owner): raise ValueError('owner cell index out of range')
    if any(not c or len(c)<4 or len(set(c))!=len(c) for c in cells): raise ValueError('cell is empty, open, or repeats a face')
    uses=[0]*len(faces)
    for ci,c in enumerate(cells):
        edges={}
        for fi in c:
            if fi<0 or fi>=len(faces): raise ValueError('cell face index out of range')
            uses[fi]+=1; f=faces[fi]; direction=1 if owner[fi]==ci else -1
            seq=f if direction==1 else list(reversed(f))
            for a,b in zip(seq,seq[1:]+seq[:1]): edges.setdefault(tuple(sorted((a,b))),[]).append((a,b))
        if any(len(v)!=2 or v[0]!=(v[1][1],v[1][0]) for v in edges.values()): raise ValueError('cell is not closed with consistent face orientation')
    if any(x>2 for x in uses): raise ValueError('face used by more than two cells')
    if any(owner[fi] not in range(len(cells)) or fi not in cells[owner[fi]] for fi in range(len(faces))): raise ValueError('owner does not match cell incidence')
    geom=[face_geometry(f,points) for f in faces]; original=set(); positive=set(); near=set(); safe_pair=set(); pairs={'positive':0,'near_coplanar':0,'safe':0}; maxima=[]
    for ci, cfaces in enumerate(cells):
        best=-2.; best_pair=None; has_original=False; has_positive=False; has_near=False
        for i,fi in enumerate(cfaces):
            fc,av=geom[fi]; n=tuple(x/max(mag(av),1e-300) for x in av)
            if owner[fi]!=ci: n=tuple(-x for x in n)
            for fj in cfaces:
                if fj==fi: continue
                p=sub(geom[fj][0],fc); d=dot(tuple(x/max(mag(p),1e-300) for x in p),n)
                if d>best: best,best_pair=d,(fi,fj)
                if d>PLANAR_COS_ANGLE: pairs['positive']+=1; has_positive=True
                elif d>=-PLANAR_COS_ANGLE: pairs['near_coplanar']+=1; has_near=True
                else: pairs['safe']+=1; safe_pair.add(ci)
                if d>-PLANAR_COS_ANGLE: has_original=True
        maxima.append({'cell':ci,'maximum_dot':best,'pair':best_pair})
        if has_original: original.add(ci)
        if has_positive: positive.add(ci)
        elif has_near: near.add(ci)
    safe=set(range(len(cells)))-positive-near
    return {'original_flagged':original,'positive_side':positive,'near_coplanar_only':near,'safe':safe,'cells_with_safe_pair':safe_pair,'pair_counts':pairs,'maxima':maxima}

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''): h.update(chunk)
    return h.hexdigest()

def run_case(case: Path, expected_set: Path|None=None):
    points=read_vectors(case/'points'); faces=read_faces(case/'faces'); owner=read_labels(case/'owner'); neighbour=read_labels(case/'neighbour')
    if len(neighbour)>len(faces): raise ValueError('neighbour count exceeds faces')
    if len(owner)!=len(faces): raise ValueError('owner/face count mismatch')
    if not owner or min(owner+neighbour)<0: raise ValueError('empty mesh or negative cell index')
    n=max(owner+neighbour)+1
    if n>len(faces): raise ValueError('cell index exceeds possible closed mesh size')
    cells=[[] for _ in range(n)]
    for fi,ci in enumerate(owner):
        if ci<0 or ci>=n: raise ValueError('owner cell index out of range')
        cells[ci].append(fi)
    for fi,ci in enumerate(neighbour):
        if ci<0 or ci>=n or owner[fi]==ci: raise ValueError('invalid neighbour relation')
        cells[ci].append(fi)
    result=classify(cells,points,faces,owner)
    actual=set(read_labels(expected_set)) if expected_set else None
    if actual is not None and any(x<0 or x>=n for x in actual): raise ValueError('expected set cell index out of range')
    files={name:sha256(case/name) for name in ('points','faces','owner','neighbour')}
    out={'format':'openfoam-face-plane-audit-v1','scope':'face-plane classification only; not overall mesh PASS','case':str(case),'input_files_sha256':files,'threshold':PLANAR_COS_ANGLE,'cell_count':n,'face_count':len(faces),'point_count':len(points),'original_flagged':len(result['original_flagged']),'positive_side':len(result['positive_side']),'near_coplanar_only':len(result['near_coplanar_only']),'safe':len(result['safe']),'cells_with_safe_pair':len(result['cells_with_safe_pair']),'pair_counts':result['pair_counts'],'maximum_dot':max((x['maximum_dot'] for x in result['maxima']),default=None)}
    if actual is None: out.update({'expected_set':None,'expected_set_matches':None,'actual_set_count':None,'actual_set_only':None,'missing_set':None})
    else: out.update({'expected_set':str(expected_set),'expected_set_matches':actual==result['original_flagged'],'actual_set_count':len(actual),'actual_set_only':sorted(actual-result['original_flagged']),'missing_set':sorted(result['original_flagged']-actual)})
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('case',type=Path,help='OpenFOAM constant/polyMesh directory'); p.add_argument('--output',type=Path); p.add_argument('--expected-set',type=Path); a=p.parse_args()
    out=run_case(a.case,a.expected_set); text=json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+'\n'; print(text,end='')
    if a.output: a.output.write_text(text)
    if out['expected_set_matches'] is False: raise SystemExit(1)
if __name__=='__main__': main()
