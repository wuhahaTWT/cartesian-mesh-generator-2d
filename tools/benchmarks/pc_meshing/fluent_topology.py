"""Exact same-cell identity across Fluent's parallel cell renumbering.

Uses the accepted case's node coordinates and complete cell-face incidence.
No nearest-neighbour lookup, no wall movement and no tolerance-based matching.
"""
from pathlib import Path
import h5py

def case_cell_map(mesh,path):
    native_nodes={tuple(p):i for i,p in enumerate(mesh.vertices)}
    if len(native_nodes)!=len(mesh.vertices):raise ValueError('native vertices are not geometrically unique')
    signatures={tuple(sorted(set(c.vertices))):c.id for c in mesh.cells}
    if len(signatures)!=len(mesh.cells):raise ValueError('native cell signatures are not unique')
    with h5py.File(path,'r') as f:
        base=f['meshes/1']
        coords=base['nodes/coords/1'][:]
        if coords.shape!=(len(mesh.vertices),2):raise ValueError('Fluent node count or dimension changed')
        try:nodes=[native_nodes[tuple(p)] for p in coords]
        except KeyError as error:raise ValueError('Fluent changed a node coordinate; exact identity unavailable') from error
        if len(set(nodes))!=len(nodes):raise ValueError('Fluent vertex mapping is not bijective')
        sizes=base['faces/nodes/1/nnodes'][:]
        if any(sizes!=2):raise ValueError('Expected genuine 2D segment faces')
        connectivity=base['faces/nodes/1/nodes'][:].reshape(-1,2)
        if len(connectivity)!=len(mesh.edges):raise ValueError('Fluent changed the face count')
        c0=base['faces/c0/1'][:]
        c1=base['faces/c1/1'][:]
        begin0=int(base['faces/c0/1'].attrs.get('minId',[1])[0])-1
        begin1=int(base['faces/c1/1'].attrs.get('minId',[1])[0])-1
        cell_nodes=[set() for _ in mesh.cells]
        for owners,begin in ((c0,begin0),(c1,begin1)):
            for offset,owner in enumerate(owners):
                if int(owner)==0:continue
                i=int(owner)-1
                if not 0<=i<len(cell_nodes):raise ValueError('Fluent owner out of input cell range')
                edge=connectivity[begin+offset]
                cell_nodes[i].update(nodes[int(n)-1] for n in edge)
        result={}
        for i,vertices in enumerate(cell_nodes):
            key=tuple(sorted(vertices))
            if key not in signatures:raise ValueError('Fluent cell incidence differs from input')
            result[i+1]=signatures[key]
        if len(set(result.values()))!=len(mesh.cells):raise ValueError('Fluent cells are not bijective')
    return result
