#!/usr/bin/env python3
"""Independent dyadic coverage and polygon classification audit (no native imports)."""
import argparse
import collections
import json
import math
from pathlib import Path


def intersects(a, b, box):
    # Closed-box segment clipping, also marks exact tangencies as intersected.
    # A dyadic coordinate constructed from a translated domain can miss an input
    # vertex by a few ulps. Account only for arithmetic roundoff here, not a
    # cell-sized geometric gap (the independent audit must still catch those).
    eps = 16 * max(math.ulp(x) for x in (*a, *b, *box))
    box = [box[0]-eps, box[1]-eps, box[2]+eps, box[3]+eps]
    t0, t1 = 0., 1.
    for axis in (0, 1):
        delta = b[axis] - a[axis]
        if delta == 0:
            if not box[axis] <= a[axis] <= box[axis + 2]:
                return False
        else:
            lo = (box[axis] - a[axis]) / delta
            hi = (box[axis + 2] - a[axis]) / delta
            t0, t1 = max(t0, min(lo, hi)), min(t1, max(lo, hi))
            if t0 > t1:
                return False
    return True


def inside(point, loops):
    x, y = point
    result = False
    for loop in loops:
        for a, b in zip(loop, loop[1:] + loop[:1]):
            if (a[1] > y) != (b[1] > y):
                if x < a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]):
                    result = not result
    return result


def audit(path, check_classification=True):
    data = json.loads(Path(path).read_text())
    assert data['format'] == 'cartmesh2d-background-v1'
    assert data['solver_ready'] is False and data['retains_solid_interior'] is True
    assert data['classification_names'] == ['outside', 'inside', 'intersected']
    assert data['mode'] in ('uniform', 'adaptive')
    domain = data['domain']
    assert len(domain) == 4 and all(math.isfinite(x) for x in domain)
    assert domain[2] > domain[0] and domain[3] > domain[1]
    cells, loops = data['cells'], data['boundary_loops']
    assert cells and loops and all(len(loop) >= 3 for loop in loops)
    max_level = max(c['level'] for c in cells)
    assert 0 <= max_level <= 28
    leaves, keys = set(), set()
    counts = collections.Counter()
    area_units = 0
    for index, cell in enumerate(cells):
        level, ix, iy = cell['level'], cell['ix'], cell['iy']
        assert all(type(x) is int for x in (level, ix, iy))
        assert 0 <= level <= max_level and 0 <= ix < 2**level and 0 <= iy < 2**level
        assert cell['id'] == index
        key = cell['key']
        assert isinstance(key, str) and key.isdecimal() and key not in keys
        keys.add(key)
        node = (level, ix, iy)
        assert node not in leaves, 'duplicate leaf'
        leaves.add(node)
        box = cell['bounds']
        expected = [domain[0] + (domain[2]-domain[0])*ix/2**level,
                    domain[1] + (domain[3]-domain[1])*iy/2**level,
                    domain[0] + (domain[2]-domain[0])*(ix+1)/2**level,
                    domain[1] + (domain[3]-domain[1])*(iy+1)/2**level]
        assert len(box) == 4 and all(math.isfinite(x) for x in box)
        assert all(math.isclose(x, y, rel_tol=2e-14, abs_tol=1e-14*max(domain[2]-domain[0], domain[3]-domain[1])) for x, y in zip(box, expected))
        area_units += 4**(max_level-level)
        kind = cell['classification']
        assert type(kind) is int and kind in (0, 1, 2)
        counts[kind] += 1
        if check_classification:
            cut = any(intersects(a, b, box) for loop in loops for a, b in zip(loop, loop[1:] + loop[:1]))
            expected_kind = 2 if cut else int(inside(((box[0]+box[2])/2, (box[1]+box[3])/2), loops))
            assert kind == expected_kind, f'classification mismatch in cell {index}'
    # Dyadic squares only overlap if one is an ancestor of the other.
    for level, ix, iy in leaves:
        for ancestor in range(level):
            shift = level - ancestor
            assert (ancestor, ix >> shift, iy >> shift) not in leaves, 'overlapping leaves'
    assert area_units == 4**max_level, 'domain not completely covered'
    assert [counts[i] for i in range(3)] == data['counts']
    if data['mode'] == 'uniform':
        assert all(c['level'] == max_level for c in cells)
        assert len(cells) == 4**max_level
    # Independent face-neighbour sweep in integer finest-grid coordinates.
    faces = [collections.defaultdict(list) for _ in range(4)]
    for level, ix, iy in leaves:
        scale = 2**(max_level-level)
        x0, x1, y0, y1 = ix*scale, (ix+1)*scale, iy*scale, (iy+1)*scale
        for axis, at, lo, hi in ((0,x0,y0,y1),(1,x1,y0,y1),(2,y0,x0,x1),(3,y1,x0,x1)):
            faces[axis][at].append((lo,hi,level))
    for left, right in ((0,1),(2,3)):
        for at in faces[left].keys() & faces[right].keys():
            a, b = sorted(faces[left][at]), sorted(faces[right][at])
            i = j = 0
            while i < len(a) and j < len(b):
                if max(a[i][0], b[j][0]) < min(a[i][1], b[j][1]):
                    assert abs(a[i][2]-b[j][2]) <= 1, '2:1 violation'
                if a[i][1] <= b[j][1]: i += 1
                else: j += 1
    return {'cells': len(cells), 'counts': data['counts'], 'full_domain_coverage': True,
            'classification_checked': check_classification, 'balanced': True, 'solver_ready': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path')
    args = parser.parse_args()
    print(json.dumps(audit(args.path)))
