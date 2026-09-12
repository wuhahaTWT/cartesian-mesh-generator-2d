#!/usr/bin/env python3
"""Read-only OpenFOAM 2606 determinant/face-plane audit for a foam case."""
from __future__ import annotations
import argparse, json, math, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/verification"))
from check_face_planes import read_vectors, read_faces, read_labels, face_geometry, classify, sha256

THRESHOLD = 1.0e-3
NORMALIZATION = 0.037037037037037  # OpenFOAM primitiveMeshGeometry: 1/27

def foam_set(path: Path) -> list[int]:
    s = path.read_text()
    m = re.search(r"\n\s*(\d+)\s*\((.*?)\n\)", s, re.S)
    if not m:
        raise ValueError(f"invalid OpenFOAM set: {path}")
    vals = [int(x) for x in m.group(2).split()]
    if len(vals) != int(m.group(1)):
        raise ValueError(f"set count mismatch: {path}")
    return vals

def determinant_all_faces(cell_faces, celli, owner, geometry):
    # OpenFOAM 2606 primitiveMeshGeometry::checkCellDeterminant:
    # areaSum += Sf*(Sf/mag(Sf)); scaledDet =
    # det(areaSum/magAreaSum) / (1/27).  The cell's outward sign cancels.
    area_sum = [[0.0] * 3 for _ in range(3)]
    magnitude_sum = 0.0
    for fi in cell_faces:
        av = geometry[fi][1]
        length = math.sqrt(sum(x*x for x in av))
        magnitude_sum += length
        for i in range(3):
            for j in range(3):
                area_sum[i][j] += av[i] * av[j] / length
    t = [[x / magnitude_sum for x in row] for row in area_sum]
    det = (t[0][0]*(t[1][1]*t[2][2]-t[1][2]*t[2][1])
           - t[0][1]*(t[1][0]*t[2][2]-t[1][2]*t[2][0])
           + t[0][2]*(t[1][0]*t[2][1]-t[1][1]*t[2][0]))
    return det / NORMALIZATION

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("poly_mesh", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    p = args.poly_mesh
    points, faces = read_vectors(p/"points"), read_faces(p/"faces")
    owner, neighbour = read_labels(p/"owner"), read_labels(p/"neighbour")
    n = max(owner + neighbour) + 1
    cells = [[] for _ in range(n)]
    for fi, ci in enumerate(owner): cells[ci].append(fi)
    for fi, ci in enumerate(neighbour): cells[ci].append(fi)
    geometry = [face_geometry(f, points) for f in faces]
    determinants = [determinant_all_faces(c, ci, owner, geometry) for ci, c in enumerate(cells)]
    low = [ci for ci, d in enumerate(determinants) if d < THRESHOLD]

    low_detail = []
    for ci in low:
        fs = cells[ci]
        verts = sorted({v for fi in fs for v in faces[fi]})
        xy = sorted(set((points[v][0], points[v][1]) for v in verts))
        cx, cy = (sum(x for x, _ in xy)/len(xy), sum(y for _, y in xy)/len(xy))
        polygon = sorted(set(xy), key=lambda q: math.atan2(q[1]-cy, q[0]-cx))
        adjacent = sorted({neighbour[fi] if owner[fi] == ci else owner[fi]
                           for fi in fs if fi < len(neighbour) and
                           (owner[fi] == ci or neighbour[fi] == ci)})
        low_detail.append({"cell": ci, "all_face_determinant": determinants[ci],
                           "faces": fs,
                           "polygon_xy": polygon, "neighbor_cells": adjacent})
    expected = p/"sets/concaveCells"
    plane = classify(cells, points, faces, owner)
    out = {"format": "cartmesh2d-external-quality-audit-v1",
           "scope": "independent extruded determinant and face-plane diagnostics; not overall mesh PASS",
           "input_files_sha256": {name: sha256(p/name) for name in ("points", "faces", "owner", "neighbour")},
           "formula": {"areaSum": "sum(Sf outer Sf / mag(Sf)) over all cell faces",
                       "magnitudeSum": "sum(mag(Sf)) over all cell faces",
                       "scaledDeterminant": "det(areaSum/magnitudeSum)/0.037037037037037",
                       "normalization": NORMALIZATION, "threshold": THRESHOLD},
           "mesh": {"cells": n, "faces": len(faces), "points": len(points)},
           "all_face_determinant": {"low_cell_count": len(low), "low_cells": low_detail,
                                    "minimum": min(determinants)},
           "face_plane": {"original_flagged": len(plane["original_flagged"]),
                          "positive_side": len(plane["positive_side"]),
                          "near_coplanar_only": len(plane["near_coplanar_only"]),
                          "maximum_dot": max(x["maximum_dot"] for x in plane["maxima"]),
                          "expected_set_matches": (set(foam_set(expected)) == plane["original_flagged"])
                          if expected.exists() else None}}
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
