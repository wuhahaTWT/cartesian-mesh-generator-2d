#!/usr/bin/env python3
"""Independent directional-connectivity check for uncoupled 2D CM2D exports.

This evaluates the OpenFOAM 2606 cellDeterminant formula for a planar mesh
extruded with uniform thickness and empty front/back patches. Boundary faces
are excluded. Coupled/periodic boundaries are outside this reader's scope.
It is a diagnostic, not a substitute for topology checks or actual checkMesh.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from check_mesh_resolution import read_cm2d

MINIMUM_DETERMINANT = 0.001


def determinant(internal_edge_vectors: list[tuple[float, float]]) -> float:
    """det(sum(n n^T))/8 with normals normalized by mean internal edge length.

    Uniform extrusion thickness cancels. Rotating tangents to normals leaves
    this determinant unchanged. The pairwise cross-product identity avoids
    subtracting almost equal tensor products for nearly parallel directions.
    """
    if any(not math.isfinite(v) for edge in internal_edge_vectors for v in edge):
        raise ValueError("non-finite internal edge vector")
    if not internal_edge_vectors:
        return 0.0
    lengths = [math.hypot(*edge) for edge in internal_edge_vectors]
    if any(length <= 0 or not math.isfinite(length) for length in lengths):
        raise ValueError("degenerate internal edge")
    mean_length = math.fsum(length / len(lengths) for length in lengths)
    normalized = [(x / mean_length, y / mean_length) for x, y in internal_edge_vectors]
    return math.fsum(
        (a[0] * b[1] - a[1] * b[0]) ** 2
        for i, a in enumerate(normalized) for b in normalized[i + 1:]
    ) / 8.0


def measure(cm2d: Path) -> dict:
    vertices, edges, cells = read_cm2d(cm2d)
    directions: list[list[tuple[float, float]]] = [[] for _ in cells]
    internal_edge_ids: list[list[int]] = [[] for _ in cells]
    for edge in edges:
        if edge["neighbour"] < 0:
            continue
        if edge["owner"] == edge["neighbour"]:
            raise ValueError("internal face has the same owner and neighbour")
        a, b = vertices[edge["v0"]], vertices[edge["v1"]]
        vector = (b[0] - a[0], b[1] - a[1])
        for cell_id in (edge["owner"], edge["neighbour"]):
            if not 0 <= cell_id < len(cells):
                raise ValueError("internal face references an invalid cell")
            if cells[cell_id]["edges"].count(edge["id"]) != 1:
                raise ValueError("internal face incidence disagrees with cell record")
            directions[cell_id].append(vector)
            internal_edge_ids[cell_id].append(edge["id"])
    values = [determinant(vectors) for vectors in directions]
    failures = [{"cell_id": i, "determinant": value,
                 "internal_edge_ids": internal_edge_ids[i]}
                for i, value in enumerate(values) if value < MINIMUM_DETERMINANT]
    return {
        "format": "cartmesh2d-directional-connectivity-v1",
        "scope": "uncoupled planar mesh; uniform extrusion; empty front/back",
        "external_checkMesh": "not_run_by_this_reader",
        "cell_count": len(cells),
        "minimum_required": MINIMUM_DETERMINANT,
        "minimum_measured": min(values) if values else None,
        "failed_cell_count": len(failures),
        "failed_cells": failures,
        "valid": bool(cells) and not failures,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cm2d", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = measure(args.cm2d)
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
