#!/usr/bin/env python3
"""Independent H4-3 pure Cut-cell fallback VTK/JSON checker."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path


class CheckError(ValueError):
    pass


def signed_area(poly: list[tuple[float, float]]) -> float:
    return 0.5 * sum(
        poly[i][0] * poly[(i + 1) % len(poly)][1]
        - poly[(i + 1) % len(poly)][0] * poly[i][1]
        for i in range(len(poly))
    )


def orient(a, b, c) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def proper_cross(a, b, c, d, eps: float) -> bool:
    ab_c, ab_d = orient(a, b, c), orient(a, b, d)
    cd_a, cd_b = orient(c, d, a), orient(c, d, b)
    return (((ab_c > eps and ab_d < -eps) or (ab_c < -eps and ab_d > eps))
            and ((cd_a > eps and cd_b < -eps) or (cd_a < -eps and cd_b > eps)))


def read_vtk(path: Path):
    tokens = path.read_text(encoding="utf-8").split()
    p_at = tokens.index("POINTS")
    n_points = int(tokens[p_at + 1])
    start = p_at + 3
    points = [(float(tokens[start + 3*i]), float(tokens[start + 3*i + 1]))
              for i in range(n_points)]
    c_at = tokens.index("CELLS", start + 3*n_points)
    n_cells = int(tokens[c_at + 1])
    cursor = c_at + 3
    cells = []
    for _ in range(n_cells):
        count = int(tokens[cursor]); cursor += 1
        cell = tuple(int(v) for v in tokens[cursor:cursor + count]); cursor += count
        if count < 3 or len(set(cell)) != count or any(v < 0 or v >= n_points for v in cell):
            raise CheckError("invalid polygon vertex IDs")
        cells.append(cell)
    return points, cells


def check(vtk_path: Path, report_path: Path) -> dict[str, object]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("h4_status") != "success" or report.get("mesh_mode") != "pure_cutcell_fallback":
        raise CheckError("report does not identify a successful pure Cut-cell fallback")
    if report.get("layer_status") != "failed":
        raise CheckError("fallback is incorrectly reported as an enabled layer")
    points, cells = read_vtk(vtk_path)
    scale = max(1.0, *(abs(value) for point in points for value in point))
    eps = 1.0e-11 * scale
    owners: dict[tuple[int, int], list[int]] = defaultdict(list)
    areas: list[float] = []
    boxes = []
    for cell_id, cell in enumerate(cells):
        poly = [points[vertex] for vertex in cell]
        area = signed_area(poly)
        if not math.isfinite(area) or area <= eps * eps:
            raise CheckError(f"cell {cell_id} has non-positive signed area")
        areas.append(area)
        boxes.append((min(p[0] for p in poly), min(p[1] for p in poly),
                      max(p[0] for p in poly), max(p[1] for p in poly)))
        for i, first in enumerate(cell):
            owners[tuple(sorted((first, cell[(i + 1) % len(cell)])))].append(cell_id)
    if any(len(incidence) not in (1, 2) for incidence in owners.values()):
        raise CheckError("duplicate or non-manifold edge incidence")
    for lhs in range(len(cells)):
        lhs_edges = [(cells[lhs][i], cells[lhs][(i + 1) % len(cells[lhs])])
                     for i in range(len(cells[lhs]))]
        for rhs in range(lhs + 1, len(cells)):
            if (boxes[lhs][2] < boxes[rhs][0] - eps or boxes[rhs][2] < boxes[lhs][0] - eps
                    or boxes[lhs][3] < boxes[rhs][1] - eps or boxes[rhs][3] < boxes[lhs][1] - eps):
                continue
            rhs_edges = [(cells[rhs][i], cells[rhs][(i + 1) % len(cells[rhs])])
                         for i in range(len(cells[rhs]))]
            if any(proper_cross(points[a], points[b], points[c], points[d], eps)
                   for a, b in lhs_edges for c, d in rhs_edges):
                raise CheckError(f"cell interiors cross: {lhs}, {rhs}")
    total = sum(areas)
    actual = float(report["actual_fluid_area"])
    tolerance = 1.0e-10 * max(1.0, abs(total), abs(actual))
    if abs(total - actual) > tolerance:
        raise CheckError("VTK cell area disagrees with the H4 report")
    if int(report.get("cell_count", -1)) != len(cells):
        raise CheckError("VTK cell count disagrees with the H4 report")
    return {
        "valid": True,
        "vertex_count": len(points),
        "cell_count": len(cells),
        "boundary_edge_count": sum(len(value) == 1 for value in owners.values()),
        "internal_edge_count": sum(len(value) == 2 for value in owners.values()),
        "minimum_signed_area": min(areas),
        "area_sum": total,
        "area_error": total - float(report["expected_fluid_area"]),
        "overlap_count": 0,
        "non_manifold_edge_count": 0,
        "issues": [],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vtk", type=Path)
    parser.add_argument("report", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = check(args.vtk, args.report)
    except (CheckError, OSError, ValueError) as exc:
        result = {"valid": False, "issues": [str(exc)]}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
