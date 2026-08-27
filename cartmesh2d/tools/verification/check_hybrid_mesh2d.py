#!/usr/bin/env python3
"""Independent H4-2 VTK/JSON conformal-topology checker."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path


class CheckError(ValueError):
    pass


def area(poly: list[tuple[float, float]]) -> float:
    return 0.5 * sum(
        poly[i][0] * poly[(i + 1) % len(poly)][1]
        - poly[(i + 1) % len(poly)][0] * poly[i][1]
        for i in range(len(poly))
    )


def orient(a, b, c) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(p, a, b, eps: float) -> bool:
    return (abs(orient(a, b, p)) <= eps
            and min(a[0], b[0]) - eps <= p[0] <= max(a[0], b[0]) + eps
            and min(a[1], b[1]) - eps <= p[1] <= max(a[1], b[1]) + eps)


def proper_cross(a, b, c, d, eps: float) -> bool:
    ab_c, ab_d = orient(a, b, c), orient(a, b, d)
    cd_a, cd_b = orient(c, d, a), orient(c, d, b)
    return (((ab_c > eps and ab_d < -eps) or (ab_c < -eps and ab_d > eps))
            and ((cd_a > eps and cd_b < -eps) or (cd_a < -eps and cd_b > eps)))


def point_inside(p, poly, eps: float) -> bool:
    inside = False
    for i, a in enumerate(poly):
        b = poly[(i + 1) % len(poly)]
        if on_segment(p, a, b, eps):
            return False
        if (a[1] > p[1]) != (b[1] > p[1]):
            crossing_x = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            if crossing_x > p[0] + eps:
                inside = not inside
    return inside


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
    data_at = tokens.index("CELL_DATA", cursor)
    scalar_at = tokens.index("hybrid_kind", data_at)
    lookup_at = tokens.index("LOOKUP_TABLE", scalar_at)
    kinds = [int(v) for v in tokens[lookup_at + 2:lookup_at + 2 + n_cells]]
    if len(kinds) != n_cells or any(v not in (0, 1, 2) for v in kinds):
        raise CheckError("invalid hybrid_kind cell data")
    return points, cells, kinds


def check(vtk_path: Path, report_path: Path) -> dict[str, object]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("hybrid_status") != "success":
        raise CheckError("hybrid report is not successful")
    points, cells, kinds = read_vtk(vtk_path)
    scale = max(1.0, *(abs(v) for p in points for v in p))
    eps = 1.0e-11 * scale
    owners: dict[tuple[int, int], list[int]] = defaultdict(list)
    areas = []
    boxes = []
    for cid, cell in enumerate(cells):
        poly = [points[v] for v in cell]
        signed = area(poly)
        if not math.isfinite(signed) or signed <= eps * eps:
            raise CheckError(f"cell {cid} has non-positive signed area")
        areas.append(signed)
        boxes.append((min(p[0] for p in poly), min(p[1] for p in poly),
                      max(p[0] for p in poly), max(p[1] for p in poly)))
        for i, first in enumerate(cell):
            edge = tuple(sorted((first, cell[(i + 1) % len(cell)])))
            owners[edge].append(cid)
    if any(len(value) not in (1, 2) for value in owners.values()):
        raise CheckError("dangling duplicate or non-manifold edge incidence")

    interface = {edge: own for edge, own in owners.items()
                 if len(own) == 2 and ((kinds[own[0]] == 0) != (kinds[own[1]] == 0))}
    if any(sorted(kinds[c] == 0 for c in own) != [False, True] for own in interface.values()):
        raise CheckError("interface edge is not a layer/remainder pair")
    degrees: dict[int, int] = defaultdict(int)
    interface_length = 0.0
    for first, second in interface:
        degrees[first] += 1; degrees[second] += 1
        interface_length += math.dist(points[first], points[second])
    if not interface or any(value != 2 for value in degrees.values()):
        raise CheckError("outer-envelope interface is empty or not two-valent")

    # Broad-phase pair scan plus strict containment/proper crossing rejects overlap.
    for lhs in range(len(cells)):
        lp = [points[v] for v in cells[lhs]]
        le = [(cells[lhs][i], cells[lhs][(i + 1) % len(cells[lhs])])
              for i in range(len(cells[lhs]))]
        for rhs in range(lhs + 1, len(cells)):
            if (boxes[lhs][2] < boxes[rhs][0] - eps or boxes[rhs][2] < boxes[lhs][0] - eps
                    or boxes[lhs][3] < boxes[rhs][1] - eps or boxes[rhs][3] < boxes[lhs][1] - eps):
                continue
            re = [(cells[rhs][i], cells[rhs][(i + 1) % len(cells[rhs])])
                  for i in range(len(cells[rhs]))]
            for a_id, b_id in le:
                for c_id, d_id in re:
                    if proper_cross(points[a_id], points[b_id], points[c_id], points[d_id], eps):
                        raise CheckError(f"cell interiors cross: {lhs}, {rhs}")
            rp = [points[v] for v in cells[rhs]]
            if any(v not in cells[rhs] and point_inside(points[v], rp, eps) for v in cells[lhs]):
                raise CheckError(f"cell interior overlap: {lhs}, {rhs}")
            if any(v not in cells[lhs] and point_inside(points[v], lp, eps) for v in cells[rhs]):
                raise CheckError(f"cell interior overlap: {lhs}, {rhs}")

    total = sum(areas)
    expected = float(report["expected_fluid_area"])
    tolerance = 1.0e-10 * max(1.0, abs(total), abs(expected))
    checks = {
        "vertex_count": len(points), "cell_count": len(cells),
        "interface_edge_count": len(interface), "interface_vertex_count": len(degrees),
    }
    for key, actual in checks.items():
        if int(report.get(key, -1)) != actual:
            raise CheckError(f"JSON {key} disagrees with independent VTK result")
    if abs(total - expected) > tolerance:
        raise CheckError("positive cell areas do not close the expected fluid domain")
    claimed_length = float(report["actual_interface_length"])
    if abs(interface_length - claimed_length) > 1.0e-10 * max(1.0, interface_length):
        raise CheckError("JSON interface length disagrees with independent VTK result")
    return {
        "valid": True, **checks,
        "boundary_edge_count": sum(len(value) == 1 for value in owners.values()),
        "internal_edge_count": sum(len(value) == 2 for value in owners.values()),
        "layer_cell_count": sum(kind == 0 for kind in kinds),
        "remainder_cell_count": sum(kind != 0 for kind in kinds),
        "minimum_signed_area": min(areas), "area_sum": total,
        "area_error": total - expected, "interface_length": interface_length,
        "overlap_count": 0, "non_manifold_edge_count": 0, "issues": [],
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
