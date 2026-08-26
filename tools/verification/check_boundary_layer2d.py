#!/usr/bin/env python3
"""Independent H4-1 legacy-VTK/JSON topology checker.

This script does not import or link cartmesh2d. It re-parses the debug artifact,
checks positive quad area and edge incidence, and cross-checks summary metrics.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path


class CheckError(ValueError):
    pass


def signed_area(points: list[tuple[float, float]]) -> float:
    return 0.5 * sum(
        points[i][0] * points[(i + 1) % len(points)][1]
        - points[(i + 1) % len(points)][0] * points[i][1]
        for i in range(len(points))
    )


def orient(a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(
    point: tuple[float, float],
    segment: tuple[tuple[float, float], tuple[float, float]],
    epsilon: float,
) -> bool:
    a, b = segment
    return (
        abs(orient(a, b, point)) <= epsilon
        and min(a[0], b[0]) - epsilon <= point[0] <= max(a[0], b[0]) + epsilon
        and min(a[1], b[1]) - epsilon <= point[1] <= max(a[1], b[1]) + epsilon
    )


def segments_intersect(
    lhs: tuple[tuple[float, float], tuple[float, float]],
    rhs: tuple[tuple[float, float], tuple[float, float]],
    epsilon: float,
) -> bool:
    a, b = lhs
    c, d = rhs
    values = (orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b))
    if ((values[0] > epsilon and values[1] < -epsilon) or
        (values[0] < -epsilon and values[1] > epsilon)) and (
        (values[2] > epsilon and values[3] < -epsilon) or
        (values[2] < -epsilon and values[3] > epsilon)
    ):
        return True
    return (
        (abs(values[0]) <= epsilon and on_segment(c, lhs, epsilon))
        or (abs(values[1]) <= epsilon and on_segment(d, lhs, epsilon))
        or (abs(values[2]) <= epsilon and on_segment(a, rhs, epsilon))
        or (abs(values[3]) <= epsilon and on_segment(b, rhs, epsilon))
    )


def point_inside(point: tuple[float, float], polygon: list[tuple[float, float]], epsilon: float) -> bool:
    winding = 0
    for i, a in enumerate(polygon):
        b = polygon[(i + 1) % len(polygon)]
        if on_segment(point, (a, b), epsilon):
            return False
        if a[1] <= point[1] < b[1] and orient(a, b, point) > epsilon:
            winding += 1
        elif b[1] <= point[1] < a[1] and orient(a, b, point) < -epsilon:
            winding -= 1
    return winding != 0


def read_vtk(path: Path) -> tuple[list[tuple[float, float]], list[tuple[int, ...]]]:
    tokens = path.read_text(encoding="utf-8").split()
    try:
        points_at = tokens.index("POINTS")
        point_count = int(tokens[points_at + 1])
        point_start = points_at + 3
        raw_points = tokens[point_start : point_start + 3 * point_count]
        if len(raw_points) != 3 * point_count:
            raise CheckError("truncated POINTS section")
        points = [
            (float(raw_points[3 * i]), float(raw_points[3 * i + 1]))
            for i in range(point_count)
        ]
        if any(not math.isfinite(x) or not math.isfinite(y) for x, y in points):
            raise CheckError("non-finite point coordinate")

        cells_at = tokens.index("CELLS", point_start + 3 * point_count)
        cell_count = int(tokens[cells_at + 1])
        cell_start = cells_at + 3
        cells: list[tuple[int, ...]] = []
        cursor = cell_start
        for _ in range(cell_count):
            count = int(tokens[cursor])
            cursor += 1
            if count != 4:
                raise CheckError("H4-1 artifact contains a non-quad cell")
            cell = tuple(int(value) for value in tokens[cursor : cursor + count])
            cursor += count
            if len(cell) != 4 or len(set(cell)) != 4:
                raise CheckError("quad contains duplicate vertex IDs")
            if any(vertex < 0 or vertex >= point_count for vertex in cell):
                raise CheckError("quad references an invalid vertex ID")
            cells.append(cell)

        types_at = tokens.index("CELL_TYPES", cursor)
        type_count = int(tokens[types_at + 1])
        if type_count != cell_count:
            raise CheckError("CELL_TYPES count does not match CELLS")
        cell_types = [int(value) for value in tokens[types_at + 2 : types_at + 2 + type_count]]
        if len(cell_types) != cell_count or any(value != 9 for value in cell_types):
            raise CheckError("all H4-1 cells must use VTK_QUAD type 9")
    except (ValueError, IndexError) as exc:
        if isinstance(exc, CheckError):
            raise
        raise CheckError(f"invalid legacy VTK structure: {exc}") from exc
    return points, cells


def check(vtk_path: Path, report_path: Path) -> dict[str, object]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("layer_status") != "success":
        raise CheckError("JSON report is not a successful layer result")
    points, cells = read_vtk(vtk_path)
    coordinate_scale = max(
        1.0,
        max(abs(value) for point in points for value in point),
        max(point[0] for point in points) - min(point[0] for point in points),
        max(point[1] for point in points) - min(point[1] for point in points),
    )
    epsilon = 1.0e-11 * coordinate_scale
    areas: list[float] = []
    incidence: Counter[tuple[int, int]] = Counter()
    for cell in cells:
        area = signed_area([points[vertex] for vertex in cell])
        if not math.isfinite(area) or area <= 0.0:
            raise CheckError(f"quad has non-positive signed area: {area}")
        polygon = [points[vertex] for vertex in cell]
        if segments_intersect((polygon[0], polygon[1]), (polygon[2], polygon[3]), epsilon) or \
           segments_intersect((polygon[1], polygon[2]), (polygon[3], polygon[0]), epsilon):
            raise CheckError("quad is self-intersecting")
        areas.append(area)
        for i in range(4):
            incidence[tuple(sorted((cell[i], cell[(i + 1) % 4])))] += 1
    if any(count not in (1, 2) for count in incidence.values()):
        raise CheckError("edge incidence is not one (boundary) or two (internal)")

    strip_reports = report.get("strips")
    if not isinstance(strip_reports, list) or len(strip_reports) != 1:
        raise CheckError("independent checker currently requires one strip")
    strip = strip_reports[0]
    wall_count = int(strip["wall_vertex_count"])
    layer_count = int(strip["n_layers"])
    if len(points) != wall_count * (layer_count + 1):
        raise CheckError("ring-major vertex count is inconsistent")
    closed = bool(strip["closed"])
    ring_edge_count = wall_count if closed else wall_count - 1
    for ring in range(layer_count + 1):
        offset = ring * wall_count
        for i in range(ring_edge_count):
            first = (points[offset + i], points[offset + (i + 1) % wall_count])
            for j in range(i + 1, ring_edge_count):
                adjacent = abs(i - j) == 1 or (closed and {i, j} == {0, ring_edge_count - 1})
                if adjacent:
                    continue
                second = (points[offset + j], points[offset + (j + 1) % wall_count])
                if segments_intersect(first, second, epsilon):
                    raise CheckError(f"ring {ring} self-intersects")
    for vertex in range(wall_count):
        distances = []
        wall = points[vertex]
        for ring in range(layer_count + 1):
            point = points[ring * wall_count + vertex]
            distances.append(math.hypot(point[0] - wall[0], point[1] - wall[1]))
        if any(distances[i + 1] <= distances[i] + epsilon for i in range(layer_count)):
            raise CheckError("hair-edge distance is not strictly increasing")

    # Exact edge intersections plus strict vertex containment reject overlapping
    # quads that would pass a mere edge-incidence count.
    for lhs_id, lhs in enumerate(cells):
        lhs_edges = [tuple(sorted((lhs[i], lhs[(i + 1) % 4]))) for i in range(4)]
        lhs_polygon = [points[index] for index in lhs]
        for rhs in cells[lhs_id + 1 :]:
            rhs_edges = [tuple(sorted((rhs[i], rhs[(i + 1) % 4]))) for i in range(4)]
            shared_edges = set(lhs_edges) & set(rhs_edges)
            for lhs_edge in lhs_edges:
                for rhs_edge in rhs_edges:
                    if lhs_edge == rhs_edge:
                        continue
                    shared_vertices = set(lhs_edge) & set(rhs_edge)
                    left_segment = (points[lhs_edge[0]], points[lhs_edge[1]])
                    right_segment = (points[rhs_edge[0]], points[rhs_edge[1]])
                    if not segments_intersect(left_segment, right_segment, epsilon):
                        continue
                    if shared_vertices:
                        shared_point = points[next(iter(shared_vertices))]
                        other_lhs = points[next(value for value in lhs_edge if value not in shared_vertices)]
                        other_rhs = points[next(value for value in rhs_edge if value not in shared_vertices)]
                        if (not on_segment(other_lhs, right_segment, epsilon) and
                            not on_segment(other_rhs, left_segment, epsilon) and
                            on_segment(shared_point, left_segment, epsilon) and
                            on_segment(shared_point, right_segment, epsilon)):
                            continue
                    raise CheckError("non-shared cell edges intersect")
            if shared_edges:
                continue
            rhs_polygon = [points[index] for index in rhs]
            if any(point_inside(points[index], rhs_polygon, epsilon) for index in lhs) or \
               any(point_inside(points[index], lhs_polygon, epsilon) for index in rhs):
                raise CheckError("quad interiors overlap")
    boundary_edges = sum(count == 1 for count in incidence.values())
    internal_edges = sum(count == 2 for count in incidence.values())
    if int(report.get("layer_vertex_count", -1)) != len(points):
        raise CheckError("JSON layer_vertex_count disagrees with VTK")
    if int(report.get("layer_cell_count", -1)) != len(cells):
        raise CheckError("JSON layer_cell_count disagrees with VTK")
    for key, actual in (("min_cell_area", min(areas)), ("max_cell_area", max(areas))):
        claimed = float(report[key])
        tolerance = 1.0e-12 * max(1.0, abs(actual), abs(claimed))
        if abs(actual - claimed) > tolerance:
            raise CheckError(f"JSON {key} disagrees with independently computed VTK value")
    return {
        "valid": True,
        "point_count": len(points),
        "quad_count": len(cells),
        "boundary_edge_count": boundary_edges,
        "internal_edge_count": internal_edges,
        "min_signed_area": min(areas),
        "max_signed_area": max(areas),
        "issues": [],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vtk", type=Path)
    parser.add_argument("json", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    try:
        result = check(args.vtk, args.json)
    except (OSError, json.JSONDecodeError, CheckError) as exc:
        print(f"boundary-layer independent check: FAIL: {exc}")
        return 1
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        "boundary-layer independent check: PASS "
        f"points={result['point_count']} quads={result['quad_count']} "
        f"min_area={result['min_signed_area']:.17g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
