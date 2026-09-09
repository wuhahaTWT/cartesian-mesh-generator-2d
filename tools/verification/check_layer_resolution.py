#!/usr/bin/env python3
"""Independently verify hybrid boundary-layer resolution claims.

The resolution JSON is an untrusted claim.  This checker cross-references the
final solver CM2D polygons with ``hybrid_kind=0`` source VTK polygons, checks
column adjacency and wall ownership, and recomputes wall coverage and first
layer normal heights.  It is a topology/geometry checker only; it does not
claim Solver or Q1 acceptance.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_mesh_resolution import read_cm2d  # noqa: E402


ABS_TOL = 2.0e-10
REL_TOL = 2.0e-9
COORD_TOL = 2.0e-10


class CheckError(ValueError):
    pass


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=REL_TOL, abs_tol=ABS_TOL)


def vtk_arrays(path: Path):
    tokens = path.read_text(encoding="utf-8").split()
    try:
        p_at = tokens.index("POINTS")
        n_points = int(tokens[p_at + 1])
        start = p_at + 3
        points = [(float(tokens[start + 3 * i]), float(tokens[start + 3 * i + 1]))
                  for i in range(n_points)]
        c_at = tokens.index("CELLS", start + 3 * n_points)
        n_cells = int(tokens[c_at + 1])
        cursor = c_at + 3
        cells = []
        for _ in range(n_cells):
            count = int(tokens[cursor]); cursor += 1
            cells.append(tuple(int(v) for v in tokens[cursor:cursor + count]))
            cursor += count
        data_at = tokens.index("CELL_DATA", cursor)
        if int(tokens[data_at + 1]) != n_cells:
            raise CheckError("VTK CELL_DATA count differs from CELLS")
    except (ValueError, IndexError) as exc:
        raise CheckError(f"malformed VTK: {exc}") from exc

    arrays = {}
    cursor = data_at + 2
    while cursor < len(tokens):
        if tokens[cursor] != "SCALARS":
            cursor += 1
            continue
        if cursor + 3 >= len(tokens):
            raise CheckError("truncated VTK SCALARS declaration")
        name = tokens[cursor + 1]
        cursor += 4
        if tokens[cursor] == "LOOKUP_TABLE":
            cursor += 2
        values = tokens[cursor:cursor + n_cells]
        if len(values) != n_cells:
            raise CheckError(f"VTK scalar {name} count mismatch")
        arrays[name] = values
        cursor += n_cells
    if "hybrid_kind" not in arrays or "layer_index" not in arrays:
        raise CheckError("VTK must contain hybrid_kind and layer_index")
    kinds = [int(v) for v in arrays["hybrid_kind"]]
    layers = [int(v) for v in arrays["layer_index"]]
    if any(kind not in (0, 1, 2, 3, 4) for kind in kinds):
        raise CheckError("VTK contains invalid hybrid_kind")
    return points, cells, kinds, layers


def quantized(point: tuple[float, float], eps: float) -> tuple[int, int]:
    return (round(point[0] / eps), round(point[1] / eps))


def source_id_mapper(points: list[tuple[float, float]], eps: float):
    buckets: dict[tuple[int, int], list[int]] = defaultdict(list)
    for index, point in enumerate(points):
        buckets[quantized(point, eps)].append(index)

    def map_point(point: tuple[float, float]) -> int:
        key = quantized(point, eps)
        candidates = []
        for ix in range(key[0] - 1, key[0] + 2):
            for iy in range(key[1] - 1, key[1] + 2):
                candidates.extend(buckets.get((ix, iy), []))
        candidates = sorted(set(candidates),
                            key=lambda i: math.dist(point, points[i]))
        if not candidates or math.dist(point, points[candidates[0]]) > eps:
            raise CheckError(f"final vertex {point} is absent from source VTK coordinates")
        if len(candidates) > 1 and math.isclose(
                math.dist(point, points[candidates[0]]),
                math.dist(point, points[candidates[1]]),
                rel_tol=0.0, abs_tol=eps * 1.0e-3):
            raise CheckError(f"source coordinate identity is ambiguous for {point}")
        return candidates[0]

    return map_point


def cross(a, b, c) -> float:
    return ((b[0] - a[0]) * (c[1] - a[1])
            - (b[1] - a[1]) * (c[0] - a[0]))


def between(a, b, c, eps: float) -> bool:
    edge_length = math.dist(a, b)
    if edge_length <= eps:
        return math.dist(a, c) <= eps
    return (abs(cross(a, b, c)) / edge_length <= eps
            and min(a[0], b[0]) - eps <= c[0] <= max(a[0], b[0]) + eps
            and min(a[1], b[1]) - eps <= c[1] <= max(a[1], b[1]) + eps)


def remove_collinear(ids: list[int], points: list[tuple[float, float]], eps: float) -> list[int]:
    result = list(ids)
    changed = True
    while changed and len(result) > 3:
        changed = False
        for index in range(len(result)):
            a = points[result[index - 1]]
            b = points[result[index]]
            c = points[result[(index + 1) % len(result)]]
            if between(a, c, b, eps):
                result.pop(index)
                changed = True
                break
    if len(result) < 3:
        raise CheckError("polygon collapsed while removing collinear vertices")
    return result


def cycle_key(ids: list[int]) -> tuple[int, ...]:
    def rotations(values):
        return [tuple(values[i:] + values[:i]) for i in range(len(values))]
    forward = rotations(ids)
    backward = rotations(list(reversed(ids)))
    return min(forward + backward)


def polygon_key(points, polygon, mapper, eps, canonical_points=None):
    ids = [mapper(points[index]) for index in polygon]
    return cycle_key(remove_collinear(ids, canonical_points or points, eps))


def vector(a, b):
    return (b[0] - a[0], b[1] - a[1])


def normal_extent(cell_points, wall_a, wall_b) -> float:
    edge = vector(wall_a, wall_b)
    length = math.hypot(*edge)
    if length <= 0.0:
        raise CheckError("zero-length wall segment")
    normal = (-edge[1] / length, edge[0] / length)
    projections = [(p[0] - wall_a[0]) * normal[0]
                   + (p[1] - wall_a[1]) * normal[1] for p in cell_points]
    return max(projections) - min(projections)


def measure(source_vtk: Path, final_cm2d: Path, report_path: Path) -> dict:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    boundary = report.get("boundary_layers")
    if not isinstance(boundary, dict) or boundary.get("status") != "evaluated":
        raise CheckError("resolution report lacks evaluated boundary_layers")
    columns = boundary.get("columns")
    if not isinstance(columns, list) or not columns:
        raise CheckError("resolution report has no boundary-layer columns")
    reference = float(report["reference_length"])
    if not math.isfinite(reference) or reference <= 0.0:
        raise CheckError("reference_length must be finite and positive")
    requested_height = report.get("requested", {}).get("first_layer_h_over_reference")
    if (requested_height is None or not math.isfinite(float(requested_height))
            or float(requested_height) <= 0.0):
        raise CheckError("positive requested first-layer height is required")
    requested_height = float(requested_height)

    source_points_raw, source_cells, kinds, layer_indices = vtk_arrays(source_vtk)
    solver_points_raw, solver_edges, solver_cells = read_cm2d(final_cm2d)
    origin = source_points_raw[0]
    source_points = [((x - origin[0]) / reference, (y - origin[1]) / reference)
                     for x, y in source_points_raw]
    solver_points = [((x - origin[0]) / reference, (y - origin[1]) / reference)
                     for x, y in solver_points_raw]
    eps = COORD_TOL
    mapper = source_id_mapper(source_points, eps)

    source_keys = {}
    layer_source = set()
    for cid, (cell, kind) in enumerate(zip(source_cells, kinds)):
        if kind == 0:
            key = polygon_key(source_points, cell, mapper, eps)
            if key in source_keys:
                raise CheckError(f"duplicate source layer polygon key for cells {source_keys[key]} and {cid}")
            source_keys[key] = cid
            layer_source.add(cid)

    solver_key_by_id = {}

    def solver_key(cell_id: int):
        if cell_id not in solver_key_by_id:
            solver_key_by_id[cell_id] = polygon_key(
                solver_points, solver_cells[cell_id]["vertices"], mapper, eps, source_points)
        return solver_key_by_id[cell_id]

    # Build solver edge incidence independently from the CM2D AUDIT record.
    pair_edges: dict[tuple[int, int], list[int]] = defaultdict(list)
    cell_edges: dict[int, list[int]] = defaultdict(list)
    for eid, edge in enumerate(solver_edges):
        owner = edge["owner"]
        neighbour = edge["neighbour"]
        if owner < 0 or owner >= len(solver_cells):
            raise CheckError(f"edge {eid} has invalid owner")
        cell_edges[owner].append(eid)
        if neighbour >= 0:
            if neighbour >= len(solver_cells):
                raise CheckError(f"edge {eid} has invalid neighbour")
            cell_edges[neighbour].append(eid)
            pair_edges[tuple(sorted((owner, neighbour)))].append(eid)

    # Source physical wall edges are the boundary edges of kind=0 polygons
    # after all source cells are included; source coordinates remain canonical.
    source_edge_owners: dict[tuple[int, int], list[int]] = defaultdict(list)
    for cid, cell in enumerate(source_cells):
        ids = [mapper(source_points[index]) for index in cell]
        for index, first in enumerate(ids):
            edge = tuple(sorted((first, ids[(index + 1) % len(ids)])))
            source_edge_owners[edge].append(cid)
    wall_edges = {}
    for edge, owners in source_edge_owners.items():
        if len(owners) != 1 or kinds[owners[0]] != 0:
            continue
        a, b = source_points[edge[0]], source_points[edge[1]]
        wall_edges[edge] = math.dist(a, b)
    if not wall_edges:
        raise CheckError("no source boundary-layer physical wall edges found")
    wall_total = math.fsum(wall_edges.values())

    claimed_solver_ids: set[int] = set()
    claimed_source_ids: set[int] = set()
    covered_first_edges: set[tuple[int, int]] = set()
    covered_full_edges: set[tuple[int, int]] = set()
    declared_wall_edges: set[tuple[int, int]] = set()
    max_first_height = 0.0
    height_exceeded_length = 0.0
    height_by_wall_edge: dict[tuple[int, int], float] = {}
    column_results = []
    for column in columns:
        required = int(column["requested_layers"])
        retained = int(column["retained_layers"])
        ids = [int(value) for value in column["solver_cell_ids"]]
        if required < 1 or retained != len(ids) or not 0 <= retained <= required:
            raise CheckError("invalid requested/retained column layer counts")
        if len(set(ids)) != len(ids) or any(cid < 0 or cid >= len(solver_cells) for cid in ids):
            raise CheckError("column contains duplicate or invalid solver cell ID")
        endpoint_values = column.get("wall_segment_endpoints")
        if not isinstance(endpoint_values, list) or len(endpoint_values) != 2:
            raise CheckError("column lacks two wall_segment_endpoints")
        endpoint_points = []
        for endpoint in endpoint_values:
            raw = tuple(map(float, endpoint))
            if len(raw) != 2 or not all(math.isfinite(value) for value in raw):
                raise CheckError("column wall endpoint is not finite 2D data")
            endpoint_points.append(((raw[0] - origin[0]) / reference,
                                    (raw[1] - origin[1]) / reference))
        endpoint_ids = tuple(sorted((mapper(endpoint_points[0]),
                                     mapper(endpoint_points[1]))))
        if endpoint_ids not in wall_edges:
            raise CheckError("column wall endpoints do not identify a source wall edge")
        if endpoint_ids in declared_wall_edges:
            raise CheckError("multiple columns declare the same wall segment")
        declared_wall_edges.add(endpoint_ids)
        wall_length = wall_edges[endpoint_ids]
        if retained:
            covered_first_edges.add(endpoint_ids)
        if retained == required:
            covered_full_edges.add(endpoint_ids)
        heights = []
        source_column_ids = []
        for layer, solver_id in enumerate(ids):
            source_id = source_keys.get(solver_key(solver_id))
            if source_id not in layer_source or layer_indices[source_id] != layer:
                raise CheckError("solver column polygon does not match ordered source layer polygon")
            if source_id in claimed_source_ids or solver_id in claimed_solver_ids:
                raise CheckError("layer polygon is claimed by multiple columns")
            claimed_source_ids.add(source_id); claimed_solver_ids.add(solver_id)
            source_column_ids.append(source_id)
            if layer == 0:
                solver_cell = solver_cells[solver_id]
                wall_intervals = []
                wall_a, wall_b = source_points[endpoint_ids[0]], source_points[endpoint_ids[1]]
                wall_vector = vector(wall_a, wall_b)
                wall_norm = math.hypot(*wall_vector)
                for edge_id in cell_edges[solver_id]:
                    edge = solver_edges[edge_id]
                    if edge["neighbour"] < 0 and edge["patch"] == 1:
                        actual_a = solver_points[edge["v0"]]
                        actual_b = solver_points[edge["v1"]]
                        if wall_norm <= eps:
                            raise CheckError("declared wall segment has zero length")
                        if (abs(cross(wall_a, wall_b, actual_a)) / wall_norm > eps or
                                abs(cross(wall_a, wall_b, actual_b)) / wall_norm > eps):
                            raise CheckError("first-layer embedded edge is outside declared wall segment")
                        def projection(point):
                            return ((point[0] - wall_a[0]) * wall_vector[0]
                                    + (point[1] - wall_a[1]) * wall_vector[1]) / (wall_norm * wall_norm)
                        lo = min(projection(actual_a), projection(actual_b)) * wall_norm
                        hi = max(projection(actual_a), projection(actual_b)) * wall_norm
                        if lo < -eps or hi > wall_length + eps:
                            raise CheckError("first-layer embedded edge does not lie within declared wall segment")
                        wall_intervals.append((max(0.0, lo), min(wall_length, hi)))
                if not wall_intervals:
                    raise CheckError("first-layer solver cell lacks an embedded wall edge")
                wall_intervals.sort()
                covered_to = 0.0
                for lo, hi in wall_intervals:
                    if lo > covered_to + eps:
                        raise CheckError("first-layer embedded edges leave a gap on declared wall segment")
                    covered_to = max(covered_to, hi)
                if covered_to < wall_length - eps:
                    raise CheckError("first-layer embedded edges do not cover declared wall segment")
                points = [solver_points[index] for index in solver_cell["vertices"]]
                height = normal_extent(points, wall_a, wall_b)
                heights.append(height)
                max_first_height = max(max_first_height, height)
                height_by_wall_edge[endpoint_ids] = height
        for first, second in zip(ids, ids[1:]):
            if len(pair_edges.get(tuple(sorted((first, second))), [])) != 1:
                raise CheckError("column layers are not connected by one internal owner/neighbour edge")
        claimed_min = float(column["first_layer_normal_height_min_over_reference"])
        claimed_max = float(column["first_layer_normal_height_max_over_reference"])
        if heights and (not close(min(heights), claimed_min) or not close(max(heights), claimed_max)):
            raise CheckError("column first-layer normal height differs from independent measurement")
        column_results.append({"strip_id": int(column["strip_id"]),
                               "wall_segment": int(column["wall_segment"]),
                               "requested_layers": required,
                               "retained_layers": retained,
                               "wall_length": wall_length,
                               "source_cell_ids": source_column_ids,
                               "solver_cell_ids": ids,
                               "first_layer_height_min": min(heights) if heights else None,
                               "first_layer_height_max": max(heights) if heights else None})

    if len(layer_source) != int(boundary["constructed_cells"]):
        raise CheckError("source layer polygon count differs from constructed_cells")
    if len(claimed_source_ids) != int(boundary["retained_cells"]):
        raise CheckError("claimed layer polygon count differs from retained_cells")
    if sum(int(column["requested_layers"]) for column in columns) != int(boundary["requested_cells"]):
        raise CheckError("column requested layer count differs from requested_cells")
    if int(boundary.get("source_mismatches", -1)) != 0:
        raise CheckError("native report claims source polygon mismatches")
    if int(boundary.get("continuity_mismatches", -1)) != 0:
        raise CheckError("native report claims column continuity mismatches")
    first_fraction = math.fsum(wall_edges[e] for e in covered_first_edges) / wall_total
    full_fraction = math.fsum(wall_edges[e] for e in covered_full_edges) / wall_total
    for edge in covered_first_edges:
        height = height_by_wall_edge[edge]
        if height > requested_height and not close(height, requested_height):
            height_exceeded_length += wall_edges[edge]
    exceeded_fraction = height_exceeded_length / wall_total
    if not close(first_fraction, float(boundary["first_layer_wall_length_fraction"])):
        raise CheckError("first-layer wall coverage differs from report")
    if not close(full_fraction, float(boundary["full_requested_layers_wall_length_fraction"])):
        raise CheckError("full-layer wall coverage differs from report")
    if not close(exceeded_fraction, float(boundary["first_layer_height_exceedance_wall_length_fraction"])):
        raise CheckError("first-layer height exceedance differs from report")
    return {"valid": True, "source_vtk": str(source_vtk), "final_cm2d": str(final_cm2d),
            "report": str(report_path),
            "source_layer_polygon_count": len(layer_source),
            "claimed_layer_polygon_count": len(claimed_source_ids),
            "source_wall_length": wall_total,
            "first_layer_wall_length_fraction": first_fraction,
            "full_requested_layers_wall_length_fraction": full_fraction,
            "first_layer_height_exceedance_wall_length_fraction": exceeded_fraction,
            "first_layer_normal_height_max_over_reference": max_first_height,
            "columns": column_results, "issues": []}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cm2d", type=Path)
    parser.add_argument("--source-vtk", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = measure(args.source_vtk, args.cm2d, args.report)
    except (CheckError, OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
        result = {"valid": False, "source_vtk": str(args.source_vtk),
                  "final_cm2d": str(args.cm2d), "report": str(args.report),
                  "issues": [str(exc)]}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
