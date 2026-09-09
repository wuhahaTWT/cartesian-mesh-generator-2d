#!/usr/bin/env python3
"""Independently verify a final CM2D resolution report.

The native report is intentionally treated as an untrusted claim.  This reader
recomputes cell areas and embedded-wall owner extents from the serialized solver
topology, then compares the measured distributions and exceedance fraction with
the companion ``.resolution.json``.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


ABS_TOL = 2.0e-12
REL_TOL = 2.0e-10
NATIVE_ABS_TOL = 1.0e-12
NATIVE_REL_TOL = 1.0e-10


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=REL_TOL, abs_tol=ABS_TOL)


def native_nearly_equal(a: float, b: float) -> bool:
    """Mirror TolerancePolicy2D::nearlyEqual(a, b) with default magnitude."""
    difference = abs(a - b)
    return difference <= NATIVE_ABS_TOL + NATIVE_REL_TOL * max(1.0, difference)


def percentile(values: list[float], p: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    index = p * (len(ordered) - 1)
    lo = math.floor(index)
    hi = math.ceil(index)
    return ordered[lo] + (index - lo) * (ordered[hi] - ordered[lo])


def read_cm2d(path: Path) -> tuple[list[tuple[float, float]], list[dict], list[dict]]:
    tokens = path.read_text(encoding="utf-8-sig").split()
    pos = 0

    def take() -> str:
        nonlocal pos
        if pos >= len(tokens):
            raise ValueError(f"{path}: truncated CM2D record")
        value = tokens[pos]
        pos += 1
        return value

    if (take(), take(), take()) != ("CM2D", "1", "VERTICES"):
        raise ValueError(f"{path}: unsupported CM2D header")
    vertices = []
    for expected in range(int(take())):
        if int(take()) != expected:
            raise ValueError(f"{path}: non-contiguous vertex ids")
        vertices.append((float(take()), float(take())))
    if take() != "EDGES":
        raise ValueError(f"{path}: missing EDGES")
    edges = []
    for expected in range(int(take())):
        values = [int(take()) for _ in range(6)]
        if values[0] != expected:
            raise ValueError(f"{path}: non-contiguous edge ids")
        edges.append({"id": values[0], "v0": values[1], "v1": values[2],
                      "owner": values[3], "neighbour": values[4],
                      "patch": values[5]})
    for edge in edges:
        if not (0 <= edge["v0"] < len(vertices) and 0 <= edge["v1"] < len(vertices)):
            raise ValueError(f"{path}: edge {edge['id']} references invalid vertex")
    if take() != "CELLS":
        raise ValueError(f"{path}: missing CELLS")
    cells = []
    for expected in range(int(take())):
        cid = int(take())
        source_id, source_key = int(take()), int(take())
        declared_area, nvertices = float(take()), int(take())
        vertices_ids = [int(take()) for _ in range(nvertices)]
        nedges = int(take())
        edge_ids = [int(take()) for _ in range(nedges)]
        if cid != expected or nvertices < 3 or nedges != nvertices:
            raise ValueError(f"{path}: invalid cell record {cid}")
        if any(index < 0 or index >= len(vertices) for index in vertices_ids):
            raise ValueError(f"{path}: cell {cid} references invalid vertex")
        if any(index < 0 or index >= len(edges) for index in edge_ids):
            raise ValueError(f"{path}: cell {cid} references invalid edge")
        cells.append({"id": cid, "source_id": source_id, "source_key": source_key,
                      "declared_area": declared_area, "vertices": vertices_ids,
                      "edges": edge_ids})
    if take() != "AUDIT":
        raise ValueError(f"{path}: missing AUDIT")
    audit = [int(take()) for _ in range(7)]
    if take() != "END" or pos != len(tokens):
        raise ValueError(f"{path}: trailing or missing CM2D data")
    # The AUDIT record is metadata, not an acceptance oracle.  All connectivity,
    # geometry and statistics used below are reconstructed independently.
    return vertices, edges, cells


def polygon_area(points: list[tuple[float, float]]) -> float:
    origin = points[0]
    shifted = [(x - origin[0], y - origin[1]) for x, y in points]
    return abs(math.fsum(a[0] * b[1] - b[0] * a[1]
                         for a, b in zip(shifted, shifted[1:] + shifted[:1]))) * 0.5


def distribution(values: list[float]) -> dict[str, float | None]:
    return {"min": min(values, default=None), "p50": percentile(values, 0.50),
            "p95": percentile(values, 0.95), "max": max(values, default=None)}


def measure(cm2d: Path, report_path: Path) -> dict:
    vertices, edges, cells = read_cm2d(cm2d)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    reference = float(report["reference_length"])
    if not math.isfinite(reference) or reference <= 0.0:
        raise ValueError("report reference_length must be finite and positive")

    areas = []
    for cell in cells:
        points = [vertices[index] for index in cell["vertices"]]
        area = polygon_area(points)
        if not math.isfinite(area) or area <= 0.0:
            raise ValueError(f"cell {cell['id']} has non-positive area")
        areas.append(math.sqrt(area) / reference)

    wall_edges = [edge for edge in edges
                  if edge["neighbour"] < 0 and edge["patch"] == 1]
    edge_lengths, tangential, normal = [], [], []
    wall_length = exceeded_length = 0.0
    requested = report.get("requested", {}).get("wall_h_over_reference")
    for edge in wall_edges:
        a, b = vertices[edge["v0"]], vertices[edge["v1"]]
        dx, dy = b[0] - a[0], b[1] - a[1]
        length = math.hypot(dx, dy)
        if not math.isfinite(length) or length <= 0.0:
            raise ValueError(f"embedded edge {edge['id']} has invalid length")
        tx, ty = dx / length, dy / length
        nx, ny = -ty, tx
        projections = []
        for vertex in cells[edge["owner"]]["vertices"]:
            x, y = vertices[vertex]
            rx, ry = x - a[0], y - a[1]
            projections.append((rx * tx + ry * ty, rx * nx + ry * ny))
        t_extent = max(p[0] for p in projections) - min(p[0] for p in projections)
        n_extent = max(p[1] for p in projections) - min(p[1] for p in projections)
        edge_lengths.append(length / reference)
        tangential.append(t_extent / reference)
        normal.append(n_extent / reference)
        wall_length += length
        if requested is not None:
            extent = t_extent / reference
            target = float(requested)
            if extent > target and not native_nearly_equal(extent, target):
                exceeded_length += length

    actual = {"sqrt_area_over_reference": distribution(areas),
              "wall_edge_length_over_reference": distribution(edge_lengths),
              "wall_owner_tangential_extent_over_reference": distribution(tangential),
              "wall_owner_normal_extent_over_reference": distribution(normal)}
    measured_fraction = ((exceeded_length / wall_length)
                        if requested is not None and wall_length else None)
    measured = {"cell_count": len(cells), "wall_edge_count": len(wall_edges),
                "wall_length_over_reference": wall_length / reference,
                "wall_owner_tangential_exceedance_length_fraction": measured_fraction,
                "actual": actual}
    issues = []
    if report.get("cell_count") != measured["cell_count"]:
        issues.append("cell_count differs from native report")
    for name, values in actual.items():
        claimed = report.get("actual", {}).get(name)
        if not isinstance(claimed, dict):
            issues.append(f"missing actual distribution: {name}")
            continue
        for key, value in values.items():
            if value is None:
                if claimed.get(key) is not None:
                    issues.append(f"actual.{name}.{key} differs from native report")
            elif claimed.get(key) is None or not close(float(value), float(claimed.get(key))):
                issues.append(f"actual.{name}.{key} differs from native report")
    for name, value in measured.items():
        if name in ("actual", "cell_count", "wall_edge_count"):
            continue
        claimed = report.get(name)
        if value is None:
            if claimed is not None:
                issues.append(f"{name} expected null")
        elif claimed is None or not close(float(value), float(claimed)):
            issues.append(f"{name} differs from native report")
    return {"valid": not issues, "cm2d": str(cm2d), "report": str(report_path),
            "reference_length": reference, "measured": measured,
            "tolerance": {"absolute": ABS_TOL, "relative": REL_TOL},
            "issues": issues}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cm2d", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = measure(args.cm2d, args.report)
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
        result = {"valid": False, "cm2d": str(args.cm2d),
                  "report": str(args.report), "issues": [str(exc)]}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
