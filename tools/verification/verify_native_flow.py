#!/usr/bin/env python3
"""Generate and independently verify the native steady 2-D flow milestone.

The native JSON is treated as a report, never as the acceptance oracle.  This
tool reads the final CM2D polygons, cell field CSV and owner-oriented face flux
CSV independently.  It reconstructs incidence and finite-volume continuity,
then applies case-specific checks:

* channel: plane Poiseuille velocity, pressure gradient and flow rate;
* cavity: Re=100 centreline values from Ghia, Ghia & Shin (1982);
* external: finite conservative steady output and recorded force values.

The external case deliberately does not compare against the DFG cylinder drag:
the bundled circle uses a symmetric slip far field, whereas the DFG 2D-1 case
uses an offset cylinder in a no-slip channel.  Similar Reynolds number alone
does not make their force coefficients interchangeable.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = REPO / "outputs/native-flow"
DEFAULT_CIRCLE = REPO / "examples/acceptance/circle.xy"

# Primary references used for the acceptance definitions.  The Ghia table is
# Table I/II of the cited paper.  MIT gives the plane-Poiseuille closed form.
# The DFG page records the Re=20 force convention and, just as importantly, the
# geometry/BC differences that prevent using its drag as a gate here.
REFERENCES = {
    "poiseuille": "https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/resources/mit2_25f13_couet_and_pois/",
    "cavity": "https://doi.org/10.1016/0021-9991(82)90058-4",
    "dfg_re20": "https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html",
    "open_cylinder_re20": "https://doi.org/10.1017/S0022112070001428",
}

# Re=100, unit square, unit lid speed.  Coordinates and velocities are
# non-dimensional.  Wall endpoints are retained as provenance but omitted from
# interpolation error because cell-centred output does not sample the wall.
GHIA_U = (
    (1.0000, 1.00000), (0.9766, 0.84123), (0.9688, 0.78871),
    (0.9609, 0.73722), (0.9531, 0.68717), (0.8516, 0.23151),
    (0.7344, 0.00332), (0.6172, -0.13641), (0.5000, -0.20581),
    (0.4531, -0.21090), (0.2813, -0.15662), (0.1719, -0.10150),
    (0.1016, -0.06434), (0.0703, -0.04775), (0.0625, -0.04192),
    (0.0547, -0.03717), (0.0000, 0.00000),
)
GHIA_V = (
    (1.0000, 0.00000), (0.9688, -0.05906), (0.9609, -0.07391),
    (0.9531, -0.08864), (0.9453, -0.10313), (0.9063, -0.16914),
    (0.8594, -0.22445), (0.8047, -0.24533), (0.5000, 0.05454),
    (0.2344, 0.17527), (0.2266, 0.17507), (0.1563, 0.16077),
    (0.0938, 0.12317), (0.0781, 0.10890), (0.0703, 0.10091),
    (0.0625, 0.09233), (0.0000, 0.00000),
)


class VerificationError(ValueError):
    pass


@dataclass(frozen=True)
class Edge:
    id: int
    v0: int
    v1: int
    owner: int
    neighbour: int
    patch: int


@dataclass(frozen=True)
class Cell:
    id: int
    stored_area: float
    vertices: tuple[int, ...]
    edges: tuple[int, ...]


@dataclass(frozen=True)
class Mesh:
    path: Path
    vertices: tuple[tuple[float, float], ...]
    edges: tuple[Edge, ...]
    cells: tuple[Cell, ...]
    audit: tuple[int, ...]


@dataclass(frozen=True)
class Measurement:
    areas: tuple[float, ...]
    centroids: tuple[tuple[float, float], ...]
    bounds: tuple[float, float, float, float]
    total_area: float
    characteristic_h: float
    issues: tuple[str, ...]


@dataclass(frozen=True)
class Request:
    label: str
    case: str
    level: int | None
    padding: float | None
    nu: float
    speed: float


GENERATED = (
    Request("channel-l5", "channel", 5, 1.0 / 30.0, 0.01, 1.0),
    Request("channel-l6", "channel", 6, 1.0 / 62.0, 0.01, 1.0),
    Request("cavity-l4", "cavity", 4, 1.0 / 14.0, 0.01, 1.0),
    Request("cavity-l5", "cavity", 5, 1.0 / 30.0, 0.01, 1.0),
    Request("circle-re20", "external", None, None, 0.1, 1.0),
)


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise VerificationError(f"{label} is boolean, not numeric")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise VerificationError(f"{label} is not numeric: {value!r}") from exc
    if not math.isfinite(number):
        raise VerificationError(f"{label} is not finite: {number!r}")
    return number


def integer(value: Any, label: str) -> int:
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise VerificationError(f"{label} is not an integer: {value!r}") from exc
    if str(value).strip() not in {str(number), f"+{number}"}:
        raise VerificationError(f"{label} is not an exact integer: {value!r}")
    return number


def close(a: float, b: float, absolute: float, relative: float) -> bool:
    return abs(a - b) <= absolute + relative * max(abs(a), abs(b))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_final_solver_path(path: Path) -> None:
    text = str(path)
    if not text.endswith(".solver.cm2d") or text.endswith(".failed.solver.cm2d"):
        raise VerificationError(f"requires final *.solver.cm2d, got {path}")


def read_cm2d(path: Path) -> Mesh:
    tokens = path.read_text(encoding="utf-8-sig").split()
    position = 0

    def take() -> str:
        nonlocal position
        if position >= len(tokens):
            raise VerificationError(f"{path}: truncated CM2D record")
        token = tokens[position]
        position += 1
        return token

    if (take(), take(), take()) != ("CM2D", "1", "VERTICES"):
        raise VerificationError(f"{path}: unsupported CM2D header")
    vertex_count = integer(take(), "vertex count")
    vertices: list[tuple[float, float]] = []
    for expected in range(vertex_count):
        vertex_id = integer(take(), f"vertex {expected} id")
        if vertex_id != expected:
            raise VerificationError(f"{path}: non-contiguous vertex ids")
        vertices.append((finite(take(), f"vertex {expected} x"),
                         finite(take(), f"vertex {expected} y")))
    if take() != "EDGES":
        raise VerificationError(f"{path}: missing EDGES")
    edge_count = integer(take(), "edge count")
    edges: list[Edge] = []
    for expected in range(edge_count):
        values = [integer(take(), f"edge {expected} field") for _ in range(6)]
        edge = Edge(*values)
        if edge.id != expected or not (0 <= edge.v0 < vertex_count and 0 <= edge.v1 < vertex_count):
            raise VerificationError(f"{path}: invalid edge {expected}")
        edges.append(edge)
    if take() != "CELLS":
        raise VerificationError(f"{path}: missing CELLS")
    cell_count = integer(take(), "cell count")
    cells: list[Cell] = []
    for expected in range(cell_count):
        cell_id = integer(take(), f"cell {expected} id")
        _source_id, _source_key = integer(take(), "source id"), integer(take(), "source key")
        area = finite(take(), f"cell {expected} area")
        nv = integer(take(), f"cell {expected} vertex count")
        vertex_ids = tuple(integer(take(), "cell vertex") for _ in range(nv))
        ne = integer(take(), f"cell {expected} edge count")
        edge_ids = tuple(integer(take(), "cell edge") for _ in range(ne))
        if cell_id != expected or nv < 3 or ne != nv:
            raise VerificationError(f"{path}: invalid cell {expected} loop")
        if any(v < 0 or v >= vertex_count for v in vertex_ids) or any(e < 0 or e >= edge_count for e in edge_ids):
            raise VerificationError(f"{path}: invalid reference in cell {expected}")
        cells.append(Cell(cell_id, area, vertex_ids, edge_ids))
    if take() != "AUDIT":
        raise VerificationError(f"{path}: missing AUDIT")
    audit = tuple(integer(take(), "audit value") for _ in range(7))
    if take() != "END" or position != len(tokens):
        raise VerificationError(f"{path}: trailing or missing CM2D data")
    return Mesh(path.resolve(), tuple(vertices), tuple(edges), tuple(cells), audit)


def polygon(points: list[tuple[float, float]]) -> tuple[float, tuple[float, float]]:
    pairs = list(zip(points, points[1:] + points[:1]))
    cross = [a[0] * b[1] - b[0] * a[1] for a, b in pairs]
    twice = math.fsum(cross)
    if not math.isfinite(twice) or twice <= 0.0:
        raise VerificationError("cell polygon is not finite counter-clockwise positive area")
    cx = math.fsum((a[0] + b[0]) * q for (a, b), q in zip(pairs, cross)) / (3.0 * twice)
    cy = math.fsum((a[1] + b[1]) * q for (a, b), q in zip(pairs, cross)) / (3.0 * twice)
    return 0.5 * twice, (cx, cy)


def measure(mesh: Mesh, absolute: float, relative: float) -> Measurement:
    issues: list[str] = []
    incidence: list[list[int]] = [[] for _ in mesh.edges]
    areas: list[float] = []
    centres: list[tuple[float, float]] = []
    for cell in mesh.cells:
        try:
            area, centre = polygon([mesh.vertices[v] for v in cell.vertices])
        except VerificationError as exc:
            issues.append(f"cell {cell.id}: {exc}")
            area, centre = math.nan, (math.nan, math.nan)
        areas.append(area)
        centres.append(centre)
        if not (cell.stored_area > 0.0 and math.isfinite(cell.stored_area)):
            issues.append(f"cell {cell.id}: stored area is not positive finite")
        elif math.isfinite(area) and not close(area, cell.stored_area, absolute, relative):
            issues.append(f"cell {cell.id}: stored and reconstructed areas differ")
        if len(set(cell.vertices)) != len(cell.vertices) or len(set(cell.edges)) != len(cell.edges):
            issues.append(f"cell {cell.id}: repeated loop entity")
        for local, edge_id in enumerate(cell.edges):
            edge = mesh.edges[edge_id]
            segment = {cell.vertices[local], cell.vertices[(local + 1) % len(cell.vertices)]}
            if segment != {edge.v0, edge.v1}:
                issues.append(f"cell {cell.id}: edge {edge_id} does not match loop")
            incidence[edge_id].append(cell.id)
    for edge, adjacent in zip(mesh.edges, incidence):
        expected = [edge.owner] + ([] if edge.neighbour < 0 else [edge.neighbour])
        if len(adjacent) != len(set(adjacent)) or sorted(adjacent) != sorted(expected):
            issues.append(f"edge {edge.id}: incidence {adjacent} differs from owner/neighbour {expected}")
        if edge.neighbour < 0 and edge.patch not in (1, 2):
            issues.append(f"edge {edge.id}: boundary patch is {edge.patch}")
        if edge.neighbour >= 0 and edge.patch != 0:
            issues.append(f"edge {edge.id}: internal patch is {edge.patch}")
    if any(mesh.audit):
        issues.append(f"CM2D audit is nonzero: {list(mesh.audit)}")
    xs, ys = [p[0] for p in mesh.vertices], [p[1] for p in mesh.vertices]
    bounds = (min(xs), min(ys), max(xs), max(ys)) if xs else (math.nan,) * 4
    total = math.fsum(a for a in areas if math.isfinite(a))
    h = math.sqrt(total / len(areas)) if areas and total > 0.0 else math.nan
    return Measurement(tuple(areas), tuple(centres), bounds, total, h, tuple(issues))


def load_csv(path: Path, columns: Iterable[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        missing = [name for name in columns if name not in fields]
        if missing:
            raise VerificationError(f"{path}: missing columns {missing}; got {fields}")
        rows = list(reader)
    if not rows:
        raise VerificationError(f"{path}: no data rows")
    return rows


def read_cells(path: Path, mesh: Mesh, measured: Measurement) -> list[dict[str, float]]:
    rows = load_csv(path, ("cell", "x", "y", "area", "u", "v", "p", "speed"))
    values: dict[int, dict[str, float]] = {}
    for line, row in enumerate(rows, 2):
        cell = integer(row["cell"], f"{path}:{line} cell")
        if cell in values or not (0 <= cell < len(mesh.cells)):
            raise VerificationError(f"{path}:{line}: duplicate/out-of-range cell {cell}")
        item = {name: finite(row[name], f"{path}:{line} {name}") for name in ("x", "y", "area", "u", "v", "p", "speed")}
        if not close(item["x"], measured.centroids[cell][0], 1e-11, 1e-9) or not close(item["y"], measured.centroids[cell][1], 1e-11, 1e-9):
            raise VerificationError(f"{path}:{line}: cell centre differs from CM2D")
        if not close(item["area"], measured.areas[cell], 1e-11, 1e-9):
            raise VerificationError(f"{path}:{line}: cell area differs from CM2D")
        if not close(item["speed"], math.hypot(item["u"], item["v"]), 1e-12, 1e-9):
            raise VerificationError(f"{path}:{line}: speed differs from hypot(u,v)")
        values[cell] = item
    if len(values) != len(mesh.cells):
        raise VerificationError(f"{path}: got {len(values)} cells, expected {len(mesh.cells)}")
    return [values[i] for i in range(len(mesh.cells))]


def read_faces(path: Path, mesh: Mesh) -> list[float]:
    rows = load_csv(path, ("face", "owner", "neighbour", "flux"))
    fluxes: dict[int, float] = {}
    for line, row in enumerate(rows, 2):
        face = integer(row["face"], f"{path}:{line} face")
        if face in fluxes or not (0 <= face < len(mesh.edges)):
            raise VerificationError(f"{path}:{line}: duplicate/out-of-range face {face}")
        edge = mesh.edges[face]
        if integer(row["owner"], "owner") != edge.owner or integer(row["neighbour"], "neighbour") != edge.neighbour:
            raise VerificationError(f"{path}:{line}: owner/neighbour differs from CM2D")
        fluxes[face] = finite(row["flux"], f"{path}:{line} flux")
    if len(fluxes) != len(mesh.edges):
        raise VerificationError(f"{path}: got {len(fluxes)} faces, expected {len(mesh.edges)}")
    return [fluxes[i] for i in range(len(mesh.edges))]


def continuity(mesh: Mesh, measured: Measurement, fluxes: list[float], speed: float, case: str,
               absolute: float, relative: float) -> dict[str, Any]:
    balances = [0.0] * len(mesh.cells)
    boundary = 0.0
    for edge, flux in zip(mesh.edges, fluxes):
        balances[edge.owner] += flux
        if edge.neighbour >= 0:
            balances[edge.neighbour] -= flux
        else:
            boundary += flux
    local_limits = [absolute + relative * speed * math.sqrt(a) for a in measured.areas]
    ratios = [abs(value) / limit for value, limit in zip(balances, local_limits)]
    native_continuity = max((abs(value) / (speed * math.sqrt(area))
                             for value, area in zip(balances, measured.areas)), default=0.0)
    global_limit = absolute + relative * speed * math.sqrt(measured.total_area)
    inflow = math.fsum(max(0.0, -flux) for edge, flux in zip(mesh.edges, fluxes)
                       if edge.neighbour < 0)
    flow_scale = speed * (measured.bounds[3] - measured.bounds[1]) if case == "cavity" else inflow
    if not (math.isfinite(flow_scale) and flow_scale > 0.0):
        raise VerificationError("independent continuity has no positive reference throughput")
    global_relative = abs(boundary) / flow_scale
    return {
        "maxCellImbalance": max(map(abs, balances), default=0.0),
        "l2CellImbalance": math.hypot(*balances),
        "maxScaledCellImbalance": max(ratios, default=0.0),
        "nativeDefinitionContinuity": native_continuity,
        "boundaryFluxSum": boundary,
        "referenceThroughput": flow_scale,
        "globalRelativeImbalance": global_relative,
        "cellValid": all(r <= 1.0 for r in ratios),
        "globalValid": abs(boundary) <= global_limit and global_relative < 1e-8,
        "globalLimit": global_limit,
    }


def read_residuals(path: Path) -> dict[str, Any]:
    required = ("iteration", "momentumResidual", "continuity", "velocityChange", "pressureChange")
    rows = load_csv(path, required)
    parsed: list[dict[str, float]] = []
    last_iteration = 0
    for line, row in enumerate(rows, 2):
        iteration = integer(row["iteration"], f"{path}:{line} iteration")
        if iteration <= last_iteration:
            raise VerificationError(f"{path}:{line}: iterations are not strictly increasing")
        last_iteration = iteration
        item = {name: finite(row[name], f"{path}:{line} {name}") for name in required[1:]}
        if any(value < 0.0 for value in item.values()):
            raise VerificationError(f"{path}:{line}: residual/change is negative")
        item["iteration"] = float(iteration)
        parsed.append(item)
    if not parsed:
        raise VerificationError(f"{path}: residual history is empty")
    return {"rows": len(parsed), "first": parsed[0], "last": parsed[-1]}


def weighted_l2(errors: Iterable[float], areas: Iterable[float]) -> float:
    pairs = list(zip(errors, areas))
    total = math.fsum(a for _, a in pairs)
    return math.sqrt(math.fsum(a * e * e for e, a in pairs) / total)


def face_centre(mesh: Mesh, edge: Edge) -> tuple[float, float]:
    a, b = mesh.vertices[edge.v0], mesh.vertices[edge.v1]
    return ((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5)


def channel_checks(mesh: Mesh, measured: Measurement, cells: list[dict[str, float]],
                   fluxes: list[float], nu: float, speed: float, args: argparse.Namespace) -> dict[str, Any]:
    xmin, ymin, xmax, ymax = measured.bounds
    height, length = ymax - ymin, xmax - xmin
    gradient = -8.0 * nu * speed / (height * height)
    expected_q = 2.0 * speed * height / 3.0
    u_errors, v_errors, p_errors = [], [], []
    for row in cells:
        eta = (row["y"] - ymin) / height
        u_exact = 4.0 * speed * eta * (1.0 - eta)
        p_exact = -gradient * (xmax - row["x"])
        u_errors.append(row["u"] - u_exact)
        v_errors.append(row["v"])
        p_errors.append(row["p"] - p_exact)
    pressure_scale = max(abs(gradient) * length, 1e-300)
    u_l2 = weighted_l2(u_errors, measured.areas) / speed
    v_l2 = weighted_l2(v_errors, measured.areas) / speed
    p_l2 = weighted_l2(p_errors, measured.areas) / pressure_scale
    mean_x = math.fsum(row["x"] for row in cells) / len(cells)
    mean_p = math.fsum(row["p"] for row in cells) / len(cells)
    denom = math.fsum((row["x"] - mean_x) ** 2 for row in cells)
    fitted_gradient = math.fsum((row["x"] - mean_x) * (row["p"] - mean_p) for row in cells) / denom
    gradient_error = abs(fitted_gradient - gradient) / abs(gradient)
    span_tol = 1e-9 * max(length, height) + 1e-12
    boundary = {"left": 0.0, "right": 0.0, "bottom": 0.0, "top": 0.0}
    for edge, flux in zip(mesh.edges, fluxes):
        if edge.neighbour >= 0:
            continue
        x, y = face_centre(mesh, edge)
        if abs(x - xmin) <= span_tol:
            boundary["left"] += flux
        elif abs(x - xmax) <= span_tol:
            boundary["right"] += flux
        elif abs(y - ymin) <= span_tol:
            boundary["bottom"] += flux
        elif abs(y - ymax) <= span_tol:
            boundary["top"] += flux
    flux_error = max(abs(boundary["left"] + expected_q), abs(boundary["right"] - expected_q)) / expected_q
    wall_leak = max(abs(boundary["bottom"]), abs(boundary["top"])) / expected_q
    checks = {
        "velocityL2Relative": u_l2,
        "transverseVelocityL2Relative": v_l2,
        "pressureL2Relative": p_l2,
        "expectedPressureGradient": gradient,
        "fittedPressureGradient": fitted_gradient,
        "pressureGradientRelativeError": gradient_error,
        "expectedFlowRatePerDepth": expected_q,
        "boundaryFluxes": boundary,
        "flowRateRelativeError": flux_error,
        "wallFluxRelative": wall_leak,
    }
    checks["valid"] = (u_l2 <= args.channel_velocity_l2 and v_l2 <= args.channel_transverse_l2
                       and p_l2 <= args.channel_pressure_l2 and gradient_error <= args.channel_gradient_error
                       and flux_error <= args.channel_flux_error and wall_leak <= args.channel_flux_error)
    return checks


def idw(cells: list[dict[str, float]], x: float, y: float, field: str, count: int = 8,
        boundary: Iterable[dict[str, float]] = ()) -> float:
    # A cell-centred field has no sample on a Dirichlet wall.  Include the
    # exact wall value instead of extrapolating the outermost cell value to a
    # Ghia point close to that wall.
    nearest = sorted(cells, key=lambda row: (row["x"] - x) ** 2 + (row["y"] - y) ** 2)[:count]
    nearest.extend(boundary)
    weights = [1.0 / max((row["x"] - x) ** 2 + (row["y"] - y) ** 2, 1e-30) for row in nearest]
    return math.fsum(w * row[field] for w, row in zip(weights, nearest)) / math.fsum(weights)


def cavity_checks(measured: Measurement, cells: list[dict[str, float]], nu: float,
                  speed: float, args: argparse.Namespace) -> dict[str, Any]:
    xmin, ymin, xmax, ymax = measured.bounds
    width, height = xmax - xmin, ymax - ymin
    re = speed * width / nu
    samples: list[dict[str, float | str]] = []
    errors: list[float] = []
    if abs(re - 100.0) <= 1e-8 and abs(width / height - 1.0) <= 1e-8:
        for coordinate, reference in GHIA_U[1:-1]:
            x = xmin + 0.5 * width
            actual = idw(cells, x, ymin + coordinate * height, "u", boundary=(
                {"x": x, "y": ymin, "u": 0.0}, {"x": x, "y": ymax, "u": speed})) / speed
            samples.append({"field": "u", "coordinate": coordinate, "reference": reference, "actual": actual})
            errors.append(actual - reference)
        for coordinate, reference in GHIA_V[1:-1]:
            y = ymin + 0.5 * height
            actual = idw(cells, xmin + coordinate * width, y, "v", boundary=(
                {"x": xmin, "y": y, "v": 0.0}, {"x": xmax, "y": y, "v": 0.0})) / speed
            samples.append({"field": "v", "coordinate": coordinate, "reference": reference, "actual": actual})
            errors.append(actual - reference)
    rmse = math.sqrt(math.fsum(e * e for e in errors) / len(errors)) if errors else math.inf
    maximum = max(map(abs, errors), default=math.inf)
    centre_u = idw(cells, xmin + 0.5 * width, ymin + 0.5 * height, "u") / speed
    left_v = idw(cells, xmin + 0.25 * width, ymin + 0.5 * height, "v") / speed
    right_v = idw(cells, xmin + 0.80 * width, ymin + 0.5 * height, "v") / speed
    max_speed = max(row["speed"] for row in cells) / speed
    checks = {
        "reynolds": re,
        "reference": "Ghia, Ghia & Shin 1982, Re=100 centreline tables",
        "samplingMethod": "inverse-distance cell-centre interpolation with exact Dirichlet wall anchors",
        "samples": samples,
        "centrelineRmse": rmse,
        "centrelineMaxError": maximum,
        "centreU": centre_u,
        "leftV": left_v,
        "rightV": right_v,
        "maximumSpeedRatio": max_speed,
    }
    checks["valid"] = (bool(errors) and rmse <= args.cavity_rmse and maximum <= args.cavity_max_error
                       and centre_u < -0.03 and left_v > 0.03 and right_v < -0.03
                       and max_speed <= args.max_speed_ratio)
    return checks


def flattened_numbers(value: Any, prefix: str = "") -> dict[str, float]:
    result: dict[str, float] = {}
    if isinstance(value, dict):
        for key, child in value.items():
            result.update(flattened_numbers(child, f"{prefix}.{key}" if prefix else str(key)))
    elif isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)):
        result[prefix] = float(value)
    return result


def force_values(payload: dict[str, Any]) -> dict[str, float]:
    aliases = {
        "dragCoefficient": {"dragcoefficient", "cd"}, "liftCoefficient": {"liftcoefficient", "cl"},
        "dragForce": {"dragforce", "drag", "forcex"},
        "liftForce": {"liftforce", "lift", "forcey"},
    }
    found: dict[str, float] = {}
    for path, value in flattened_numbers(payload).items():
        leaf = "".join(ch for ch in path.rsplit(".", 1)[-1].lower() if ch.isalnum())
        for canonical, names in aliases.items():
            if leaf in names and canonical not in found:
                found[canonical] = value
    return found


def external_checks(mesh: Mesh, measured: Measurement, cells: list[dict[str, float]],
                    payload: dict[str, Any], nu: float, speed: float,
                    args: argparse.Namespace) -> dict[str, Any]:
    wall_vertices = {vertex for edge in mesh.edges if edge.neighbour < 0 and edge.patch == 1
                     for vertex in (edge.v0, edge.v1)}
    if not wall_vertices:
        return {"valid": False, "issues": ["external mesh has no embedded wall"]}
    xs = [mesh.vertices[v][0] for v in wall_vertices]
    ys = [mesh.vertices[v][1] for v in wall_vertices]
    diameter = max(max(xs) - min(xs), max(ys) - min(ys))
    obstacle_xmax = max(xs)
    re = speed * diameter / nu
    forces = force_values(payload)
    drag = forces.get("dragCoefficient", forces.get("dragForce"))
    lift = forces.get("liftCoefficient", forces.get("liftForce"))
    drag_coefficient = 2.0 * drag / (speed * speed * diameter) if drag is not None else None
    lift_coefficient = 2.0 * lift / (speed * speed * diameter) if lift is not None else None
    reference_drag = 2.045  # Dennis & Chang (1970), open circular cylinder, Re=20.
    wake = [row["u"] / speed for row in cells if row["x"] > obstacle_xmax and min(ys) < row["y"] < max(ys)]
    max_speed = max(row["speed"] for row in cells) / speed
    issues: list[str] = []
    if drag is None or lift is None:
        issues.append("native JSON does not expose both drag and lift force/coefficient")
    elif not (drag > 0.0):
        issues.append("drag is not positive")
    elif abs(lift) > args.external_lift_drag_ratio * abs(drag):
        issues.append("symmetric-circle lift/drag ratio is too large")
    if max_speed > args.max_speed_ratio:
        issues.append("external maximum speed ratio exceeds stability limit")
    return {
        "valid": not issues, "issues": issues, "reynoldsByMaximumSpeedAndDiameter": re,
        "diameter": diameter, "forces": forces, "maximumSpeedRatio": max_speed,
        "derivedBodyCoefficients": {"drag": drag_coefficient, "lift": lift_coefficient,
                                    "definition": "2 F / (U^2 D), unit density and depth"},
        "openCylinderReference": {"source": "Dennis & Chang 1970", "reynolds": 20.0,
                                  "dragCoefficient": reference_drag,
                                  "relativeDifference": (abs(drag_coefficient - reference_drag) / reference_drag
                                                         if drag_coefficient is not None else None),
                                  "acceptanceGate": False},
        "minimumWakeUOverSpeed": min(wake) if wake else None,
        "benchmarkCaveat": "Reference drag is context only: this 32-gon, finite slip-domain, coarse upwind result is not an accuracy certification; DFG 2D-1 has different geometry and boundary conditions.",
    }


def verify_case(mesh_path: Path, prefix: Path, case: str, nu: float, speed: float,
                args: argparse.Namespace) -> dict[str, Any]:
    issues: list[str] = []
    require_final_solver_path(mesh_path)
    mesh = read_cm2d(mesh_path)
    measured = measure(mesh, args.geometry_absolute_tolerance, args.geometry_relative_tolerance)
    issues.extend(measured.issues)
    if issues:
        return {"valid": False, "case": case, "issues": issues, "mesh": str(mesh.path)}
    cells_path, faces_path = Path(str(prefix) + ".cells.csv"), Path(str(prefix) + ".faces.csv")
    residual_path, json_path = Path(str(prefix) + ".residuals.csv"), Path(str(prefix) + ".json")
    missing = [str(path) for path in (cells_path, faces_path, residual_path, json_path) if not path.is_file()]
    if missing:
        return {"valid": False, "case": case, "issues": issues + [f"missing artifact: {p}" for p in missing]}
    cells = read_cells(cells_path, mesh, measured)
    fluxes = read_faces(faces_path, mesh)
    residuals = read_residuals(residual_path)
    with json_path.open(encoding="utf-8-sig") as stream:
        payload = json.load(stream, parse_constant=lambda token: (_ for _ in ()).throw(
            VerificationError(f"native JSON contains non-finite token {token}")))
    if not isinstance(payload, dict):
        raise VerificationError("native JSON root is not an object")
    required_json = ("format", "case", "status", "converged", "cells", "iterations", "nu", "speed",
                     "continuity", "globalImbalance", "momentumResidual", "velocityChange",
                     "pressureChange", "forceX", "forceY", "globalRelativeImbalance", "domainHeight",
                     "tolerance", "units", "pressureReference", "method", "scope")
    for field in required_json:
        if field not in payload:
            issues.append(f"native JSON missing {field}")
    if payload.get("format") != "cartmesh2d-flow-summary-v1":
        issues.append("native JSON format is not cartmesh2d-flow-summary-v1")
    if payload.get("case") != case:
        issues.append("native JSON case differs from requested case")
    if payload.get("status") != "converged":
        issues.append("native JSON status is not converged")
    if payload.get("converged") is not True:
        issues.append("native solver did not report convergence")
    try:
        if integer(payload.get("cells"), "native cells") != len(mesh.cells):
            issues.append("native JSON cell count differs from CM2D")
        native_iterations = integer(payload.get("iterations"), "native iterations")
        if native_iterations != int(residuals["last"]["iteration"]) or native_iterations > args.max_iterations:
            issues.append("native iteration count differs from residual history or exceeds limit")
        native_continuity = finite(payload.get("continuity"), "native continuity")
        if native_continuity > args.max_reported_continuity:
            issues.append("native reported continuity exceeds limit")
        if not close(finite(payload.get("nu"), "native nu"), nu, 1e-15, 1e-12):
            issues.append("native JSON viscosity differs from invocation")
        if not close(finite(payload.get("speed"), "native speed"), speed, 1e-15, 1e-12):
            issues.append("native JSON speed differs from invocation")
        for field in ("globalImbalance", "globalRelativeImbalance", "momentumResidual", "velocityChange",
                      "pressureChange", "forceX", "forceY", "domainHeight", "tolerance"):
            finite(payload.get(field), f"native {field}")
    except VerificationError as exc:
        issues.append(str(exc))
    independent = continuity(mesh, measured, fluxes, speed, case,
                             args.continuity_absolute_tolerance, args.continuity_relative_tolerance)
    if not independent["cellValid"]:
        issues.append("independent per-cell face-flux continuity failed")
    if not independent["globalValid"]:
        issues.append("independent global boundary-flux balance failed")
    try:
        if not close(finite(payload.get("continuity"), "native continuity"),
                     independent["nativeDefinitionContinuity"], 1e-14, 1e-9):
            issues.append("native continuity differs from independent face-flux reconstruction")
        if not close(finite(payload.get("globalImbalance"), "native globalImbalance"),
                     independent["boundaryFluxSum"], 1e-14, 1e-9):
            issues.append("native global imbalance differs from independent boundary-flux sum")
        if not close(finite(payload.get("globalRelativeImbalance"), "native globalRelativeImbalance"),
                     independent["globalRelativeImbalance"], 1e-14, 1e-9):
            issues.append("native global relative imbalance differs from independent reconstruction")
        domain_height = measured.bounds[3] - measured.bounds[1]
        if not close(finite(payload.get("domainHeight"), "native domainHeight"), domain_height, 1e-12, 1e-10):
            issues.append("native domain height differs from final CM2D bounds")
        expected_pressure_reference = ("cell 0, kinematic pressure zero" if case == "cavity"
                                       else "right outlet faces, kinematic pressure zero")
        if payload.get("pressureReference") != expected_pressure_reference:
            issues.append("native pressure reference differs from case boundary conditions")
        for native_name, csv_name in (("momentumResidual", "momentumResidual"),
                                      ("continuity", "continuity"),
                                      ("velocityChange", "velocityChange"),
                                      ("pressureChange", "pressureChange")):
            if not close(finite(payload.get(native_name), f"native {native_name}"),
                         residuals["last"][csv_name], 1e-14, 1e-9):
                issues.append(f"native {native_name} differs from final residual CSV row")
    except VerificationError as exc:
        issues.append(str(exc))
    if case == "channel":
        benchmark = channel_checks(mesh, measured, cells, fluxes, nu, speed, args)
    elif case == "cavity":
        benchmark = cavity_checks(measured, cells, nu, speed, args)
    elif case == "external":
        benchmark = external_checks(mesh, measured, cells, payload, nu, speed, args)
    else:
        raise VerificationError(f"unsupported case {case}")
    if not benchmark.get("valid"):
        issues.append(f"{case} benchmark checks failed")
    return {
        "valid": not issues, "case": case, "issues": issues,
        "mesh": str(mesh.path), "meshSha256": sha256_file(mesh.path),
        "prefix": str(prefix.resolve()), "nu": nu, "speed": speed,
        "counts": {"cells": len(mesh.cells), "faces": len(mesh.edges)},
        "meshMeasurement": {"area": measured.total_area, "characteristicH": measured.characteristic_h,
                            "bounds": measured.bounds},
        "independentContinuity": independent, "residualHistory": residuals,
        "benchmark": benchmark, "native": payload,
        "artifactSha256": {str(path): sha256_file(path) for path in
                           (cells_path, faces_path, residual_path, json_path)},
    }


def run(command: list[str], log: Path, timeout: int) -> dict[str, Any]:
    log.parent.mkdir(parents=True, exist_ok=True)
    record: dict[str, Any] = {"command": command, "commandText": shlex.join(command), "returncode": None,
                              "stdout": str(log.with_suffix(".stdout.log")),
                              "stderr": str(log.with_suffix(".stderr.log")), "timedOut": False}
    try:
        completed = subprocess.run(command, cwd=REPO, text=True, capture_output=True,
                                   check=False, timeout=timeout)
        Path(record["stdout"]).write_text(completed.stdout, encoding="utf-8")
        Path(record["stderr"]).write_text(completed.stderr, encoding="utf-8")
        record.update(returncode=completed.returncode,
                      status="passed" if completed.returncode == 0 else "failed")
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        Path(record["stdout"]).write_text(stdout, encoding="utf-8")
        Path(record["stderr"]).write_text(stderr, encoding="utf-8")
        record.update(status="failed", timedOut=True, error=f"timed out after {timeout}s")
    except OSError as exc:
        Path(record["stdout"]).write_text("", encoding="utf-8")
        Path(record["stderr"]).write_text(str(exc) + "\n", encoding="utf-8")
        record.update(status="failed", error=str(exc))
    return record


def geometry_path(root: Path, case: str) -> Path:
    path = root / "geometry" / f"{case}.xy"
    path.parent.mkdir(parents=True, exist_ok=True)
    points = ((0.0, 0.0), (4.0, 0.0), (4.0, 1.0), (0.0, 1.0)) if case == "channel" else (
        (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
    path.write_text("".join(f"{x:.17g} {y:.17g}\n" for x, y in points), encoding="utf-8")
    return path


def mesh_command(mesh_cli: Path, output: Path, request: Request, root: Path) -> list[str]:
    # Supplying a case directory makes the mesher construct and quality-check
    # the final solver partition and emit *.solver.cm2d.  It only writes case
    # files; this verifier never launches OpenFOAM.
    solver_case = output.parent / "openfoam"
    if request.case == "external":
        return [str(mesh_cli), str(DEFAULT_CIRCLE), str(output), "8", "0.25", "0.1", "exterior", str(solver_case), "0", "0",
                "--size-field", "--reference-length", "2", "--wall-relative-size", "0.0625",
                "--background-relative-size", "0.5", "--far-field-spans", "10", "--cells-per-level", "3"]
    assert request.level is not None and request.padding is not None
    return [str(mesh_cli), str(geometry_path(root, request.case)), str(output), str(request.level),
            f"{request.padding:.17g}", "0.1", "interior", str(solver_case), str(request.level), "0"]


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def sequence_checks(cases: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {"valid": True, "issues": [], "series": {}}
    for case_name, metric in (("channel", "velocityL2Relative"), ("cavity", "centrelineRmse")):
        selected = [item for item in cases if item.get("case") == case_name and item.get("benchmark", {}).get(metric) is not None]
        selected.sort(key=lambda item: item["meshMeasurement"]["characteristicH"], reverse=True)
        values = [{"label": item.get("label"), "h": item["meshMeasurement"]["characteristicH"],
                   "error": item["benchmark"][metric]} for item in selected]
        result["series"][case_name] = values
        for coarse, fine in zip(values, values[1:]):
            if not (fine["h"] < coarse["h"] and fine["error"] <= 1.10 * coarse["error"]):
                result["issues"].append(f"{case_name} finer mesh did not preserve/improve {metric}")
    result["valid"] = not result["issues"]
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", action="store_true", help="generate the default five-case suite")
    parser.add_argument("--verify-only", action="store_true",
                        help="reuse existing output-root/runs artifacts without launching the flow CLI")
    parser.add_argument("--mesh", action="append", nargs=3, metavar=("CASE", "LABEL", "PATH"),
                        help="verify an existing mesh; repeat for multiple cases")
    parser.add_argument("--cases", nargs="+", choices=("channel", "cavity", "external"),
                        default=("channel", "cavity", "external"))
    parser.add_argument("--mesh-cli", type=Path, default=REPO / "build/cartmesh2d_cli")
    parser.add_argument("--flow-cli", type=Path, default=REPO / "build/cartmesh2d_flow_cli")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--channel-nu", type=float, default=0.01)
    parser.add_argument("--cavity-nu", type=float, default=0.01)
    parser.add_argument("--external-nu", type=float, default=0.1)
    parser.add_argument("--speed", type=float, default=1.0)
    parser.add_argument("--max-iterations", type=int, default=1500)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--geometry-absolute-tolerance", type=float, default=1e-11)
    parser.add_argument("--geometry-relative-tolerance", type=float, default=1e-9)
    parser.add_argument("--continuity-absolute-tolerance", type=float, default=1e-10)
    parser.add_argument("--continuity-relative-tolerance", type=float, default=1e-7)
    parser.add_argument("--max-reported-continuity", type=float, default=1e-8)
    parser.add_argument("--channel-velocity-l2", type=float, default=0.08)
    parser.add_argument("--channel-transverse-l2", type=float, default=0.05)
    parser.add_argument("--channel-pressure-l2", type=float, default=0.12)
    parser.add_argument("--channel-gradient-error", type=float, default=0.12)
    parser.add_argument("--channel-flux-error", type=float, default=0.08)
    parser.add_argument("--cavity-rmse", type=float, default=0.12)
    parser.add_argument("--cavity-max-error", type=float, default=0.30)
    parser.add_argument("--max-speed-ratio", type=float, default=4.0)
    parser.add_argument("--external-lift-drag-ratio", type=float, default=0.30)
    args = parser.parse_args()

    output_root = args.output_root.resolve()
    summary_path = args.summary.resolve() if args.summary else output_root / "summary.json"
    summary: dict[str, Any] = {"format": "cartmesh2d-native-flow-verification-v1", "valid": False,
                              "references": REFERENCES, "meshGeneration": [], "runs": [], "cases": [],
                              "sequenceChecks": None, "issues": []}
    try:
        if args.generate and args.mesh:
            raise VerificationError("--generate and --mesh are mutually exclusive")
        if args.verify_only and not args.mesh:
            raise VerificationError("--verify-only requires explicit --mesh entries")
        positive = ("geometry_absolute_tolerance", "geometry_relative_tolerance",
                    "continuity_absolute_tolerance", "continuity_relative_tolerance",
                    "max_reported_continuity", "channel_velocity_l2", "channel_transverse_l2",
                    "channel_pressure_l2", "channel_gradient_error", "channel_flux_error",
                    "cavity_rmse", "cavity_max_error", "max_speed_ratio", "external_lift_drag_ratio")
        for name in positive:
            if not math.isfinite(getattr(args, name)) or getattr(args, name) < 0.0:
                raise VerificationError(f"--{name.replace('_', '-')} must be finite non-negative")
        if args.timeout <= 0 or args.max_iterations <= 0 or not (math.isfinite(args.speed) and args.speed > 0.0):
            raise VerificationError("timeout, iteration limit and speed must be positive")
        nu_by_case = {"channel": args.channel_nu, "cavity": args.cavity_nu, "external": args.external_nu}
        if any(not math.isfinite(value) or value <= 0.0 for value in nu_by_case.values()):
            raise VerificationError("all viscosities must be positive finite")
        mesh_cli, flow_cli = args.mesh_cli.resolve(), args.flow_cli.resolve()
        summary["executables"] = {
            "meshCli": {"path": str(mesh_cli), "sha256": sha256_file(mesh_cli)},
            "flowCli": {"path": str(flow_cli), "sha256": sha256_file(flow_cli)},
        }
        summary["acceptanceScope"] = {
            "statement": "Project milestone regression gates for these cases and meshes; not universal CFD quality standards.",
            "maxIterations": args.max_iterations,
            "maxReportedContinuity": args.max_reported_continuity,
            "channel": {"velocityL2": args.channel_velocity_l2,
                        "transverseL2": args.channel_transverse_l2,
                        "pressureL2": args.channel_pressure_l2,
                        "gradientError": args.channel_gradient_error,
                        "fluxError": args.channel_flux_error},
            "cavity": {"centrelineRmse": args.cavity_rmse,
                       "centrelineMaxError": args.cavity_max_error},
            "external": {"maximumSpeedRatio": args.max_speed_ratio,
                         "maximumLiftDragRatio": args.external_lift_drag_ratio,
                         "referenceDragIsGate": False},
        }

        entries: list[tuple[str, str, Path, float, float]] = []
        if args.mesh:
            labels: set[str] = set()
            for case, label, raw_path in args.mesh:
                if case not in args.cases or case not in nu_by_case:
                    raise VerificationError(f"unsupported or filtered case {case}")
                if label in labels:
                    raise VerificationError(f"duplicate label {label}")
                labels.add(label)
                mesh_path = Path(raw_path).resolve()
                require_final_solver_path(mesh_path)
                entries.append((label, case, mesh_path, nu_by_case[case], args.speed))
        else:
            for request in GENERATED:
                if request.case not in args.cases:
                    continue
                adjusted = Request(request.label, request.case, request.level, request.padding,
                                   nu_by_case[request.case], args.speed)
                prefix = output_root / "meshes" / adjusted.label / adjusted.label
                stage = run(mesh_command(args.mesh_cli.resolve(), prefix, adjusted, output_root),
                            output_root / "logs" / f"mesh-{adjusted.label}", args.timeout)
                stage["label"], stage["case"] = adjusted.label, adjusted.case
                summary["meshGeneration"].append(stage)
                entries.append((adjusted.label, adjusted.case, Path(str(prefix) + ".solver.cm2d"),
                                adjusted.nu, adjusted.speed))

        for label, case, mesh_path, nu, speed in entries:
            prefix = output_root / "runs" / label / case
            generation = next((stage for stage in summary["meshGeneration"] if stage["label"] == label), None)
            if generation is not None and generation.get("status") != "passed":
                summary["cases"].append({"label": label, "case": case, "valid": False,
                                         "issues": ["mesh generation failed"]})
                continue
            command = [str(args.flow_cli.resolve()), "--mesh", str(mesh_path), "--output", str(prefix),
                       "--case", case, "--nu", f"{nu:.17g}", "--speed", f"{speed:.17g}",
                       "--max-iterations", str(args.max_iterations)]
            if args.verify_only:
                log = output_root / "logs" / f"flow-{label}"
                stage = {"command": command, "commandText": shlex.join(command), "returncode": None,
                         "stdout": str(log.with_suffix(".stdout.log")),
                         "stderr": str(log.with_suffix(".stderr.log")), "timedOut": False,
                         "status": "passed", "execution": "reused-existing-artifacts"}
            else:
                stage = run(command, output_root / "logs" / f"flow-{label}", args.timeout)
            stage["label"], stage["case"] = label, case
            summary["runs"].append(stage)
            if stage.get("status") != "passed":
                item = {"label": label, "case": case, "valid": False,
                        "issues": ["native flow command failed"], "stage": stage}
            else:
                try:
                    item = verify_case(mesh_path, prefix, case, nu, speed, args)
                except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
                    item = {"case": case, "valid": False, "issues": [str(exc)]}
                item["label"], item["stage"] = label, stage
            summary["cases"].append(item)

        summary["sequenceChecks"] = sequence_checks(summary["cases"])
        summary["valid"] = (len(summary["cases"]) == len(entries) and bool(entries)
                            and all(item.get("valid") is True for item in summary["cases"])
                            and summary["sequenceChecks"]["valid"] is True)
        if not summary["valid"]:
            summary["issues"].append("one or more native flow acceptance gates failed")
    except (OSError, ValueError) as exc:
        summary["issues"].append(str(exc))
    write_json(summary_path, summary)
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0 if summary["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
