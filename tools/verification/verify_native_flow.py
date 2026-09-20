#!/usr/bin/env python3
"""Generate and independently verify the native steady 2-D flow milestone.

The native JSON is treated as a report, never as the acceptance oracle.  This
tool reads the final CM2D polygons, cell field CSV and owner-oriented face flux
CSV independently.  It reconstructs incidence and finite-volume continuity,
then applies case-specific checks:

* channel: plane Poiseuille velocity, pressure gradient and flow rate;
* cavity: Re=100 centreline values from Ghia, Ghia & Shin (1982);
* external: finite conservative steady output and recorded force values.
* manufactured: independent smooth unit-square vortex/source reconstruction,
  closed-wall/source momentum balance, and refinement error reporting.

The external case deliberately does not compare against the DFG cylinder drag:
the bundled circle uses a symmetric slip far field, whereas the DFG 2D-1 case
uses an offset cylinder in a no-slip channel.  Similar Reynolds number alone
does not make their force coefficients interchangeable.
"""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
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
class FaceGeometry:
    centre: tuple[float, float]
    area_vector: tuple[float, float]
    correction: tuple[float, float]
    transmissibility: float
    neighbour_weight: float


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
    # Filtered out by default --cases; explicit --cases manufactured enables
    # these two square MMS meshes without changing the default five-case suite.
    Request("manufactured-l4", "manufactured", 4, 1.0 / 14.0, 0.1, 1.0),
    Request("manufactured-l5", "manufactured", 5, 1.0 / 30.0, 0.1, 1.0),
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
    with path.open("r", encoding="utf-8-sig") as stream:
        tokens = (token for line in stream for token in line.split())

        def take() -> str:
            try:
                return next(tokens)
            except StopIteration as exc:
                raise VerificationError(f"{path}: truncated CM2D record") from exc

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
        if take() != "END":
            raise VerificationError(f"{path}: trailing or missing CM2D data")
        try:
            next(tokens)
        except StopIteration:
            return Mesh(path.resolve(), tuple(vertices), tuple(edges), tuple(cells), audit)
        raise VerificationError(f"{path}: trailing or missing CM2D data")


def polygon(points: list[tuple[float, float]]) -> tuple[float, tuple[float, float]]:
    if not points:
        raise VerificationError("cell polygon is not finite counter-clockwise positive area")
    # Translate before cross products: small Cut-cells far from the global
    # origin otherwise subtract nearly equal products, corrupting their area
    # and centroid even when the final sums use fsum. This remains a direct
    # measurement of actual vertices, independent of exported cell geometry.
    ox, oy = points[0]
    local = [(x - ox, y - oy) for x, y in points]
    pairs = list(zip(local, local[1:] + local[:1]))
    cross = [a[0] * b[1] - b[0] * a[1] for a, b in pairs]
    twice = math.fsum(cross)
    if not math.isfinite(twice) or twice <= 0.0:
        raise VerificationError("cell polygon is not finite counter-clockwise positive area")
    # This is the independent reference measurement, not the native mesher.
    # Exact binary64 vertices -> rational polygon moments -> one final rounding.
    # Selective high precision misses neighbours just below an aspect trigger;
    # their one-ulp centroid drift can create a spurious transverse gradient.
    # Finite Decimal arithmetic can also perturb exact halfway rounding ties.
    origin=tuple(Fraction.from_float(x) for x in points[0])
    vertices=[tuple(Fraction.from_float(p[k])-origin[k] for k in (0,1)) for p in points]
    edges=list(zip(vertices,vertices[1:]+vertices[:1]))
    products=[a[0]*b[1]-b[0]*a[1] for a,b in edges]
    total=sum(products)
    if total<=0:raise VerificationError("cell polygon is not finite counter-clockwise positive area")
    centre=tuple(float(origin[k]+sum((a[k]+b[k])*q for (a,b),q in zip(edges,products))/(3*total)) for k in (0,1))
    return float(total/2),centre


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


def csv_fields(path: Path) -> tuple[str, ...]:
    """Return the header without accepting a malformed/empty CSV as legacy."""
    with path.open(newline="", encoding="utf-8-sig") as stream:
        fields = csv.DictReader(stream).fieldnames or []
    return tuple(fields)


MANUFACTURED_CELL_COLUMNS = ("sourceX", "sourceY", "exactU", "exactV", "exactP")


def read_cells(path: Path, mesh: Mesh, measured: Measurement,
               case: str | None = None) -> list[dict[str, float]]:
    fields = csv_fields(path)
    manufactured_fields = [name for name in MANUFACTURED_CELL_COLUMNS if name in fields]
    if case in ("manufactured", "counterflow"):
        if tuple(manufactured_fields) != MANUFACTURED_CELL_COLUMNS:
            missing = [name for name in MANUFACTURED_CELL_COLUMNS if name not in fields]
            raise VerificationError(f"{path}: manufactured CSV missing columns {missing}")
    elif manufactured_fields:
        raise VerificationError(f"{path}: source/exact columns are only valid for manufactured/counterflow cases")
    rows = load_csv(path, ("cell", "x", "y", "area", "u", "v", "p", "speed"))
    values: dict[int, dict[str, float]] = {}
    for line, row in enumerate(rows, 2):
        cell = integer(row["cell"], f"{path}:{line} cell")
        if cell in values or not (0 <= cell < len(mesh.cells)):
            raise VerificationError(f"{path}:{line}: duplicate/out-of-range cell {cell}")
        item = {name: finite(row[name], f"{path}:{line} {name}") for name in ("x", "y", "area", "u", "v", "p", "speed")}
        if case in ("manufactured", "counterflow"):
            item.update({name: finite(row[name], f"{path}:{line} {name}")
                         for name in MANUFACTURED_CELL_COLUMNS})
        if all(name in fields for name in ("previousU", "previousV", "temporalX", "temporalY")):
            item.update({name: finite(row[name], f"{path}:{line} {name}")
                         for name in ("previousU", "previousV", "temporalX", "temporalY")})
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


FACE_MOMENTUM_COLUMNS = ("pressure", "advectionX", "advectionY", "diffusionX", "diffusionY")


def read_faces(path: Path, mesh: Mesh) -> tuple[list[dict[str, float]], bool]:
    rows = load_csv(path, ("face", "owner", "neighbour", "flux"))
    fields = csv_fields(path)
    present = [name for name in FACE_MOMENTUM_COLUMNS if name in fields]
    if present and len(present) != len(FACE_MOMENTUM_COLUMNS):
        missing = [name for name in FACE_MOMENTUM_COLUMNS if name not in fields]
        raise VerificationError(f"{path}: partial momentum face schema; missing {missing}")
    has_momentum = bool(present)
    has_wall = "wall" in fields
    records: dict[int, dict[str, float]] = {}
    for line, row in enumerate(rows, 2):
        face = integer(row["face"], f"{path}:{line} face")
        if face in records or not (0 <= face < len(mesh.edges)):
            raise VerificationError(f"{path}:{line}: duplicate/out-of-range face {face}")
        edge = mesh.edges[face]
        if integer(row["owner"], "owner") != edge.owner or integer(row["neighbour"], "neighbour") != edge.neighbour:
            raise VerificationError(f"{path}:{line}: owner/neighbour differs from CM2D")
        item = {"flux": finite(row["flux"], f"{path}:{line} flux")}
        if has_momentum:
            item.update({name: finite(row[name], f"{path}:{line} {name}")
                         for name in FACE_MOMENTUM_COLUMNS})
        if "viscosity" in fields:
            item["viscosity"] = finite(row["viscosity"], "face viscosity")
        if has_wall:
            wall = integer(row["wall"], f"{path}:{line} wall")
            if wall not in (0, 1):
                raise VerificationError(f"{path}:{line}: wall must be 0 or 1")
            item["wall"] = float(wall)
        records[face] = item
    if len(records) != len(mesh.edges):
        raise VerificationError(f"{path}: got {len(records)} faces, expected {len(mesh.edges)}")
    return [records[i] for i in range(len(mesh.edges))], has_momentum


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
    if case in ("cavity", "manufactured", "taylor-green"):
        # Both are closed domains, so inlet flux is identically zero.  Keep a
        # positive scale for the reported relative global balance while the
        # actual gate remains the closed-boundary flux sum.
        flow_scale = speed * math.sqrt(measured.total_area)
    else:
        flow_scale = inflow
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


def manufactured_sample(x: float, y: float, speed: float, nu: float,
                       pressure_slope: float = 0.0, viscosity_slope: float = 0.0,
                       symmetric: bool = True) -> dict[str, float]:
    """Evaluate the smooth unit-square manufactured field independently.

    The expressions below are an analytic differentiation of the stream
    function and pressure definition.  They intentionally do not consume any
    exported exact/source columns; those columns are checked against this
    independent reconstruction later.
    """
    pi = math.pi
    sx, sy = math.sin(pi * x), math.sin(pi * y)
    cx, cy = math.cos(pi * x), math.cos(pi * y)
    s2x, s2y = math.sin(2.0 * pi * x), math.sin(2.0 * pi * y)
    u = speed * sx * sx * s2y
    v = -speed * s2x * sy * sy
    # Simplified convective products and Laplacians are written out here so
    # this oracle has a visibly independent form from the native implementation.
    convective_x = 2.0 * speed * speed * pi * sx * sx * s2x * sy * sy
    convective_y = 2.0 * speed * speed * pi * sx * sx * sy * sy * s2y
    lap_u = 2.0 * speed * pi * pi * s2y * (1.0 - 4.0 * sx * sx)
    lap_v = -2.0 * speed * pi * pi * s2x * (1.0 - 4.0 * sy * sy)
    pressure = speed * speed * (cx * cy + pressure_slope * (x + y))
    pressure_x = speed * speed * (-pi * sx * cy + pressure_slope)
    pressure_y = speed * speed * (-pi * cx * sy + pressure_slope)
    local_nu = nu*(1+viscosity_slope*x)
    ux = speed*pi*s2x*s2y
    uy = 2*speed*pi*sx*sx*(1-2*sy*sy)
    vx = -2*speed*pi*(1-2*sx*sx)*sy*sy
    source_x = convective_x + pressure_x - local_nu * lap_u - nu*viscosity_slope*(2*ux if symmetric else ux)
    source_y = convective_y + pressure_y - local_nu * lap_v - nu*viscosity_slope*(vx+uy if symmetric else vx)
    return {"u": u, "v": v, "p": pressure, "sourceX": source_x, "sourceY": source_y}


def manufactured_source_integral(mesh: Mesh, measured: Measurement, speed: float,
                                 nu: float, pressure_slope: float = 0.0, viscosity_slope: float = 0.0,
                                 symmetric: bool = True) -> list[tuple[float, float]]:
    return [
        (area * sample["sourceX"], area * sample["sourceY"])
        for area, (x, y) in zip(measured.areas, measured.centroids)
        for sample in (manufactured_sample(x, y, speed, nu, pressure_slope, viscosity_slope, symmetric),)
    ]


def face_centre(mesh: Mesh, edge: Edge) -> tuple[float, float]:
    a, b = mesh.vertices[edge.v0], mesh.vertices[edge.v1]
    return ((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5)


def face_fluxes(records: list[dict[str, float]]) -> list[float]:
    return [record["flux"] for record in records]


def face_geometry(mesh: Mesh, measured: Measurement) -> list[FaceGeometry]:
    """Rebuild FvMesh2D's owner-normal coefficients from CM2D polygons."""
    result: list[FaceGeometry] = []
    for edge_id, edge in enumerate(mesh.edges):
        cell = mesh.cells[edge.owner]
        try:
            local = cell.edges.index(edge_id)
        except ValueError as exc:
            raise VerificationError(f"face {edge_id}: owner does not reference face") from exc
        a = mesh.vertices[cell.vertices[local]]
        b = mesh.vertices[cell.vertices[(local + 1) % len(cell.vertices)]]
        # CM2D polygons are positive CCW.  (dy,-dx) is outward for owner.
        s = (b[1] - a[1], a[0] - b[0])
        fc = ((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5)
        centre = measured.centroids[edge.owner]
        other = measured.centroids[edge.neighbour] if edge.neighbour >= 0 else fc
        d = (other[0] - centre[0], other[1] - centre[1])
        sd = s[0] * d[0] + s[1] * d[1]
        s2 = s[0] * s[0] + s[1] * s[1]
        if not (sd > 0.0 and s2 > 0.0):
            raise VerificationError(f"face {edge_id}: nonpositive normal-centre distance")
        transmissibility = s2 / sd
        correction = (s[0] - transmissibility * d[0], s[1] - transmissibility * d[1])
        weight = ((s[0] * (fc[0] - centre[0]) + s[1] * (fc[1] - centre[1])) / sd
                  if edge.neighbour >= 0 else 0.0)
        result.append(FaceGeometry(fc, s, correction, transmissibility, weight))
    return result


def outlet_backflow_mode(payload: dict[str, Any]) -> str:
    mode = payload.get("outletBackflow", "reject")
    if mode not in ("reject", "normal-inlet"):
        raise VerificationError(f"native outletBackflow is unsupported: {mode!r}")
    return mode


def audit_explicit_boundaries(prefix: Path, mesh: Mesh, measured: Measurement,
                              payload: dict[str, Any]) -> list[dict[str, Any]]:
    """Read the separate input snapshot and bind every record to CM2D geometry."""
    if payload.get('boundaryFileSuffix') != '.boundaries' or payload.get('referenceSpeedRole') != 'normalization-only':
        raise VerificationError('custom boundary input metadata is missing or unsupported')
    path = Path(str(prefix) + '.boundaries')
    lines = [shlex.split(line) for line in path.read_text(encoding='utf-8').splitlines() if line.strip()]
    nb = sum(edge.neighbour < 0 for edge in mesh.edges)
    if (len(lines) != nb + 3 or lines[0] != ['CARTMESH2D_FLOW_BOUNDARIES', '1'] or
            lines[1] != ['COUNTS', str(len(mesh.cells)), str(len(mesh.edges)), str(nb)] or
            lines[-1] != ['END']):
        raise VerificationError('custom boundary file header/counts/termination differ from mesh')
    geometries = face_geometry(mesh, measured)
    xs, ys = zip(*(g.centre for g in geometries))
    position_tolerance = 1e-12 + 1e-10 * max(max(xs)-min(xs), max(ys)-min(ys))
    records = []
    for tokens in lines[2:-1]:
        if len(tokens) != 12 or tokens[0] != 'BOUNDARY':
            raise VerificationError('malformed custom boundary record')
        face, owner = integer(tokens[1], 'boundary face'), integer(tokens[2], 'boundary owner')
        if not 0 <= face < len(mesh.edges) or mesh.edges[face].neighbour >= 0 or mesh.edges[face].owner != owner:
            raise VerificationError('custom boundary face/owner differs from mesh')
        geom = geometries[face]
        values = [finite(value, 'boundary geometry') for value in tokens[3:7]]
        vector_tolerance = 1e-12 + 1e-10 * math.hypot(*geom.area_vector)
        if any(abs(a-b) > tolerance for a, b, tolerance in zip(values, (*geom.centre, *geom.area_vector),
                (position_tolerance, position_tolerance, vector_tolerance, vector_tolerance))):
            raise VerificationError('custom boundary geometry differs from final CM2D')
        records.append(dict(face=face, type=tokens[7], name=tokens[8],
                            **{k: finite(v, 'boundary '+k) for k, v in zip(('u','v','p'), tokens[9:])}))
    # Independently validate coverage and physical roles before using either
    # the input or summary in momentum reconstruction.
    flow_boundaries(mesh, measured, 'custom', finite(payload.get('speed'), 'custom speed'),
                    outlet_backflow_mode(payload), explicit_boundaries=records)
    metadata = payload.get('boundaryConditions')
    if not isinstance(metadata, list) or len(metadata) != nb:
        raise VerificationError('custom summary boundary count differs from input')
    for entry in metadata:
        if not isinstance(entry, dict) or set(entry) != {'face','type','name','u','v','p'}:
            raise VerificationError('custom summary boundary record is malformed')
        integer(entry['face'], 'summary boundary face')
        for key in ('u','v','p'):
            finite(entry[key], 'summary boundary '+key)
    records.sort(key=lambda entry: entry['face'])
    if sorted(metadata, key=lambda entry: entry['face']) != records:
        raise VerificationError('custom summary boundary values differ from input snapshot')
    return records


def closed_flow_case(case: str, payload: dict[str, Any]) -> bool:
    return case in ('cavity', 'manufactured', 'taylor-green') or (
        case == 'custom' and not any(b['type'] == 'pressure-outlet' for b in payload['boundaryConditions']))


def pressure_reference(case: str, payload: dict[str, Any]) -> str:
    if closed_flow_case(case, payload):
        return 'cell 0, kinematic pressure zero'
    return ('explicit pressure outlet faces, prescribed kinematic pressure' if case == 'custom'
            else 'right outlet faces, kinematic pressure zero')


def flow_boundaries(mesh: Mesh, measured: Measurement, case: str, speed: float,
                    outlet_backflow: str = "reject",
                    fluxes: list[float] | None = None,
                    flat_plate_leading_edge: float | None = None,
                    flat_plate_top: str = "pressure-farfield",
                    explicit_boundaries: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    """Mirror the solver's explicit role/fixed-value classification."""
    outlet_backflow_mode({"outletBackflow": outlet_backflow})
    if outlet_backflow == "normal-inlet" and fluxes is None:
        raise VerificationError("normal-inlet boundary reconstruction requires final face fluxes")
    if fluxes is not None and (len(fluxes) != len(mesh.edges) or
                              any(not math.isfinite(q) for q in fluxes)):
        raise VerificationError("outlet boundary reconstruction requires finite flux per face")
    xmin, ymin, xmax, ymax = measured.bounds
    eps = 1e-10 * max(xmax - xmin, ymax - ymin) + 1e-12
    if case == 'flatplate':
        if (flat_plate_leading_edge is None or not math.isfinite(flat_plate_leading_edge) or
                not xmin <= flat_plate_leading_edge < xmax or
                flat_plate_top not in ('pressure-farfield', 'symmetry') or fluxes is None):
            raise VerificationError('invalid flat plate boundary configuration')
    roles: list[str] = ["internal"] * len(mesh.edges)
    fixed_u = [False] * len(mesh.edges)
    fixed_v = [False] * len(mesh.edges)
    fixed_p = [False] * len(mesh.edges)
    constant_u = [False] * len(mesh.edges)
    constant_v = [False] * len(mesh.edges)
    bc_u = [0.0] * len(mesh.edges)
    bc_v = [0.0] * len(mesh.edges)
    bc_p = [0.0] * len(mesh.edges)
    result = {"roles": roles, "fixedU": fixed_u, "fixedV": fixed_v, "fixedP": fixed_p,
              "constantU": constant_u, "constantV": constant_v, "u": bc_u, "v": bc_v, "p": bc_p}
    if case == 'custom':
        if outlet_backflow != 'reject' or not isinstance(explicit_boundaries, list):
            raise VerificationError('custom requires explicit boundaries and reject backflow')
        geometries = face_geometry(mesh, measured)
        seen, names = set(), {}
        for entry in explicit_boundaries:
            face = integer(entry.get('face'), 'custom face')
            if not 0 <= face < len(mesh.edges) or face in seen or mesh.edges[face].neighbour >= 0:
                raise VerificationError('custom duplicate, internal or out-of-range face')
            seen.add(face)
            kind, name = entry.get('type'), entry.get('name')
            if (not isinstance(name, str) or not name or len(name.encode('utf-8')) > 128 or
                    any(ord(c) < 32 or ord(c) == 127 or c in ',"' for c in name) or
                    (name in names and names[name] != kind)):
                raise VerificationError('custom invalid name or inconsistent named patch type')
            names[name] = kind
            u, v, p = [finite(entry.get(k), 'custom '+k) for k in ('u','v','p')]
            sx, sy = geometries[face].area_vector
            q = u*sx + v*sy
            if kind == 'velocity-inlet':
                if p != 0 or not q < 0:
                    raise VerificationError('custom velocity inlet must point inward and omit pressure')
                roles[face] = 'inlet'
                fixed_u[face] = fixed_v[face] = True
            elif kind == 'pressure-outlet':
                if u != 0 or v != 0:
                    raise VerificationError('custom pressure outlet cannot prescribe velocity')
                roles[face] = 'outlet'
                fixed_p[face] = True
            elif kind in ('wall', 'moving-wall'):
                tolerance = (1e-12 + 1e-10*max(speed, math.hypot(u,v))) * math.hypot(sx,sy)
                if p != 0 or (kind == 'wall' and (u != 0 or v != 0)) or abs(q) > tolerance:
                    raise VerificationError('custom wall has pressure or penetrating velocity')
                roles[face] = 'wall' if kind == 'wall' else 'lid'
                fixed_u[face] = fixed_v[face] = constant_u[face] = constant_v[face] = True
            else:
                raise VerificationError('custom unsupported boundary type')
            bc_u[face], bc_v[face], bc_p[face] = u, v, p
        if seen != {e.id for e in mesh.edges if e.neighbour < 0}:
            raise VerificationError('custom missing boundary faces')
        if ('inlet' in roles) != ('outlet' in roles):
            raise VerificationError('custom requires both inlet and outlet, or only walls')
        return result
    height = ymax - ymin
    for edge in mesh.edges:
        if edge.neighbour >= 0:
            continue
        x, y = face_centre(mesh, edge)
        left, right = abs(x - xmin) <= eps, abs(x - xmax) <= eps
        top, bottom = abs(y - ymax) <= eps, abs(y - ymin) <= eps
        if case == "manufactured":
            # The MMS is a closed unit-square problem: all four sides are
            # stationary no-slip walls.  In particular, there is no lid,
            # inlet, outlet, or pressure boundary.
            role = "wall"
        elif case == "taylor-green":
            # Free slip on all four axis-aligned sides: tangential velocity is
            # unconstrained, normal velocity is fixed to zero.
            role = "slip"
        elif case == "external" and edge.patch == 1:
            role = "wall"
        elif case == "cavity":
            role = "lid" if top else "wall"
        elif case == "duct":
            a, z = mesh.vertices[edge.v0], mesh.vertices[edge.v1]
            at_inlet = abs(a[0] - xmin) <= eps and abs(z[0] - xmin) <= eps
            at_outlet = abs(a[0] - xmax) <= eps and abs(z[0] - xmax) <= eps
            role = "inlet" if at_inlet else "outlet" if at_outlet else "wall"
        elif left:
            role = "inlet"
        elif right:
            role = "outlet"
        elif case == 'flatplate':
            if not (top or bottom):
                raise VerificationError('flat plate requires rectangular boundary')
            if top:
                role = 'farfield' if flat_plate_top == 'pressure-farfield' else 'slip'
            else:
                lo, hi = sorted((mesh.vertices[edge.v0][0], mesh.vertices[edge.v1][0]))
                if lo < flat_plate_leading_edge - eps and hi > flat_plate_leading_edge + eps:
                    raise VerificationError('flat plate leading edge bisects a mesh face')
                role = 'wall' if x >= flat_plate_leading_edge else 'slip'
        else:
            role = "slip" if case in ("external", "counterflow") else "wall"
        roles[edge.id] = role
        constant_u[edge.id] = role in ("wall", "lid")
        constant_v[edge.id] = constant_u[edge.id] or role == "slip"
        if role == "wall":
            fixed_u[edge.id] = fixed_v[edge.id] = True
        elif role == "lid":
            fixed_u[edge.id] = fixed_v[edge.id] = True
            bc_u[edge.id] = speed
        elif role == "inlet":
            fixed_u[edge.id] = fixed_v[edge.id] = True
            bc_u[edge.id] = (4.0 * speed * (y - ymin) * (ymax - y) / (height * height)
                             if case == "channel" else speed)
            if case == "counterflow":
                bc_u[edge.id] = speed * (1.0 + 2.0 * math.cos(2.0 * math.pi * y))
        elif role == "outlet":
            fixed_p[edge.id] = True
            if outlet_backflow == "normal-inlet" and fluxes[edge.id] < 0.0:
                fixed_v[edge.id] = constant_v[edge.id] = True
        elif role == 'farfield':
            fixed_p[edge.id] = True
            bc_u[edge.id] = speed
            fixed_u[edge.id] = constant_u[edge.id] = fluxes[edge.id] < 0
        elif role == "slip":
            if case == "taylor-green":
                fixed_u[edge.id] = left or right
                fixed_v[edge.id] = top or bottom
                # The solver's constant-trace optimization applies only to
                # the constrained normal component on each slip wall.
                constant_u[edge.id] = fixed_u[edge.id]
                constant_v[edge.id] = fixed_v[edge.id]
            else:
                fixed_v[edge.id] = True
    if case == "duct" and not all(role in roles for role in ("inlet", "outlet", "wall")):
        raise VerificationError("duct requires vertical end openings and stationary walls")
    return {"roles": roles, "fixedU": fixed_u, "fixedV": fixed_v, "fixedP": fixed_p,
            "constantU": constant_u, "constantV": constant_v,
            "u": bc_u, "v": bc_v, "p": bc_p}


def reconstruct_gradient(mesh: Mesh, measured: Measurement, geometries: list[FaceGeometry],
                         values: list[float], boundary: list[float], fixed: list[bool],
                         skip_unknown_boundary: bool = False,
                         boundary_second_ring: bool = True,
                         boundary_adaptive: bool = False) -> list[tuple[float, float]]:
    result: list[tuple[float, float]] = []
    for cell in mesh.cells:
        xx = xy = yy = bx = by = 0.0
        omitted_boundary = False
        ci = measured.centroids[cell.id]
        for edge_id in cell.edges:
            edge = mesh.edges[edge_id]
            if edge.neighbour < 0 and skip_unknown_boundary and not fixed[edge_id]:
                # One-sided pressure reconstruction uses only neighbouring
                # cell values; an unknown pressure boundary contributes no
                # artificial zero-normal row to the least-squares stencil.
                omitted_boundary = True
                continue
            if edge.neighbour >= 0:
                other = edge.neighbour
                if edge.owner == cell.id:
                    other = edge.neighbour
                else:
                    other = edge.owner
                d = (measured.centroids[other][0] - ci[0], measured.centroids[other][1] - ci[1])
                delta = values[other] - values[cell.id]
            elif fixed[edge_id]:
                fc = geometries[edge_id].centre
                d = (fc[0] - ci[0], fc[1] - ci[1])
                delta = boundary[edge_id] - values[cell.id]
            else:
                s = geometries[edge_id].area_vector
                length = math.hypot(*s)
                d = (s[0] / length, s[1] / length)
                delta = 0.0
            length = math.hypot(*d)
            if not (length > 0.0):
                raise VerificationError(f"cell {cell.id}: degenerate gradient stencil")
            if edge.neighbour >= 0 or fixed[edge_id]:
                d = (d[0] / length, d[1] / length)
                delta /= length
            xx += d[0] * d[0]
            xy += d[0] * d[1]
            yy += d[1] * d[1]
            bx += d[0] * delta
            by += d[1] * delta
        if skip_unknown_boundary:
            determinant = xx * yy - xy * xy
            rank_limit = 64.0 * 2.220446049250313e-16 * (xx + yy) * (xx + yy)
            if (boundary_second_ring and omitted_boundary) or not determinant > rank_limit:
                # Expand complete graph shells using geometry only. Legacy
                # summaries retain the old two-ring operator; the new name
                # selects a condition-based expansion up to six rings.
                reached, shell = {cell.id}, {cell.id}
                for ring in range(1, 7 if boundary_adaptive else 3):
                    fresh = set()
                    for neighbour in shell:
                        for edge_id in mesh.cells[neighbour].edges:
                            edge = mesh.edges[edge_id]
                            if edge.neighbour < 0:
                                continue
                            other = edge.neighbour if edge.owner == neighbour else edge.owner
                            if other not in reached:
                                fresh.add(other)
                    if not fresh:
                        break
                    if ring > 1:
                        for other in sorted(fresh):
                            d = (measured.centroids[other][0] - ci[0],
                                 measured.centroids[other][1] - ci[1])
                            length = math.hypot(*d)
                            if not length > 0:
                                raise VerificationError(f"cell {cell.id}: degenerate extended gradient stencil")
                            unit = (d[0] / length, d[1] / length)
                            delta = (values[other] - values[cell.id]) / length
                            xx += unit[0] * unit[0]; xy += unit[0] * unit[1]; yy += unit[1] * unit[1]
                            bx += unit[0] * delta; by += unit[1] * delta
                    reached.update(fresh); shell = fresh
                    maximum = .5 * (xx + yy + math.hypot(xx-yy, 2*xy))
                    if ring >= 2 and (not boundary_adaptive or xx*yy-xy*xy >= maximum*maximum/16.):
                        break
        det = xx * yy - xy * xy
        if not (det > 64.0 * 2.220446049250313e-16 * (xx + yy) * (xx + yy)):
            raise VerificationError(f"cell {cell.id}: rank-deficient gradient stencil")
        result.append(((yy * bx - xy * by) / det, (xx * by - xy * bx) / det))
    return result


def face_limiter(mesh: Mesh, measured: Measurement, values: list[float],
                 gradients: list[tuple[float, float]], boundary: list[float],
                 fixed: list[bool]) -> list[float]:
    limiter = [1.0] * len(mesh.cells)
    for cell in mesh.cells:
        lo = hi = values[cell.id]
        for edge_id in cell.edges:
            edge = mesh.edges[edge_id]
            if edge.neighbour >= 0:
                other = edge.neighbour if edge.owner == cell.id else edge.owner
                lo, hi = min(lo, values[other]), max(hi, values[other])
            elif fixed[edge_id]:
                lo, hi = min(lo, boundary[edge_id]), max(hi, boundary[edge_id])
        ci = measured.centroids[cell.id]
        # The field reconstruction is evaluated at every actual face centre.
        for edge_id in cell.edges:
            fc = face_centre(mesh, mesh.edges[edge_id])
            d = (fc[0] - ci[0], fc[1] - ci[1])
            delta = gradients[cell.id][0] * d[0] + gradients[cell.id][1] * d[1]
            if delta > 0.0:
                limiter[cell.id] = min(limiter[cell.id], (hi - values[cell.id]) / delta)
            elif delta < 0.0:
                limiter[cell.id] = min(limiter[cell.id], (lo - values[cell.id]) / delta)
        limiter[cell.id] = min(1.0, max(0.0, limiter[cell.id]))
    return limiter


def _interpolated_gradient(edge: Edge, gradients: list[tuple[float, float]], weight: float) -> tuple[float, float]:
    if edge.neighbour < 0:
        return gradients[edge.owner]
    owner, neighbour = gradients[edge.owner], gradients[edge.neighbour]
    return ((1.0 - weight) * owner[0] + weight * neighbour[0],
            (1.0 - weight) * owner[1] + weight * neighbour[1])


def face_frame_velocity(mesh, measured, geometries, edge, flux, u, v, gu, gv, boundaries):
    """Reconstruct in the geometric face frame using only input fields/stencils."""
    up = edge.neighbour if edge.neighbour >= 0 and flux < 0 else edge.owner
    geom = geometries[edge.id]
    length = math.hypot(*geom.area_vector)
    nx, ny = (value/length for value in geom.area_vector)
    result = [u[up],v[up]]
    for ax, ay in ((nx,ny),(-ny,nx)):
        samples = [0.]
        for fid in mesh.cells[up].edges:
            face = mesh.edges[fid]
            if face.neighbour >= 0:
                other = face.neighbour if face.owner == up else face.owner
                du,dv = u[other]-u[up],v[other]-v[up]
            elif boundaries['fixedU'][fid] or boundaries['fixedV'][fid]:
                du = boundaries['u'][fid]-u[up] if boundaries['fixedU'][fid] else 0.
                dv = boundaries['v'][fid]-v[up] if boundaries['fixedV'][fid] else 0.
            else:
                continue
            samples.append(ax*du+ay*dv)
        lo,hi = min(samples),max(samples)
        # Project the gradient tensor into this velocity direction first.
        gx,gy = ax*gu[up][0]+ay*gv[up][0],ax*gu[up][1]+ay*gv[up][1]
        changes = [gx*(geometries[fid].centre[0]-measured.centroids[up][0]) +
                   gy*(geometries[fid].centre[1]-measured.centroids[up][1]) for fid in mesh.cells[up].edges]
        phi = max(0.,min([1.]+[hi/d if d>0 else lo/d for d in changes if d!=0]))
        increment = phi*(gx*(geom.centre[0]-measured.centroids[up][0])+gy*(geom.centre[1]-measured.centroids[up][1]))
        result[0] += ax*increment;result[1] += ay*increment
    return result


def _viscous_face_gradient(mesh: Mesh, measured: Measurement, edge: Edge, geometry: FaceGeometry,
                           values: list[float], gradients: list[tuple[float, float]],
                           boundary: list[float], fixed: list[bool],
                           constant_trace: list[bool]) -> tuple[float, float]:
    """Rebuild the core's normal-corrected gradient used by symmetric stress."""
    i = edge.owner
    sx, sy = geometry.area_vector
    area = math.hypot(sx, sy)
    normal = (sx / area, sy / area)
    g = gradients[i]
    if edge.neighbour >= 0:
        other = gradients[edge.neighbour]
        w = geometry.neighbour_weight
        g = ((1.0 - w) * g[0] + w * other[0],
             (1.0 - w) * g[1] + w * other[1])
    elif constant_trace[edge.id]:
        g = (0.0, 0.0)
    gn = g[0] * normal[0] + g[1] * normal[1]
    if edge.neighbour < 0 and not fixed[edge.id]:
        return (g[0] - normal[0] * gn, g[1] - normal[1] * gn)
    if edge.neighbour >= 0:
        other_centre = measured.centroids[edge.neighbour]
        d = (other_centre[0] - measured.centroids[i][0],
             other_centre[1] - measured.centroids[i][1])
    else:
        d = (geometry.centre[0] - measured.centroids[i][0],
             geometry.centre[1] - measured.centroids[i][1])
    dn = d[0] * normal[0] + d[1] * normal[1]
    if not dn > 0.0:
        raise VerificationError(f"face {edge.id}: non-positive viscous normal distance")
    other = values[edge.neighbour] if edge.neighbour >= 0 else boundary[edge.id]
    correction = (other - values[i] - g[0] * d[0] - g[1] * d[1]) / dn
    return (g[0] + normal[0] * correction, g[1] + normal[1] * correction)


def _symmetric_viscous_correction(mesh: Mesh, measured: Measurement,
                                  geometry: FaceGeometry, edge: Edge, u: list[float], v: list[float],
                                  gu: list[tuple[float, float]], gv: list[tuple[float, float]],
                                  boundaries: dict[str, Any], nu: float) -> tuple[float, float]:
    """Return the core's added -nu*(gradU+gradU^T).S face flux."""
    au = _viscous_face_gradient(mesh, measured, edge, geometry,
                                u, gu, boundaries["u"], boundaries["fixedU"],
                                boundaries["constantU"])
    av = _viscous_face_gradient(mesh, measured, edge, geometry,
                                v, gv, boundaries["v"], boundaries["fixedV"],
                                boundaries["constantV"])
    sx, sy = geometry.area_vector
    result = (-nu * (au[0] * sx + av[0] * sy),
              -nu * (au[1] * sx + av[1] * sy))
    if edge.neighbour < 0:
        correction = geometry.correction
        if boundaries["constantU"][edge.id] and boundaries["fixedU"][edge.id]:
            result = (result[0] + nu * (gu[edge.owner][0] * correction[0] + gu[edge.owner][1] * correction[1]), result[1])
        if boundaries["constantV"][edge.id] and boundaries["fixedV"][edge.id]:
            result = (result[0], result[1] + nu * (gv[edge.owner][0] * correction[0] + gv[edge.owner][1] * correction[1]))
    return result


def _advective_value(mesh: Mesh, measured: Measurement, edge: Edge, geometry: FaceGeometry,
                     flux: float, values: list[float], gradients: list[tuple[float, float]],
                     limiter: list[float] | None, fixed: list[bool], boundary: list[float],
                     normal_inlet: bool = False) -> float:
    if edge.neighbour < 0 and fixed[edge.id]:
        return boundary[edge.id]
    if normal_inlet and edge.neighbour < 0 and flux < 0.0:
        return values[edge.owner]
    up = edge.neighbour if edge.neighbour >= 0 and flux < 0.0 else edge.owner
    if limiter is None:
        return values[up]
    ci = measured.centroids[up]
    d = (geometry.centre[0] - ci[0], geometry.centre[1] - ci[1])
    return values[up] + limiter[up] * (gradients[up][0] * d[0] + gradients[up][1] * d[1])


def _deviation(actual: float, expected: float) -> dict[str, float]:
    return {"actual": actual, "expected": expected, "absolute": abs(actual - expected),
            "relative": abs(actual - expected) / max(1.0, abs(actual), abs(expected))}


def prescribed_face_viscosity(mesh, geometry, records, nu, payload):
    """Derive coefficients from analytic definition or original file, not output."""
    model = payload.get('viscosityModel', 'uniform')
    alpha = finite(payload.get('manufacturedViscositySlope', 0.), 'viscosity slope')
    if model == 'uniform':
        if alpha != 0 or payload.get('viscosityFile') or any('viscosity' in r for r in records):
            raise VerificationError('uniform viscosity metadata conflicts with face field')
        return [nu]*len(mesh.edges), {'model':'uniform'}
    if model != 'face-values':raise VerificationError('unknown viscosity model')
    source = {'model':model}
    if payload.get('case') == 'manufactured':
        if payload.get('viscosityFile') or alpha <= -1:
            raise VerificationError('invalid manufactured viscosity definition')
        values = [nu*(1+alpha*g.centre[0]) for g in geometry]
        source['manufacturedViscositySlope'] = alpha
    else:
        if alpha != 0 or payload.get('case') in ('counterflow','taylor-green'):
            raise VerificationError('unsupported variable viscosity verification case')
        path = Path(payload.get('viscosityFile',''))
        if not path.is_file():raise VerificationError('missing viscosity input file')
        lines=path.read_text().splitlines()
        if not lines or not all(lines):raise VerificationError('empty viscosity CSV row')
        reader=csv.DictReader(lines)
        if reader.fieldnames != ['face','viscosity']:raise VerificationError('invalid viscosity CSV header')
        values=[None]*len(mesh.edges)
        for row in reader:
            if set(row) != {'face','viscosity'}:raise VerificationError('invalid viscosity CSV row')
            idx=integer(row['face'],'viscosity face')
            if idx<0 or idx>=len(values) or values[idx] is not None:raise VerificationError('duplicate/out of range viscosity face')
            values[idx]=finite(row['viscosity'],'input face viscosity')
        source.update(file=str(path), sha256=sha256_file(path))
    for value,record in zip(values,records):
        if value is None or not math.isfinite(value) or value<=0:raise VerificationError('missing or nonpositive face viscosity')
        if 'viscosity' not in record or not close(record['viscosity'],value,0.,1e-13):
            raise VerificationError('exported viscosity differs from independent input')
    return values,source


def audit_named_wall_loads(mesh, measured, faces, payload):
    """Integrate exported traction on actual edges; momentum auditing separately
    reconstructs that traction from cell fields. Legacy results may omit loads.
    """
    keys = ('namedWallLoads', 'wallLoadReference', 'wallLoadDefinition')
    if not any(k in payload for k in keys):
        return {'status': 'unavailable'}
    definition = 'fluid-on-wall / density / depth; shared-face pressure and selected viscous flux; torque positive counterclockwise'
    if (payload.get('case') != 'custom' or payload.get('wallLoadReference') != [0,0]
            or payload.get('wallLoadDefinition') != definition or not isinstance(payload.get('namedWallLoads'), list)):
        raise VerificationError('invalid named wall load definition/reference')
    groups = {}
    geometry = face_geometry(mesh, measured)
    for bc in payload['boundaryConditions']:
        if bc['type'] not in ('wall','moving-wall'): continue
        edge = mesh.edges[bc['face']]; f = geometry[edge.id]; row = faces[edge.id]
        x,y = f.centre; sx,sy = f.area_vector
        px,py = row['pressure']*sx,row['pressure']*sy
        vx,vy = row['diffusionX'],row['diffusionY']
        groups.setdefault(bc['name'], []).append(dict(length=math.hypot(sx,sy),
            pressureForceX=px,pressureForceY=py,viscousForceX=vx,viscousForceY=vy,
            forceX=px+vx,forceY=py+vy,pressureTorque=x*py-y*px,viscousTorque=x*vy-y*vx,
            torque=x*(py+vy)-y*(px+vx)))
    seen = set()
    for actual in payload['namedWallLoads']:
        if not isinstance(actual,dict) or actual.get('name') not in groups or actual['name'] in seen:
            raise VerificationError('missing, duplicate or unknown named wall load')
        name=actual['name'];seen.add(name);samples=groups[name]
        if set(actual) != {'name','faces',*samples[0]} or integer(actual['faces'],'wall face count') != len(samples):
            raise VerificationError('named wall load count/fields disagree with boundary input')
        for key in samples[0]:
            expected=math.fsum(s[key] for s in samples)
            if not close(finite(actual[key],'wall load '+key),expected,5e-10,1e-10):
                raise VerificationError('named wall load differs from face integration: '+name+' '+key)
    if seen != set(groups):
        raise VerificationError('named wall load does not cover every no-slip patch')
    return {'status':'verified','patches':len(groups),'reference':[0,0]}


def reconstruct_momentum_audit(mesh: Mesh, measured: Measurement, cells: list[dict[str, float]],
                               face_records: list[dict[str, float]], nu: float, speed: float,
                               case: str, payload: dict[str, Any],
                               manufactured_pressure_slope: float = 0.0,
                               pressure_boundary_reconstruction: str = "zero-normal",
                               time_step: float | None = None) -> dict[str, Any]:
    """Rebuild face momentum terms, equation residuals and embedded-wall forces.

    All expected values here come from CM2D geometry, exported cell fields and
    the explicit scenario boundary definitions. Actual final CSV flux selects
    upwind direction and the opt-in outlet mask; exported momentum terms never
    define the reconstructed stress, pressure or advected velocity.
    """
    geometries = face_geometry(mesh, measured)
    viscosities, viscosity_source = prescribed_face_viscosity(mesh, geometries, face_records, nu, payload)
    viscosity_slope = payload.get("manufacturedViscositySlope", 0.)
    fluxes = face_fluxes(face_records)
    outlet_backflow = outlet_backflow_mode(payload)
    boundaries = flow_boundaries(mesh, measured, case, speed, outlet_backflow, fluxes,
                                payload.get('flatPlateLeadingEdge'), payload.get('flatPlateTop', 'pressure-farfield'),
                                payload.get('boundaryConditions'))
    backflow_faces = [edge.id for edge in mesh.edges
                      if boundaries["roles"][edge.id] == "outlet" and fluxes[edge.id] < 0.0]
    outlet_inflow = math.fsum(-fluxes[i] for i in backflow_faces)
    if "outletBackflowFaces" in payload and integer(payload["outletBackflowFaces"], "native outletBackflowFaces") != len(backflow_faces):
        raise VerificationError("native outletBackflowFaces differs from actual negative outlet flux count")
    if "outletInflow" in payload and not close(finite(payload["outletInflow"], "native outletInflow"),
                                               outlet_inflow, 1e-14, 1e-9):
        raise VerificationError("native outletInflow differs from actual inward outlet volume flux")
    if outlet_backflow == "reject":
        for fid in backflow_faces:
            if fluxes[fid] < -1e-12 * speed * math.hypot(*geometries[fid].area_vector):
                raise VerificationError(f"face {fid}: outlet backflow violates reject boundary mode")
    viscous_stress = payload.get("viscousStress", "laplacian")
    if viscous_stress not in ("symmetric", "laplacian"):
        raise VerificationError(f"native viscousStress is unsupported: {viscous_stress!r}")
    wall_schema = all("wall" in record for record in face_records)
    expected_wall_flags = [float(edge.neighbour < 0 and boundaries["roles"][edge.id] in ("wall", "lid"))
                           for edge in mesh.edges]
    if viscous_stress == "symmetric" and not wall_schema:
        raise VerificationError("symmetric viscousStress requires face wall flags")
    if wall_schema:
        for edge, record, expected in zip(mesh.edges, face_records, expected_wall_flags):
            if record["wall"] != expected:
                raise VerificationError(f"face {edge.id}: wall flag differs from boundary role")
    u = [row["u"] for row in cells]
    v = [row["v"] for row in cells]
    p = [row["p"] for row in cells]
    gu = reconstruct_gradient(mesh, measured, geometries, u, boundaries["u"], boundaries["fixedU"])
    gv = reconstruct_gradient(mesh, measured, geometries, v, boundaries["v"], boundaries["fixedV"])
    gp = reconstruct_gradient(mesh, measured, geometries, p, boundaries["p"], boundaries["fixedP"],
                             skip_unknown_boundary=pressure_boundary_reconstruction in ("one-sided-linear", "one-sided-linear-2ring", "one-sided-linear-adaptive"),
                             boundary_second_ring=(pressure_boundary_reconstruction != "one-sided-linear"),
                             boundary_adaptive=(pressure_boundary_reconstruction == "one-sided-linear-adaptive"))
    pressure_faces: list[float] = []
    for edge, geom in zip(mesh.edges, geometries):
        i = edge.owner
        if edge.neighbour >= 0:
            j = edge.neighbour
            w = geom.neighbour_weight
            point = ((1.0 - w) * measured.centroids[i][0] + w * measured.centroids[j][0],
                     (1.0 - w) * measured.centroids[i][1] + w * measured.centroids[j][1])
            g = ((1.0 - w) * gp[i][0] + w * gp[j][0],
                 (1.0 - w) * gp[i][1] + w * gp[j][1])
            d = (geom.centre[0] - point[0], geom.centre[1] - point[1])
            pressure_faces.append((1.0 - w) * p[i] + w * p[j] + g[0] * d[0] + g[1] * d[1])
        elif boundaries["fixedP"][edge.id]:
            pressure_faces.append(boundaries["p"][edge.id])
        else:
            d = (geom.centre[0] - measured.centroids[i][0], geom.centre[1] - measured.centroids[i][1])
            pressure_faces.append(p[i] + gp[i][0] * d[0] + gp[i][1] * d[1])
    convection = payload.get("convection")
    if convection not in ("upwind", "limited-linear", "face-limited-linear"):
        raise VerificationError(f"native convection is unsupported: {convection!r}")
    lu = (face_limiter(mesh, measured, u, gu, boundaries["u"], boundaries["fixedU"])
          if convection == "limited-linear" else None)
    lv = (face_limiter(mesh, measured, v, gv, boundaries["v"], boundaries["fixedV"])
          if convection == "limited-linear" else None)
    cell_residuals = [(0.0, 0.0) for _ in mesh.cells]
    source_integrals = (manufactured_source_integral(mesh, measured, speed, nu,
                                                     manufactured_pressure_slope, viscosity_slope, viscous_stress=="symmetric")
                        if case == "manufactured" else [(0.0, 0.0)] * len(mesh.cells))
    if case == "counterflow":
        source_integrals = [(area * counterflow_sample(y, speed, nu)["sourceX"], 0.0)
                            for area, (_, y) in zip(measured.areas, measured.centroids)]
    source_sum = (math.fsum(x for x, _ in source_integrals),
                  math.fsum(y for _, y in source_integrals))
    boundary_vector = [0.0, 0.0]
    boundary_scale = 0.0
    max_deviation = {name: 0.0 for name in FACE_MOMENTUM_COLUMNS}
    deviations = {name: {"actual": 0.0, "expected": 0.0, "absolute": 0.0, "relative": 0.0}
                  for name in FACE_MOMENTUM_COLUMNS}
    pressure_force = [0.0, 0.0]
    discrete_force = [0.0, 0.0]
    wall_force = [0.0, 0.0]
    wall_viscous_force = [0.0, 0.0]
    for edge, geom, record, pf, flux in zip(mesh.edges, geometries, face_records, pressure_faces, fluxes):
        i = edge.owner
        face_nu = viscosities[edge.id]
        normal_inlet = boundaries["roles"][edge.id] == 'farfield' or (
            outlet_backflow == "normal-inlet" and boundaries["roles"][edge.id] == "outlet")
        if edge.neighbour >= 0:
            other_u, other_v = u[edge.neighbour], v[edge.neighbour]
        else:
            other_u, other_v = boundaries["u"][edge.id], boundaries["v"][edge.id]
        av_u = flux * _advective_value(mesh, measured, edge, geom, flux, u, gu, lu,
                                       boundaries["fixedU"], boundaries["u"], normal_inlet)
        av_v = flux * _advective_value(mesh, measured, edge, geom, flux, v, gv, lv,
                                       boundaries["fixedV"], boundaries["v"], normal_inlet)
        if convection == 'face-limited-linear':
            vector = face_frame_velocity(mesh, measured, geometries, edge, flux, u, v, gu, gv, boundaries)
            for component, values, fixed, bc in ((0,u,boundaries['fixedU'],boundaries['u']),
                                                 (1,v,boundaries['fixedV'],boundaries['v'])):
                if edge.neighbour < 0:
                    if fixed[edge.id]: vector[component] = bc[edge.id]
                    elif flux < 0 and normal_inlet: vector[component] = values[i]
            av_u,av_v = flux*vector[0],flux*vector[1]
        if edge.neighbour >= 0 or boundaries["fixedU"][edge.id]:
            gi = _interpolated_gradient(edge, gu, geom.neighbour_weight)
            dx = -face_nu * (geom.transmissibility * (other_u - u[i]) +
                        gi[0] * geom.correction[0] + gi[1] * geom.correction[1])
        else:
            dx = 0.0
        if edge.neighbour >= 0 or boundaries["fixedV"][edge.id]:
            gi = _interpolated_gradient(edge, gv, geom.neighbour_weight)
            dy = -face_nu * (geom.transmissibility * (other_v - v[i]) +
                        gi[0] * geom.correction[0] + gi[1] * geom.correction[1])
        else:
            dy = 0.0
        if viscous_stress == "symmetric":
            sx, sy = _symmetric_viscous_correction(
                mesh, measured, geom, edge, u, v, gu, gv, boundaries, face_nu)
            dx += sx
            dy += sy
        expected = {"pressure": pf, "advectionX": av_u, "advectionY": av_v,
                    "diffusionX": dx, "diffusionY": dy}
        for name in FACE_MOMENTUM_COLUMNS:
            item = _deviation(record[name], expected[name])
            deviations[name] = item if item["absolute"] > deviations[name]["absolute"] else deviations[name]
            max_deviation[name] = max(max_deviation[name], item["absolute"])
        sx, sy = geom.area_vector
        rx = pf * sx + av_u + dx
        ry = pf * sy + av_v + dy
        old = cell_residuals[i]
        cell_residuals[i] = (old[0] + rx, old[1] + ry)
        if edge.neighbour >= 0:
            old = cell_residuals[edge.neighbour]
            cell_residuals[edge.neighbour] = (old[0] - rx, old[1] - ry)
        else:
            boundary_vector[0] += rx
            boundary_vector[1] += ry
            boundary_scale += abs(pf * sx) + abs(pf * sy) + abs(av_u) + abs(av_v) + abs(dx) + abs(dy)
        if edge.patch == 1 and boundaries["roles"][edge.id] == "wall":
            px, py = pf * sx, pf * sy
            pressure_force[0] += px
            pressure_force[1] += py
            discrete_force[0] += px + dx
            discrete_force[1] += py + dy
        if boundaries["roles"][edge.id] in ("wall", "lid") and edge.neighbour < 0:
            wall_force[0] += pf * sx + dx
            wall_force[1] += pf * sy + dy
            wall_viscous_force[0] += dx
            wall_viscous_force[1] += dy
    if case in ("manufactured", "counterflow"):
        # Subtract the independent cell-volume forcing before measuring the
        # equation residual.  Exported source columns are checked separately
        # and cannot influence this audit.
        cell_residuals = [(rx - source[0], ry - source[1])
                          for (rx, ry), source in zip(cell_residuals, source_integrals)]
    temporal = [(0.0, 0.0)] * len(mesh.cells)
    if time_step is not None:
        if not time_step > 0.0:
            raise VerificationError("transient audit dt must be positive")
        for i, row in enumerate(cells):
            for name in ("previousU", "previousV", "temporalX", "temporalY"):
                if name not in row:
                    raise VerificationError(f"transient cell schema missing {name}")
            tx = measured.areas[i] * (row["u"] - row["previousU"]) / time_step
            ty = measured.areas[i] * (row["v"] - row["previousV"]) / time_step
            if not close(row["temporalX"], tx, 1e-12, 1e-9) or not close(row["temporalY"], ty, 1e-12, 1e-9):
                raise VerificationError(f"cell {i}: exported temporal integral differs from independent area/dt value")
            temporal[i] = (tx, ty)
            rx, ry = cell_residuals[i]
            cell_residuals[i] = (rx + tx, ry + ty)
    diagonal_u = [0.0] * len(mesh.cells)
    diagonal_v = [0.0] * len(mesh.cells)
    for edge, flux in zip(mesh.edges, fluxes):
        q = flux
        d = viscosities[edge.id] * geometries[edge.id].transmissibility
        if edge.neighbour >= 0:
            j = edge.neighbour
            diagonal_u[edge.owner] += d + max(q, 0.0)
            diagonal_u[j] += d + max(-q, 0.0)
            diagonal_v[edge.owner] += d + max(q, 0.0)
            diagonal_v[j] += d + max(-q, 0.0)
        else:
            normal_inflow = q < 0 and (boundaries['roles'][edge.id] == 'farfield' or
                (outlet_backflow == 'normal-inlet' and boundaries['roles'][edge.id] == 'outlet'))
            if boundaries["fixedU"][edge.id]:
                diagonal_u[edge.owner] += d
            elif not normal_inflow:
                diagonal_u[edge.owner] += q
            if boundaries["fixedV"][edge.id]:
                diagonal_v[edge.owner] += d
            elif not normal_inflow:
                diagonal_v[edge.owner] += q
    if time_step is not None:
        for i, area in enumerate(measured.areas):
            diagonal_u[i] += area / time_step
            diagonal_v[i] += area / time_step
    denominators = [(du + dv) * speed for du, dv in zip(diagonal_u, diagonal_v)]
    if any(not math.isfinite(value) or value <= 0.0 for value in denominators):
        raise VerificationError("independent momentum diagonal has no positive finite scale")
    normalized = [math.hypot(rx, ry) / denominator
                  for (rx, ry), denominator in zip(cell_residuals, denominators)]
    momentum_residual = max(normalized, default=0.0)
    worst_index = max(range(len(normalized)), key=normalized.__getitem__, default=None)
    if worst_index is None:
        worst_normalized_x = worst_normalized_y = 0.0
        worst_centre = None
    else:
        worst_denominator = denominators[worst_index]
        worst_normalized_x = cell_residuals[worst_index][0] / worst_denominator
        worst_normalized_y = cell_residuals[worst_index][1] / worst_denominator
        worst_centre = {"x": measured.centroids[worst_index][0],
                        "y": measured.centroids[worst_index][1]}
    final_native = finite(payload.get("momentumResidual"), "native momentumResidual")
    summary_deviations = {
        "momentumResidual": _deviation(final_native, momentum_residual),
        "pressureForceX": _deviation(finite(payload.get("pressureForceX"), "native pressureForceX"), pressure_force[0]),
        "pressureForceY": _deviation(finite(payload.get("pressureForceY"), "native pressureForceY"), pressure_force[1]),
        "discreteForceX": _deviation(finite(payload.get("discreteForceX"), "native discreteForceX"), discrete_force[0]),
        "discreteForceY": _deviation(finite(payload.get("discreteForceY"), "native discreteForceY"), discrete_force[1]),
    }
    reconstructed_force = [0.0, 0.0]
    for edge, geom, pf in zip(mesh.edges, geometries, pressure_faces):
        if edge.patch == 1 and boundaries["roles"][edge.id] == "wall":
            i = edge.owner
            sx, sy = geom.area_vector
            reconstructed_force[0] += pf * sx - viscosities[edge.id] * (2 * gu[i][0] * sx + (gu[i][1] + gv[i][0]) * sy)
            reconstructed_force[1] += pf * sy - viscosities[edge.id] * ((gu[i][1] + gv[i][0]) * sx + 2 * gv[i][1] * sy)
    summary_deviations.update({
        "reconstructedForceX": _deviation(finite(payload.get("reconstructedForceX", payload.get("forceX")), "native reconstructedForceX"), reconstructed_force[0]),
        "reconstructedForceY": _deviation(finite(payload.get("reconstructedForceY", payload.get("forceY")), "native reconstructedForceY"), reconstructed_force[1]),
        "forceX": _deviation(finite(payload.get("forceX"), "native forceX"), discrete_force[0] if viscous_stress == "symmetric" else reconstructed_force[0]),
        "forceY": _deviation(finite(payload.get("forceY"), "native forceY"), discrete_force[1] if viscous_stress == "symmetric" else reconstructed_force[1]),
    })
    for name, actual, expected in (
        ("wallForceX", payload.get("wallForceX"), wall_force[0]),
        ("wallForceY", payload.get("wallForceY"), wall_force[1]),
        ("wallViscousForceX", payload.get("wallViscousForceX"), wall_viscous_force[0]),
        ("wallViscousForceY", payload.get("wallViscousForceY"), wall_viscous_force[1]),
    ):
        if actual is not None:
            summary_deviations[name] = _deviation(finite(actual, f"native {name}"), expected)
    cell_sum = (math.fsum(rx for rx, _ in cell_residuals),
                math.fsum(ry for _, ry in cell_residuals))
    temporal_sum = (math.fsum(x for x, _ in temporal), math.fsum(y for _, y in temporal))
    expected_cell_sum = (boundary_vector[0] - source_sum[0] + temporal_sum[0],
                         boundary_vector[1] - source_sum[1] + temporal_sum[1])
    conservation_difference = (cell_sum[0] - expected_cell_sum[0],
                               cell_sum[1] - expected_cell_sum[1])
    source_adjusted_boundary = (boundary_vector[0] - source_sum[0],
                                boundary_vector[1] - source_sum[1])
    boundary_relative = (math.hypot(*boundary_vector) / boundary_scale
                         if boundary_scale > 0.0 else 0.0)
    source_definition = ("independent analytic source integrated at CM2D cell centroids; subtracted from face balance"
                         if case in ("manufactured", "counterflow") else "zero/unforced non-manufactured case")
    return {
        "status": "available", "valid": True, "viscosityInput": viscosity_source, "convection": convection,
        "namedWallLoads": audit_named_wall_loads(mesh, measured, face_records, payload),
        "pressureDiscretization": payload.get("pressureDiscretization"),
        "pressureBoundaryReconstruction": pressure_boundary_reconstruction,
        "outletBackflow": outlet_backflow,
        "outletBackflowAudit": {"faceIds": backflow_faces, "faceCount": len(backflow_faces),
                                "inwardVolumeFlux": outlet_inflow},
        "faceCount": len(mesh.edges), "maxFaceDeviation": max_deviation,
        "faceConsistencyTolerance": {"absolute": 5e-10, "relative": 0.0,
                                      "meaning": "CSV reconstruction comparison only; pressure in m2/s2 and momentum flux in m3/s2; not a CFD accuracy gate"},
        "faceDeviation": deviations, "cellResidual": {
            "maxNormalized": momentum_residual,
            "worstCellIndex": worst_index,
            "worstNormalizedX": worst_normalized_x,
            "worstNormalizedY": worst_normalized_y,
            "worstCentre": worst_centre,
            "l2": math.sqrt(math.fsum(x * x + y * y for x, y in cell_residuals)),
            "maxCellVector": max((math.hypot(x, y) for x, y in cell_residuals), default=0.0),
            "denominatorDefinition": "(uDiagonal+vDiagonal)*speed; nu*T plus upwind flux; fixed component uses nu*T; incoming normal-inlet u flux is on RHS, not diagonal; transient adds measuredArea/dt per component",
        },
        "summaryDeviation": summary_deviations,
        "sourceForcing": {
            "case": case, "integralX": source_sum[0], "integralY": source_sum[1],
            "definition": source_definition,
        },
        "temporal": {"enabled": time_step is not None,
                     "dt": time_step,
                     "integralX": temporal_sum[0], "integralY": temporal_sum[1],
                     "maxAbsX": max((abs(x) for x, _ in temporal), default=0.0),
                     "maxAbsY": max((abs(y) for _, y in temporal), default=0.0)},
        "pressureForce": {"x": pressure_force[0], "y": pressure_force[1]},
        "discreteForce": {"x": discrete_force[0], "y": discrete_force[1]},
        "globalBoundaryMomentum": {"x": boundary_vector[0], "y": boundary_vector[1],
                                    "magnitude": math.hypot(*boundary_vector),
                                    "boundaryTermScale": boundary_scale,
                                    "relativeToBoundary": boundary_relative,
                                    "denominatorDefinition": "sum over boundary faces of |p*Sx|+|p*Sy|+|advectionX|+|advectionY|+|diffusionX|+|diffusionY|; zero if all terms zero"},
        "globalSourceBalance": {"boundaryMinusSourceX": source_adjusted_boundary[0],
                                 "boundaryMinusSourceY": source_adjusted_boundary[1],
                                 "magnitude": math.hypot(*source_adjusted_boundary),
                                 "sourceIntegralX": source_sum[0], "sourceIntegralY": source_sum[1],
                                 "definition": "global boundary momentum minus independent volume source integral"},
        "cellResidualSum": {"x": cell_sum[0], "y": cell_sum[1]},
        "expectedCellResidualSum": {"x": expected_cell_sum[0], "y": expected_cell_sum[1]},
        "cellBoundaryConservationDifference": {"x": conservation_difference[0],
                                                 "y": conservation_difference[1],
                                                 "magnitude": math.hypot(*conservation_difference)},
        "viscousStress": viscous_stress,
        "forceDefinition": payload.get("forceDefinition", "reconstructed-newtonian-traction"),
        "reconstructedForce": {"x": reconstructed_force[0], "y": reconstructed_force[1]},
        "wallForce": {"x": wall_force[0], "y": wall_force[1]},
        "wallViscousForce": {"x": wall_viscous_force[0], "y": wall_viscous_force[1]},
        "wallFaceFlags": "checked" if wall_schema else "legacy-unavailable",
        "fullNewtonianForceAudit": "forceX/forceY selected by viscousStress; reconstructedForceX/Y remain the cell-gradient diagnostic",
    }


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


def affine_sample(cells: list[dict[str, float]], x: float, y: float, field: str,
                  count: int = 8, boundary: Iterable[dict[str, float]] = ()) -> float:
    """Local weighted plane fit, exact on affine fields including near walls.

    Coordinates are relative to the query and scaled by the local stencil.
    Unlike an inverse-distance average, the fit reproduces a gradient.  Wall
    samples are constraints with the same distance weights, not a global blend.
    Rank-deficient stencils fail explicitly rather than returning biased data.
    This is a point-value diagnostic, not a conservative field remapping.
    """
    if count < 3 or not math.isfinite(x) or not math.isfinite(y):
        raise VerificationError("invalid affine sampling controls")
    rows = sorted([*cells, *boundary],
                  key=lambda row: (math.hypot(row["x"] - x, row["y"] - y),
                                   row["x"], row["y"]))[:count]
    if len(rows) < 3 or any(not math.isfinite(row[key]) for row in rows for key in ("x", "y", field)):
        raise VerificationError("affine sampling needs finite samples")
    distances = [math.hypot(row["x"] - x, row["y"] - y) for row in rows]
    exact = [row[field] for row, distance in zip(rows, distances) if distance == 0]
    if exact:
        if any(value != exact[0] for value in exact):
            raise VerificationError("conflicting coincident samples")
        return exact[0]
    scale = max(distances)
    minimum = min(distances)
    # Normalize weights to avoid overflow on physically small geometries.
    weights = [(minimum / distance) ** 2 for distance in distances]
    weight_sum = math.fsum(weights)
    dx = [(row["x"] - x) / scale for row in rows]
    dy = [(row["y"] - y) / scale for row in rows]
    values = [row[field] for row in rows]
    def average(array: list[float]) -> float:
        return math.fsum(w * value for w, value in zip(weights, array)) / weight_sum
    mx, my, mv = average(dx), average(dy), average(values)
    xx = average([(v - mx) ** 2 for v in dx])
    yy = average([(v - my) ** 2 for v in dy])
    xy = average([(a - mx) * (b - my) for a, b in zip(dx, dy)])
    xv = average([(a - mx) * (v - mv) for a, v in zip(dx, values)])
    yv = average([(b - my) * (v - mv) for b, v in zip(dy, values)])
    determinant = xx * yy - xy * xy
    if not determinant > 1e-12 * (xx + yy) ** 2:
        raise VerificationError("rank-deficient affine sampling stencil")
    gx = (yy * xv - xy * yv) / determinant
    gy = (xx * yv - xy * xv) / determinant
    result = mv - gx * mx - gy * my
    if not math.isfinite(result):
        raise VerificationError("non-finite affine sample")
    return result


def cavity_checks(measured: Measurement, cells: list[dict[str, float]], nu: float,
                  speed: float, args: argparse.Namespace, *, legacy: bool = False) -> dict[str, Any]:
    sample = idw if legacy else affine_sample
    xmin, ymin, xmax, ymax = measured.bounds
    width, height = xmax - xmin, ymax - ymin
    re = speed * width / nu
    samples: list[dict[str, float | str]] = []
    errors: list[float] = []
    if abs(re - 100.0) <= 1e-8 and abs(width / height - 1.0) <= 1e-8:
        for coordinate, reference in GHIA_U[1:-1]:
            x = xmin + 0.5 * width
            actual = sample(cells, x, ymin + coordinate * height, "u", boundary=(
                {"x": x, "y": ymin, "u": 0.0}, {"x": x, "y": ymax, "u": speed})) / speed
            samples.append({"field": "u", "coordinate": coordinate, "reference": reference, "actual": actual})
            errors.append(actual - reference)
        for coordinate, reference in GHIA_V[1:-1]:
            y = ymin + 0.5 * height
            actual = sample(cells, xmin + coordinate * width, y, "v", boundary=(
                {"x": xmin, "y": y, "v": 0.0}, {"x": xmax, "y": y, "v": 0.0})) / speed
            samples.append({"field": "v", "coordinate": coordinate, "reference": reference, "actual": actual})
            errors.append(actual - reference)
    rmse = math.sqrt(math.fsum(e * e for e in errors) / len(errors)) if errors else math.inf
    maximum = max(map(abs, errors), default=math.inf)
    centre_u = sample(cells, xmin + 0.5 * width, ymin + 0.5 * height, "u") / speed
    left_v = sample(cells, xmin + 0.25 * width, ymin + 0.5 * height, "v") / speed
    right_v = sample(cells, xmin + 0.80 * width, ymin + 0.5 * height, "v") / speed
    max_speed = max(row["speed"] for row in cells) / speed
    checks = {
        "reynolds": re,
        "reference": "Ghia, Ghia & Shin 1982, Re=100 centreline tables",
        "samplingMethod": ("legacy IDW8 with Dirichlet wall anchors" if legacy else
                           "distance-weighted affine fit with Dirichlet wall anchors; v2"),
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
    if not legacy:
        checks["legacyIdw"] = cavity_checks(measured, cells, nu, speed, args, legacy=True)
    return checks


def counterflow_sample(y: float, speed: float, nu: float) -> dict[str, float]:
    cosine = math.cos(2.0 * math.pi * y)
    return {"u": speed * (1.0 + 2.0 * cosine), "v": 0.0, "p": 0.0,
            "sourceX": 8.0 * math.pi**2 * nu * speed * cosine, "sourceY": 0.0}


def counterflow_checks(mesh: Mesh, measured: Measurement, cells: list[dict[str, float]],
                       fluxes: list[float], nu: float, speed: float,
                       payload: dict[str, Any]) -> dict[str, Any]:
    issues = []
    if (any(not close(x, y, 1e-12, 0.0) for x, y in zip(measured.bounds, (0., 0., 1., 1.)))
            or not close(measured.total_area, 1.0, 1e-12, 0.0)):
        issues.append("counterflow requires the complete unit square")
    mode = outlet_backflow_mode(payload)
    boundaries = flow_boundaries(mesh, measured, "counterflow", speed, mode, fluxes)
    geometries = face_geometry(mesh, measured)
    outgoing, incoming = [], []
    for edge, geom in zip(mesh.edges, geometries):
        if edge.neighbour >= 0:
            continue
        role = boundaries["roles"][edge.id]
        q = fluxes[edge.id]
        if role == "inlet" and not close(q, boundaries["u"][edge.id]*geom.area_vector[0], 1e-12, 1e-9):
            issues.append(f"counterflow face {edge.id}: prescribed inlet flux differs")
        elif role == "slip" and not close(q, 0.0, 1e-12, 0.0):
            issues.append(f"counterflow face {edge.id}: slip wall leaks")
        elif role == "outlet":
            (incoming if q < 0.0 else outgoing).append(q)
    if not incoming or not outgoing:
        issues.append("counterflow must exercise simultaneous right-outlet inflow and outflow")
    errors = {key: [] for key in ("u", "v", "p")}
    column_errors = {key: 0.0 for key in MANUFACTURED_CELL_COLUMNS}
    for row, area, (_, y) in zip(cells, measured.areas, measured.centroids):
        exact = counterflow_sample(y, speed, nu)
        for key in errors:
            errors[key].append(row[key]-exact[key])
        for column, key in (("exactU", "u"), ("exactV", "v"), ("exactP", "p"),
                            ("sourceX", "sourceX"), ("sourceY", "sourceY")):
            expected = exact[key] * (area if column.startswith("source") else 1.0)
            column_errors[column] = max(column_errors[column], abs(row[column]-expected))
    if any(error > 2e-12 for error in column_errors.values()):
        issues.append("counterflow CSV exact/source columns differ from independent analytic formula")
    return {"valid": not issues, "issues": issues,
            "definition": "u=U*(1+2*cos(2*pi*y)), v=p=0; source=(8*pi^2*nu*U*cos(2*pi*y),0)",
            "sourceSampling": "analytic source density at actual cell centroid times measured area",
            "columnMaxDeviation": column_errors, "columnAbsoluteTolerance": 2e-12,
            "analyticErrors": {key+"L2": weighted_l2(values, measured.areas) for key, values in errors.items()},
            "velocityL2OverSpeed": math.hypot(weighted_l2(errors["u"], measured.areas),
                                               weighted_l2(errors["v"], measured.areas))/speed,
            "outlet": {"incomingFaces": len(incoming), "outgoingFaces": len(outgoing),
                       "inwardVolumeFlux": -math.fsum(incoming), "outwardVolumeFlux": math.fsum(outgoing)},
            "scope": "analytic errors are diagnostics; valid verifies source/BC consistency and exercises reverse flow, not accuracy certification"}


def manufactured_wall_traction(mesh, measured, face_records, nu, speed, viscosity_slope, symmetric):
    """Exact face-centre Newtonian traction of the verification velocity field."""
    geometry=face_geometry(mesh,measured);squared=[];lengths=[];errors=[]
    for edge,g,row in zip(mesh.edges,geometry,face_records):
        if edge.neighbour>=0:continue
        x,y=g.centre;sx,sy=g.area_vector;length=math.hypot(sx,sy);pi=math.pi
        ux=speed*pi*math.sin(2*pi*x)*math.sin(2*pi*y)
        uy=2*speed*pi*math.sin(pi*x)**2*math.cos(2*pi*y)
        vx=-2*speed*pi*math.cos(2*pi*x)*math.sin(pi*y)**2;vy=-ux
        material=nu*(1+viscosity_slope*x)
        if symmetric:exact=(-material*(2*ux*sx+(uy+vx)*sy),-material*((vx+uy)*sx+2*vy*sy))
        else:exact=(-material*(ux*sx+uy*sy),-material*(vx*sx+vy*sy))
        error=math.hypot(row['diffusionX']-exact[0],row['diffusionY']-exact[1])/length
        squared.append(length*error*error);lengths.append(length);errors.append(error)
    return {'l2':math.sqrt(math.fsum(squared)/math.fsum(lengths)), 'linf':max(errors),
            'definition':'Viscous traction / density at boundary face centre; edge-length weighted vector RMS. Diagnostic, not a universal tolerance.'}


def manufactured_checks(mesh: Mesh, measured: Measurement, cells: list[dict[str, float]],
                        nu: float, speed: float, pressure_slope: float = 0.0,
                        viscosity_slope: float = 0.0, symmetric: bool = True) -> dict[str, Any]:
    """Audit the analytic fields and source columns of the closed MMS case."""
    xmin, ymin, xmax, ymax = measured.bounds
    geometry_ok = (close(xmin, 0.0, 1e-12, 0.0) and close(ymin, 0.0, 1e-12, 0.0)
                   and close(xmax, 1.0, 1e-12, 0.0) and close(ymax, 1.0, 1e-12, 0.0))
    boundaries = flow_boundaries(mesh, measured, "manufactured", speed)
    boundary_edges = [edge for edge in mesh.edges if edge.neighbour < 0]
    roles = {role: boundaries["roles"].count(role) for role in set(boundaries["roles"])}
    boundary_ok = (bool(boundary_edges) and all(boundaries["roles"][edge.id] == "wall"
                                                for edge in boundary_edges)
                   and all(boundaries["fixedU"][edge.id] and boundaries["fixedV"][edge.id]
                           and not boundaries["fixedP"][edge.id] for edge in boundary_edges))

    first = manufactured_sample(*measured.centroids[0], speed, nu, pressure_slope, viscosity_slope, symmetric)
    gauge = first["p"]
    velocity_errors: list[float] = []
    pressure_errors: list[float] = []
    source_errors: list[float] = []
    exact_errors = {name: [] for name in ("exactU", "exactV", "exactP")}
    expected_sources: list[tuple[float, float]] = []
    for row, area, (x, y) in zip(cells, measured.areas, measured.centroids):
        sample = manufactured_sample(x, y, speed, nu, pressure_slope, viscosity_slope, symmetric)
        expected_sources.append((area * sample["sourceX"], area * sample["sourceY"]))
        exact_u, exact_v, exact_p = sample["u"], sample["v"], sample["p"] - gauge
        velocity_errors.extend((row["u"] - exact_u, row["v"] - exact_v))
        pressure_errors.append(row["p"] - exact_p)
        exact_errors["exactU"].append(row["exactU"] - exact_u)
        exact_errors["exactV"].append(row["exactV"] - exact_v)
        exact_errors["exactP"].append(row["exactP"] - exact_p)
        source_errors.extend((row["sourceX"] - area * sample["sourceX"],
                              row["sourceY"] - area * sample["sourceY"]))
    velocity_l2 = math.sqrt(
        weighted_l2([row["u"] - manufactured_sample(x, y, speed, nu, pressure_slope, viscosity_slope, symmetric)["u"]
                     for row, (x, y) in zip(cells, measured.centroids)], measured.areas) ** 2
        + weighted_l2([row["v"] - manufactured_sample(x, y, speed, nu, pressure_slope, viscosity_slope, symmetric)["v"]
                       for row, (x, y) in zip(cells, measured.centroids)], measured.areas) ** 2
    ) / speed
    pressure_l2 = weighted_l2(pressure_errors, measured.areas) / max(speed * speed, 1e-300)
    velocity_linf = max((abs(value) for value in velocity_errors), default=0.0) / speed
    pressure_linf = max((abs(value) for value in pressure_errors), default=0.0) / max(speed * speed, 1e-300)
    exact_column_max = {name: max((abs(value) for value in values), default=0.0)
                        for name, values in exact_errors.items()}
    source_column_max = max((abs(value) for value in source_errors), default=0.0)
    schema_tolerance = {"absolute": 2e-12, "relative": 0.0,
                        "meaning": "CSV exact/source consistency only; not a physical accuracy gate"}
    gauge_ok = close(cells[0]["p"], 0.0, schema_tolerance["absolute"], 0.0)
    columns_ok = (all(close(value, 0.0, schema_tolerance["absolute"], schema_tolerance["relative"])
                      for value in exact_column_max.values())
                  and close(source_column_max, 0.0, schema_tolerance["absolute"], schema_tolerance["relative"]))
    return {
        "valid": geometry_ok and boundary_ok and columns_ok and gauge_ok,
        "issues": (["manufactured domain is not exactly [0,1]^2"] if not geometry_ok else [])
                   + ([] if boundary_ok else ["manufactured boundary is not closed stationary no-slip"])
                   + ([] if columns_ok else ["manufactured CSV exact/source columns differ from independent analytic reconstruction"])
                   + ([] if gauge_ok else ["manufactured cell 0 pressure is not the prescribed zero gauge"]),
        "domain": {"bounds": measured.bounds, "expected": [0.0, 0.0, 1.0, 1.0]},
        "boundaryAudit": {"roles": roles, "boundaryEdges": len(boundary_edges),
                           "allStationaryNoSlipWalls": boundary_ok,
                           "pressureBoundaryEdges": sum(boundaries["fixedP"])},
        "pressureGauge": {"definition": "raw analytic p minus raw p at CM2D cell 0 centroid",
                           "cell0RawPressure": gauge, "cell0ComputedPressure": cells[0]["p"],
                           "valid": gauge_ok},
        "pressureSlope": pressure_slope,
        "velocityL2Relative": velocity_l2, "velocityLinfRelative": velocity_linf,
        "pressureL2Relative": pressure_l2, "pressureLinfRelative": pressure_linf,
        "exactColumnMaxAbsolute": exact_column_max,
        "sourceColumnMaxAbsolute": source_column_max,
        "sourceIntegralExpected": {"x": math.fsum(x for x, _ in expected_sources),
                                    "y": math.fsum(y for _, y in expected_sources)},
        "csvConsistencyTolerance": schema_tolerance,
    }


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
                                  "reynoldsMatches": math.isclose(re,20.0,rel_tol=1e-10),
                                  "relativeDifference": (abs(drag_coefficient - reference_drag) / reference_drag
                                                         if drag_coefficient is not None and math.isclose(re,20.0,rel_tol=1e-10) else None),
                                  "acceptanceGate": False},
        "minimumWakeUOverSpeed": min(wake) if wake else None,
        "benchmarkCaveat": "Reference drag is context only: this 32-gon, finite slip-domain result is not an accuracy certification; DFG 2D-1 has different geometry and boundary conditions.",
    }


def verify_case(mesh_path: Path, prefix: Path, case: str, nu: float, speed: float,
                args: argparse.Namespace) -> dict[str, Any]:
    issues: list[str] = []
    manufactured_pressure_slope = float(getattr(args, "manufactured_pressure_slope", 0.0))
    if not math.isfinite(manufactured_pressure_slope):
        raise VerificationError("manufactured pressure slope must be finite")
    if case != "manufactured" and manufactured_pressure_slope != 0.0:
        raise VerificationError("nonzero manufactured pressure slope is only valid for case manufactured")
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
    cells = read_cells(cells_path, mesh, measured, case)
    face_records, face_momentum_available = read_faces(faces_path, mesh)
    fluxes = face_fluxes(face_records)
    residuals = read_residuals(residual_path)
    with json_path.open(encoding="utf-8-sig") as stream:
        payload = json.load(stream, parse_constant=lambda token: (_ for _ in ()).throw(
            VerificationError(f"native JSON contains non-finite token {token}")))
    if not isinstance(payload, dict):
        raise VerificationError("native JSON root is not an object")
    if case == 'custom':
        audit_explicit_boundaries(prefix, mesh, measured, payload)
    expected_viscosity_slope = getattr(args, 'manufactured_viscosity_slope', None)
    if expected_viscosity_slope is not None:
        actual_slope = finite(payload.get('manufacturedViscositySlope',0.), 'viscosity slope')
        if not math.isfinite(expected_viscosity_slope) or actual_slope != expected_viscosity_slope:
            issues.append('manufactured viscosity slope differs from requested configuration')
    # The steady verifier must never silently reinterpret a transient field as
    # a steady result.  Transient artifacts have both explicit metadata and
    # previous-state columns; reject either marker here.
    if any(name in payload for name in ("temporalDiscretization", "dt", "acceptedTime",
                                        "requestedSteps", "completedSteps")):
        issues.append("steady verification received transient summary metadata")
    cell_fields = set(csv_fields(cells_path))
    transient_cell_fields = {"previousU", "previousV", "temporalX", "temporalY"}
    if cell_fields & transient_cell_fields:
        issues.append("steady verification received transient cell fields")
    required_json = ("format", "case", "status", "converged", "cells", "iterations", "nu", "speed",
                     "continuity", "globalImbalance", "momentumResidual", "velocityChange",
                     "pressureChange", "forceX", "forceY", "globalRelativeImbalance", "domainHeight",
                     "tolerance", "units", "pressureReference", "method", "scope")
    for field in required_json:
        if field not in payload:
            issues.append(f"native JSON missing {field}")
    summary_momentum_columns = ("pressureForceX", "pressureForceY", "discreteForceX", "discreteForceY",
                                "pressureDiscretization", "convection")
    summary_present = [name for name in summary_momentum_columns if name in payload]
    summary_momentum_available = bool(summary_present)
    if summary_present and len(summary_present) != len(summary_momentum_columns):
        missing = [name for name in summary_momentum_columns if name not in payload]
        issues.append(f"native JSON partial momentum schema; missing {missing}")
    if face_momentum_available != summary_momentum_available:
        issues.append("face and summary momentum schemas are inconsistent")
    if payload.get("format") != "cartmesh2d-flow-summary-v1":
        issues.append("native JSON format is not cartmesh2d-flow-summary-v1")
    if payload.get("case") != case:
        issues.append("native JSON case differs from requested case")
    try:
        outlet_backflow_mode(payload)
    except VerificationError as exc:
        issues.append(str(exc))
    viscous_stress = payload.get("viscousStress", "laplacian")
    if viscous_stress not in ("symmetric", "laplacian"):
        issues.append(f"native viscousStress is unsupported: {viscous_stress!r}")
    if "viscousStress" in payload:
        expected_definition = ("shared-face-newtonian-traction" if viscous_stress == "symmetric"
                               else "reconstructed-newtonian-traction")
        if payload.get("forceDefinition") != expected_definition:
            issues.append("native forceDefinition does not match viscousStress")
        if viscous_stress == "symmetric":
            for field in ("reconstructedForceX", "reconstructedForceY", "wallForceX", "wallForceY",
                          "wallViscousForceX", "wallViscousForceY", "wallForceDefinition"):
                if field not in payload:
                    issues.append(f"native JSON missing {field} for symmetric viscous stress")
    if case == "manufactured" and not isinstance(payload.get("manufacturedDefinition"), str):
        issues.append("native JSON missing manufacturedDefinition")
    if case == "manufactured":
        if "manufacturedPressureSlope" not in payload:
            if manufactured_pressure_slope != 0.0:
                issues.append("native JSON missing manufacturedPressureSlope for nonzero expected slope")
        else:
            try:
                if not close(finite(payload.get("manufacturedPressureSlope"),
                                    "native manufacturedPressureSlope"),
                             manufactured_pressure_slope, 1e-15, 1e-12):
                    issues.append("native manufacturedPressureSlope differs from expected invocation")
            except VerificationError as exc:
                issues.append(str(exc))
    elif "manufacturedPressureSlope" in payload:
        try:
            if finite(payload.get("manufacturedPressureSlope"), "native manufacturedPressureSlope") != 0.0:
                issues.append("nonzero manufacturedPressureSlope is invalid for ordinary cases")
        except VerificationError as exc:
            issues.append(str(exc))
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
    independent = continuity(mesh, measured, fluxes, speed, 'cavity' if closed_flow_case(case, payload) else case,
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
        expected_pressure_reference = pressure_reference(case, payload)
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
    if payload.get('viscosityModel') == 'face-values' and case in ('channel','cavity','external'):
        benchmark = {'valid':True, 'status':'not-applicable',
                     'scope':'Prescribed variable viscosity: constant-property reference is inapplicable; only geometry/conservation/constitutive audit is performed, no accuracy qualification.'}
    elif case in ("duct", "custom"):
        benchmark = {'valid': True, 'status': 'not-qualified',
                     'scope': 'Arbitrary duct: geometry, constitutive and conservation audits only; no physical reference or accuracy qualification.'}
    elif case == "channel":
        benchmark = channel_checks(mesh, measured, cells, fluxes, nu, speed, args)
    elif case == "cavity":
        benchmark = cavity_checks(measured, cells, nu, speed, args)
    elif case == "external":
        benchmark = external_checks(mesh, measured, cells, payload, nu, speed, args)
    elif case == "manufactured":
        benchmark = manufactured_checks(mesh, measured, cells, nu, speed, manufactured_pressure_slope,
                                        payload.get("manufacturedViscositySlope",0.), payload.get("viscousStress")=="symmetric")
    elif case == "counterflow":
        benchmark = counterflow_checks(mesh, measured, cells, fluxes, nu, speed, payload)
    else:
        raise VerificationError(f"unsupported case {case}")
    if case == 'manufactured' and payload.get('viscosityModel') == 'face-values' and face_momentum_available:
        benchmark['wallTraction'] = manufactured_wall_traction(mesh,measured,face_records,nu,speed,
            payload.get('manufacturedViscositySlope',0.),payload.get('viscousStress')=='symmetric')
    if not benchmark.get("valid"):
        issues.append(f"{case} benchmark checks failed")
    pressure_boundary_reconstruction = payload.get("pressureBoundaryReconstruction", "zero-normal")
    if pressure_boundary_reconstruction not in ("zero-normal", "one-sided-linear", "one-sided-linear-2ring", "one-sided-linear-adaptive"):
        issues.append(f"native pressureBoundaryReconstruction is unsupported: {pressure_boundary_reconstruction!r}")
    if face_momentum_available and summary_momentum_available:
        try:
            momentum_audit = reconstruct_momentum_audit(
                mesh, measured, cells, face_records, nu, speed, case, payload,
                manufactured_pressure_slope, pressure_boundary_reconstruction
            )
            if pressure_boundary_reconstruction not in ("zero-normal", "one-sided-linear", "one-sided-linear-2ring", "one-sided-linear-adaptive"):
                momentum_audit["valid"] = False
                momentum_audit.setdefault("issues", []).append(
                    "unsupported pressureBoundaryReconstruction metadata")
            if payload.get("pressureDiscretization") != "shared-face-gauss":
                momentum_audit["valid"] = False
                momentum_audit.setdefault("issues", []).append(
                    "pressureDiscretization is not shared-face-gauss")
            # Metadata is part of the schema contract; the field-value audit
            # must not silently accept a different convection operator.
            if payload.get("convection") not in ("upwind", "limited-linear", "face-limited-linear"):
                momentum_audit["valid"] = False
                momentum_audit.setdefault('issues', []).append('unsupported momentum convection metadata')
            face_tol = 5e-10
            summary_tol = 5e-10
            if any(value > face_tol for value in momentum_audit["maxFaceDeviation"].values()):
                momentum_audit["valid"] = False
                momentum_audit.setdefault("issues", []).append("independent face momentum reconstruction differs")
            if momentum_audit["summaryDeviation"]["momentumResidual"]["absolute"] > summary_tol:
                momentum_audit["valid"] = False
                momentum_audit.setdefault("issues", []).append("independent momentum residual differs from summary")
            if case == "counterflow" and momentum_audit["cellResidual"]["maxNormalized"] >= finite(payload.get("tolerance"), "native tolerance"):
                momentum_audit["valid"] = False
                momentum_audit.setdefault("issues", []).append("independent counterflow momentum residual is not below tolerance")
            if any(momentum_audit["summaryDeviation"][name]["absolute"] > summary_tol
                   for name in ("pressureForceX", "pressureForceY", "discreteForceX", "discreteForceY")):
                momentum_audit["valid"] = False
                momentum_audit.setdefault("issues", []).append("independent force reconstruction differs from summary")
            force_names = ("reconstructedForceX", "reconstructedForceY", "forceX", "forceY")
            if any(momentum_audit["summaryDeviation"][name]["absolute"] > summary_tol
                   for name in force_names):
                momentum_audit["valid"] = False
                momentum_audit.setdefault("issues", []).append("independent selected/reconstructed force differs from summary")
            wall_names = ("wallForceX", "wallForceY", "wallViscousForceX", "wallViscousForceY")
            if viscous_stress == "symmetric":
                if any(name not in momentum_audit["summaryDeviation"] for name in wall_names):
                    momentum_audit["valid"] = False
                    momentum_audit.setdefault("issues", []).append("symmetric wall force summary is unavailable")
                elif any(momentum_audit["summaryDeviation"][name]["absolute"] > summary_tol for name in wall_names):
                    momentum_audit["valid"] = False
                    momentum_audit.setdefault("issues", []).append("independent wall force reconstruction differs from summary")
            if not momentum_audit["valid"]:
                issues.extend(momentum_audit.get("issues", ["momentum audit failed"]))
        except (VerificationError, ArithmeticError) as exc:
            momentum_audit = {"status": "available", "valid": False, "issues": [str(exc)]}
            issues.append(f"momentum audit failed: {exc}")
    else:
        momentum_audit = {"status": "unavailable", "valid": False,
                          "reason": "legacy faces.csv/summary schema has no momentum face terms"}
    return {
        "valid": not issues, "case": case, "issues": issues,
        "mesh": str(mesh.path), "meshSha256": sha256_file(mesh.path),
        "prefix": str(prefix.resolve()), "nu": nu, "speed": speed,
        "manufacturedPressureSlope": manufactured_pressure_slope if case == "manufactured" else 0.0,
        "pressureBoundaryReconstruction": pressure_boundary_reconstruction,
        "outletBackflow": payload.get("outletBackflow", "reject"),
        "counts": {"cells": len(mesh.cells), "faces": len(mesh.edges)},
        "meshMeasurement": {"area": measured.total_area, "characteristicH": measured.characteristic_h,
                            "bounds": measured.bounds},
        "independentContinuity": independent, "residualHistory": residuals,
        "momentumAudit": momentum_audit,
        "benchmark": benchmark, "native": payload,
        "artifactSha256": {str(path): sha256_file(path) for path in
                           (cells_path, faces_path, residual_path, json_path,
                            *((Path(str(prefix)+'.boundaries'),) if case == 'custom' else ()))},
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
    for case_name, metric in (("channel", "velocityL2Relative"),
                              ("cavity", "centrelineRmse"),
                              ("manufactured", "velocityL2Relative")):
        selected = [item for item in cases if item.get("case") == case_name and item.get("benchmark", {}).get(metric) is not None]
        if case_name == "cavity" and len({item["benchmark"].get("samplingMethod", "unrecorded")
                                         for item in selected}) > 1:
            result["issues"].append("cavity sequence mixes different sampling methods")
        selected.sort(key=lambda item: item["meshMeasurement"]["characteristicH"], reverse=True)
        values = [{"label": item.get("label"), "h": item["meshMeasurement"]["characteristicH"],
                   "error": item["benchmark"][metric]} for item in selected]
        result["series"][case_name] = values
        for coarse, fine in zip(values, values[1:]):
            if not fine["h"] < coarse["h"]:
                result["issues"].append(f"{case_name} sequence is not geometrically refined")
                continue
            if coarse["error"] > 0.0 and fine["error"] > 0.0:
                observed_order = math.log(coarse["error"] / fine["error"]) / math.log(coarse["h"] / fine["h"])
            else:
                observed_order = None
            coarse["observedOrderToNext"] = observed_order
            fine["observedOrder"] = observed_order
            if case_name == "manufactured":
                if not fine["error"] < coarse["error"]:
                    result["issues"].append(f"manufactured finer mesh did not reduce velocity error")
            elif not fine["error"] <= 1.10 * coarse["error"]:
                result["issues"].append(f"{case_name} finer mesh did not preserve/improve {metric}")
    result["valid"] = not result["issues"]
    return result


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", action="store_true", help="generate the default five-case suite")
    parser.add_argument("--verify-only", action="store_true",
                        help="reuse existing output-root/runs artifacts without launching the flow CLI")
    parser.add_argument("--mesh", action="append", nargs=3, metavar=("CASE", "LABEL", "PATH"),
                        help="verify an existing mesh; repeat for multiple cases")
    parser.add_argument("--cases", nargs="+", choices=("channel", "cavity", "external", "manufactured", "counterflow"),
                        default=("channel", "cavity", "external"))
    parser.add_argument("--mesh-cli", type=Path, default=REPO / "build/cartmesh2d_cli")
    parser.add_argument("--flow-cli", type=Path, default=REPO / "build/cartmesh2d_flow_cli")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--channel-nu", type=float, default=0.01)
    parser.add_argument("--duct-nu", type=float, default=0.1,
                        help="viscosity for explicit --mesh duct LABEL PATH; conservation audit only")
    parser.add_argument("--cavity-nu", type=float, default=0.01)
    parser.add_argument("--external-nu", type=float, default=0.1)
    parser.add_argument("--manufactured-nu", type=float, default=0.1)
    parser.add_argument("--counterflow-nu", type=float, default=0.1)
    parser.add_argument("--outlet-backflow", choices=("reject", "normal-inlet"), default="reject")
    parser.add_argument("--manufactured-pressure-slope", type=float, default=0.0)
    parser.add_argument("--manufactured-viscosity-slope",type=float,default=None,
                        help="manufactured nu(x)=nu*(1+slope*x); omitted audit uses recorded material definition, generation uses zero")
    parser.add_argument("--pressure-preconditioner", choices=("ic0", "jacobi", "aggregation"), default="ic0",
                        help="native pressure preconditioner; aggregation remains experimental")
    parser.add_argument("--viscous-stress", choices=("symmetric", "laplacian"), default="symmetric",
                        help="native viscous stress mode for generated runs; legacy metadata remains laplacian")
    parser.add_argument("--convection", choices=("upwind", "limited-linear", "face-limited-linear"), default="upwind",
                        help="native convection mode for generated runs")
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
    return parser


def main() -> int:
    args = argument_parser().parse_args()

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
        if not math.isfinite(args.manufactured_pressure_slope):
            raise VerificationError("--manufactured-pressure-slope must be finite")
        if args.manufactured_pressure_slope != 0.0 and any(case != "manufactured" for case in args.cases):
            raise VerificationError("nonzero --manufactured-pressure-slope is only valid for --cases manufactured")
        nu_by_case = {"channel": args.channel_nu, "duct": args.duct_nu, "cavity": args.cavity_nu,
                      "external": args.external_nu, "manufactured": args.manufactured_nu,
                      "counterflow": args.counterflow_nu}
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
            "viscousStress": args.viscous_stress,
            "convection": args.convection,
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
            "manufactured": {"nu": args.manufactured_nu,
                             "pressureSlope": args.manufactured_pressure_slope,
                             "definition": "analytic stream-function vortex on [0,1]^2; source forcing is verification-only",
                             "velocityErrorGate": "reported with refinement decrease and observed order; no universal physical threshold"},
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
                       "--max-iterations", str(args.max_iterations),
                       "--viscous-stress", args.viscous_stress,
                       "--pressure-preconditioner", args.pressure_preconditioner,
                       "--convection", args.convection]
            if args.outlet_backflow != "reject":
                command.extend(["--outlet-backflow", args.outlet_backflow])
            if case == "manufactured":
                command.extend(["--manufactured-pressure-slope",
                                f"{args.manufactured_pressure_slope:.17g}"])
            if case == 'manufactured' and args.manufactured_viscosity_slope is not None:
                command.extend(['--manufactured-viscosity-slope',str(args.manufactured_viscosity_slope)])
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
