"""Read native mesh/field data and compute geometry for tools and plots."""
from __future__ import annotations
import csv
import hashlib
import json
import math
from fractions import Fraction
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

class MeshFormatError(ValueError):
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


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise MeshFormatError(f"{label} is boolean, not numeric")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise MeshFormatError(f"{label} is not numeric: {value!r}") from exc
    if not math.isfinite(number):
        raise MeshFormatError(f"{label} is not finite: {number!r}")
    return number


def integer(value: Any, label: str) -> int:
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise MeshFormatError(f"{label} is not an integer: {value!r}") from exc
    if str(value).strip() not in {str(number), f"+{number}"}:
        raise MeshFormatError(f"{label} is not an exact integer: {value!r}")
    return number


def read_cm2d(path: Path) -> Mesh:
    with path.open("r", encoding="utf-8-sig") as stream:
        tokens = (token for line in stream for token in line.split())

        def take() -> str:
            try:
                return next(tokens)
            except StopIteration as exc:
                raise MeshFormatError(f"{path}: truncated CM2D record") from exc

        if (take(), take(), take()) != ("CM2D", "1", "VERTICES"):
            raise MeshFormatError(f"{path}: unsupported CM2D header")
        vertex_count = integer(take(), "vertex count")
        vertices: list[tuple[float, float]] = []
        for expected in range(vertex_count):
            vertex_id = integer(take(), f"vertex {expected} id")
            if vertex_id != expected:
                raise MeshFormatError(f"{path}: non-contiguous vertex ids")
            vertices.append((finite(take(), f"vertex {expected} x"),
                             finite(take(), f"vertex {expected} y")))
        if take() != "EDGES":
            raise MeshFormatError(f"{path}: missing EDGES")
        edge_count = integer(take(), "edge count")
        edges: list[Edge] = []
        for expected in range(edge_count):
            values = [integer(take(), f"edge {expected} field") for _ in range(6)]
            edge = Edge(*values)
            if edge.id != expected or not (0 <= edge.v0 < vertex_count and 0 <= edge.v1 < vertex_count):
                raise MeshFormatError(f"{path}: invalid edge {expected}")
            edges.append(edge)
        if take() != "CELLS":
            raise MeshFormatError(f"{path}: missing CELLS")
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
                raise MeshFormatError(f"{path}: invalid cell {expected} loop")
            if any(v < 0 or v >= vertex_count for v in vertex_ids) or any(e < 0 or e >= edge_count for e in edge_ids):
                raise MeshFormatError(f"{path}: invalid reference in cell {expected}")
            cells.append(Cell(cell_id, area, vertex_ids, edge_ids))
        if take() != "AUDIT":
            raise MeshFormatError(f"{path}: missing AUDIT")
        audit = tuple(integer(take(), "audit value") for _ in range(7))
        if take() != "END":
            raise MeshFormatError(f"{path}: trailing or missing CM2D data")
        try:
            next(tokens)
        except StopIteration:
            return Mesh(path.resolve(), tuple(vertices), tuple(edges), tuple(cells), audit)
        raise MeshFormatError(f"{path}: trailing or missing CM2D data")


def polygon(points: list[tuple[float, float]]) -> tuple[float, tuple[float, float]]:
    if not points:
        raise MeshFormatError("cell polygon is not finite counter-clockwise positive area")
    ox, oy = points[0]
    local = [(x - ox, y - oy) for x, y in points]
    pairs = list(zip(local, local[1:] + local[:1]))
    cross = [a[0] * b[1] - b[0] * a[1] for a, b in pairs]
    twice = math.fsum(cross)
    if not math.isfinite(twice) or twice <= 0.0:
        raise MeshFormatError("cell polygon is not finite counter-clockwise positive area")
    origin=tuple(Fraction.from_float(x) for x in points[0])
    vertices=[tuple(Fraction.from_float(p[k])-origin[k] for k in (0,1)) for p in points]
    edges=list(zip(vertices,vertices[1:]+vertices[:1]))
    products=[a[0]*b[1]-b[0]*a[1] for a,b in edges]
    total=sum(products)
    if total<=0:raise MeshFormatError("cell polygon is not finite counter-clockwise positive area")
    centre=tuple(float(origin[k]+sum((a[k]+b[k])*q for (a,b),q in zip(edges,products))/(3*total)) for k in (0,1))
    return float(total/2),centre


def measure(mesh):
    geometry = [polygon([mesh.vertices[v] for v in cell.vertices]) for cell in mesh.cells]
    areas, centres = zip(*geometry)
    xs, ys = zip(*mesh.vertices)
    total = math.fsum(areas)
    return Measurement(tuple(areas), tuple(centres), (min(xs),min(ys),max(xs),max(ys)),
                       total, math.sqrt(total/len(areas)))


def read_cells(path):
    with Path(path).open(newline='',encoding='utf-8-sig') as stream:
        rows = list(csv.DictReader(stream))
    return [{key:float(value) for key,value in row.items()} for row in sorted(rows,key=lambda r:int(r['cell']))]


def face_geometry(mesh: Mesh, measured: Measurement) -> list[FaceGeometry]:
    """Rebuild FvMesh2D's owner-normal coefficients from CM2D polygons."""
    result: list[FaceGeometry] = []
    for edge_id, edge in enumerate(mesh.edges):
        cell = mesh.cells[edge.owner]
        try:
            local = cell.edges.index(edge_id)
        except ValueError as exc:
            raise MeshFormatError(f"face {edge_id}: owner does not reference face") from exc
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
            raise MeshFormatError(f"face {edge_id}: nonpositive normal-centre distance")
        transmissibility = s2 / sd
        correction = (s[0] - transmissibility * d[0], s[1] - transmissibility * d[1])
        weight = ((s[0] * (fc[0] - centre[0]) + s[1] * (fc[1] - centre[1])) / sd
                  if edge.neighbour >= 0 else 0.0)
        result.append(FaceGeometry(fc, s, correction, transmissibility, weight))
    return result


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
        raise MeshFormatError("invalid affine sampling controls")
    rows = sorted([*cells, *boundary],
                  key=lambda row: (math.hypot(row["x"] - x, row["y"] - y),
                                   row["x"], row["y"]))[:count]
    if len(rows) < 3 or any(not math.isfinite(row[key]) for row in rows for key in ("x", "y", field)):
        raise MeshFormatError("affine sampling needs finite samples")
    distances = [math.hypot(row["x"] - x, row["y"] - y) for row in rows]
    exact = [row[field] for row, distance in zip(rows, distances) if distance == 0]
    if exact:
        if any(value != exact[0] for value in exact):
            raise MeshFormatError("conflicting coincident samples")
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
        raise MeshFormatError("rank-deficient affine sampling stencil")
    gx = (yy * xv - xy * yv) / determinant
    gy = (xx * yv - xy * xv) / determinant
    result = mv - gx * mx - gy * my
    if not math.isfinite(result):
        raise MeshFormatError("non-finite affine sample")
    return result


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


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")

