#!/usr/bin/env python3
"""Render and summarize a completed 2D OpenFOAM airfoil RANS trial.

The script reads the actual final fields and mesh.  It reports measurements;
it deliberately does not decide whether the steady calculation converged.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import Normalize, TwoSlopeNorm  # noqa: E402
from matplotlib.path import Path as MplPath  # noqa: E402
from PIL import Image, ImageDraw  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "verification"))
from openfoam_nozzle_flow import (  # noqa: E402
    boundary_type,
    boundary_values,
    face_area,
    read_boundary,
    read_faces,
    read_labels,
    read_points,
    scalar_field,
    vector_field,
)


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
FIELD_NAMES = ("U", "p", "k", "omega", "nut", "phi")
FLUX_PATCHES = ("left", "right", "top", "bottom")


def finite_range(values) -> dict[str, float]:
    data = np.asarray(values, dtype=float)
    if data.size == 0 or not np.all(np.isfinite(data)):
        raise ValueError("cannot report an empty or non-finite field range")
    return {"min": float(np.min(data)), "max": float(np.max(data))}


def expand(values, count: int, source: Path):
    if len(values) == count:
        return values
    if len(values) == 1:
        return values * count
    raise ValueError(f"{source}: internal field length {len(values)} != mesh cell count {count}")


def numeric_directories(case: Path):
    result = []
    for path in case.iterdir():
        if path.is_dir() and re.fullmatch(FLOAT, path.name) and float(path.name) > 0:
            result.append(path)
    return sorted(result, key=lambda path: float(path.name))


def latest_complete_time(case: Path) -> Path:
    candidates = [path for path in numeric_directories(case)
                  if all((path / name).is_file() for name in FIELD_NAMES)]
    if not candidates:
        raise ValueError(f"{case}: no positive time contains {', '.join(FIELD_NAMES)}")
    return candidates[-1]


def polygon_centroid(polygon):
    x = np.asarray([p[0] for p in polygon], dtype=float)
    y = np.asarray([p[1] for p in polygon], dtype=float)
    cross = x * np.roll(y, -1) - np.roll(x, -1) * y
    twice_area = float(np.sum(cross))
    if abs(twice_area) <= 1.0e-30:
        return float(np.mean(x)), float(np.mean(y))
    return (float(np.sum((x + np.roll(x, -1)) * cross) / (3.0 * twice_area)),
            float(np.sum((y + np.roll(y, -1)) * cross) / (3.0 * twice_area)))


def load_mesh(case: Path, wall_name: str):
    mesh = case / "constant" / "polyMesh"
    points = read_points(mesh / "points")
    faces = read_faces(mesh / "faces")
    owners = read_labels(mesh / "owner")
    neighbours = read_labels(mesh / "neighbour")
    patches = read_boundary(mesh / "boundary")
    if len(owners) != len(faces):
        raise ValueError("polyMesh owner length must equal face length")
    if len(neighbours) > len(faces):
        raise ValueError("polyMesh neighbour length exceeds face length")
    cell_count = 1 + max(owners + neighbours, default=-1)
    if cell_count <= 0 or any(value < 0 or value >= cell_count for value in owners + neighbours):
        raise ValueError("polyMesh owner/neighbour contains an invalid cell label")
    front = next((patch for patch in patches if patch["name"] == "frontAndBack"), None)
    wall = next((patch for patch in patches if patch["name"] == wall_name), None)
    if front is None or front["type"] != "empty":
        raise ValueError("frontAndBack empty patch is required to recover 2D cell polygons")
    if wall is None or wall["type"] != "wall":
        raise ValueError(f"{wall_name}: wall patch is missing")

    candidates: dict[int, tuple[float, list[tuple[float, float]]]] = {}
    front_ids = range(front["startFace"], front["startFace"] + front["nFaces"])
    for face_id in front_ids:
        if not 0 <= face_id < len(faces):
            raise ValueError("frontAndBack face range exceeds polyMesh faces")
        point_ids = faces[face_id]
        polygon = [(points[index][0], points[index][1]) for index in point_ids]
        if len(set(polygon)) < 3:
            raise ValueError(f"front/back face {face_id} does not define a 2D polygon")
        cell = owners[face_id]
        mean_z = math.fsum(points[index][2] for index in point_ids) / len(point_ids)
        if cell not in candidates or mean_z < candidates[cell][0]:
            candidates[cell] = (mean_z, polygon)
    if set(candidates) != set(range(cell_count)):
        missing = cell_count - len(candidates)
        raise ValueError(f"frontAndBack/owner mapping does not cover every cell ({missing} missing)")
    polygons = [candidates[cell][1] for cell in range(cell_count)]

    graph = defaultdict(list)
    wall_face_ids = list(range(wall["startFace"], wall["startFace"] + wall["nFaces"]))
    for face_id in wall_face_ids:
        xy = list(dict.fromkeys((points[index][0], points[index][1]) for index in faces[face_id]))
        if len(xy) != 2:
            raise ValueError(f"wall face {face_id} is not one extruded 2D segment")
        a, b = xy
        graph[a].append((b, face_id))
        graph[b].append((a, face_id))
    if not graph or any(len(edges) != 2 for edges in graph.values()):
        raise ValueError(f"{wall_name}: boundary segments do not form one degree-two loop")
    start = min(graph, key=lambda point: (point[0], point[1]))
    loop = [start]
    edge_ids = []
    previous = None
    current = start
    for _ in range(len(graph)):
        choices = [entry for entry in graph[current] if entry[0] != previous]
        if not choices:
            raise ValueError(f"{wall_name}: open boundary loop")
        nxt, face_id = choices[0]
        edge_ids.append(face_id)
        previous, current = current, nxt
        if current == start:
            break
        loop.append(current)
    if current != start or len(loop) != len(graph) or len(set(edge_ids)) != len(wall_face_ids):
        raise ValueError(f"{wall_name}: wall faces are disconnected, duplicated, or multi-loop")
    return {
        "mesh": mesh, "points": points, "faces": faces, "owners": owners,
        "neighbours": neighbours, "patches": patches, "cell_count": cell_count,
        "polygons": polygons, "wall": wall, "wall_loop": loop,
        "wall_edge_face_ids": edge_ids,
    }


def cell_centres(time: Path, polygons, cell_count: int):
    if (time / "C").is_file():
        centres = expand(vector_field(time / "C"), cell_count, time / "C")
        return np.asarray([(row[0], row[1]) for row in centres]), "final C field"
    if (time / "Cx").is_file() and (time / "Cy").is_file():
        cx = expand(scalar_field(time / "Cx"), cell_count, time / "Cx")
        cy = expand(scalar_field(time / "Cy"), cell_count, time / "Cy")
        return np.column_stack((cx, cy)), "final Cx/Cy fields"
    return np.asarray([polygon_centroid(polygon) for polygon in polygons]), "2D polygon area centroids"


def read_fields(time: Path, cell_count: int):
    u = np.asarray(expand(vector_field(time / "U"), cell_count, time / "U"), dtype=float)
    fields = {name: np.asarray(expand(scalar_field(time / name), cell_count, time / name), dtype=float)
              for name in ("p", "k", "omega", "nut")}
    if u.shape != (cell_count, 3) or any(values.shape != (cell_count,) for values in fields.values()):
        raise ValueError("final field arrays do not align one-to-one with owner cell labels")
    if not np.all(np.isfinite(u)) or any(not np.all(np.isfinite(values)) for values in fields.values()):
        raise ValueError("final fields contain non-finite values")
    fields["U"] = u
    fields["speed"] = np.linalg.norm(u, axis=1)
    return fields


def parse_log(path: Path, stage: str):
    if not path.is_file():
        return []
    text = path.read_text(encoding="utf-8", errors="replace")
    matches = list(re.finditer(rf"^Time\s*=\s*({FLOAT})\s*$", text, re.M))
    rows = []
    for index, match in enumerate(matches):
        body = text[match.end():matches[index + 1].start() if index + 1 < len(matches) else len(text)]
        residuals = {}
        for field in ("Ux", "Uy", "Uz", "p", "k", "omega"):
            found = re.search(rf"Solving for\s+{field},\s+Initial residual\s*=\s*({FLOAT})", body)
            if found:
                residuals[field] = float(found.group(1))
        continuity = re.search(
            rf"time step continuity errors\s*:\s*sum local\s*=\s*({FLOAT}),\s*"
            rf"global\s*=\s*({FLOAT}),\s*cumulative\s*=\s*({FLOAT})", body)
        rows.append({"iteration": float(match.group(1)), "stage": stage,
                     "initialResiduals": residuals,
                     "continuity": ({"sumLocal": float(continuity.group(1)),
                                      "global": float(continuity.group(2)),
                                      "cumulative": float(continuity.group(3))}
                                     if continuity else None)})
    return rows


def read_solver_history(root: Path):
    log_paths = sorted((root / "stages").glob("*/stdout.log"),
                       key=lambda path: (path.stat().st_mtime_ns, str(path)))
    by_iteration = {}
    for path in log_paths:
        stage = path.parent.name
        for row in parse_log(path, stage):
            if row["initialResiduals"]:
                # A later-executed stage wins when restart logs contain the
                # same iteration.  File mtime establishes execution order.
                by_iteration[row["iteration"]] = row
    if not by_iteration:
        raise ValueError(f"{root / 'stages'}: no solver Time blocks found")
    return [by_iteration[key] for key in sorted(by_iteration)]


def solver_log_sources(root: Path):
    result = []
    for path in sorted((root / "stages").glob("*/stdout.log")):
        rows = [row for row in parse_log(path, path.parent.name) if row["initialResiduals"]]
        if rows:
            result.append({"stage": path.parent.name, "path": str(path),
                           "firstIteration": rows[0]["iteration"],
                           "lastObservedIteration": rows[-1]["iteration"],
                           "lastObservedBlockFields": sorted(rows[-1]["initialResiduals"])})
    return result


def parse_openfoam_tables(paths, key_columns=1):
    records = {}
    for path in paths:
        header = None
        for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
            line = raw.strip()
            if line.startswith("# Time"):
                header = line[1:].split()
                continue
            if not line or line.startswith("#") or header is None:
                continue
            parts = line.split()
            if len(parts) != len(header):
                raise ValueError(f"{path}: table row has {len(parts)} columns, expected {len(header)}")
            values = []
            for part in parts:
                try:
                    values.append(float(part))
                except ValueError:
                    values.append(part)
            key = tuple(values[:key_columns])
            records[key] = dict(zip(header, values))
    return [records[key] for key in sorted(records, key=lambda item: tuple(str(v) if isinstance(v, str) else v for v in item))]


def table_files(case: Path, function: str, filename: str):
    base = case / "postProcessing" / function
    return sorted(base.glob(f"*/{filename}"), key=lambda path: float(path.parent.name)) if base.is_dir() else []


def force_history(case: Path):
    paths = table_files(case, "forces", "coefficient.dat")
    if not paths:
        raise ValueError("forceCoeffs coefficient.dat is missing")
    rows = parse_openfoam_tables(paths)
    required = {"Time", "Cd", "Cl", "CmPitch"}
    if not rows or not required.issubset(rows[-1]):
        raise ValueError("force coefficient header lacks Time/Cd/Cl/CmPitch")
    return sorted(rows, key=lambda row: row["Time"])


def tail_change(rows, field: str):
    count = min(len(rows), max(5, math.ceil(0.1 * len(rows))))
    tail = rows[-count:]
    values = np.asarray([row[field] for row in tail], dtype=float)
    mean = float(np.mean(values))
    return {"sampleCount": count, "firstIteration": tail[0]["Time"],
            "lastIteration": tail[-1]["Time"], "last": float(values[-1]),
            "mean": mean, "min": float(np.min(values)), "max": float(np.max(values)),
            "range": float(np.ptp(values)), "lastMinusFirst": float(values[-1] - values[0]),
            "rangeOverAbsMean": (float(np.ptp(values) / abs(mean)) if mean != 0 else None)}


def yplus_summary(case: Path, time: Path, mesh):
    paths = table_files(case, "yPlus", "yPlus.dat")
    rows = parse_openfoam_tables(paths, key_columns=2) if paths else []
    rows = [row for row in rows if row.get("patch") == mesh["wall"]["name"]
            and row["Time"] <= float(time.name)]
    result = {"functionObjectLatest": None, "independentWallField": None}
    if rows:
        row = max(rows, key=lambda item: item["Time"])
        result["functionObjectLatest"] = {
            "iteration": row["Time"], "min": row["min"], "max": row["max"],
            "arithmeticAverage": row["average"],
            "averageDefinition": "OpenFOAM yPlus.dat function-object average; no area weighting inferred",
        }
    y_times = [path for path in numeric_directories(case) if (path / "yPlus").is_file()
               and float(path.name) <= float(time.name)]
    if y_times:
        y_time = y_times[-1]
        values_by_patch = boundary_values(y_time / "yPlus")
        values = values_by_patch.get(mesh["wall"]["name"])
        count = mesh["wall"]["nFaces"]
        if values is not None and len(values) in (1, count):
            values = values * count if len(values) == 1 else values
            start = mesh["wall"]["startFace"]
            areas = []
            for face_id in range(start, start + count):
                vector = face_area(mesh["faces"][face_id], mesh["points"])
                areas.append(math.sqrt(math.fsum(value * value for value in vector)))
            if any(area <= 0 or not math.isfinite(area) for area in areas):
                raise ValueError("wall contains non-positive face area while weighting yPlus")
            result["independentWallField"] = {
                "iteration": float(y_time.name), "faceCount": count,
                "min": min(values), "max": max(values),
                "arithmeticAverage": math.fsum(values) / count,
                "areaWeightedAverage": math.fsum(v * a for v, a in zip(values, areas)) / math.fsum(areas),
                "source": str(y_time / "yPlus"),
            }
    return result


def final_flux_summary(time: Path, mesh):
    patch_values = boundary_values(time / "phi")
    patches = {patch["name"]: patch for patch in mesh["patches"]}
    fluxes = {}
    values_expanded = {}
    for name in FLUX_PATCHES:
        patch = patches.get(name)
        values = patch_values.get(name)
        if patch is None or values is None or len(values) not in (1, patch["nFaces"]):
            raise ValueError(f"{time / 'phi'}: missing or mismatched {name} boundary flux")
        values = values * patch["nFaces"] if len(values) == 1 else values
        values_expanded[name] = values
        fluxes[name] = math.fsum(values)
    incoming = -math.fsum(min(value, 0.0) for values in values_expanded.values() for value in values)
    outgoing = math.fsum(max(value, 0.0) for values in values_expanded.values() for value in values)
    signed_net = math.fsum(fluxes.values())
    right = values_expanded["right"]
    reverse = max(0.0, -math.fsum(min(value, 0.0) for value in right))
    right_positive = math.fsum(max(value, 0.0) for value in right)
    return {
        "iteration": float(time.name), "patchSignedFlux": fluxes,
        "signedNetFlux": signed_net, "totalIncomingMagnitude": incoming,
        "totalOutgoingMagnitude": outgoing,
        "relativeImbalance": abs(signed_net) / max(incoming, outgoing, 1.0e-300),
        "rightBackflowFaceCount": sum(value < 0.0 for value in right),
        "rightFaceCount": len(right), "rightBackflowFluxMagnitude": reverse,
        "rightBackflowFractionOfGrossPositiveOutflow": reverse / max(right_positive, 1.0e-300),
        "signConvention": "OpenFOAM outward-positive phi; negative right-patch faces are counted as backflow",
        "note": "No reverse-flow velocity value is inferred from pressureInletOutletVelocity.",
    }


def discretization_summary(root: Path, case: Path):
    candidates = list((root / "stages").glob("*/fvSchemes"))
    candidates.append(case / "system" / "fvSchemes")
    candidates = [candidate for candidate in candidates if candidate.is_file()]
    path = max(candidates, key=lambda candidate: candidate.stat().st_mtime_ns) if candidates else None
    if path is None:
        return {"source": None, "divPhiU": None, "divPhiK": None, "divPhiOmega": None}
    text = path.read_text(encoding="utf-8", errors="replace")
    def scheme(field):
        match = re.search(rf"div\(phi,{field}\)\s+([^;]+);", text)
        return " ".join(match.group(1).split()) if match else None
    return {"source": str(path), "divPhiU": scheme("U"), "divPhiK": scheme("k"),
            "divPhiOmega": scheme("omega"),
            "description": "mixed convection discretization; velocity and turbulence schemes are reported separately"}


def plot_collection(ax, polygons, values, cmap, norm, title, wall_loop):
    collection = PolyCollection(polygons, array=np.asarray(values), cmap=cmap, norm=norm,
                                edgecolors="none", rasterized=True)
    ax.add_collection(collection)
    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="box")
    outline = np.asarray(wall_loop + [wall_loop[0]])
    ax.plot(outline[:, 0], outline[:, 1], color="black", lw=0.7)
    ax.set_title(title, loc="left", fontweight="semibold")
    ax.set_xlabel("x / c")
    ax.set_ylabel("y / c")
    return collection


def pressure_reference_label(value: float):
    return ("p_ref = outlet gauge pressure (0 m²/s²)" if value == 0.0
            else f"p_ref = {value:g} m²/s² (user supplied)")


def render_fields(output: Path, mesh, fields, cp, final_time: str, config, p_reference: float):
    speed_norm = Normalize(vmin=float(np.min(fields["speed"])), vmax=float(np.max(fields["speed"])))
    cp_min, cp_max = float(np.min(cp)), float(np.max(cp))
    cp_norm = (TwoSlopeNorm(vmin=cp_min, vcenter=0.0, vmax=cp_max)
               if cp_min < 0.0 < cp_max else Normalize(vmin=cp_min, vmax=cp_max))
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.6), constrained_layout=True)
    speed = plot_collection(axes[0], mesh["polygons"], fields["speed"], "viridis", speed_norm,
                            f"Velocity magnitude |U| at iteration {final_time}", mesh["wall_loop"])
    pressure = plot_collection(axes[1], mesh["polygons"], cp, "coolwarm", cp_norm,
                               f"Pressure coefficient Cp at iteration {final_time}", mesh["wall_loop"])
    fig.colorbar(speed, ax=axes[0], shrink=0.82, label="|U| [m/s]")
    fig.colorbar(pressure, ax=axes[1], shrink=0.82, label="Cp = (p - p_ref)/(0.5 U_ref²)")
    fig.suptitle(
        f"NACA 2412 — Re={float(config.get('Re', 1.0e6)):.3g} | "
        f"{mesh['cell_count']:,} cells | U∞={float(config.get('velocity', 15.0)):.3g} m/s\n"
        "Confined-domain k-ω SST trial | " + pressure_reference_label(p_reference), fontsize=12)
    fig.savefig(output, dpi=190)
    plt.close(fig)


def dilate(mask):
    expanded = mask.copy()
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            shifted = np.zeros_like(mask)
            i_src = slice(max(0, -di), mask.shape[0] - max(0, di))
            j_src = slice(max(0, -dj), mask.shape[1] - max(0, dj))
            i_dst = slice(max(0, di), mask.shape[0] - max(0, -di))
            j_dst = slice(max(0, dj), mask.shape[1] - max(0, -dj))
            shifted[i_dst, j_dst] = mask[i_src, j_src]
            expanded |= shifted
    return expanded


def render_streamlines(output: Path, centres, polygons, u, wall_loop, final_time: str):
    wall = np.asarray(wall_loop)
    chord = float(np.max(wall[:, 0]) - np.min(wall[:, 0]))
    xlim = (float(np.min(wall[:, 0]) - 0.18 * chord), float(np.max(wall[:, 0]) + 0.35 * chord))
    ylim = (float(np.min(wall[:, 1]) - 0.25 * chord), float(np.max(wall[:, 1]) + 0.25 * chord))
    margin = 0.08 * chord
    selection = ((centres[:, 0] >= xlim[0] - margin) & (centres[:, 0] <= xlim[1] + margin) &
                 (centres[:, 1] >= ylim[0] - margin) & (centres[:, 1] <= ylim[1] + margin))
    selected = np.flatnonzero(selection)
    if len(selected) < 3:
        raise ValueError("not enough cell centres in streamline window")
    nx, ny = 620, 260
    gx = np.linspace(*xlim, nx)
    gy = np.linspace(*ylim, ny)
    xx, yy = np.meshgrid(gx, gy)
    # Rasterize the actual owner-indexed cell polygons.  This avoids a Delaunay
    # bridge across the airfoil hole and keeps the mapping to OpenFOAM cells
    # explicit even on highly structured/collinear Cartesian centres.
    index_image = Image.new("I", (nx, ny), 0)
    draw = ImageDraw.Draw(index_image)
    xscale = (nx - 1) / (xlim[1] - xlim[0])
    yscale = (ny - 1) / (ylim[1] - ylim[0])
    for cell in selected:
        pixel_polygon = [((point[0] - xlim[0]) * xscale,
                          (ylim[1] - point[1]) * yscale) for point in polygons[cell]]
        draw.polygon(pixel_polygon, fill=int(cell) + 1)
    cell_index = np.flipud(np.asarray(index_image, dtype=np.int64)) - 1
    solid = MplPath(np.vstack((wall, wall[0]))).contains_points(
        np.column_stack((xx.ravel(), yy.ravel()))).reshape(xx.shape)
    solid = dilate(solid)
    # Polygon rasterization can leave one-pixel cracks on shared edges.  Fill
    # only unassigned, non-solid pixels from an assigned direct neighbour.
    for _ in range(3):
        missing = (cell_index < 0) & ~solid
        if not np.any(missing):
            break
        updated = cell_index.copy()
        for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            shifted = np.full_like(cell_index, -1)
            i_src = slice(max(0, -di), cell_index.shape[0] - max(0, di))
            j_src = slice(max(0, -dj), cell_index.shape[1] - max(0, dj))
            i_dst = slice(max(0, di), cell_index.shape[0] - max(0, -di))
            j_dst = slice(max(0, dj), cell_index.shape[1] - max(0, -dj))
            shifted[i_dst, j_dst] = cell_index[i_src, j_src]
            take = missing & (updated < 0) & (shifted >= 0)
            updated[take] = shifted[take]
        cell_index = updated
    covered = cell_index >= 0
    ui = np.zeros_like(xx)
    vi = np.zeros_like(xx)
    ui[covered] = u[cell_index[covered], 0]
    vi[covered] = u[cell_index[covered], 1]
    # A few Jacobi averaging passes make the piecewise-constant cell raster
    # usable by streamplot's continuous integrator.  Solid pixels never
    # participate, so smoothing cannot create a path through the body.
    valid = covered & ~solid
    for _ in range(8):
        u_sum, v_sum = ui.copy(), vi.copy()
        count = valid.astype(float)
        for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            shifted_u = np.zeros_like(ui)
            shifted_v = np.zeros_like(vi)
            shifted_valid = np.zeros_like(valid)
            i_src = slice(max(0, -di), ui.shape[0] - max(0, di))
            j_src = slice(max(0, -dj), ui.shape[1] - max(0, dj))
            i_dst = slice(max(0, di), ui.shape[0] - max(0, -di))
            j_dst = slice(max(0, dj), ui.shape[1] - max(0, -dj))
            shifted_u[i_dst, j_dst] = ui[i_src, j_src]
            shifted_v[i_dst, j_dst] = vi[i_src, j_src]
            shifted_valid[i_dst, j_dst] = valid[i_src, j_src]
            u_sum += shifted_u * shifted_valid
            v_sum += shifted_v * shifted_valid
            count += shifted_valid
        ui[valid] = u_sum[valid] / count[valid]
        vi[valid] = v_sum[valid] / count[valid]
    invalid = solid | ~covered
    ui = np.ma.array(ui, mask=invalid)
    vi = np.ma.array(vi, mask=invalid)
    speed = np.ma.sqrt(ui * ui + vi * vi)
    true_speed = np.linalg.norm(u, axis=1)
    norm = Normalize(vmin=float(np.min(true_speed)), vmax=float(np.max(true_speed)))
    fig, ax = plt.subplots(figsize=(13.0, 5.2), constrained_layout=True)
    seeds = np.column_stack((np.full(55, gx[2]), np.linspace(gy[2], gy[-3], 55)))
    stream = ax.streamplot(gx, gy, ui, vi, color=speed, cmap="viridis", norm=norm,
                           start_points=seeds, integration_direction="forward",
                           broken_streamlines=False, density=2.0, linewidth=0.7, arrowsize=0.7)
    closed = np.vstack((wall, wall[0]))
    ax.fill(closed[:, 0], closed[:, 1], color="white", zorder=5)
    ax.plot(closed[:, 0], closed[:, 1], color="black", lw=1.0, zorder=6)
    ax.set(xlim=xlim, ylim=ylim, aspect="equal", xlabel="x / c", ylabel="y / c")
    ax.set_title(f"Near-airfoil streamlines at iteration {final_time}", loc="left", pad=8)
    fig.suptitle("Confined-domain SST trial; solid mask blocks streamline integration", fontsize=11)
    fig.colorbar(stream.lines, ax=ax, pad=0.02, label="locally interpolated |U| [m/s]")
    fig.savefig(output, dpi=190)
    plt.close(fig)


def surface_cp(mesh, time: Path, p_internal, p_reference: float, dynamic_pressure: float):
    wall = mesh["wall"]
    values = boundary_values(time / "p").get(wall["name"])
    if values is None and boundary_type(time / "p", wall["name"]) == "zeroGradient":
        values = [p_internal[mesh["owners"][face_id]]
                  for face_id in range(wall["startFace"], wall["startFace"] + wall["nFaces"])]
    if values is None or len(values) not in (1, wall["nFaces"]):
        raise ValueError(f"{time / 'p'}: wall pressure values do not match wall face count")
    values = values * wall["nFaces"] if len(values) == 1 else values
    by_face = {}
    for offset, face_id in enumerate(range(wall["startFace"], wall["startFace"] + wall["nFaces"])):
        xy = list(dict.fromkeys((mesh["points"][index][0], mesh["points"][index][1])
                                for index in mesh["faces"][face_id]))
        by_face[face_id] = {"x": 0.5 * (xy[0][0] + xy[1][0]),
                            "y": 0.5 * (xy[0][1] + xy[1][1]),
                            "Cp": (values[offset] - p_reference) / dynamic_pressure}
    loop = mesh["wall_loop"]
    leading = min(range(len(loop)), key=lambda i: loop[i][0])
    rotated_edges = mesh["wall_edge_face_ids"][leading:] + mesh["wall_edge_face_ids"][:leading]
    rotated_vertices = loop[leading:] + loop[:leading]
    trailing = max(range(len(rotated_vertices)), key=lambda i: rotated_vertices[i][0])
    branch_a = rotated_edges[:trailing]
    branch_b = rotated_edges[trailing:]
    if not branch_a or not branch_b:
        raise ValueError("wall loop cannot be split into upper/lower airfoil surfaces")
    if np.mean([by_face[face]["y"] for face in branch_a]) >= np.mean([by_face[face]["y"] for face in branch_b]):
        upper, lower = branch_a, branch_b
    else:
        upper, lower = branch_b, branch_a
    return ([by_face[face] for face in upper], [by_face[face] for face in lower])


def render_surface_cp(output: Path, surfaces, final_time: str, chord: float, p_reference: float):
    fig, ax = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
    for rows, label, color in ((surfaces[0], "upper surface", "#2563eb"),
                               (surfaces[1], "lower surface", "#dc2626")):
        rows = sorted(rows, key=lambda row: row["x"])
        ax.plot([row["x"] / chord for row in rows], [row["Cp"] for row in rows],
                color=color, lw=1.25, label=label)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.25)
    ax.set(xlabel="x / c", ylabel="Cp", title=f"Wall pressure coefficient at iteration {final_time}\n"
           + pressure_reference_label(p_reference) + "\nConfined-domain trial")
    ax.legend(frameon=False)
    fig.savefig(output, dpi=190)
    plt.close(fig)


def render_histories(output: Path, solver, forces, flux, yplus, discretization):
    fig, axes = plt.subplots(3, 1, figsize=(10.8, 10.5), constrained_layout=True)
    colors = {"Ux": "#2563eb", "Uy": "#0f766e", "p": "#dc2626", "k": "#7c3aed", "omega": "#d97706"}
    for field, color in colors.items():
        x = [row["iteration"] for row in solver if field in row["initialResiduals"]]
        y = [row["initialResiduals"][field] for row in solver if field in row["initialResiduals"]]
        if x:
            axes[0].semilogy(x, y, color=color, lw=1.0, label=field)
    axes[0].set(ylabel="first initial residual", title="Solver residual history")
    axes[0].grid(True, which="both", alpha=0.22)
    axes[0].legend(ncol=5, frameon=False)
    previous_stage = solver[0]["stage"]
    for row in solver[1:]:
        if row["stage"] != previous_stage:
            axes[0].axvline(row["iteration"], color="black", ls="--", lw=0.8, alpha=0.7)
            axes[0].text(row["iteration"], axes[0].get_ylim()[1], f" {row['stage']}",
                         va="top", fontsize=8)
            previous_stage = row["stage"]

    x = [row["Time"] for row in forces]
    axes[1].plot(x, [row["Cd"] for row in forces], label="Cd", color="#dc2626")
    axes[1].plot(x, [row["Cl"] for row in forces], label="Cl", color="#2563eb")
    axes[1].plot(x, [row["CmPitch"] for row in forces], label="CmPitch", color="#7c3aed")
    axes[1].set(ylabel="coefficient", title="Aerodynamic coefficient history")
    axes[1].grid(True, alpha=0.22)
    axes[1].legend(ncol=3, frameon=False)

    continuity = [row for row in solver if row["continuity"]]
    axes[2].semilogy([row["iteration"] for row in continuity],
                     [abs(row["continuity"]["global"]) for row in continuity],
                     label="|global continuity error|", color="#334155")
    ax2 = axes[2].twinx()
    # yPlus history is supplied separately below when available.
    axes[2].set(xlabel="Iteration (steady SIMPLE count, not seconds)", ylabel="continuity error",
                title="Continuity and wall y+")
    axes[2].grid(True, which="both", alpha=0.22)
    handles, labels = axes[2].get_legend_handles_labels()
    if yplus.get("history"):
        history = yplus["history"]
        ax2.plot([row["Time"] for row in history], [row["average"] for row in history],
                 color="#d97706", lw=1.1, label="y+ arithmetic average")
        ax2.fill_between([row["Time"] for row in history], [row["min"] for row in history],
                         [row["max"] for row in history], color="#f59e0b", alpha=0.16,
                         label="y+ min–max")
        ax2.set_ylabel("y+")
        h2, l2 = ax2.get_legend_handles_labels()
        handles += h2
        labels += l2
    axes[2].legend(handles, labels, frameon=False, loc="best")
    fig.suptitle("Close-domain RANS trial histories\nU: linearUpwind | k, omega: upwind", fontsize=12)
    fig.savefig(output, dpi=190)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="run root containing case/ and stages/, or the case directory")
    parser.add_argument("--output-dir", type=Path, help="default: RUN_ROOT/visualization")
    parser.add_argument("--wall", default="wall_0")
    parser.add_argument("--p-reference", type=float, default=0.0,
                        help="kinematic gauge-pressure reference for Cp (default: 0)")
    args = parser.parse_args()
    case = args.root if (args.root / "constant" / "polyMesh").is_dir() else args.root / "case"
    root = args.root.parent if case == args.root else args.root
    output = args.output_dir or root / "visualization"
    config_path = root / "configuration.json"
    config = json.loads(config_path.read_text(encoding="utf-8")) if config_path.is_file() else {}
    u_reference = float(config.get("velocity", 15.0))
    chord = float(config.get("referenceChord", 1.0))
    if not math.isfinite(u_reference) or u_reference <= 0 or not math.isfinite(chord) or chord <= 0:
        raise ValueError("configuration velocity/referenceChord must be finite and positive")
    dynamic_pressure = 0.5 * u_reference * u_reference
    time = latest_complete_time(case)
    mesh = load_mesh(case, args.wall)
    fields = read_fields(time, mesh["cell_count"])
    centres, centre_source = cell_centres(time, mesh["polygons"], mesh["cell_count"])
    if centres.shape != (mesh["cell_count"], 2) or not np.all(np.isfinite(centres)):
        raise ValueError("cell-centre array does not align with owner labels")
    cp = (fields["p"] - args.p_reference) / dynamic_pressure
    solver = [row for row in read_solver_history(root) if row["iteration"] <= float(time.name)]
    forces = [row for row in force_history(case) if row["Time"] <= float(time.name)]
    if not solver or not forces:
        raise ValueError("no solver/force history aligns at or before the final field iteration")
    flux = final_flux_summary(time, mesh)
    yplus = yplus_summary(case, time, mesh)
    ypaths = table_files(case, "yPlus", "yPlus.dat")
    yplus["history"] = sorted(
        [row for row in parse_openfoam_tables(ypaths, key_columns=2)
         if row.get("patch") == args.wall and row["Time"] <= float(time.name)],
        key=lambda row: row["Time"]) if ypaths else []
    discretization = discretization_summary(root, case)
    surfaces = surface_cp(mesh, time, fields["p"], args.p_reference, dynamic_pressure)
    output.mkdir(parents=True, exist_ok=True)
    paths = {
        "fields": output / "airfoil_fields.png",
        "streamlines": output / "airfoil_streamlines.png",
        "surfaceCp": output / "airfoil_surface_cp.png",
        "histories": output / "airfoil_histories.png",
        "summary": output / "airfoil_rans_summary.json",
    }
    scope = config.get("scope", "Close-domain steady RANS trial; not free-air or grid-independence validation.")
    render_fields(paths["fields"], mesh, fields, cp, time.name, config, args.p_reference)
    render_streamlines(paths["streamlines"], centres, mesh["polygons"], fields["U"],
                       mesh["wall_loop"], time.name)
    render_surface_cp(paths["surfaceCp"], surfaces, time.name, chord, args.p_reference)
    render_histories(paths["histories"], solver, forces, flux, yplus, discretization)

    required_residuals = {"Ux", "Uy", "p", "k", "omega"}
    last_solver = next((row for row in reversed(solver)
                        if required_residuals.issubset(row["initialResiduals"])), None)
    if last_solver is None:
        raise ValueError("no complete Time block contains first Initial residuals for Ux, Uy, p, k, omega")
    last_continuity = last_solver["continuity"]
    field_ranges = {
        "U": {"magnitude": finite_range(fields["speed"]),
              "x": finite_range(fields["U"][:, 0]), "y": finite_range(fields["U"][:, 1]),
              "z": finite_range(fields["U"][:, 2])},
        "p": finite_range(fields["p"]), "Cp": finite_range(cp),
        "k": finite_range(fields["k"]), "omega": finite_range(fields["omega"]),
        "nut": finite_range(fields["nut"]),
    }
    summary = {
        "measurementOnly": True,
        "convergenceDecision": "not made by this tool",
        "scope": scope,
        "case": str(case.resolve()), "finalFieldIteration": float(time.name),
        "iterationUnit": "steady SIMPLE iteration count, not seconds",
        "mesh": {"cells": mesh["cell_count"], "faces": len(mesh["faces"]),
                 "internalFaces": len(mesh["neighbours"]), "points": len(mesh["points"]),
                 "wallFaces": mesh["wall"]["nFaces"], "cellCentreSource": centre_source,
                 "ownerFieldAlignmentChecked": True},
        "reference": {"U": u_reference, "chord": chord, "pKinematic": args.p_reference,
                      "CpDefinition": "(p - p_reference)/(0.5*U_reference^2); p is kinematic pressure"},
        "units": {"U": "m/s", "p": "m2/s2 (kinematic)", "k": "m2/s2",
                  "omega": "1/s", "nut": "m2/s", "phi": "m3/s", "Cp": "dimensionless",
                  "forceCoefficients": "dimensionless", "yPlus": "dimensionless"},
        "fieldRanges": field_ranges,
        "fieldSignDiagnostics": {name: {"negativeCellCount": int(np.count_nonzero(fields[name] < 0.0))}
                                 for name in ("k", "omega", "nut")},
        "lastInitialResiduals": {"iteration": last_solver["iteration"], "stage": last_solver["stage"],
                                 "values": last_solver["initialResiduals"],
                                 "selection": "last complete Time block; first Initial residual for each of Ux, Uy, p, k, omega"},
        "lastContinuity": last_continuity,
        "solverLogSources": solver_log_sources(root),
        "forceCoefficients": {
            "coordinateHeaders": {
                str(path): [line for line in path.read_text().splitlines()
                            if re.match(r"#\s+(dragDir|liftDir|pitchAxis|CofR|magUInf|lRef|Aref)\s+:", line)]
                for path in table_files(case, "forces", "coefficient.dat")
            },
            "lastSampleIteration": forces[-1]["Time"],
            "lastHeaderParsedRow": {name: forces[-1][name] for name in forces[-1] if name != "Time"},
            "tailChange": {name: tail_change(forces, name) for name in ("Cd", "Cl", "CmPitch")},
        },
        "flow": flux, "yPlus": {key: value for key, value in yplus.items() if key != "history"},
        "discretization": discretization,
        "meshChecksFromLogs": {},
        "outputs": {name: str(path.resolve()) for name, path in paths.items() if name != "summary"},
        "notes": [
            "Color limits use the full finite final-field minima and maxima; values were not percentile-clipped.",
            "Streamlines use owner-indexed cell-polygon rasterization, local fluid-only interpolation, and a dilated solid mask.",
            "Force tail statistics describe variation only and are not a convergence verdict.",
            "This is a close-domain RANS trial with mixed convection discretization, not an all-second-order claim.",
        ],
    }
    standard_log = root / "stages/check-standard/stdout.log"
    expanded_log = root / "stages/check-expanded/stdout.log"
    if standard_log.is_file():
        text = standard_log.read_text(encoding="utf-8", errors="replace")
        summary["meshChecksFromLogs"]["standardMeshOK"] = "Mesh OK." in text and "Failed " not in text
    else:
        summary["meshChecksFromLogs"]["standardMeshOK"] = None
    if expanded_log.is_file():
        text = expanded_log.read_text(encoding="utf-8", errors="replace")
        failed = re.search(r"Failed (\d+) mesh checks", text)
        concave = re.search(r"Writing (\d+) concave cells to set concaveCells", text)
        summary["meshChecksFromLogs"].update({
            "expandedFailedChecks": int(failed.group(1)) if failed else 0,
            "expandedConcaveCells": int(concave.group(1)) if concave else 0,
            "expandedLog": str(expanded_log),
        })
    else:
        summary["meshChecksFromLogs"].update({"expandedFailedChecks": None,
                                               "expandedConcaveCells": None,
                                               "expandedLog": None})
    paths["summary"].write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                                encoding="utf-8")
    print(json.dumps({"finalFieldIteration": float(time.name), "cells": mesh["cell_count"],
                      "outputs": {key: str(value) for key, value in paths.items()}}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
