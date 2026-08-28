#!/usr/bin/env python3
"""Generate deterministic Q0 quality baselines without changing meshing algorithms."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Any, Iterable


FORMAT_VERSION = "cartmesh2d-q0-quality-v1"
KIND_NAMES = {0: "boundary_layer", 1: "remainder_cut",
              2: "remainder_cartesian", 3: "termination", 4: "transition"}
CASES = {
    "circle": ("examples/acceptance/circle.xy", [6, 3, 6, 4, 0.02, 1.2, 1.0]),
    "superellipse": ("examples/complex/superellipse_24.xy", [6, 3, 6, 3, 0.015, 1.15, 1.0]),
    "concave_l": ("examples/h4_3/concave_l.xy", [8, 3, 8, 4, 0.012, 1.15, 1.0]),
    "narrow_gap": ("examples/h4_3/narrow_gap.xy", [8, 3, 8, 4, 0.012, 1.15, 1.0]),
    "sharp_trailing_edge": ("examples/h4_3/sharp_trailing_edge.xy", [8, 3, 8, 4, 0.012, 1.15, 1.0]),
}


@dataclass
class Edge:
    id: int
    v0: int
    v1: int
    owner: int
    neighbour: int | None
    patch: int


@dataclass
class Cell:
    id: int
    source_id: int
    source_key: int
    geometry_area: float
    vertices: list[int]
    edges: list[int]


@dataclass
class Mesh:
    vertices: list[tuple[float, float]]
    edges: list[Edge]
    cells: list[Cell]


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command: list[str], cwd: pathlib.Path, capture: bool = False) -> str:
    completed = subprocess.run(command, cwd=cwd, check=True, text=True,
                               stdout=subprocess.PIPE if capture else None,
                               stderr=subprocess.STDOUT if capture else None)
    return completed.stdout or ""


def read_cm2d(path: pathlib.Path) -> Mesh:
    tokens = path.read_text(encoding="utf-8").split()
    cursor = 0

    def take() -> str:
        nonlocal cursor
        value = tokens[cursor]
        cursor += 1
        return value

    if take() != "CM2D" or take() != "1" or take() != "VERTICES":
        raise ValueError(f"unsupported CM2D file: {path}")
    vertices: list[tuple[float, float]] = []
    for expected in range(int(take())):
        if int(take()) != expected:
            raise ValueError("non-deterministic vertex ids")
        vertices.append((float(take()), float(take())))
    if take() != "EDGES":
        raise ValueError("missing CM2D EDGES")
    edges: list[Edge] = []
    for expected in range(int(take())):
        values = [int(take()) for _ in range(6)]
        if values[0] != expected:
            raise ValueError("non-deterministic edge ids")
        edges.append(Edge(values[0], values[1], values[2], values[3],
                          None if values[4] < 0 else values[4], values[5]))
    if take() != "CELLS":
        raise ValueError("missing CM2D CELLS")
    cells: list[Cell] = []
    for expected in range(int(take())):
        cell_id, source_id, source_key = int(take()), int(take()), int(take())
        area, vertex_count = float(take()), int(take())
        cell_vertices = [int(take()) for _ in range(vertex_count)]
        edge_count = int(take())
        cell_edges = [int(take()) for _ in range(edge_count)]
        if cell_id != expected:
            raise ValueError("non-deterministic cell ids")
        cells.append(Cell(cell_id, source_id, source_key, area,
                          cell_vertices, cell_edges))
    return Mesh(vertices, edges, cells)


def read_hybrid_kinds(path: pathlib.Path, cell_count: int) -> list[str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    start = lines.index("SCALARS hybrid_kind int 1") + 2
    values = [KIND_NAMES[int(line)] for line in lines[start:start + cell_count]]
    if len(values) != cell_count:
        raise ValueError("hybrid kind count does not match topology")
    return values


def polygon(mesh: Mesh, cell: Cell) -> list[tuple[float, float]]:
    return [mesh.vertices[index] for index in cell.vertices]


def polygon_geometry(points: list[tuple[float, float]]) -> tuple[float, tuple[float, float]]:
    twice_area = 0.0
    cx = 0.0
    cy = 0.0
    for a, b in zip(points, points[1:] + points[:1]):
        cross = a[0] * b[1] - b[0] * a[1]
        twice_area += cross
        cx += (a[0] + b[0]) * cross
        cy += (a[1] + b[1]) * cross
    if twice_area <= 0.0:
        raise ValueError("quality input contains a non-positive cell")
    return twice_area / 2.0, (cx / (3.0 * twice_area), cy / (3.0 * twice_area))


def length(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(b[0] - a[0], b[1] - a[1])


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def cell_context(mesh: Mesh, cell: Cell, source_type: str) -> dict[str, Any]:
    area, centre = polygon_geometry(polygon(mesh, cell))
    return {"cell_id": cell.id, "edge_id": None, "coordinates": list(centre),
            "cell_source_type": source_type, "source_id": cell.source_id,
            "source_key": cell.source_key, "owner": cell.id, "neighbour": None,
            "local_h": math.sqrt(area)}


def edge_context(mesh: Mesh, edge: Edge, kinds: list[str]) -> dict[str, Any]:
    a, b = mesh.vertices[edge.v0], mesh.vertices[edge.v1]
    owner_area, _ = polygon_geometry(polygon(mesh, mesh.cells[edge.owner]))
    areas = [owner_area]
    neighbour_kind = None
    if edge.neighbour is not None:
        neighbour_area, _ = polygon_geometry(polygon(mesh, mesh.cells[edge.neighbour]))
        areas.append(neighbour_area)
        neighbour_kind = kinds[edge.neighbour]
    source_type = kinds[edge.owner]
    if neighbour_kind is not None and neighbour_kind != source_type:
        source_type = f"{source_type}|{neighbour_kind}"
    return {"cell_id": edge.owner, "edge_id": edge.id,
            "coordinates": [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0],
            "cell_source_type": source_type,
            "source_id": mesh.cells[edge.owner].source_id,
            "source_key": mesh.cells[edge.owner].source_key,
            "owner": edge.owner, "neighbour": edge.neighbour,
            "local_h": sum(math.sqrt(value) for value in areas) / len(areas)}


def distribution(samples: list[tuple[float, dict[str, Any]]], worst: str) -> dict[str, Any]:
    if not samples:
        return {"count": 0, "p50": None, "p95": None, "p99": None,
                "worst": None, "worst_entity": None}
    values = [sample[0] for sample in samples]
    selected = (min if worst == "min" else max)(samples, key=lambda sample: sample[0])
    return {"count": len(values), "p50": percentile(values, 0.50),
            "p95": percentile(values, 0.95), "p99": percentile(values, 0.99),
            "worst": selected[0], "worst_direction": worst,
            "worst_entity": selected[1]}


def quality_metrics(mesh: Mesh, kinds: list[str], solver: bool) -> dict[str, Any]:
    centres: list[tuple[float, float]] = []
    areas: list[float] = []
    construction: dict[str, list[tuple[float, dict[str, Any]]]] = {
        "cell_edge_length_ratio": [], "centroid_vertex_mean_offset_normalized": [],
        "cell_area": [], "edge_length": []}
    solver_metrics: dict[str, list[tuple[float, dict[str, Any]]]] = {
        "cell_hydraulic_aspect": [], "cell_concavity_deg": [],
        "cell_min_interior_angle_deg": [], "cell_compactness": [],
        "face_length": [], "internal_non_orthogonality_deg": [],
        "internal_skewness": [], "boundary_skewness": [],
        "internal_face_weight": [], "internal_volume_ratio": []}

    for cell, kind in zip(mesh.cells, kinds):
        points = polygon(mesh, cell)
        area, centre = polygon_geometry(points)
        areas.append(area)
        centres.append(centre)
        edge_lengths = [length(a, b) for a, b in zip(points, points[1:] + points[:1])]
        context = cell_context(mesh, cell, kind)
        mean = (sum(p[0] for p in points) / len(points),
                sum(p[1] for p in points) / len(points))
        construction["cell_edge_length_ratio"].append(
            (max(edge_lengths) / min(edge_lengths), context))
        construction["centroid_vertex_mean_offset_normalized"].append(
            (length(centre, mean) / math.sqrt(area), context))
        construction["cell_area"].append((area, context))

        perimeter = sum(edge_lengths)
        angles: list[float] = []
        concavity = 0.0
        for index, current in enumerate(points):
            previous, following = points[index - 1], points[(index + 1) % len(points)]
            a = (previous[0] - current[0], previous[1] - current[1])
            b = (following[0] - current[0], following[1] - current[1])
            base = math.degrees(math.acos(max(-1.0, min(1.0,
                (a[0] * b[0] + a[1] * b[1]) /
                (math.hypot(*a) * math.hypot(*b))))))
            orientation = ((current[0] - previous[0]) * (following[1] - current[1]) -
                           (current[1] - previous[1]) * (following[0] - current[0]))
            angle = base if orientation >= 0.0 else 360.0 - base
            angles.append(angle)
            concavity = max(concavity, angle - 180.0)
        solver_metrics["cell_hydraulic_aspect"].append(
            (max(edge_lengths) / (2.0 * area / perimeter), context))
        solver_metrics["cell_concavity_deg"].append((concavity, context))
        solver_metrics["cell_min_interior_angle_deg"].append((min(angles), context))
        solver_metrics["cell_compactness"].append(
            (4.0 * math.pi * area / (perimeter * perimeter), context))

    for edge in mesh.edges:
        a, b = mesh.vertices[edge.v0], mesh.vertices[edge.v1]
        edge_len = length(a, b)
        context = edge_context(mesh, edge, kinds)
        construction["edge_length"].append((edge_len, context))
        solver_metrics["face_length"].append((edge_len, context))
        owner = centres[edge.owner]
        face = ((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)
        e = (b[0] - a[0], b[1] - a[1])
        normal = (e[1], -e[0])
        if edge.neighbour is None:
            owner_face = (face[0] - owner[0], face[1] - owner[1])
            normal_mag = math.hypot(*normal)
            unit = (normal[0] / normal_mag, normal[1] / normal_mag)
            projection_scale = owner_face[0] * unit[0] + owner_face[1] * unit[1]
            projection = (unit[0] * projection_scale, unit[1] * projection_scale)
            tangent = (owner_face[0] - projection[0], owner_face[1] - projection[1])
            skew_mag = math.hypot(*tangent)
            face_distance = 0.0
            if skew_mag > 0.0:
                direction = (tangent[0] / skew_mag, tangent[1] / skew_mag)
                face_distance = max(abs(direction[0] * (a[0] - face[0]) +
                                        direction[1] * (a[1] - face[1])),
                                    abs(direction[0] * (b[0] - face[0]) +
                                        direction[1] * (b[1] - face[1])))
            normalization = max(0.4 * math.hypot(*projection), face_distance)
            solver_metrics["boundary_skewness"].append(
                (skew_mag / max(normalization, sys.float_info.min), context))
            continue

        neighbour = centres[edge.neighbour]
        d = (neighbour[0] - owner[0], neighbour[1] - owner[1])
        dot_dn = d[0] * normal[0] + d[1] * normal[1]
        non_orth = math.degrees(math.acos(max(0.0, min(1.0,
            abs(dot_dn) / (math.hypot(*d) * math.hypot(*normal))))))
        solver_metrics["internal_non_orthogonality_deg"].append((non_orth, context))
        line_denom = e[0] * d[1] - e[1] * d[0]
        if abs(line_denom) > 1.0e-30:
            owner_a = (owner[0] - a[0], owner[1] - a[1])
            t = (owner_a[0] * d[1] - owner_a[1] * d[0]) / line_denom
            intersection = (a[0] + e[0] * t, a[1] + e[1] * t)
            skew = length(intersection, face) / max(0.2 * math.hypot(*d), 0.5 * edge_len)
        else:
            skew = float("inf")
        solver_metrics["internal_skewness"].append((skew, context))
        own_distance = abs(normal[0] * (face[0] - owner[0]) +
                           normal[1] * (face[1] - owner[1]))
        nei_distance = abs(normal[0] * (neighbour[0] - face[0]) +
                           normal[1] * (neighbour[1] - face[1]))
        solver_metrics["internal_face_weight"].append(
            (min(own_distance, nei_distance) / (own_distance + nei_distance), context))
        solver_metrics["internal_volume_ratio"].append(
            (min(areas[edge.owner], areas[edge.neighbour]) /
             max(areas[edge.owner], areas[edge.neighbour]), context))

    selected = solver_metrics if solver else construction
    low_worst = {"cell_area", "edge_length", "cell_min_interior_angle_deg",
                 "cell_compactness", "face_length", "internal_face_weight",
                 "internal_volume_ratio"}
    return {name: distribution(samples, "min" if name in low_worst else "max")
            for name, samples in selected.items()}


def read_foam_scalar_field(path: pathlib.Path, count: int) -> list[float]:
    text = path.read_text(encoding="utf-8")
    nonuniform = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
        text, re.DOTALL)
    if nonuniform:
        declared = int(nonuniform.group(1))
        values = [float(value) for value in nonuniform.group(2).split()]
        if declared != count or len(values) != count:
            raise ValueError(f"OpenFOAM field size mismatch: {path}")
        return values
    uniform = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", text)
    if uniform:
        return [float(uniform.group(1))] * count
    raise ValueError(f"cannot parse OpenFOAM scalar field: {path}")


def parse_checkmesh(log: str, foam_case: pathlib.Path, mesh: Mesh,
                    kinds: list[str]) -> dict[str, Any]:
    numeric = r"([-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)"

    def number(pattern: str) -> float | None:
        match = re.search(pattern, log, re.IGNORECASE)
        return float(match.group(1)) if match else None

    version_match = re.search(r"Version:\s*([^\s]+)", log)
    summary_metrics = {
        "max_aspect_ratio": number(r"Max aspect ratio\s*=\s*" + numeric),
        "max_non_orthogonality_deg": number(r"non-orthogonality Max:\s*" + numeric),
        "max_skewness": number(r"Max skewness\s*=\s*" + numeric),
        "min_cell_volume": number(r"Min volume\s*=\s*" + numeric),
    }
    fields = {
        "aspect_ratio": ("aspectRatio", "max"),
        "approximate_cell_aspect_ratio": ("cellAspectRatio", "max"),
        "non_orthogonality_deg": ("nonOrthoAngle", "max"),
        "skewness": ("skewness", "max"),
        "face_weight": ("faceWeight", "min"),
        "cell_volume_ratio": ("cellVolumeRatio", "min"),
        "cell_determinant": ("cellDeterminant", "min"),
        "cell_volume": ("cellVolume", "min"),
    }
    external: dict[str, Any] = {}
    for metric_name, (field_name, direction) in fields.items():
        values = read_foam_scalar_field(foam_case / "0" / field_name, len(mesh.cells))
        samples = [(value, cell_context(mesh, cell, kind))
                   for value, cell, kind in zip(values, mesh.cells, kinds)]
        external[metric_name] = distribution(samples, direction)
        external[metric_name]["openfoam_field"] = field_name
    return {"quality_class": "openfoam_quality", "source": "external_checkMesh",
            "openfoam_version": version_match.group(1) if version_match else "v2606",
            "status": "pass" if "Mesh OK" in log else "fail",
            "summary_extrema": summary_metrics,
            "distribution_source": "checkMesh -writeAllFields volScalarField",
            "metrics": external}


def normalized_log_hash(log: str) -> str:
    lines = [line for line in log.splitlines()
             if not line.startswith(("Time =", "ExecutionTime =", "ClockTime ="))]
    return hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest()


def generate_case(repo: pathlib.Path, executable: pathlib.Path, evidence: pathlib.Path,
                  name: str, image: str, skip_openfoam: bool) -> pathlib.Path:
    boundary, config = CASES[name]
    prefix = evidence / name
    foam_case = evidence / f"{name}-case"
    command = [str(executable), str(repo / boundary), str(prefix),
               *[str(value) for value in config], str(foam_case), "0.01"]
    run(command, repo)
    run([sys.executable, str(repo / "tools/verification/check_openfoam2d.py"),
         str(foam_case), "--report", str(foam_case / "independent_check.json")], repo)
    log_path = evidence / f"{name}.checkMesh.log"
    if not skip_openfoam:
        container_root = "/home/openfoam/workingDir"
        relative_case = foam_case.relative_to(repo)
        log = run(["docker", "run", "--rm", "-v", f"{repo}:{container_root}",
                   "-w", container_root, image, "checkMesh", "-case",
                   f"{container_root}/{relative_case}", "-writeAllFields"],
                  repo, capture=True)
        log_path.write_text(log, encoding="utf-8")
        if "Mesh OK" not in log:
            raise RuntimeError(f"OpenFOAM rejected {name}")
    elif not log_path.exists():
        raise RuntimeError(f"--skip-openfoam requires {log_path}")
    return log_path


def build_report(repo: pathlib.Path, evidence: pathlib.Path, name: str,
                 source_commit: str, log_path: pathlib.Path) -> dict[str, Any]:
    boundary, config = CASES[name]
    construction_path = evidence / f"{name}.hybrid.cm2d"
    solver_path = evidence / f"{name}.hybrid.solver.cm2d"
    vtk_path = evidence / f"{name}.hybrid.vtk"
    construction_mesh = read_cm2d(construction_path)
    construction_kinds = read_hybrid_kinds(vtk_path, len(construction_mesh.cells))
    source_kinds: dict[int, str] = {}
    for cell, kind in zip(construction_mesh.cells, construction_kinds):
        previous = source_kinds.get(cell.source_id)
        source_kinds[cell.source_id] = kind if previous in (None, kind) else "mixed"
    solver_mesh = read_cm2d(solver_path)
    solver_kinds = [source_kinds.get(cell.source_id, "solver_repartition")
                    for cell in solver_mesh.cells]
    log = log_path.read_text(encoding="utf-8")
    command = ["cartmesh2d_hybrid_cli", boundary, f"$EVIDENCE/{name}",
               *[str(value) for value in config], f"$EVIDENCE/{name}-case", "0.01"]
    return {
        "format_version": FORMAT_VERSION,
        "case": name,
        "generator_commit": source_commit,
        "generation": {"command": command,
                       "config": {"max_level": config[0], "minimum_level": config[1],
                                  "boundary_level": config[2], "layers": config[3],
                                  "first_thickness": config[4], "growth_ratio": config[5],
                                  "domain_padding": config[6], "extrusion_thickness": 0.01}},
        "input": {"path": boundary, "sha256": sha256(repo / boundary)},
        "construction_quality": {
            "quality_class": "construction_quality",
            "source": "native_unified_2d_topology",
            "metric_semantics": {
                "cell_edge_length_ratio": "max edge / min edge within one polygon",
                "centroid_vertex_mean_offset_normalized":
                    "distance(area centroid, vertex mean) / sqrt(area)"},
            "metrics": quality_metrics(construction_mesh, construction_kinds, False)},
        "solver_quality": {
            "quality_class": "solver_quality",
            "source": "native_solver_repaired_topology",
            "metric_semantics":
                "native 2D solver metrics; not external OpenFOAM evidence",
            "metrics": quality_metrics(solver_mesh, solver_kinds, True)},
        "openfoam_quality": parse_checkmesh(
            log, evidence / f"{name}-case", solver_mesh, solver_kinds),
        "evidence": {
            "construction_topology_sha256": sha256(construction_path),
            "solver_topology_sha256": sha256(solver_path),
            "openfoam_poly_mesh": {
                filename: sha256(evidence / f"{name}-case/constant/polyMesh/{filename}")
                for filename in ("points", "faces", "owner", "neighbour", "boundary")},
            "openfoam_quality_fields": {
                filename: sha256(evidence / f"{name}-case/0/{filename}")
                for filename in ("aspectRatio", "cellAspectRatio", "nonOrthoAngle",
                                 "skewness", "faceWeight", "cellVolumeRatio",
                                 "cellDeterminant", "cellVolume")},
            "normalized_checkmesh_log_sha256": normalized_log_hash(log)},
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--build-dir", type=pathlib.Path, default=pathlib.Path("build/q0"))
    parser.add_argument("--evidence-dir", type=pathlib.Path,
                        default=pathlib.Path("build/q0_evidence"))
    parser.add_argument("--output-dir", type=pathlib.Path,
                        default=pathlib.Path("artifacts/q0"))
    parser.add_argument("--source-commit")
    parser.add_argument("--openfoam-image", default="opencfd/openfoam-run:2606")
    parser.add_argument("--skip-openfoam", action="store_true")
    parser.add_argument("--collect-only", action="store_true",
                        help="reuse existing mesh/checkMesh evidence without regeneration")
    args = parser.parse_args()
    repo = args.repo.resolve()
    build_dir = (repo / args.build_dir).resolve() if not args.build_dir.is_absolute() else args.build_dir
    evidence = (repo / args.evidence_dir).resolve() if not args.evidence_dir.is_absolute() else args.evidence_dir
    output = (repo / args.output_dir).resolve() if not args.output_dir.is_absolute() else args.output_dir
    source_commit = args.source_commit or run(
        ["git", "rev-parse", "HEAD"], repo, capture=True).strip()
    executable = build_dir / "cartmesh2d_hybrid_cli"
    if not executable.exists():
        raise SystemExit(f"missing executable: {executable}")
    if evidence.exists() and not args.skip_openfoam and not args.collect_only:
        shutil.rmtree(evidence)
    evidence.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)

    reports: dict[str, dict[str, Any]] = {}
    for name in CASES:
        if args.collect_only:
            log_path = evidence / f"{name}.checkMesh.log"
            if not log_path.exists():
                raise SystemExit(f"missing collected checkMesh log: {log_path}")
        else:
            log_path = generate_case(repo, executable, evidence, name,
                                     args.openfoam_image, args.skip_openfoam)
        report = build_report(repo, evidence, name, source_commit, log_path)
        report_path = output / f"{name}.quality-baseline.json"
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                               encoding="utf-8")
        reports[name] = report

    manifest = {
        "format_version": FORMAT_VERSION,
        "generator_commit": source_commit,
        "openfoam_image": args.openfoam_image,
        "openfoam_version": reports["circle"]["openfoam_quality"]["openfoam_version"],
        "cases": {
            name: {"input_sha256": report["input"]["sha256"],
                   "report_sha256": sha256(output / f"{name}.quality-baseline.json"),
                   "evidence": report["evidence"]}
            for name, report in reports.items()}}
    (output / "provenance-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"generated {len(reports)} Q0 baselines in {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
