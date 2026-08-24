#!/usr/bin/env python3
"""Independent structural/geometry check for an extruded cartmesh2d polyMesh.

This deliberately does not import cartmesh2d.  It is a product reader and is
not a substitute for OpenFOAM checkMesh, which remains the V1 acceptance gate.
"""

import argparse
import json
import math
import re
import re
from pathlib import Path


def payload(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    end = text.find("}", text.find("FoamFile"))
    if end < 0:
        raise ValueError(f"{path.name}: missing FoamFile header")
    return text[end + 1 :]


def counted_body(path: Path):
    text = payload(path)
    match = re.search(r"\b(\d+)\s*\n\s*\((.*)\)\s*$", text, re.S)
    if not match:
        raise ValueError(f"{path.name}: malformed counted list")
    return int(match.group(1)), match.group(2)


def read_points(path: Path):
    count, body = counted_body(path)
    points = [tuple(map(float, row.split())) for row in re.findall(r"\(([^()]*)\)", body)]
    if len(points) != count or any(len(p) != 3 for p in points):
        raise ValueError("points: count or vector width mismatch")
    return points


def read_faces(path: Path):
    count, body = counted_body(path)
    faces = []
    for declared, row in re.findall(r"(\d+)\s*\(([^()]*)\)", body):
        vertices = [int(value) for value in row.split()]
        if len(vertices) != int(declared):
            raise ValueError("faces: declared face size mismatch")
        faces.append(vertices)
    if len(faces) != count:
        raise ValueError("faces: count mismatch")
    return faces


def read_labels(path: Path):
    count, body = counted_body(path)
    labels = [int(value) for value in re.findall(r"\b\d+\b", body)]
    if len(labels) != count:
        raise ValueError(f"{path.name}: count mismatch")
    return labels


def read_boundary(path: Path):
    count, body = counted_body(path)
    patches = []
    pattern = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)\s*\{(.*?)\}", re.S)
    for name, fields in pattern.findall(body):
        field = dict(re.findall(r"(type|nFaces|startFace)\s+([^;]+);", fields))
        if set(field) != {"type", "nFaces", "startFace"}:
            raise ValueError(f"boundary: incomplete patch {name}")
        patches.append({"name": name, "type": field["type"].strip(),
                        "nFaces": int(field["nFaces"]),
                        "startFace": int(field["startFace"])})
    if len(patches) != count:
        raise ValueError("boundary: patch count mismatch")
    return patches


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(a, factor):
    return tuple(x * factor for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def face_geometry(face, points):
    origin = points[face[0]]
    area = (0.0, 0.0, 0.0)
    for index in range(1, len(face) - 1):
        area = add(area, scale(cross(sub(points[face[index]], origin),
                                     sub(points[face[index + 1]], origin)), 0.5))
    centroid = tuple(sum(points[index][axis] for index in face) / len(face)
                     for axis in range(3))
    return area, centroid


def check(case: Path):
    mesh = case / "constant" / "polyMesh"
    points = read_points(mesh / "points")
    faces = read_faces(mesh / "faces")
    owner = read_labels(mesh / "owner")
    neighbour = read_labels(mesh / "neighbour")
    patches = read_boundary(mesh / "boundary")
    issues = []
    if len(owner) != len(faces):
        issues.append("owner count differs from face count")
    if any(len(face) < 3 or any(v < 0 or v >= len(points) for v in face) for face in faces):
        issues.append("face has invalid arity or vertex index")
    cell_count = 1 + max(owner + neighbour, default=-1)
    if any(owner[index] >= neighbour[index] for index in range(len(neighbour))):
        issues.append("internal face violates owner < neighbour ordering")
    internal_pairs = [(owner[index], neighbour[index]) for index in range(len(neighbour))]
    if internal_pairs != sorted(internal_pairs):
        issues.append("internal faces are not in OpenFOAM upper-triangular order")
    next_face = len(neighbour)
    for patch in patches:
        if patch["startFace"] != next_face:
            issues.append(f"patch {patch['name']} is not contiguous")
        next_face += patch["nFaces"]
    if next_face != len(faces):
        issues.append("boundary patches do not cover all boundary faces")
    empty = [p for p in patches if p["name"] == "frontAndBack"]
    if len(empty) != 1 or empty[0]["type"] != "empty" or empty[0]["nFaces"] != 2 * cell_count:
        issues.append("frontAndBack is not an empty two-face-per-cell patch")

    closure = [(0.0, 0.0, 0.0) for _ in range(cell_count)]
    volumes = [0.0 for _ in range(cell_count)]
    used = [False for _ in points]
    for face_id, face in enumerate(faces):
        for vertex in face:
            if 0 <= vertex < len(used):
                used[vertex] = True
        if any(vertex < 0 or vertex >= len(points) for vertex in face) or len(face) < 3:
            continue
        area, centroid = face_geometry(face, points)
        own = owner[face_id]
        closure[own] = add(closure[own], area)
        volumes[own] += dot(centroid, area) / 3.0
        if face_id < len(neighbour):
            nei = neighbour[face_id]
            closure[nei] = add(closure[nei], scale(area, -1.0))
            volumes[nei] -= dot(centroid, area) / 3.0
    scale_length = max((abs(value) for point in points for value in point), default=1.0)
    closure_tolerance = max(1.0, scale_length * scale_length) * 1.0e-10
    bad_closure = [i for i, value in enumerate(closure)
                   if math.sqrt(dot(value, value)) > closure_tolerance]
    bad_volume = [i for i, value in enumerate(volumes) if not value > 0.0]
    if bad_closure:
        issues.append(f"{len(bad_closure)} cells have non-closed oriented faces")
    if bad_volume:
        issues.append(f"{len(bad_volume)} cells have non-positive volume")
    if not all(used):
        issues.append(f"{used.count(False)} points are unreferenced")
    constant_flux = (0.7, -0.3, 0.0)
    divergence_residual = [abs(dot(value, constant_flux)) / volume
                           for value, volume in zip(closure, volumes) if volume > 0.0]
    physical_boundary_area = (0.0, 0.0, 0.0)
    front = next((patch for patch in patches if patch["name"] == "frontAndBack"), None)
    front_faces = (set(range(front["startFace"], front["startFace"] + front["nFaces"]))
                   if front else set())
    for face_id in range(len(neighbour), len(faces)):
        if face_id not in front_faces:
            area, _ = face_geometry(faces[face_id], points)
            physical_boundary_area = add(physical_boundary_area, area)
    return {"valid": not issues, "point_count": len(points), "face_count": len(faces),
            "internal_face_count": len(neighbour), "cell_count": cell_count,
            "min_volume": min(volumes, default=0.0),
            "max_closure_residual": max((math.sqrt(dot(v, v)) for v in closure), default=0.0),
            "max_constant_flux_divergence_residual": max(divergence_residual, default=0.0),
            "physical_boundary_area_vector": physical_boundary_area,
            "global_constant_flux_balance": dot(physical_boundary_area, constant_flux),
            "patches": patches, "issues": issues}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--require-patch", action="append", default=[],
                        help="require an exact boundary patch as name:type")
    parser.add_argument("--require-field-boundary", action="append", default=[],
                        help="require name:type in both 0/U and 0/p")
    args = parser.parse_args()
    try:
        report = check(args.case)
    except (OSError, ValueError) as exc:
        report = {"valid": False, "issues": [str(exc)]}
    if report.get("valid"):
        available = {f'{patch["name"]}:{patch["type"]}' for patch in report["patches"]}
        missing = sorted(set(args.require_patch) - available)
        if missing:
            report["valid"] = False
            report["issues"].append("missing required patches: " + ", ".join(missing))
        for required in args.require_field_boundary:
            if ":" not in required:
                report["valid"] = False
                report["issues"].append("invalid required field boundary: " + required)
                continue
            name, boundary_type = required.split(":", 1)
            pattern = re.compile(rf"\b{re.escape(name)}\s*\{{[^}}]*\btype\s+"
                                 rf"{re.escape(boundary_type)}\s*;", re.DOTALL)
            for relative in (Path("0/U"), Path("0/p")):
                try:
                    content = (args.case / relative).read_text(encoding="utf-8")
                except OSError as exc:
                    report["valid"] = False
                    report["issues"].append(str(exc))
                    continue
                if not pattern.search(content):
                    report["valid"] = False
                    report["issues"].append(
                        f"missing {required} in {relative.as_posix()}")
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    raise SystemExit(0 if report["valid"] else 1)


if __name__ == "__main__":
    main()
