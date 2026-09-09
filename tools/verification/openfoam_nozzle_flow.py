#!/usr/bin/env python3
"""Prepare and evaluate a fixed low-Re nozzle-flow OpenFOAM case.

``configure`` only copies and edits cases; it never starts an OpenFOAM solver.
The inlet/outlet segment IDs and coordinates are explicit properties of the
checked nozzle input, rather than an x-extreme heuristic.  ``evaluate`` reads
real final-time ``phi``, ``p`` and ``U`` boundary values after a later solve.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
INPUT = REPO / "examples/complex/nozzle_profile.xy"
INPUT_SHA256 = "8a7ee17ec10a01ea6f628839ae9c28162e9ab2de80f0b799001f6c450692a18e"
THICKNESS = 0.02
U_IN = (0.1, 0.0, 0.0)
NU = 0.01
P_OUT = 0.0
# Explicit source loop segment IDs: line 40 is x=+3, line 81 closes x=-3.
OUTLET_SEGMENT = 40
INLET_SEGMENT = 81
INLET_ENDPOINTS = ((-3.0, 0.98), (-3.0, -0.98))
OUTLET_ENDPOINTS = ((3.0, -0.98), (3.0, 0.98))
EPS = 1.0e-9


def payload(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    at = text.find("FoamFile")
    end = text.find("}", at)
    if at < 0 or end < 0:
        raise ValueError(f"{path}: missing FoamFile header")
    return text[end + 1:]


def counted(path: Path) -> tuple[int, str]:
    match = re.search(r"\b(\d+)\s*\n\s*\((.*)\)\s*$", payload(path), re.S)
    if not match:
        raise ValueError(f"{path}: malformed counted list")
    return int(match.group(1)), match.group(2)


def read_points(path: Path):
    count, body = counted(path)
    rows = [tuple(map(float, row.split())) for row in re.findall(r"\(([^()]*)\)", body)]
    if len(rows) != count or any(len(row) != 3 for row in rows):
        raise ValueError(f"{path}: point count/vector mismatch")
    return rows


def read_faces(path: Path):
    count, body = counted(path)
    faces = []
    for declared, row in re.findall(r"(\d+)\s*\(([^()]*)\)", body):
        face = [int(value) for value in row.split()]
        if len(face) != int(declared):
            raise ValueError(f"{path}: face arity mismatch")
        faces.append(face)
    if len(faces) != count:
        raise ValueError(f"{path}: face count mismatch")
    return faces


def read_labels(path: Path):
    count, body = counted(path)
    values = [int(value) for value in re.findall(r"\b-?\d+\b", body)]
    if len(values) != count:
        raise ValueError(f"{path}: label count mismatch")
    return values


def read_boundary(path: Path):
    count, body = counted(path)
    patches = []
    for name, fields in re.findall(r"([A-Za-z_][A-Za-z0-9_]*)\s*\{(.*?)\}", body, re.S):
        values = dict(re.findall(r"(type|nFaces|startFace)\s+([^;]+);", fields))
        if set(values) != {"type", "nFaces", "startFace"}:
            raise ValueError(f"{path}: incomplete patch {name}")
        patches.append({"name": name, "type": values["type"].strip(),
                        "nFaces": int(values["nFaces"]),
                        "startFace": int(values["startFace"])})
    if len(patches) != count:
        raise ValueError(f"{path}: patch count mismatch")
    return patches


def header(cls: str, location: str, obj: str) -> str:
    return ("FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            f"    class {cls};\n    location \"{location}\";\n    object {obj};\n}}\n\n")


def write_faces(path: Path, faces) -> None:
    text = [header("faceList", "constant/polyMesh", "faces"), str(len(faces)), "\n(\n"]
    text.extend(f"{len(face)}({' '.join(map(str, face))})\n" for face in faces)
    text.append(")\n")
    path.write_text("".join(text), encoding="utf-8")


def write_labels(path: Path, cls: str, obj: str, values) -> None:
    path.write_text(header(cls, "constant/polyMesh", obj) + str(len(values)) +
                    "\n(\n" + "".join(f"{value}\n" for value in values) + ")\n",
                    encoding="utf-8")


def write_boundary(path: Path, patches) -> None:
    text = [header("polyBoundaryMesh", "constant/polyMesh", "boundary"),
            str(len(patches)), "\n(\n"]
    for patch in patches:
        text.append(f"{patch['name']}\n{{\n    type {patch['type']};\n"
                    f"    nFaces {patch['nFaces']};\n    startFace {patch['startFace']};\n}}\n")
    text.append(")\n")
    path.write_text("".join(text), encoding="utf-8")


def source_points() -> list[tuple[float, float]]:
    digest = hashlib.sha256(INPUT.read_bytes()).hexdigest()
    if digest != INPUT_SHA256:
        raise ValueError("nozzle input SHA-256 differs from fixed verified source")
    rows = []
    for line in INPUT.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 2 and not fields[0].startswith("#"):
            rows.append((float(fields[0]), float(fields[1])))
    if len(rows) != 82:
        raise ValueError("fixed nozzle source must contain exactly 82 vertices")
    for segment, expected in ((INLET_SEGMENT, INLET_ENDPOINTS),
                              (OUTLET_SEGMENT, OUTLET_ENDPOINTS)):
        actual = (rows[segment], rows[(segment + 1) % len(rows)])
        if actual != expected:
            raise ValueError(f"fixed nozzle segment {segment} endpoints changed")
    return rows


def on_segment(point, a, b, eps=EPS) -> bool:
    cross = (b[0] - a[0]) * (point[1] - a[1]) - (b[1] - a[1]) * (point[0] - a[0])
    return (abs(cross) <= eps and min(a[0], b[0]) - eps <= point[0] <= max(a[0], b[0]) + eps
            and min(a[1], b[1]) - eps <= point[1] <= max(a[1], b[1]) + eps)


def patch_face_groups(case: Path, geometry: list[tuple[float, float]]):
    mesh = case / "constant/polyMesh"
    points = read_points(mesh / "points")
    z_planes = sorted({point[2] for point in points})
    if len(z_planes) != 2 or not all(math.isfinite(z) for z in z_planes) or z_planes[0] != 0.0 or z_planes[1] <= 0.0:
        raise ValueError("source case must have exactly two finite extrusion planes z=0 and z>0")
    faces = read_faces(mesh / "faces")
    owners = read_labels(mesh / "owner")
    neighbours = read_labels(mesh / "neighbour")
    patches = read_boundary(mesh / "boundary")
    if len(owners) != len(faces) or len(owners) < len(neighbours):
        raise ValueError("mesh owner/face/neighbour counts disagree")
    wall = next((p for p in patches if p["name"] == "wall_0"), None)
    if wall is None or wall["type"] != "wall":
        raise ValueError("source case must contain wall_0 wall patch")
    inlet = (geometry[INLET_SEGMENT], geometry[(INLET_SEGMENT + 1) % len(geometry)])
    outlet = (geometry[OUTLET_SEGMENT], geometry[(OUTLET_SEGMENT + 1) % len(geometry)])
    groups = {"inlet": [], "outlet": [], "wall_0": []}
    for face_id in range(wall["startFace"], wall["startFace"] + wall["nFaces"]):
        xy = {(round(points[index][0], 12), round(points[index][1], 12)) for index in faces[face_id]}
        if not xy:
            raise ValueError("empty wall face")
        matches = []
        for name, segment in (("inlet", inlet), ("outlet", outlet)):
            if all(on_segment(point, segment[0], segment[1]) for point in xy):
                matches.append(name)
        if len(matches) > 1:
            raise ValueError(f"wall face {face_id} has ambiguous inlet/outlet classification")
        groups[matches[0] if matches else "wall_0"].append(face_id)
    if not groups["inlet"] or not groups["outlet"]:
        raise ValueError("explicit inlet/outlet source segments matched no wall faces")
    expected_wall = wall["nFaces"]
    if sum(len(value) for value in groups.values()) != expected_wall:
        raise ValueError("wall face classification does not conserve wall face count")
    order = groups["inlet"] + groups["outlet"] + groups["wall_0"]
    boundary_start = len(neighbours)
    boundary_ids = list(range(boundary_start, len(faces)))
    if set(order) != set(boundary_ids) - set(range(wall["startFace"] + wall["nFaces"], boundary_start)):
        # This source exporter has wall_0 immediately after internal faces;
        # reject unexpected boundary patches instead of silently dropping them.
        if any(p["name"] not in ("wall_0", "frontAndBack") for p in patches):
            raise ValueError("unexpected source boundary patch")
    front = next((p for p in patches if p["name"] == "frontAndBack"), None)
    if front is None or front["type"] != "empty":
        raise ValueError("source case must contain frontAndBack empty patch")
    front_ids = list(range(front["startFace"], front["startFace"] + front["nFaces"]))
    new_boundary = order + front_ids
    if len(new_boundary) != len(set(new_boundary)) or set(new_boundary) != set(boundary_ids):
        raise ValueError("boundary face partition is incomplete or duplicated")
    new_order = list(range(boundary_start)) + new_boundary
    write_faces(mesh / "faces", [faces[index] for index in new_order])
    write_labels(mesh / "owner", "labelList", "owner", [owners[index] for index in new_order])
    new_patches = []
    cursor = boundary_start
    for name, typ, ids in (("inlet", "patch", groups["inlet"]),
                           ("outlet", "patch", groups["outlet"]),
                           ("wall", "wall", groups["wall_0"]),
                           ("frontAndBack", "empty", front_ids)):
        new_patches.append({"name": name, "type": typ, "nFaces": len(ids), "startFace": cursor})
        cursor += len(ids)
    write_boundary(mesh / "boundary", new_patches)
    # Rescale only the extrusion coordinate in the copied case.
    points_text = [header("vectorField", "constant/polyMesh", "points"), str(len(points)), "\n(\n"]
    for x, y, z in points:
        mapped_z = 0.0 if z == z_planes[0] else THICKNESS
        points_text.append(f"({x:.17g} {y:.17g} {mapped_z:.17g})\n")
    points_text.append(")\n")
    (mesh / "points").write_text("".join(points_text), encoding="utf-8")
    return {"point_count": len(points), "face_count": len(faces),
            "internal_face_count": len(neighbours), "wall_source_faces": expected_wall,
            "wall_faces": len(groups["wall_0"]),
            "inlet_faces": len(groups["inlet"]), "outlet_faces": len(groups["outlet"]),
            "front_back_faces": len(front_ids)}


def configure_fields(case: Path) -> None:
    (case / "constant/transportProperties").write_text(
        header("dictionary", "constant", "transportProperties") +
        "transportModel Newtonian;\nnu [0 2 -1 0 0 0 0] 0.01;\n", encoding="utf-8")
    (case / "constant/turbulenceProperties").write_text(
        header("dictionary", "constant", "turbulenceProperties") +
        "simulationType laminar;\n", encoding="utf-8")
    (case / "system/controlDict").write_text(
        header("dictionary", "system", "controlDict") +
        "application simpleFoam;\nstartFrom startTime;\nstartTime 0;\n"
        "stopAt endTime;\nendTime 2000;\ndeltaT 1;\nwriteControl timeStep;\n"
        "writeInterval 100;\nwriteFormat ascii;\nwritePrecision 12;\n"
        "writeCompression off;\ntimeFormat general;\ntimePrecision 8;\n"
        "runTimeModifiable false;\nfunctions\n{\n"
        "    inletFlux { type surfaceFieldValue; libs (\"libfieldFunctionObjects.so\"); "
        "regionType patch; name inlet; operation sum; fields (phi); writeFields false; writeControl timeStep; writeInterval 10; }\n"
        "    outletFlux { type surfaceFieldValue; libs (\"libfieldFunctionObjects.so\"); "
        "regionType patch; name outlet; operation sum; fields (phi); writeFields false; writeControl timeStep; writeInterval 10; }\n"
        "    inletPressure { type surfaceFieldValue; libs (\"libfieldFunctionObjects.so\"); "
        "regionType patch; name inlet; operation areaAverage; fields (p); writeFields false; writeControl timeStep; writeInterval 10; }\n"
        "    outletPressure { type surfaceFieldValue; libs (\"libfieldFunctionObjects.so\"); "
        "regionType patch; name outlet; operation areaAverage; fields (p); writeFields false; writeControl timeStep; writeInterval 10; }\n"
        "    outletProfile { type sets; libs (\"libsampling.so\"); setFormat raw; fields (U p); "
        "interpolationScheme cellPoint; sets { outlet { type uniform; axis y; start (3 -0.98 0.01); "
        "end (3 0.98 0.01); nPoints 200; } } writeFields false; writeControl timeStep; writeInterval 10; }\n}\n",
        encoding="utf-8")
    (case / "system/fvSchemes").write_text(
        header("dictionary", "system", "fvSchemes") +
        "ddtSchemes { default steadyState; }\n"
        "gradSchemes { default cellLimited Gauss linear 1; }\n"
        "divSchemes { default none; div(phi,U) bounded Gauss upwind; "
        "div((nuEff*dev2(T(grad(U))))) Gauss linear; }\n"
        "laplacianSchemes { default Gauss linear corrected; }\n"
        "interpolationSchemes { default linear; }\n"
        "snGradSchemes { default corrected; }\nwallDist { method meshWave; }\n", encoding="utf-8")
    (case / "system/fvSolution").write_text(
        header("dictionary", "system", "fvSolution") +
        "solvers { p { solver GAMG; tolerance 1e-10; relTol 0.02; smoother GaussSeidel; } "
        "U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.02; } }\n"
        "SIMPLE { nNonOrthogonalCorrectors 1; residualControl { p 1e-8; U 1e-8; } }\n"
        "relaxationFactors { fields { p 0.3; } equations { U 0.6; } }\n", encoding="utf-8")
    (case / "0/U").write_text(
        header("volVectorField", "0", "U") +
        "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0.1 0 0);\n"
        "boundaryField { inlet { type fixedValue; value uniform (0.1 0 0); } "
        "outlet { type zeroGradient; } wall { type noSlip; } frontAndBack { type empty; } }\n",
        encoding="utf-8")
    (case / "0/p").write_text(
        header("volScalarField", "0", "p") +
        "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\n"
        "boundaryField { inlet { type zeroGradient; } outlet { type fixedValue; value uniform 0; } "
        "wall { type zeroGradient; } frontAndBack { type empty; } }\n", encoding="utf-8")


def configure(args) -> int:
    geometry = source_points()
    output_root = args.output_root
    output_root.mkdir(parents=True, exist_ok=True)
    results = []
    for name, source in args.case:
        if not source.is_dir():
            raise ValueError(f"source case does not exist: {source}")
        destination = output_root / name
        if destination.exists():
            raise ValueError(f"refusing to overwrite existing configured case: {destination}")
        shutil.copytree(source, destination)
        mesh_result = patch_face_groups(destination, geometry)
        configure_fields(destination)
        zs = {point[2] for point in read_points(destination / "constant/polyMesh/points")}
        if zs != {0.0, THICKNESS}:
            raise ValueError(f"{destination}: configured extrusion thickness is not {THICKNESS}")
        results.append({"name": name, "source": str(source), "case": str(destination),
                        "source_input_sha256": INPUT_SHA256, "thickness": THICKNESS,
                        "nu": NU, "U_in": U_IN, "p_out": P_OUT, **mesh_result})
    report = {"valid": True, "mode": "configure", "input": str(INPUT),
              "input_sha256": INPUT_SHA256, "thickness": THICKNESS, "cases": results,
              "notes": ["configuration only; no OpenFOAM solver was run",
                        "patch segmentation is a fixed verification fixture, not product UI support"]}
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


def scalar_field(path: Path):
    text = payload(path)
    uniform = re.search(r"internalField\s+uniform\s+([^;\s]+)", text)
    if uniform:
        value = float(uniform.group(1))
        if not math.isfinite(value): raise ValueError(f"{path}: non-finite scalar")
        return [value]
    match = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not match:
        raise ValueError(f"{path}: unsupported scalar field")
    values = [float(v) for v in match.group(2).split()]
    if len(values) != int(match.group(1)) or any(not math.isfinite(v) for v in values):
        raise ValueError(f"{path}: scalar count or finiteness mismatch")
    return values


def vector_field(path: Path):
    text = payload(path)
    uniform = re.search(r"internalField\s+uniform\s*\(([^()]*)\)", text)
    if uniform:
        row = tuple(float(v) for v in uniform.group(1).split())
        if len(row) != 3 or any(not math.isfinite(v) for v in row): raise ValueError(f"{path}: invalid internal vector")
        return [row]
    match = re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not match: raise ValueError(f"{path}: unsupported vector field")
    rows = [tuple(float(v) for v in row.split()) for row in re.findall(r"\(([^()]*)\)", match.group(2))]
    if len(rows) != int(match.group(1)) or any(len(row) != 3 or any(not math.isfinite(v) for v in row) for row in rows):
        raise ValueError(f"{path}: vector count or finiteness mismatch")
    return rows


def boundary_type(path: Path, name: str):
    text = payload(path); at = text.find("boundaryField")
    match = re.search(rf"\b{re.escape(name)}\s*\{{(.*?)\n\s*\}}", text[at:], re.S)
    if not match: raise ValueError(f"{path}: missing boundary {name}")
    kind = re.search(r"\btype\s+([^;\s]+)", match.group(1))
    return kind.group(1) if kind else None


def boundary_values(path: Path, vector: bool = False):
    """Read solved boundary values, preserving uniform/nonuniform cardinality."""
    text = payload(path); at = text.find("boundaryField")
    if at < 0:
        raise ValueError(f"{path}: missing boundaryField")
    values = {}
    for name, block in re.findall(r"\n\s*([A-Za-z_][A-Za-z0-9_]*)\s*\{(.*?)\n\s*\}", text[at:], re.S):
        if vector:
            uni = re.search(r"\bvalue\s+uniform\s*\(([^()]*)\)", block)
            if uni:
                row = tuple(float(v) for v in uni.group(1).split())
                if len(row) != 3 or any(not math.isfinite(v) for v in row): raise ValueError(f"{path}: invalid vector value for {name}")
                values[name] = [row]
                continue
            non = re.search(r"\bvalue\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", block, re.S)
            if non:
                rows = [tuple(float(v) for v in row.split()) for row in re.findall(r"\(([^()]*)\)", non.group(2))]
                if len(rows) != int(non.group(1)) or any(len(row) != 3 or any(not math.isfinite(v) for v in row) for row in rows):
                    raise ValueError(f"{path}: vector value count mismatch for {name}")
                values[name] = rows
        else:
            uni = re.search(r"\bvalue\s+uniform\s+([^;\s]+)", block)
            if uni:
                value = float(uni.group(1))
                if not math.isfinite(value): raise ValueError(f"{path}: non-finite scalar value for {name}")
                values[name] = [value]
                continue
            non = re.search(r"\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", block, re.S)
            if non:
                rows = [float(v) for v in non.group(2).split()]
                if len(rows) != int(non.group(1)) or any(not math.isfinite(v) for v in rows): raise ValueError(f"{path}: scalar value count mismatch for {name}")
                values[name] = rows
    return values


def face_area(face, points):
    origin = points[face[0]]; area = [0.0, 0.0, 0.0]
    for i in range(1, len(face) - 1):
        a = points[face[i]]; b = points[face[i + 1]]
        u = [a[j] - origin[j] for j in range(3)]; v = [b[j] - origin[j] for j in range(3)]
        cross = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
        for j in range(3): area[j] += 0.5 * cross[j]
    return area


def latest_time(case: Path) -> Path:
    times = [p for p in case.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?", p.name)]
    times = [p for p in times if float(p.name) > 0.0]
    if not times: raise ValueError(f"{case}: no positive numeric solved result time")
    return max(times, key=lambda p: float(p.name))


def evaluate_case(case: Path) -> dict:
    mesh = case / "constant/polyMesh"; points = read_points(mesh / "points")
    faces = read_faces(mesh / "faces"); owners = read_labels(mesh / "owner"); patches = read_boundary(mesh / "boundary")
    time = latest_time(case); phi = boundary_values(time / "phi"); p = boundary_values(time / "p"); u = boundary_values(time / "U", vector=True)
    p_internal, u_internal = scalar_field(time / "p"), vector_field(time / "U")
    neighbours = read_labels(mesh / "neighbour")
    cell_count = 1 + max(owners + neighbours, default=-1)
    for field in ("p", "U"):
        declared = re.search(r"internalField\s+nonuniform\s+List<\w+>\s+(\d+)", payload(time / field))
        if declared and int(declared.group(1)) != cell_count:
            raise ValueError(f"{time / field}: nonuniform internal count differs from mesh")
    if len(p_internal) not in (1, cell_count) or len(u_internal) not in (1, cell_count):
        raise ValueError(f"{time}: internal p/U field count does not match cell count {cell_count}")
    patch_data = {}
    for patch in patches:
        name = patch["name"]; ids = range(patch["startFace"], patch["startFace"] + patch["nFaces"])
        if patch["type"] == "empty":
            continue
        areas = [face_area(faces[i], points) for i in ids]
        mags = [math.sqrt(sum(v * v for v in area)) for area in areas]
        flux = phi.get(name)
        if flux is None or len(flux) not in (1, patch["nFaces"]):
            raise ValueError(f"{time/'phi'}: missing or mismatched phi values for {name}")
        patch_flux = sum((flux[0] if len(flux) == 1 else flux[j]) for j in range(patch["nFaces"]))
        pressure = p.get(name)
        if pressure is None and boundary_type(time / "p", name) == "zeroGradient":
            pressure = [p_internal[owners[face_id]] if len(p_internal) > 1 else p_internal[0] for face_id in ids]
        if pressure is None or len(pressure) not in (1, patch["nFaces"]):
            raise ValueError(f"{time/'p'}: missing or mismatched p values for {name}")
        area_sum = sum(mags)
        if not math.isfinite(area_sum) or area_sum <= 0.0: raise ValueError(f"{mesh/'faces'}: non-positive area for {name}")
        pavg = sum(mags[j] * (pressure[0] if len(pressure) == 1 else pressure[j]) for j in range(patch["nFaces"])) / area_sum
        patch_data[name] = {"face_count": patch["nFaces"], "area": area_sum, "flux": patch_flux, "area_average_p": pavg}
    if "inlet" not in patch_data or "outlet" not in patch_data:
        raise ValueError("configured case lacks inlet/outlet patches")
    qin, qout = patch_data["inlet"]["flux"], patch_data["outlet"]["flux"]
    balance = abs(qin + qout) / max(abs(qin), 1.0e-30)
    if not qin < 0.0 or not qout > 0.0: raise ValueError("inlet/outlet phi signs do not match the fixed flow direction")
    outlet = next(patch for patch in patches if patch["name"] == "outlet")
    outlet_ids = list(range(outlet["startFace"], outlet["startFace"] + outlet["nFaces"]))
    outlet_u = u.get("outlet")
    if outlet_u is None and boundary_type(time / "U", "outlet") == "zeroGradient":
        outlet_u = [u_internal[owners[face_id]] if len(u_internal) > 1 else u_internal[0] for face_id in outlet_ids]
    if outlet_u is None or len(outlet_u) not in (1, outlet["nFaces"]):
        raise ValueError(f"{time/'U'}: missing or mismatched outlet velocity values")
    profile = []
    weighted_ux = 0.0; outlet_area = 0.0
    for j, face_id in enumerate(outlet_ids):
        area = face_area(faces[face_id], points); amag = math.sqrt(sum(v * v for v in area))
        if not amag: raise ValueError(f"{mesh/'faces'}: zero-area outlet face {face_id}")
        velocity = outlet_u[0] if len(outlet_u) == 1 else outlet_u[j]
        centroid_y = sum(points[index][1] for index in faces[face_id]) / len(faces[face_id])
        weighted_ux += amag * velocity[0]; outlet_area += amag
        profile.append({"y": centroid_y, "normalized_ux": velocity[0], "area": amag})
    mean_ux = weighted_ux / outlet_area
    if abs(mean_ux) <= 1.0e-30: raise ValueError("outlet area-average streamwise velocity is zero")
    for sample in profile: sample["normalized_ux"] /= mean_ux
    profile.sort(key=lambda sample: sample["y"])
    return {"valid": True, "case": str(case), "time": time.name, "patches": patch_data,
            "flow": {"inlet_flux": qin, "outlet_flux": qout, "relative_flux_imbalance": balance,
            "pressure_drop": patch_data["inlet"]["area_average_p"] - patch_data["outlet"]["area_average_p"],
                     "expected_inlet_flux": -U_IN[0] * patch_data["inlet"]["area"],
                     "inlet_flux_relative_to_fixed_U": abs(qin + U_IN[0] * patch_data["inlet"]["area"]) / abs(U_IN[0] * patch_data["inlet"]["area"]),
                     "outlet_profile": {"area_average_ux": mean_ux, "samples": profile}},
            "notes": ["flux and pressure are read from solved fields; no convergence claim is inferred"]}


def compare_profiles(profiles):
    """Use one common interval and one quadrature for all grid comparisons."""
    ordered = []
    for rows in profiles:
        rows = sorted(rows, key=lambda row: row["y"])
        if len(rows) < 2 or any(not math.isfinite(row[key]) for row in rows
                               for key in ("y", "normalized_ux")):
            raise ValueError("profile must contain at least two finite samples")
        if any(a["y"] >= b["y"] for a, b in zip(rows, rows[1:])):
            raise ValueError("profile ordinates must be unique")
        ordered.append(rows)
    lo, hi = max(rows[0]["y"] for rows in ordered), min(rows[-1]["y"] for rows in ordered)
    if not hi > lo:
        raise ValueError("outlet profiles have no common interval")
    grid = [lo + (hi - lo) * i / 63.0 for i in range(64)]
    grid[0], grid[-1] = lo, hi
    def interp(rows, y):
        for a, b in zip(rows, rows[1:]):
            if a["y"] <= y <= b["y"]:
                t = (y - a["y"]) / (b["y"] - a["y"])
                return a["normalized_ux"] + t * (b["normalized_ux"] - a["normalized_ux"])
        raise ValueError("profile interpolation would extrapolate")
    values = [[interp(rows, y) for y in grid] for rows in ordered]
    # Trapezoidal integral divided by common interval length (uniform samples).
    l2 = []
    for left, right in zip(values, values[1:]):
        squares = [(a - b) ** 2 for a, b in zip(left, right)]
        l2.append(math.sqrt((math.fsum(squares) - .5 * (squares[0] + squares[-1])) / 63.0))
    return {"common_y_interval": [lo, hi], "sample_y": grid,
            "normalized_profiles": values, "adjacent_l2": l2,
            "normalization": "sqrt(integral((Ua/meanUa-Ub/meanUb)^2)dy / interval_length)"}


def solver_convergence(log):
    blocks = re.split(r"^Time = ([0-9.eE+-]+)\s*$", log, flags=re.M)
    residuals = {}
    if len(blocks) >= 3:
        for field, initial in re.findall(r"Solving for (Ux|Uy|Uz|p), Initial residual = ([0-9.eE+-]+)", blocks[-1]):
            residuals.setdefault(field, []).append(float(initial))
    reached = ("SIMPLE solution converged" in blocks[-1] and all(k in residuals for k in ("Ux", "Uy", "p"))
               and all(math.isfinite(value) and value <= 1e-8 for values in residuals.values() for value in values))
    return {"converged_at_requested_residual": reached,
            "last_time": blocks[-2] if len(blocks) >= 3 else None,
            "last_initial_residuals": residuals, "requested_initial_residual": 1e-8}


def evaluate(args) -> int:
    results = [evaluate_case(path) for path in args.case]
    profiles = compare_profiles([item["flow"]["outlet_profile"]["samples"] for item in results])
    pressure = [item["flow"]["pressure_drop"] for item in results]
    pressure_delta = [abs(b-a)/abs(b) if b != 0 else None for a,b in zip(pressure,pressure[1:])]
    logs = getattr(args, "solve_logs", None)
    convergence = []
    if logs is not None:
        if len(logs) != len(results):
            raise ValueError("one solve log is required per evaluated case")
        for item, path in zip(results, logs):
            measured = solver_convergence(Path(path).read_text())
            measured["log"] = str(path)
            if measured["last_time"] is None or float(measured["last_time"]) != float(item["time"]):
                raise ValueError("solve log final time does not match the evaluated fields")
            convergence.append(measured)
    report = {"valid": bool(convergence) and all(c["converged_at_requested_residual"] for c in convergence),
              "fields_read_valid": True, "mode": "evaluate", "cases": results,
              "units": {"phi": "m3/s", "pressure": "m2/s2 (kinematic)", "profile": "Ux/area_average_Ux"},
              "iterative_convergence": convergence if convergence else "not_checked",
              "comparison": {"pressure_drop": pressure, "adjacent_pressure_relative_change": pressure_delta,
                              "relative_flux_imbalance": [item["flow"]["relative_flux_imbalance"] for item in results],
                              "outlet_profile": profiles},
              "notes": ["A nonconverged grid cannot establish spatial convergence.",
                        "Iteration convergence is separate from Q1, expanded checkMesh and general physical accuracy."]}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k not in ("cases", "comparison")}, indent=2))
    return 0 if report["valid"] else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)
    configure_parser = sub.add_parser("configure")
    configure_parser.add_argument("--output-root", type=Path, default=REPO / "outputs/engineering-cfd/nozzle-flow")
    configure_parser.add_argument("--report", type=Path, default=REPO / "outputs/engineering-cfd/nozzle-flow/configure.json")
    configure_parser.set_defaults(case=[("r01", REPO / "outputs/engineering-scale/nozzle/r01/openfoam"),
                                        ("r02", REPO / "outputs/engineering-scale/nozzle/r02/openfoam"),
                                        ("r03", REPO / "outputs/engineering-scale/nozzle/r03/openfoam")])
    evaluate_parser = sub.add_parser("evaluate")
    evaluate_parser.add_argument("--root", type=Path, default=REPO / "outputs/engineering-cfd/nozzle-flow")
    evaluate_parser.add_argument("--report", type=Path, default=REPO / "outputs/engineering-cfd/nozzle-flow/evaluate.json")
    evaluate_parser.add_argument("--solve-logs", nargs="+", type=Path, help="one final solve log per grid, in r01/r02/r03 order")
    evaluate_parser.set_defaults(case=None)
    args = parser.parse_args()
    if args.mode == "evaluate":
        args.case = [args.root / name for name in ("r01", "r02", "r03")]
    try:
        return configure(args) if args.mode == "configure" else evaluate(args)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        result = {"valid": False, "mode": args.mode, "issues": [str(exc)]}
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
