#!/usr/bin/env python3
"""Independent wall-traction diagnostics for already verified native-flow runs.

This is deliberately diagnostic-only: source acceptance and refinement results
are copied into the output and are never replaced by the traction measurements.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import verify_native_flow as vf


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def gauss_edge(a: tuple[float, float], b: tuple[float, float]):
    nodes, weights = np.polynomial.legendre.leggauss(16)
    length = math.hypot(b[0] - a[0], b[1] - a[1])
    for node, weight in zip(nodes, weights):
        t = (float(node) + 1.0) * 0.5
        yield ((a[0] * (1 - t) + b[0] * t,
                a[1] * (1 - t) + b[1] * t), float(weight) * length * 0.5)


def exact_channel(x: float, y: float, bounds: tuple[float, float, float, float], nu: float, speed: float):
    xmin, ymin, xmax, ymax = bounds
    h = ymax - ymin
    pressure = 8.0 * nu * speed * (xmax - x) / (h * h)
    du_dy = 4.0 * speed * (1.0 - 2.0 * (y - ymin) / h) / h
    return pressure, ((0.0, du_dy), (0.0, 0.0))


def exact_mms(x: float, y: float, speed: float, nu: float, slope: float, gauge: float):
    pi = math.pi
    pressure = speed * speed * (math.cos(pi * x) * math.cos(pi * y) + slope * (x + y)) - gauge
    ux = pi * speed * math.sin(2 * pi * x) * math.sin(2 * pi * y)
    uy = 2 * pi * speed * math.sin(pi * x) ** 2 * math.cos(2 * pi * y)
    vx = -2 * pi * speed * math.cos(2 * pi * x) * math.sin(pi * y) ** 2
    vy = -ux
    return pressure, ((ux, uy), (vx, vy))


def traction(pressure: float, grad: tuple[tuple[float, float], tuple[float, float]],
             area_vector: tuple[float, float], nu: float):
    sx, sy = area_vector
    viscous = (-nu * ((grad[0][0] + grad[0][0]) * sx + (grad[0][1] + grad[1][0]) * sy),
               -nu * ((grad[1][0] + grad[0][1]) * sx + (grad[1][1] + grad[1][1]) * sy))
    return (pressure * sx, pressure * sy), viscous


def rms(rows: list[dict[str, Any]], actual: str, exact: str) -> float:
    total = sum(row["length"] for row in rows)
    return math.sqrt(sum(row["length"] * sum((a - b) ** 2 for a, b in zip(row[actual], row[exact]))
                      for row in rows) / max(total, 1e-300))


def audit_case(v: dict[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {"label": v.get("label"), "case": v.get("case"), "valid": False,
                              "diagnosticOnly": True, "issues": []}
    declared = dict(v.get("artifactSha256", {}))
    mesh_path = Path(v["mesh"])
    checks: dict[str, Any] = {"declaredArtifactsUnchanged": True, "hashes": {}}
    for raw, expected in [(str(mesh_path), v.get("meshSha256")), *declared.items()]:
        path = Path(raw)
        if not path.is_file():
            checks["declaredArtifactsUnchanged"] = False
            result["issues"].append(f"missing declared artifact: {path}")
            continue
        actual = sha256(path)
        checks["hashes"][str(path)] = {"declared": expected, "actual": actual, "match": actual == expected}
        if expected and actual != expected:
            checks["declaredArtifactsUnchanged"] = False
            result["issues"].append(f"artifact hash changed: {path}")
    result["artifactIntegrity"] = checks
    if not checks["declaredArtifactsUnchanged"] or not v.get("valid"):
        result["issues"].append("source verification is not valid; traction audit skipped")
        return result
    try:
        mesh = vf.read_cm2d(mesh_path)
        measured = vf.measure(mesh, 1e-11, 1e-9)
        prefix = Path(v["prefix"])
        cells = vf.read_cells(Path(str(prefix) + ".cells.csv"), mesh, measured, v["case"])
        faces, has_momentum = vf.read_faces(Path(str(prefix) + ".faces.csv"), mesh)
        if not has_momentum:
            raise vf.VerificationError("face momentum columns are unavailable")
        geometries = vf.face_geometry(mesh, measured)
        boundaries = vf.flow_boundaries(mesh, measured, v["case"], float(v["speed"]))
        slope = float(v.get("manufacturedPressureSlope", 0.0))
        gauge = 0.0
        if v["case"] == "manufactured":
            gauge = vf.manufactured_sample(cells[0]["x"], cells[0]["y"], float(v["speed"]), float(v["nu"]), slope)["p"]
        rows = []
        total_actual = {"pressure": [0.0, 0.0], "viscous": [0.0, 0.0], "total": [0.0, 0.0]}
        total_exact = {"pressure": [0.0, 0.0], "viscous": [0.0, 0.0], "total": [0.0, 0.0]}
        for edge, geom, face in zip(mesh.edges, geometries, faces):
            role = boundaries["roles"][edge.id]
            if edge.neighbour >= 0 or role not in ("wall", "lid"):
                continue
            a, b = mesh.vertices[edge.v0], mesh.vertices[edge.v1]
            exact_p_force = [0.0, 0.0]; exact_v_force = [0.0, 0.0]
            length = math.hypot(*geom.area_vector)
            for point, weight in gauss_edge(a, b):
                if v["case"] == "channel":
                    p, grad = exact_channel(point[0], point[1], measured.bounds, float(v["nu"]), float(v["speed"]))
                else:
                    p, grad = exact_mms(point[0], point[1], float(v["speed"]), float(v["nu"]), slope, gauge)
                ep, ev = traction(p, grad, geom.area_vector, float(v["nu"]))
                exact_p_force[0] += ep[0] * weight / length
                exact_p_force[1] += ep[1] * weight / length
                exact_v_force[0] += ev[0] * weight / length
                exact_v_force[1] += ev[1] * weight / length
            exact_p = (exact_p_force[0] / length, exact_p_force[1] / length)
            exact_v = (exact_v_force[0] / length, exact_v_force[1] / length)
            actual_p = (face["pressure"] * geom.area_vector[0] / length,
                        face["pressure"] * geom.area_vector[1] / length)
            actual_v = (face["diffusionX"] / length, face["diffusionY"] / length)
            actual_t = (actual_p[0] + actual_v[0], actual_p[1] + actual_v[1])
            exact_t = (exact_p[0] + exact_v[0], exact_p[1] + exact_v[1])
            rows.append({"face": edge.id, "length": length, "pressureActual": actual_p, "pressureExact": tuple(exact_p),
                         "viscousActual": actual_v, "viscousExact": tuple(exact_v),
                         "totalActual": actual_t, "totalExact": exact_t})
            for kind, actual, exact in (("pressure", actual_p, exact_p), ("viscous", actual_v, exact_v), ("total", actual_t, exact_t)):
                for k in (0, 1):
                    total_actual[kind][k] += actual[k] * length
                    total_exact[kind][k] += exact[k] * length
        if not rows:
            raise vf.VerificationError("no wall/lid boundary edges")
        exact_totals = {kind: tuple(total_exact[kind]) for kind in total_exact}
        analytic_checks: dict[str, Any] = {}
        if v["case"] == "channel":
            expected = (8.0 * float(v["nu"]) * float(v["speed"]) *
                        (measured.bounds[2] - measured.bounds[0]) /
                        (measured.bounds[3] - measured.bounds[1]), 0.0)
            analytic_checks["expectedViscousTotal"] = expected
            analytic_checks["viscousTotalMatches"] = all(
                math.isclose(exact_totals["viscous"][i], expected[i], abs_tol=1e-12)
                for i in (0, 1))
        elif abs(slope - 1.0) <= 1e-14:
            analytic_checks["expectedPressureTotal"] = (1.0, 1.0)
            analytic_checks["expectedViscousTotal"] = (0.0, 0.0)
            analytic_checks["pressureTotalMatches"] = all(
                math.isclose(exact_totals["pressure"][i], (1.0, 1.0)[i], abs_tol=1e-12)
                for i in (0, 1))
            analytic_checks["viscousTotalMatches"] = all(
                math.isclose(exact_totals["viscous"][i], 0.0, abs_tol=1e-12)
                for i in (0, 1))
        if any(value is False for value in analytic_checks.values()):
            result["issues"].append("analytic traction self-check failed")
        result.update({"valid": not result["issues"], "faceCount": len(rows), "wallFaces": rows,
                       "forceTotals": {kind: {"actual": total_actual[kind], "exact": total_exact[kind],
                                                "error": [total_actual[kind][i] - total_exact[kind][i] for i in (0, 1)]}
                                        for kind in total_actual},
                       "tractionRmsPerFace": {
                           "pressure": rms(rows, "pressureActual", "pressureExact"),
                           "viscous": rms(rows, "viscousActual", "viscousExact"),
                           "total": rms(rows, "totalActual", "totalExact")},
                       "quadrature": "16-point Gauss-Legendre per straight boundary edge",
                       "definition": "fluid-on-wall pressure and -nu*(gradU+gradU^T).S",
                       "analyticSelfChecks": analytic_checks})
        if v["case"] == "channel":
            expected = 8.0 * float(v["nu"]) * float(v["speed"]) * (measured.bounds[2] - measured.bounds[0]) / (measured.bounds[3] - measured.bounds[1])
            result["channelViscousExpectedForceX"] = expected
            result["channelViscousExpectedPerFaceTractionX"] = (
                4.0 * float(v["nu"]) * float(v["speed"]) /
                (measured.bounds[3] - measured.bounds[1]))
    except (OSError, KeyError, ValueError, vf.VerificationError) as exc:
        result["issues"].append(str(exc))
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", action="append", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output: dict[str, Any] = {"format": "cartmesh2d-wall-traction-diagnostic-v1", "diagnosticOnly": True,
                              "scriptSha256": sha256(Path(__file__).resolve()),
                              "sources": [], "cases": [], "valid": True, "issues": []}
    for raw in args.summary:
        path = Path(raw).resolve()
        source = json.loads(path.read_text(encoding="utf-8"))
        output["sources"].append({"path": str(path), "sha256": sha256(path), "valid": source.get("valid"), "status": source.get("status"),
                                  "sequenceChecks": source.get("sequenceChecks"), "issues": source.get("issues", [])})
        for item in source.get("cases", []):
            verification = dict(item.get("verification", item))
            verification["label"] = item.get("label", verification.get("label"))
            verification.setdefault("scheme", verification.get("native", {}).get("convection"))
            for key in ("meshKind", "scheme", "level"):
                if key in item:
                    verification[key] = item[key]
            verification["inputSummarySha256"] = sha256(path)
            if verification.get("case") in ("channel", "manufactured"):
                audit = audit_case(verification)
                audit["sourceSummary"] = str(path)
                audit["inputSummarySha256"] = verification["inputSummarySha256"]
                audit["counts"] = verification.get("counts")
                audit["mesh"] = verification.get("mesh")
                audit["characteristicH"] = verification.get("meshMeasurement", {}).get("characteristicH")
                for key in ("label", "meshKind", "scheme", "level"):
                    if key in verification:
                        audit[key] = verification[key]
                output["cases"].append(audit)
    output["valid"] = all(case.get("valid") for case in output["cases"]) if output["cases"] else False
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2, sort_keys=True, allow_nan=False))
    return 0 if output["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
