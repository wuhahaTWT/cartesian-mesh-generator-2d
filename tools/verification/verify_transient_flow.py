#!/usr/bin/env python3
"""Independent checks for one backward-Euler native transient result.

This verifier treats the solver summary as metadata only.  It rebuilds cell
temporal integrals, face-flux continuity, Taylor--Green errors, and the time
history from the exported artifacts.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_native_flow as native  # noqa: E402


def fail(message: str) -> None:
    raise native.VerificationError(message)


def read_history(path: Path) -> list[dict[str, float]]:
    required = ("step", "time", "dt", "accepted", "innerIterations",
                "momentumResidual", "continuity", "maxCourant",
                "kineticEnergy", "forceX", "forceY")
    rows = native.load_csv(path, required)
    result = []
    last_step = 0
    last_time = -math.inf
    for line, row in enumerate(rows, 2):
        step = native.integer(row["step"], f"{path}:{line} step")
        if step != last_step + 1:
            fail(f"{path}:{line}: step sequence is not contiguous")
        item = {name: native.finite(row[name], f"{path}:{line} {name}")
                for name in required[1:]}
        item["step"] = step
        item["accepted"] = native.integer(row["accepted"], f"{path}:{line} accepted")
        if item["accepted"] not in (0, 1):
            fail(f"{path}:{line}: accepted must be 0 or 1")
        if item["time"] <= 0 or item["time"] <= last_time or item["dt"] <= 0:
            fail(f"{path}:{line}: time/dt is not strictly positive and monotonic")
        if item["accepted"] == 1 and not native.close(item["time"], last_time + item["dt"], 1e-12, 1e-10) and last_time > -math.inf:
            fail(f"{path}:{line}: accepted time does not advance by dt")
        if item["accepted"] == 0 and step != len(rows):
            fail(f"{path}:{line}: failed candidate is not the final history row")
        item["innerIterations"] = native.integer(row["innerIterations"], f"{path}:{line} innerIterations")
        if item["innerIterations"] < 1 or any(item[k] < 0 for k in
                ("momentumResidual", "continuity", "maxCourant", "kineticEnergy")):
            fail(f"{path}:{line}: negative metric or invalid inner iteration count")
        if step == 1 and item["time"] - item["dt"] < -1e-12:
            fail(f"{path}:{line}: inferred initial time is negative")
        last_step, last_time = step, item["time"]
        result.append(item)
    if not result:
        fail(f"{path}: empty time history")
    return result


def exact_tg(x: float, y: float, time: float, speed: float, nu: float) -> tuple[float, float, float]:
    a = speed * math.exp(-2.0 * nu * math.pi * math.pi * time)
    return (a * math.sin(math.pi*x) * math.cos(math.pi*y),
            -a * math.cos(math.pi*x) * math.sin(math.pi*y),
            .25 * a*a * (math.cos(2*math.pi*x) + math.cos(2*math.pi*y)))


def verify(mesh_path: Path, prefix: Path, output: Path) -> dict:
    mesh = native.read_cm2d(mesh_path)
    measured = native.measure(mesh, 1e-11, 1e-9)
    if measured.issues:
        fail("invalid CM2D geometry: " + "; ".join(measured.issues))
    cells_path = Path(str(prefix) + ".cells.csv")
    faces_path = Path(str(prefix) + ".faces.csv")
    summary_path = Path(str(prefix) + ".json")
    history_path = Path(str(prefix) + ".time-history.csv")
    residual_path = Path(str(prefix) + ".residuals.csv")
    for path in (cells_path, faces_path, summary_path, history_path, residual_path):
        if not path.is_file():
            fail(f"missing transient artifact: {path}")
    with summary_path.open(encoding="utf-8-sig") as stream:
        summary = json.load(stream, parse_constant=lambda token: fail(f"non-finite JSON token {token}"))
    if not isinstance(summary, dict) or summary.get("format") != "cartmesh2d-flow-summary-v1":
        fail("invalid native summary format")
    for field in ("nu", "speed", "tolerance"):
        if native.finite(summary.get(field), f"summary {field}") <= 0:
            fail(f"summary {field} must be positive")
    if summary.get("pressureDiscretization") != "shared-face-gauss" or summary.get("pressureBoundaryReconstruction") != "one-sided-linear":
        fail("unsupported transient pressure scheme")
    if summary.get("viscousStress") not in ("symmetric", "laplacian"):
        fail("unsupported viscous stress")
    expected_force = "shared-face-newtonian-traction" if summary["viscousStress"] == "symmetric" else "reconstructed-newtonian-traction"
    if summary.get("forceDefinition") != expected_force:
        fail("force definition disagrees with viscous stress")
    if native.integer(summary.get("cells"), "summary cells") != len(mesh.cells):
        fail("summary cell count differs from mesh")
    if summary.get("temporalDiscretization") != "backward-euler":
        fail("missing or invalid backward-Euler metadata")
    if summary.get("temporalFaceInterpolation") not in (
            "old-flux-defect-skew-corrected-v1",
            "old-and-iteration-flux-defect-skew-corrected-v2"):
        fail("missing or invalid temporal face interpolation metadata")
    relaxation = native.finite(summary.get("velocityRelaxation"), "summary velocityRelaxation")
    if not 0.0 < relaxation <= 1.0:
        fail("transient velocityRelaxation must lie in (0,1]")
    if summary.get("temporalFaceInterpolation") == "old-flux-defect-skew-corrected-v1" and relaxation != 1.0:
        fail("v1 temporal interpolation requires unrelaxed momentum")
    if summary.get("status") != "converged" or summary.get("converged") is not True:
        fail("transient final state is not converged")
    dt = native.finite(summary.get("dt"), "summary dt")
    if dt <= 0:
        fail("summary dt must be positive")
    history = read_history(history_path)
    completed = sum(row["accepted"] for row in history)
    if native.integer(summary.get("completedSteps"), "summary completedSteps") != completed:
        fail("completedSteps differs from accepted history rows")
    if native.integer(summary.get("requestedSteps"), "summary requestedSteps") != len(history):
        fail("requestedSteps differs from history length")
    if completed != len(history):
        fail("converged summary contains an unaccepted time step")
    for row in history:
        if not native.close(row["dt"], dt, 1e-14, 1e-10):
            fail("history dt differs from summary")
        if row["innerIterations"] < 10 or row["momentumResidual"] >= summary["tolerance"] or row["continuity"] >= 1e-8:
            fail("accepted history row fails native stopping conditions")
    residuals = native.read_residuals(residual_path)
    if native.integer(summary.get("iterations"), "summary iterations") != residuals["last"]["iteration"] or summary["iterations"] != history[-1]["innerIterations"]:
        fail("final inner iteration count differs between artifacts")
    for name in ("momentumResidual", "velocityChange", "pressureChange", "continuity"):
        value = native.finite(summary.get(name), f"summary {name}")
        if not native.close(value, residuals["last"][name], 1e-14, 1e-9):
            fail(f"final residual CSV differs from summary {name}")
        if value < 0 or value >= (1e-8 if name == "continuity" else summary["tolerance"]):
            fail(f"final {name} fails native stopping condition")
    accepted = [row for row in history if row["accepted"]]
    if not accepted or accepted[-1]["accepted"] != 1:
        fail("transient result has no accepted final state")
    final_time = accepted[-1]["time"]
    if not native.close(native.finite(summary.get("acceptedTime"), "summary acceptedTime"), final_time, 1e-12, 1e-10):
        fail("acceptedTime differs from history")
    if not native.close(native.finite(summary.get("time"), "summary time"), history[-1]["time"], 1e-12, 1e-10):
        fail("summary time differs from history")
    fields = set(native.csv_fields(cells_path))
    required = {"previousU", "previousV", "temporalX", "temporalY"}
    if fields & required != required:
        fail("transient cell schema is missing previous/temporal columns")
    cells_list = native.read_cells(cells_path, mesh, measured, summary.get("case", "taylor-green"))
    cells = {i: row for i, row in enumerate(cells_list)}
    for i, vals in cells.items():
        if not native.close(vals["temporalX"], vals["area"]*(vals["u"]-vals["previousU"])/dt, 1e-12, 1e-9):
            fail(f"{cells_path}:cell {i}: temporalX is not area*dU/dt")
        if not native.close(vals["temporalY"], vals["area"]*(vals["v"]-vals["previousV"])/dt, 1e-12, 1e-9):
            fail(f"{cells_path}:cell {i}: temporalY is not area*dV/dt")
    if len(cells) != len(mesh.cells):
        fail("transient cell count differs from mesh")
    if set(cells) != set(range(len(mesh.cells))):
        fail("transient cell ids are not contiguous")
    faces, has_momentum = native.read_faces(faces_path, mesh)
    if not has_momentum:
        fail("transient faces.csv lacks complete face momentum columns")
    fluxes = native.face_fluxes(faces)
    case = summary.get("case")
    if case not in ("taylor-green", "cavity", "channel", "external"):
        fail(f"unsupported transient case {case!r}")
    cont = native.continuity(mesh, measured, fluxes, native.finite(summary.get("speed"), "summary speed"), case, 1e-10, 1e-7)
    if not cont["cellValid"] or not cont["globalValid"]:
        fail("independent face-flux continuity failed")
    for name, key in (("continuity", "nativeDefinitionContinuity"), ("globalImbalance", "boundaryFluxSum"),
                      ("globalRelativeImbalance", "globalRelativeImbalance")):
        if not native.close(native.finite(summary.get(name), f"summary {name}"), cont[key], 1e-14, 1e-9):
            fail(f"summary {name} differs from independent continuity")
    expected_reference = "cell 0, kinematic pressure zero" if case in ("cavity", "taylor-green") else "right outlet faces, kinematic pressure zero"
    if summary.get("pressureReference") != expected_reference:
        fail("pressure reference disagrees with boundary conditions")
    # Recompute CFL from the exported final face fluxes and measured cell areas.
    abs_flux = [0.0] * len(mesh.cells)
    for edge, flux in zip(mesh.edges, fluxes):
        abs_flux[edge.owner] += abs(flux)
        if edge.neighbour >= 0:
            abs_flux[edge.neighbour] += abs(flux)
    computed_cfl = max((.5 * dt * q / a for q, a in zip(abs_flux, measured.areas)), default=0.0)
    if not native.close(computed_cfl, native.finite(summary.get("maxCourant"), "summary maxCourant"), 1e-12, 1e-9):
        fail("summary maxCourant differs from independent face-flux CFL")
    if not native.close(computed_cfl, history[-1]["maxCourant"], 1e-12, 1e-9):
        fail("history maxCourant differs from independent face-flux CFL")
    final_native = history[-1]
    for name in ("momentumResidual", "continuity"):
        if not native.close(native.finite(summary.get(name), f"summary {name}"), final_native[name], 1e-12, 1e-9):
            fail(f"summary {name} differs from final history")
    momentum = native.reconstruct_momentum_audit(
        mesh, measured, cells_list, faces, native.finite(summary.get("nu"), "summary nu"),
        native.finite(summary.get("speed"), "summary speed"), case, summary,
        pressure_boundary_reconstruction=summary.get("pressureBoundaryReconstruction", "zero-normal"),
        time_step=dt)
    for deviation in momentum["maxFaceDeviation"].values():
        if deviation > 5e-10:
            fail("independent transient face momentum reconstruction failed")
    if any(d["absolute"] > 5e-10 for d in momentum["summaryDeviation"].values()):
        fail("independent transient momentum/force differs from summary")
    if momentum["cellBoundaryConservationDifference"]["magnitude"] > 5e-10:
        fail("global temporal plus boundary momentum does not match cell residual sum")
    if momentum["cellResidual"]["maxNormalized"] >= native.finite(summary.get("tolerance"), "summary tolerance"):
        fail("independent transient momentum residual is not below tolerance")
    errors = None
    energy = math.fsum(.5*a*(cells[i]["u"]**2+cells[i]["v"]**2) for i, a in enumerate(measured.areas))
    for name, value in (("kineticEnergy", energy), ("forceX", summary["forceX"]), ("forceY", summary["forceY"])):
        if not native.close(history[-1][name], value, 1e-12, 1e-9):
            fail(f"final history {name} differs from final field/summary")
    if case == "taylor-green":
        if any(not native.close(a,b,1e-11,1e-9) for a,b in zip(measured.bounds, (0,0,1,1))) or not native.close(math.fsum(measured.areas),1,1e-11,1e-9):
            fail("Taylor-Green requires a complete unit square")
        if abs(cells[0]["p"]) > 1e-12:
            fail("Taylor-Green pressure gauge is not zero at cell 0")
        nu = native.finite(summary.get("nu"), "summary nu")
        speed = native.finite(summary.get("speed"), "summary speed")
        eu, ev, ep = [], [], []
        for i, cell in enumerate(mesh.cells):
            x, y = measured.centroids[i]
            u, v, p = exact_tg(x, y, final_time, speed, nu)
            eu.append(cells[i]["u"]-u); ev.append(cells[i]["v"]-v)
            # pressure is gauge-fixed at cell 0 by the native solver
            x0, y0 = measured.centroids[0]
            _, _, p0 = exact_tg(x0, y0, final_time, speed, nu)
            ep.append(cells[i]["p"]-(p-p0))
        areas = measured.areas
        errors = {"uL2": native.weighted_l2(eu, areas), "vL2": native.weighted_l2(ev, areas),
                  "pL2": native.weighted_l2(ep, areas)}
        energy = math.fsum(.5*a*(cells[i]["u"]**2+cells[i]["v"]**2) for i, a in enumerate(areas))
        exact_energy = .25*(speed*math.exp(-2*nu*math.pi*math.pi*final_time))**2
        errors["energy"] = abs(energy-exact_energy)
    result = {"valid": True, "scope": "discrete transient balance and artifact consistency; analytic errors are diagnostics, not engineering qualification",
              "sha256": {str(path): native.sha256_file(path) for path in (mesh_path, cells_path, faces_path, summary_path, history_path, residual_path)},
              "controls": {k: summary[k] for k in ("nu", "speed", "tolerance", "convection", "viscousStress", "temporalFaceInterpolation", "velocityRelaxation")}, "case": case, "time": final_time,
              "dt": dt, "history": history, "independentContinuity": cont,
              "temporalIntegral": {"maxAbsX": max(abs(c["temporalX"]) for c in cells.values()),
                                    "maxAbsY": max(abs(c["temporalY"]) for c in cells.values())},
              "analyticErrors": errors, "kineticEnergy": energy,
              "computedMaxCourant": computed_cfl, "momentumAudit": momentum,
              "counts": {"cells": len(mesh.cells), "faces": len(mesh.edges)}}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = verify(args.mesh, args.prefix, args.output)
    except (native.VerificationError, OSError, ValueError, json.JSONDecodeError) as exc:
        failure = {"valid": False, "issues": [str(exc)]}
        native.write_json(args.output, failure)
        print(json.dumps(failure, indent=2))
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
