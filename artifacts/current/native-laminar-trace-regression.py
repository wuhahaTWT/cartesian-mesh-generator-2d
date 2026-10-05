#!/usr/bin/env python3
"""Regress native prescribed-velocity trace reporting on a curved duct and cylinder.

This script only orchestrates the product mesh/flow CLIs and reads their native
summaries and output hashes.  It does not implement a flow equation or change
any mesh, solver, convergence, or quality gate.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import time


ROOT = Path("outputs/cloud-laminar/trace-regression")
DUCT_SOURCE = Path("examples/complex/nozzle_profile.xy")
CYLINDER_MESH = Path("outputs/cloud-laminar/cylinder-joint/level-1.solver.cm2d")
CYLINDER_BASELINE = Path("outputs/cloud-laminar/cylinder-joint/level-1-face-limited-linear")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command: list[str], log: Path, timeout: float) -> dict:
    started = time.monotonic()
    with log.open("w") as stream:
        try:
            result = subprocess.run(command, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=timeout)
            status = "completed"
            returncode = result.returncode
        except subprocess.TimeoutExpired:
            status = "timed-out-no-qualification"
            returncode = None
    return {"command": command, "log": str(log), "status": status,
            "returncode": returncode, "processSeconds": time.monotonic()-started}


def summary(prefix: Path) -> dict:
    data = json.loads(Path(str(prefix)+".json").read_text())
    area = pressure_integral = 0.0
    maximum_speed = 0.0
    pressure_minimum = float("inf")
    pressure_maximum = float("-inf")
    pressure_maximum_location = None
    with Path(str(prefix)+".cells.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            cell_area = float(row["area"])
            pressure = float(row["p"])
            area += cell_area
            pressure_integral += cell_area*pressure
            maximum_speed = max(maximum_speed, float(row["speed"]))
            pressure_minimum = min(pressure_minimum, pressure)
            if pressure > pressure_maximum:
                pressure_maximum = pressure
                pressure_maximum_location = [float(row["x"]), float(row["y"])]
    performance_path = Path(str(prefix)+".performance.json")
    performance = json.loads(performance_path.read_text()) if performance_path.exists() else {}
    return {
        "status": data["status"], "converged": data["converged"],
        "cells": data["cells"], "evaluations": data["iterations"],
        "steadyAcceleration": data["steadyAcceleration"],
        "convection": data["convection"],
        "continuity": data["continuity"],
        "globalRelativeImbalance": data["globalRelativeImbalance"],
        "momentumResidual": data["momentumResidual"],
        "wallForceX": data["wallForceX"],
        "wallForceY": data["wallForceY"],
        "outletBackflowFaces": data["outletBackflowFaces"],
        "wallTraceWallFaces": data["wallTraceWallFaces"],
        "wallTraceAdjacentVertices": data["wallTraceAdjacentVertices"],
        "wallTraceDiscontinuousVertices": data["wallTraceDiscontinuousVertices"],
        "wallTraceMaximumVelocityJump_m_s": data["wallTraceMaximumVelocityJump"],
        "wallTraceVelocityTolerance_m_s": data["wallTraceVelocityTolerance"],
        "velocityTracePrescribedFaces": data["velocityTracePrescribedFaces"],
        "velocityTraceAdjacentVertices": data["velocityTraceAdjacentVertices"],
        "velocityTraceDiscontinuousVertices": data["velocityTraceDiscontinuousVertices"],
        "velocityTraceMaximumVelocityJump_m_s": data["velocityTraceMaximumVelocityJump"],
        "velocityTraceMaximumJumpLocation_m": data["velocityTraceMaximumJumpLocation"],
        "velocityTraceVelocityTolerance_m_s": data["velocityTraceVelocityTolerance"],
        "maximumCellSpeed_m_s": maximum_speed,
        "pressureMinimum_m2_s2": pressure_minimum,
        "pressureMaximum_m2_s2": pressure_maximum,
        "pressureMaximumLocation_m": pressure_maximum_location,
        "areaWeightedPressureMean_m2_s2": pressure_integral/area,
        "nativeSolveSeconds": performance.get("solveSeconds"),
        "summarySha256": sha256(Path(str(prefix)+".json")),
        "cellsSha256": sha256(Path(str(prefix)+".cells.csv")),
        "facesSha256": sha256(Path(str(prefix)+".faces.csv")),
        "fieldsSha256": sha256(Path(str(prefix)+".fields.json")),
        "residualsSha256": sha256(Path(str(prefix)+".residuals.csv")),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--levels", nargs="+", type=int, default=[5, 6, 7])
    parser.add_argument("--timeout", type=float, default=900)
    args = parser.parse_args()
    ROOT.mkdir(parents=True, exist_ok=True)
    record = {
        "scope": "Native prescribed-velocity trace regression; no independent flow equation",
        "ductGeometry": str(DUCT_SOURCE),
        "ductGeometrySha256": sha256(DUCT_SOURCE),
        "cylinder": {}, "ductLevels": [],
    }

    cylinder = ROOT/"cylinder-level1"
    cylinder_run = run([
        "build/cartmesh2d_flow_cli", "--mesh", str(CYLINDER_MESH),
        "--case", "external", "--nu", ".1", "--speed", "1",
        "--convection", "face-limited-linear", "--tolerance", "1e-8",
        "--max-iterations", "2500", "--output", str(cylinder),
    ], Path(str(cylinder)+".log"), args.timeout)
    record["cylinder"] = cylinder_run
    if cylinder_run["returncode"] == 0:
        record["cylinder"]["metrics"] = summary(cylinder)
        record["cylinder"]["byteIdenticalToPriorAccepted"] = {
            suffix: sha256(Path(str(cylinder)+suffix)) == sha256(Path(str(CYLINDER_BASELINE)+suffix))
            for suffix in [".cells.csv", ".faces.csv", ".fields.json", ".residuals.csv"]
        }

    for level in args.levels:
        mesh = ROOT/f"duct-l{level}"
        row = {"level": level}
        row["mesh"] = run([
            "build/cartmesh2d_cli", str(DUCT_SOURCE), str(mesh), str(level),
            "0.03333333333333333", "0.1", "interior", str(mesh)+"-foam",
            str(level), "0",
        ], Path(str(mesh)+"-mesh.log"), args.timeout)
        if row["mesh"]["returncode"] == 0:
            row["meshSha256"] = sha256(Path(str(mesh)+".solver.cm2d"))
            flow = Path(str(mesh)+"-default")
            row["solve"] = run([
                "build/cartmesh2d_flow_cli", "--mesh", str(mesh)+".solver.cm2d",
                "--case", "duct", "--nu", ".1", "--speed", "1",
                "--tolerance", "1e-8", "--max-iterations", "2500",
                "--output", str(flow), "--profile",
            ], Path(str(flow)+".log"), args.timeout)
            if row["solve"]["returncode"] == 0:
                row["metrics"] = summary(flow)
        record["ductLevels"].append(row)
        (ROOT/"runs.json").write_text(json.dumps(record, indent=2)+"\n")
        print(json.dumps(row), flush=True)

    ok = (record["cylinder"].get("returncode") == 0 and
          all(row["mesh"].get("returncode") == 0 and
              row.get("solve", {}).get("returncode") == 0
              for row in record["ductLevels"]))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
