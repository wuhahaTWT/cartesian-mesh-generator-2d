#!/usr/bin/env python3
"""Run and summarize the native guide certificate on the finest joint annulus.

The guide changes only viscosity and its stopping tolerance.  The certified
target and independent Anderson control return to the original equation,
mesh, boundary trace and target tolerance.  No cell is removed from the field
comparison and pressure alignment removes only the incompressible gauge.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs/cloud-laminar"
CLI = ROOT / "build/cartmesh2d_flow_cli"
MESH = OUT / "annulus-joint-7.solver.cm2d"
BOUNDARY = OUT / "annulus-joint-7.boundaries"
PREFIX = OUT / "guide-annulus-joint-7"
EVIDENCE = ROOT / "artifacts/current/native-laminar-guide-crosscase.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command(output: str, nu: str, tolerance: str, *extra: str) -> list[str]:
    return [
        str(CLI), "--mesh", str(MESH), "--case", "custom", "--boundary", str(BOUNDARY),
        "--nu", nu, "--speed", ".5", "--tolerance", tolerance,
        "--max-iterations", "1800", *extra, "--output", str(PREFIX) + output,
    ]


COMMANDS = {
    "direct": command("-direct", ".1", "1e-6"),
    "guide": command("-guide", "1", "1e-3"),
    "certifiedTarget": command(
        "-certified", ".1", "1e-6", "--initial-guess", str(PREFIX) + "-guide.initial.csv"
    ),
    "anderson": command(
        "-anderson", ".1", "1e-6", "--convection", "face-limited-linear",
        "--steady-acceleration", "anderson"
    ),
}


def run(name: str) -> None:
    with (Path(str(PREFIX) + f"-{name}.log")).open("w") as log:
        subprocess.run(COMMANDS[name], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)


def make_initial_guess() -> None:
    source = Path(str(PREFIX) + "-guide.cells.csv")
    target = Path(str(PREFIX) + "-guide.initial.csv")
    with source.open(newline="") as stream, target.open("w", newline="") as output:
        reader = csv.DictReader(stream)
        writer = csv.writer(output)
        writer.writerow(("cell", "x", "y", "u", "v", "p"))
        for row in reader:
            writer.writerow((row["cell"], row["x"], row["y"], row["u"], row["v"], row["p"]))


def field_distance(left_path: Path, right_path: Path, speed: float) -> dict:
    with left_path.open(newline="") as stream:
        left = list(csv.DictReader(stream))
    with right_path.open(newline="") as stream:
        right = list(csv.DictReader(stream))
    if len(left) != len(right):
        raise RuntimeError("field comparison requires identical final cells")
    area = sum(float(row["area"]) for row in left)
    raw_pressure = [float(a["p"]) - float(b["p"]) for a, b in zip(left, right)]
    gauge = sum(float(row["area"]) * value for row, value in zip(left, raw_pressure)) / area
    velocity2 = pressure2 = maximum_velocity = maximum_pressure = 0.0
    for a, b, raw_dp in zip(left, right, raw_pressure):
        weight = float(a["area"])
        du, dv = float(a["u"]) - float(b["u"]), float(a["v"]) - float(b["v"])
        dp = raw_dp - gauge
        velocity2 += weight * (du * du + dv * dv)
        pressure2 += weight * dp * dp
        maximum_velocity = max(maximum_velocity, math.hypot(du, dv))
        maximum_pressure = max(maximum_pressure, abs(dp))
    return {
        "definition": "all-cell area-weighted; pressure difference has only its area mean removed",
        "velocityRmsOverReferenceSpeed": math.sqrt(velocity2 / area) / speed,
        "maximumVelocityDifferenceOverReferenceSpeed": maximum_velocity / speed,
        "pressureRmsOverReferenceSpeedSquared": math.sqrt(pressure2 / area) / speed**2,
        "maximumPressureDifferenceOverReferenceSpeedSquared": maximum_pressure / speed**2,
        "pressureGaugeOffset_m2_s2": gauge,
    }


def run_summary(name: str) -> dict:
    prefix = Path(str(PREFIX) + {"certifiedTarget": "-certified"}.get(name, f"-{name}"))
    raw = json.loads(prefix.with_suffix(".json").read_text())
    keep = (
        "status", "cells", "iterations", "coupledEvaluations", "nu", "tolerance", "convection",
        "steadyAcceleration", "maximumSpeedRatio", "pressureRangeRatio",
        "wallTraceDiscontinuousVertices", "wallTraceMaximumNormalVelocity", "globalRelativeImbalance",
    )
    result = {key: raw.get(key) for key in keep}
    result["evaluations"] = raw.get("coupledEvaluations") or raw["iterations"]
    result["prefix"] = str(prefix.relative_to(ROOT))
    result["sha256"] = {
        suffix: sha256(Path(str(prefix) + suffix))
        for suffix in (".json", ".cells.csv", ".faces.csv", ".residuals.csv")
    }
    return result


def write_evidence() -> None:
    annulus = {name: run_summary(name) for name in COMMANDS}
    direct_cells = Path(str(PREFIX) + "-direct.cells.csv")
    certified_cells = Path(str(PREFIX) + "-certified.cells.csv")
    anderson_cells = Path(str(PREFIX) + "-anderson.cells.csv")
    manufactured_prefix = OUT / "guide-manufactured-64"
    manufactured = json.loads(manufactured_prefix.with_suffix(".json").read_text())
    manufactured_anderson_path = OUT / "manufactured-face-limited-linear-64-anderson-tol1e6.json"
    manufactured_anderson = json.loads(manufactured_anderson_path.read_text())
    manufactured["flatAnderson"] = manufactured_anderson
    manufactured["directVsAnderson"] = field_distance(
        Path(str(manufactured_prefix) + "-direct.cells.csv"),
        OUT / "manufactured-face-limited-linear-64-anderson-tol1e6.cells.csv", 1.0,
    )
    manufactured["rawSha256"] = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (
            manufactured_prefix.with_suffix(".json"),
            Path(str(manufactured_prefix) + "-direct.cells.csv"),
            Path(str(manufactured_prefix) + "-guide.cells.csv"),
            Path(str(manufactured_prefix) + "-certified.cells.csv"),
            manufactured_anderson_path,
            OUT / "manufactured-face-limited-linear-64-anderson-tol1e6.cells.csv",
        )
    }
    evidence = {
        "schema": 1,
        "purpose": "cross-case native branch-certificate control; no physical branch selection",
        "rule": {
            "guideViscosity": "10 * target viscosity",
            "guideTolerance": "sqrt(target tolerance)",
            "certifiedTarget": "original equation and original tolerance",
        },
        "manufactured64": manufactured,
        "movingAnnulusJoint7": {
            "meshSha256": sha256(MESH),
            "boundarySha256": sha256(BOUNDARY),
            "commands": COMMANDS,
            "runs": annulus,
            "directVsCertifiedTarget": field_distance(direct_cells, certified_cells, .5),
            "directVsAnderson": field_distance(direct_cells, anderson_cells, .5),
            "guideCertificateEvaluationCostRatio": (
                annulus["guide"]["evaluations"]
                + annulus["certifiedTarget"]["evaluations"]
            ) / annulus["direct"]["evaluations"],
            "directAndAndersonEvaluationCostRatio": (
                annulus["direct"]["evaluations"]
                + annulus["anderson"]["evaluations"]
            ) / annulus["direct"]["evaluations"],
        },
        "artifacts": {
            "driver": str(Path(__file__).relative_to(ROOT)),
            "driverSha256": sha256(Path(__file__)),
            "nativeGuideSource": "artifacts/current/native-laminar-guide-certificate.cpp",
            "nativeGuideSourceSha256": sha256(ROOT / "artifacts/current/native-laminar-guide-certificate.cpp"),
            "nativeGuideBinarySha256": sha256(OUT / "guide-certificate"),
        },
        "qualification": {
            "productAlgorithmChanged": False,
            "defaultChanged": False,
            "qualityGateChanged": False,
            "relatedLinuxNativeTests": "3/3 passed: flow_boundary, solver_export, io",
            "fullNativeSuiteRerun": False,
            "frontendOrAppRerun": False,
            "conclusion": (
                "Agreement is evidence that independent paths retain the same target branch; "
                "disagreement only exposes branch risk and does not select a physical solution."
            ),
        },
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true", help="rerun all four annulus paths")
    parser.add_argument("--resume", action="store_true", help="reuse direct/guide and run target/Anderson")
    args = parser.parse_args()
    if args.run:
        run("direct")
        run("guide")
    if args.run or args.resume:
        make_initial_guess()
        run("certifiedTarget")
        run("anderson")
    write_evidence()
    print(EVIDENCE)
