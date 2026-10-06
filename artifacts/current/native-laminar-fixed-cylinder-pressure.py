#!/usr/bin/env python3
"""Refine space on one fixed 128-segment Re=20 cylinder in the 20D domain.

The driver only orchestrates the native mesh/flow CLIs and native wall-pressure
postprocessor.  It keeps the input polygon byte-identical and does not alter a
flow equation, boundary condition, mesh acceptance threshold, or result field.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time


ROOT = Path("outputs/cloud-laminar/cylinder-joint")
POLYGON = ROOT / "level-2.xy"
REFERENCE_CD = 2.045
EVIDENCE = Path("artifacts/current/native-laminar-fixed-cylinder-pressure.json")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def run(command: list[str], log: Path, timeout: float) -> dict:
    started = time.monotonic()
    with log.open("w") as stream:
        try:
            result = subprocess.run(command, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=timeout)
            returncode = result.returncode
            status = "completed"
        except subprocess.TimeoutExpired:
            returncode = None
            status = "timed-out-no-qualification"
    return {
        "command": command,
        "log": str(log),
        "status": status,
        "returncode": returncode,
        "processSeconds": time.monotonic() - started,
    }


def selected_summary(path: Path) -> dict:
    summary = json.loads(path.read_text())
    keys = [
        "status", "converged", "cells", "iterations", "coupledEvaluations",
        "coupledFailedEvaluations", "continuity", "globalImbalance",
        "globalRelativeImbalance", "momentumResidual", "velocityChange",
        "pressureChange", "wallForceX", "pressureForceX",
        "wallViscousForceX", "steadyAcceleration", "convection",
        "velocityTracePrescribedFaces", "velocityTraceAdjacentVertices",
        "velocityTraceDiscontinuousVertices",
        "velocityTraceMagnitudeDiscontinuousVertices",
        "velocityTraceEqualMagnitudeDirectionVertices",
        "velocityTraceMaximumVelocityJump", "wallTraceMaximumNormalVelocity",
    ]
    result = {key: summary[key] for key in keys if key in summary}
    result["referenceCd"] = REFERENCE_CD
    result["referenceRelativeDifference"] = (
        summary["wallForceX"] - REFERENCE_CD) / REFERENCE_CD
    return result


def field_extrema(prefix: Path) -> dict:
    with Path(str(prefix) + ".cells.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    values = {name: [float(row[name]) for row in rows]
              for name in ("u", "v", "p", "speed")}
    return {
        "minimumU_m_s": min(values["u"]),
        "maximumU_m_s": max(values["u"]),
        "minimumV_m_s": min(values["v"]),
        "maximumV_m_s": max(values["v"]),
        "minimumKinematicPressure_m2_s2": min(values["p"]),
        "maximumKinematicPressure_m2_s2": max(values["p"]),
        "maximumSpeed_m_s": max(values["speed"]),
    }


def summarized_case(name: str, mesh: Path, prefix: Path) -> dict:
    summary_path = Path(str(prefix) + ".json")
    summary = json.loads(summary_path.read_text())
    pressure = summary["pressureForceX"]
    viscous = summary["wallViscousForceX"]
    total = summary["wallForceX"]
    return {
        "name": name,
        "mesh": str(mesh),
        "meshSha256": sha256(mesh),
        "outputPrefix": str(prefix),
        "summary": selected_summary(summary_path),
        "fieldExtrema": field_extrema(prefix),
        "streamwiseWallLoadCancellationFactor":
            (abs(pressure) + abs(viscous)) / max(abs(total), 1e-300),
        "rawSha256": {
            suffix: sha256(Path(str(prefix) + suffix))
            for suffix in (".json", ".cells.csv", ".faces.csv", ".residuals.csv")
        },
    }


def analyze(profile_binary: Path) -> None:
    root = ROOT
    cases = [
        ("fixed128-grid0-fll", root / "far-20-fixed128-grid0.solver.cm2d",
         root / "far-20-fixed128-grid0-face-limited-linear"),
        ("fixed128-grid1-fll", root / "far-20-fixed128-grid1.solver.cm2d",
         root / "far-20-fixed128-grid1-face-limited-linear"),
        ("fixed128-grid1-upwind", root / "far-20-fixed128-grid1.solver.cm2d",
         root / "far-20-fixed128-grid1-upwind"),
        ("fixed128-grid1-volume075-fll",
         root / "far-20-fixed128-grid1-volume-075.solver.cm2d",
         root / "far-20-fixed128-grid1-volume-075-face-limited-linear"),
        ("fixed128-grid2-fll", root / "far-20.solver.cm2d",
         root / "far-20-face-limited-linear"),
    ]
    for _, mesh, prefix in cases:
        if not mesh.exists() or not Path(str(prefix) + ".json").exists():
            raise FileNotFoundError(f"missing accepted native case: {mesh} / {prefix}")
    profile_path = root / "cylinder-pressure-profile-fixed128-20d.json"
    profile_command = [str(profile_binary)]
    for label, mesh, prefix in (cases[0], cases[1], cases[4]):
        profile_command.extend([label, str(mesh), str(prefix) + ".faces.csv"])
    with profile_path.open("w") as output:
        subprocess.run(profile_command, stdout=output, check=True)
    repair_profile_path = root / "cylinder-pressure-profile-fixed128-grid1-repair.json"
    repair_command = [str(profile_binary)]
    for label, mesh, prefix in (cases[1], cases[3]):
        repair_command.extend([label, str(mesh), str(prefix) + ".faces.csv"])
    with repair_profile_path.open("w") as output:
        subprocess.run(repair_command, stdout=output, check=True)
    run_records = json.loads((root / "far-20-fixed128-runs.json").read_text())
    repair_record = json.loads(
        (root / "far-20-fixed128-grid1-volume-075-repair.json").read_text())
    evidence = {
        "format": "cartmesh2d-native-laminar-fixed-cylinder-pressure-v1",
        "status": "diagnostic-counterexample-not-default",
        "scope": (
            "One byte-identical 128-segment stationary-cylinder polygon and one 20D "
            "outer domain, with spatial size only refined. Native accepted fields and "
            "the existing exact-union/split repair are compared; no surrogate flow "
            "equation or relaxed acceptance threshold is used."
        ),
        "environment": {
            "system": f"{platform.system()} {platform.machine()}",
            "localHeadBeforeEvidenceCommit": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True).strip(),
            "flowCliSha256": sha256(Path("build/cartmesh2d_flow_cli")),
            "meshCliSha256": sha256(Path("build/cartmesh2d_cli")),
            "profileSourceSha256": sha256(Path(
                "artifacts/current/native-laminar-cylinder-pressure-profile.cpp")),
            "profileBinarySha256": sha256(profile_binary),
            "driverSourceSha256": sha256(Path(__file__)),
        },
        "fixedSetup": {
            "polygon": str(POLYGON),
            "polygonSha256": sha256(POLYGON),
            "segments": 128,
            "centre_m": [0.07, 0.03],
            "diameter_m": 2,
            "outerDomainDiameterSpans": 20,
            "speed_m_s": 1,
            "kinematicViscosity_m2_s": 0.1,
            "reynoldsNumber": 20,
            "solverTolerance": 1e-8,
            "unchangedSolverMinimumVolumeRatio": 0.01,
        },
        "nativeRuns": [summarized_case(*case) for case in cases],
        "newGridRunRecords": run_records,
        "fixedGeometrySpatialPressureProfile": json.loads(profile_path.read_text()),
        "repair": {
            "status": "causal-sensitivity-not-default",
            "targetQuantity": (
                "minimum internal-face area ratio min(A_owner,A_neighbour)/"
                "max(A_owner,A_neighbour), dimensionless"
            ),
            "experimentalTarget": 0.075,
            "unchangedDefaultSolverLimit": 0.01,
            "nativeRepair": repair_record,
            "pressureProfile": json.loads(repair_profile_path.read_text()),
        },
        "interpretation": {
            "finding": (
                "The fixed polygon does not yield a monotone spatial sequence. The "
                "middle grid passes every existing Solver and strict residual gate but "
                "contains a large local velocity/pressure branch and cancellation of "
                "opposite pressure and viscous wall loads. Upwind has the same failure, "
                "so the Newton default and face-limited convection are not its cause."
            ),
            "causalControl": (
                "The existing conservative exact-union/split repair changes only 24 "
                "transactions, preserves total area and every physical boundary atom, "
                "passes the unchanged Solver gate, and removes the large local branch."
            ),
            "qualification": (
                "This proves a mesh-topology sensitivity and supplies a deterministic "
                "default-path counterexample. It does not justify promoting 0.075: prior "
                "cross-case evidence shows non-monotone response and a systematic force "
                "shift. A discretization-consistency objective and cross-case regression "
                "remain required before any product repair or default change."
            ),
        },
        "productChange": False,
        "defaultChange": False,
        "solverGateChange": False,
        "rawEvidence": [
            str(root / "far-20-fixed128-runs.json"),
            str(profile_path),
            str(repair_profile_path),
            str(root / "far-20-fixed128-grid1-volume-075-repair.json"),
            str(root / "far-20-fixed128-grid0-face-limited-linear") + ".*",
            str(root / "far-20-fixed128-grid1-face-limited-linear") + ".*",
            str(root / "far-20-fixed128-grid1-upwind") + ".*",
            str(root / "far-20-fixed128-grid1-volume-075-face-limited-linear") + ".*",
            str(root / "far-20-face-limited-linear") + ".*",
        ],
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grids", nargs="+", type=int, default=[0, 1])
    parser.add_argument("--timeout", type=float, default=3600)
    parser.add_argument("--skip-runs", action="store_true")
    parser.add_argument("--profile-binary", type=Path,
                        default=Path("outputs/cloud-laminar/cylinder-pressure-profile"))
    args = parser.parse_args()
    if not POLYGON.exists():
        raise FileNotFoundError(POLYGON)
    records: list[dict] = []
    for grid in ([] if args.skip_runs else args.grids):
        if grid not in (0, 1):
            raise ValueError("only missing fixed-polygon grids 0 and 1 are generated")
        wall = 0.125 / (2**grid)
        background = 1.0 / (2**grid)
        prefix = ROOT / f"far-20-fixed128-grid{grid}"
        mesh_command = [
            "build/cartmesh2d_cli", str(POLYGON), str(prefix), "8", ".25", ".1",
            "exterior", str(prefix) + "-foam", "0", "0", "--size-field",
            "--reference-length", "2", "--wall-relative-size", str(wall),
            "--background-relative-size", str(background),
            "--far-field-spans", "20", "--cells-per-level", "3",
        ]
        record = {
            "grid": grid,
            "segments": 128,
            "wallRelativeSize": wall,
            "backgroundRelativeSize": background,
            "farFieldSpans": 20,
            "input": str(POLYGON),
            "inputSha256": sha256(POLYGON),
            "mesh": run(mesh_command, Path(str(prefix) + "-mesh.log"), args.timeout),
        }
        records.append(record)
        if record["mesh"]["returncode"] != 0:
            break
        mesh = Path(str(prefix) + ".solver.cm2d")
        record["meshSha256"] = sha256(mesh)
        quality = Path(str(prefix) + ".construction-quality.json")
        record["constructionQuality"] = json.loads(quality.read_text())
        output = Path(str(prefix) + "-face-limited-linear")
        solve_command = [
            "build/cartmesh2d_flow_cli", "--mesh", str(mesh), "--case", "external",
            "--nu", ".1", "--speed", "1", "--convection", "face-limited-linear",
            "--tolerance", "1e-8", "--max-iterations", "2500", "--output", str(output),
        ]
        record["solve"] = run(solve_command, Path(str(output) + ".log"), args.timeout)
        summary = Path(str(output) + ".json")
        if record["solve"]["returncode"] == 0 and summary.exists():
            record["summary"] = selected_summary(summary)
            record["rawSha256"] = {
                suffix: sha256(Path(str(output) + suffix))
                for suffix in (".json", ".cells.csv", ".faces.csv", ".residuals.csv")
            }
        (ROOT / "far-20-fixed128-runs.json").write_text(
            json.dumps(records, indent=2) + "\n")
        print(json.dumps(record, sort_keys=True), flush=True)
    if not args.skip_runs and not (len(records) == len(args.grids) and all(
        row["mesh"]["returncode"] == 0 and row.get("solve", {}).get("returncode") == 0
        and row.get("summary", {}).get("converged") for row in records)):
        return 1
    if args.skip_runs:
        analyze(args.profile_binary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
