#!/usr/bin/env python3
"""Complete the Re=20 cylinder refinement sequence on the same 20D domain.

This driver only orchestrates the native mesh/flow CLIs and the compiled native
wall-pressure postprocessor.  It does not implement a flow equation or change
any product acceptance threshold.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time


ROOT = Path("outputs/cloud-laminar/cylinder-joint")
REFERENCE_CD = 2.045


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--levels", nargs="+", type=int, default=[0, 1])
    parser.add_argument("--timeout", type=float, default=3600)
    args = parser.parse_args()
    ROOT.mkdir(parents=True, exist_ok=True)
    records: list[dict] = []
    for level in args.levels:
        if level not in (0, 1):
            raise ValueError("only the missing 20D levels 0 and 1 are generated")
        source = ROOT / f"level-{level}.xy"
        if not source.exists():
            raise FileNotFoundError(source)
        wall = 0.125 / (2**level)
        background = 1.0 / (2**level)
        prefix = ROOT / f"far-20-level{level}"
        mesh_command = [
            "build/cartmesh2d_cli", str(source), str(prefix), "8", ".25", ".1",
            "exterior", str(prefix) + "-foam", "0", "0", "--size-field",
            "--reference-length", "2", "--wall-relative-size", str(wall),
            "--background-relative-size", str(background),
            "--far-field-spans", "20", "--cells-per-level", "3",
        ]
        record = {
            "level": level,
            "segments": 32 * (2**level),
            "wallRelativeSize": wall,
            "backgroundRelativeSize": background,
            "farFieldSpans": 20,
            "input": str(source),
            "inputSha256": sha256(source),
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
        (ROOT / "far-20-spatial-runs.json").write_text(
            json.dumps(records, indent=2) + "\n")
        print(json.dumps(record, sort_keys=True), flush=True)
    return 0 if len(records) == len(args.levels) and all(
        row["mesh"]["returncode"] == 0 and row.get("solve", {}).get("returncode") == 0
        and row.get("summary", {}).get("converged") for row in records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
