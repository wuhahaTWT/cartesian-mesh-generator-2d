#!/usr/bin/env python3
"""Native Re=20 circular-cylinder joint geometry/grid refinement.

The script only orchestrates the product mesh and flow CLIs and reads their
native summaries.  It does not implement another flow equation or alter an
acceptance gate.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time


ROOT = Path("outputs/cloud-laminar/cylinder-joint")
REFERENCE_CD = 2.045


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    return {"command": command, "log": str(log), "status": status,
            "returncode": returncode, "processSeconds": time.monotonic()-started}


def polygon(path: Path, segments: int) -> None:
    centre = (0.07, 0.03)
    with path.open("w") as stream:
        for index in range(segments):
            angle = 2*math.pi*index/segments
            stream.write(f"{centre[0]+math.cos(angle):.17g} "
                         f"{centre[1]+math.sin(angle):.17g}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--levels", nargs="+", type=int, default=[0, 1, 2])
    parser.add_argument("--schemes", nargs="+",
                        default=["upwind", "face-limited-linear"])
    parser.add_argument("--timeout", type=float, default=900)
    args = parser.parse_args()
    ROOT.mkdir(parents=True, exist_ok=True)
    records = []
    for level in args.levels:
        segments = 32*(2**level)
        wall = 0.125/(2**level)
        background = 1.0/(2**level)
        prefix = ROOT/f"level-{level}"
        source = Path(str(prefix)+".xy")
        polygon(source, segments)
        row = {
            "level": level, "segments": segments,
            "wallRelativeSize": wall, "backgroundRelativeSize": background,
            "inputSha256": sha256(source),
            "polygonArea_m2": 0.5*segments*math.sin(2*math.pi/segments),
            "circleAreaDifference_m2": 0.5*segments*math.sin(2*math.pi/segments)-math.pi,
            "polygonPerimeter_m": 2*segments*math.sin(math.pi/segments),
        }
        mesh_command = [
            "build/cartmesh2d_cli", str(source), str(prefix), "8", ".25", ".1",
            "exterior", str(prefix)+"-foam", "0", "0", "--size-field",
            "--reference-length", "2", "--wall-relative-size", str(wall),
            "--background-relative-size", str(background), "--far-field-spans", "10",
            "--cells-per-level", "3",
        ]
        row["mesh"] = run(mesh_command, Path(str(prefix)+"-mesh.log"), args.timeout)
        row["solves"] = []
        if row["mesh"]["returncode"] == 0:
            mesh = Path(str(prefix)+".solver.cm2d")
            row["meshSha256"] = sha256(mesh)
            for scheme in args.schemes:
                output = Path(str(prefix)+"-"+scheme)
                command = [
                    "build/cartmesh2d_flow_cli", "--mesh", str(mesh), "--case", "external",
                    "--nu", ".1", "--speed", "1", "--convection", scheme,
                    "--tolerance", "1e-8", "--max-iterations", "2500",
                    "--output", str(output),
                ]
                result = run(command, Path(str(output)+".log"), args.timeout)
                result["scheme"] = scheme
                summary_path = Path(str(output)+".json")
                if result["returncode"] == 0 and summary_path.exists():
                    summary = json.loads(summary_path.read_text())
                    result["metrics"] = {
                        "cells": summary["cells"], "evaluations": summary["iterations"],
                        "converged": summary["converged"], "continuity": summary["continuity"],
                        "momentumResidual": summary["momentumResidual"],
                        "dragCoefficient": summary["wallForceX"],
                        "pressureDragCoefficient": summary["pressureForceX"],
                        "viscousDragCoefficient": summary["wallViscousForceX"],
                        "referenceDifference": summary["wallForceX"]-REFERENCE_CD,
                        "fieldSha256": sha256(Path(str(output)+".cells.csv")),
                    }
                row["solves"].append(result)
                print(json.dumps({"level": level, "segments": segments, **result}), flush=True)
        records.append(row)
        (ROOT/"runs.json").write_text(json.dumps(records, indent=2)+"\n")
    return 0 if all(x["mesh"]["returncode"] == 0 and
                    all(y["returncode"] == 0 for y in x["solves"])
                    for x in records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
