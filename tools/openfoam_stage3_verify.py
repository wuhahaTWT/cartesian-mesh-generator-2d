#!/usr/bin/env python3
"""Use an existing local OpenFOAM image to verify a complete polyMesh."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True, type=Path)
    parser.add_argument("--project-report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log-output", type=Path)
    parser.add_argument("--image", default="opencfd/openfoam-run:2606")
    parser.add_argument("--milestone", choices=("stage3", "stage4"),
                        default="stage3")
    args = parser.parse_args()

    case = args.case.resolve()
    report = json.loads(args.project_report.read_text())
    required = ["points", "faces", "owner", "neighbour", "boundary"]
    poly_mesh = case / "constant" / "polyMesh"
    missing = [name for name in required if not (poly_mesh / name).is_file()]
    if missing:
        raise SystemExit(f"missing OpenFOAM polyMesh files: {missing}")
    if not report.get("stage3GeometryTopologyComplete", False):
        raise SystemExit("project geometry/topology report did not pass")
    expected_output = report.get("completeSolverVolumeMeshOutput")
    if not expected_output:
        raise SystemExit("project report does not name a complete volume mesh")

    # checkMesh writes region sets for valid multi-region meshes. Run it on a
    # disposable local copy so the source evidence remains byte-stable.
    with tempfile.TemporaryDirectory(prefix="cartmesh-checkmesh-") as temporary:
        checker_case = Path(temporary) / "case"
        shutil.copytree(case, checker_case)
        command = [
            "docker", "run", "--rm", "--network", "none",
            "-v", f"{checker_case}:/case", args.image,
            "checkMesh", "-case", "/case", "-constant", "-allTopology",
        ]
        completed = subprocess.run(
            command, check=False, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    log = completed.stdout
    if args.log_output is not None:
        args.log_output.parent.mkdir(parents=True, exist_ok=True)
        args.log_output.write_text(log, encoding="utf-8")
    mesh_ok = "Mesh OK." in log and "Failed " not in log

    def integer(pattern: str) -> int | None:
        match = re.search(pattern, log)
        return int(match.group(1)) if match else None

    def number(pattern: str) -> float | None:
        match = re.search(pattern, log)
        return float(match.group(1)) if match else None

    numeric = r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    build_match = re.search(r"^Build\s*:\s*(\S+)", log, re.MULTILINE)
    version_match = re.search(r"Version:\s*(\S+)", log)
    result = {
        "schema": "cartmesh-stage3-openfoam-checkmesh-v1",
        "status": "pass" if mesh_ok else "fail",
        "stage3Complete": mesh_ok,
        "stage4GeometryRobustnessComplete": report.get(
            "stage4GeometryRobustnessComplete", False),
        "stage4Complete": bool(
            args.milestone == "stage4" and mesh_ok and
            report.get("stage4GeometryRobustnessComplete", False)),
        "solverReadyCutCellMesh": mesh_ok,
        "externalCfdCheckerAccepted": mesh_ok,
        "checker": "OpenFOAM checkMesh",
        "checkerImage": args.image,
        "checkerOptions": ["-constant", "-allTopology"],
        "command": command,
        "logOutput": str(args.log_output.resolve()) if args.log_output else None,
        "checkerBuild": build_match.group(1) if build_match else None,
        "checkerVersion": version_match.group(1) if version_match else None,
        "networkDisabled": True,
        "sourceCaseUnmodified": True,
        "temporaryWritableCaseCopy": True,
        "case": str(case),
        "projectResultHashFnv1a64": report.get("resultHashFnv1a64"),
        "points": integer(r"points:\s+(\d+)"),
        "faces": integer(r"faces:\s+(\d+)"),
        "internalFaces": integer(r"internal faces:\s+(\d+)"),
        "cells": integer(r"cells:\s+(\d+)"),
        "boundaryPatches": integer(r"boundary patches:\s+(\d+)"),
        "minimumVolume": number(r"Min volume =\s+" + numeric),
        "maximumVolume": number(r"Max volume =\s+" + numeric),
        "totalVolume": number(r"Total volume =\s+" + numeric),
        "maximumNonOrthogonality": number(
            r"Mesh non-orthogonality Max:\s+" + numeric),
        "maximumSkewness": number(r"Max skewness =\s+" + numeric),
        "polyMeshSha256": {name: sha256(poly_mesh / name) for name in required},
        "returnCode": completed.returncode,
        "meshOkMarker": mesh_ok,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, ensure_ascii=False))
    return 0 if mesh_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
