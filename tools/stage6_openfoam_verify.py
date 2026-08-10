#!/usr/bin/env python3
"""Capture an honest OpenFOAM checkMesh result for a Stage-6 binary polyMesh.

The checker runs on a disposable local copy with Docker networking disabled.
It separates successful format/core-topology reading from the stricter default
geometry-quality result so a readable mesh is never mislabeled as ``Mesh OK``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path


REQUIRED_POLYMESH_FILES = ("points", "faces", "owner", "neighbour", "boundary")
NUMBER = r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def integer(log: str, *patterns: str) -> int | None:
    for pattern in patterns:
        match = re.search(pattern, log, re.IGNORECASE | re.MULTILINE)
        if match:
            return int(match.group(1))
    return None


def number(log: str, *patterns: str) -> float | None:
    for pattern in patterns:
        match = re.search(pattern, log, re.IGNORECASE | re.MULTILINE)
        if match:
            return float(match.group(1))
    return None


def marker(log: str, text: str) -> bool:
    return text.lower() in log.lower()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log-output", required=True, type=Path)
    parser.add_argument("--image", default="opencfd/openfoam-run:2606")
    args = parser.parse_args()

    source_case = args.case.resolve()
    poly_mesh = source_case / "constant" / "polyMesh"
    missing = [name for name in REQUIRED_POLYMESH_FILES if not (poly_mesh / name).is_file()]
    if missing:
        raise SystemExit(f"missing OpenFOAM polyMesh files: {missing}")

    start = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="cartmesh-stage6-checkmesh-") as temporary:
        checker_case = Path(temporary) / "case"
        shutil.copytree(source_case, checker_case)
        command = [
            "docker",
            "run",
            "--rm",
            "--network",
            "none",
            "-v",
            f"{checker_case}:/case",
            args.image,
            "checkMesh",
            "-case",
            "/case",
            "-constant",
            "-allTopology",
        ]
        completed = subprocess.run(
            command,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    wall_seconds = time.perf_counter() - start
    log = completed.stdout

    mesh_ok = marker(log, "Mesh OK.") and not re.search(
        r"Failed\s+\d+\s+mesh checks", log, re.IGNORECASE
    )
    core_topology_markers = {
        "boundaryDefinitionOk": marker(log, "Boundary definition OK"),
        "cellToFaceAddressingOk": marker(log, "Cell to face addressing OK"),
        "upperTriangularOrderingOk": marker(log, "Upper triangular ordering OK"),
        "faceVerticesOk": marker(log, "Face vertices OK"),
        "cellZipUpOk": marker(log, "Cell zip-up check OK")
        or marker(log, "Topological cell zip-up check OK"),
        "singleRegionOk": bool(
            re.search(r"Number of regions:\s*1\s*\(OK\)", log, re.IGNORECASE)
        ),
    }
    format_read = all(
        integer(log, pattern) is not None
        for pattern in (
            r"points:\s+(\d+)",
            r"faces:\s+(\d+)",
            r"internal faces:\s+(\d+)",
            r"cells:\s+(\d+)",
        )
    )
    core_topology_pass = format_read and all(core_topology_markers.values())
    failed_checks = integer(log, r"Failed\s+(\d+)\s+mesh checks") or 0
    status = (
        "mesh_ok"
        if mesh_ok
        else "read_and_core_topology_pass_quality_fail"
        if core_topology_pass
        else "read_or_core_topology_fail"
    )

    build_match = re.search(r"^Build\s*:\s*(\S+)", log, re.MULTILINE)
    version_match = re.search(r"Version:\s*(\S+)", log)
    result = {
        "schema": "cartmesh-stage6-openfoam-checkmesh-v1",
        "status": status,
        "stage6Complete": False,
        "solverReadyCutCellMesh": mesh_ok,
        "externalCfdCheckerAccepted": mesh_ok,
        "checker": "OpenFOAM checkMesh",
        "checkerImage": args.image,
        "checkerOptions": ["-constant", "-allTopology"],
        "checkerBuild": build_match.group(1) if build_match else None,
        "checkerVersion": version_match.group(1) if version_match else None,
        "networkDisabled": True,
        "sourceCaseUnmodified": True,
        "temporaryWritableCaseCopy": True,
        "case": str(source_case),
        "formatRead": format_read,
        "coreTopologyPass": core_topology_pass,
        "coreTopologyMarkers": core_topology_markers,
        "meshOkMarker": mesh_ok,
        "failedMeshCheckCount": failed_checks,
        "returnCode": completed.returncode,
        "wallClockSecondsIncludingCopy": wall_seconds,
        "counts": {
            "points": integer(log, r"points:\s+(\d+)"),
            "faces": integer(log, r"faces:\s+(\d+)"),
            "internalFaces": integer(log, r"internal faces:\s+(\d+)"),
            "cells": integer(log, r"cells:\s+(\d+)"),
            "boundaryPatches": integer(log, r"boundary patches:\s+(\d+)"),
            "unusedPoints": integer(
                log,
                r"number unused by cells:\s*(\d+)",
            ) or 0,
            "wrongOrientedFacePyramids": integer(
                log,
                r"(\d+)\s+faces are incorrectly oriented",
                r"(\d+)\s+wrong-oriented face pyramids",
            ),
            "highlySkewFaces": integer(
                log,
                r"(\d+)\s+highly skew faces",
            ),
            "duplicateBaffleFaces": integer(
                log,
                r"Number of identical duplicate faces \(baffle faces\):\s*(\d+)",
            ),
            "nonStandardEdgeConnectivityFaces": integer(
                log,
                r"(\d+)\s+faces with non-standard edge connectivity",
            ),
            "neighbouringCellPairsWithMultipleFaces": integer(
                log,
                r"Found\s+(\d+)\s+neighbouring cells with multiple inbetween faces",
            ),
            "unorderedFacesWritten": integer(
                log,
                r"Writing\s+(\d+)\s+unordered faces",
            ),
            "cellsWithAtMostOneInternalFace": integer(
                log,
                r"Writing\s+(\d+)\s+cells with zero or one non-boundary face",
            ),
            "cellsWithTwoInternalFaces": integer(
                log,
                r"Writing\s+(\d+)\s+cells with two non-boundary faces",
            ),
        },
        "quality": {
            "minimumVolume": number(log, r"Min volume =\s*" + NUMBER),
            "maximumVolume": number(log, r"Max volume =\s*" + NUMBER),
            "totalVolume": number(log, r"Total volume =\s*" + NUMBER),
            "maximumNonOrthogonality": number(
                log, r"Mesh non-orthogonality Max:\s*" + NUMBER
            ),
            "maximumSkewness": number(log, r"Max skewness =\s*" + NUMBER),
            "boundaryOpennessOk": bool(
                re.search(r"Boundary openness .* OK\.", log, re.IGNORECASE)
            ),
            "faceAreasPositive": marker(log, "Face area magnitudes OK"),
            "cellVolumesPositive": marker(log, "Cell volumes OK"),
        },
        "polyMeshSha256": {
            name: sha256(poly_mesh / name) for name in REQUIRED_POLYMESH_FILES
        },
        "log": str(args.log_output),
    }

    args.log_output.parent.mkdir(parents=True, exist_ok=True)
    args.log_output.write_text(log, encoding="utf-8")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
