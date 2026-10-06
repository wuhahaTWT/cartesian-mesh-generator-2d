#!/usr/bin/env python3
"""Summarize fixed-cylinder nonlinear-branch and geometry-only diagnostics.

This postprocessor reads native CLI fields/summaries and the native topology
spectrum.  It does not solve an alternate equation or remove any cells from
the reported norms.  Paths are repository-relative so the archived raw runs
remain independently auditable.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "outputs/cloud-laminar/cylinder-joint"
OUTPUT = ROOT / "artifacts/current/native-laminar-branch-continuation.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def summary(name: str, prefix: str) -> dict:
    path = RUN / f"{prefix}.json"
    raw = json.loads(path.read_text())
    keys = (
        "case", "status", "cells", "iterations", "coupledEvaluations", "nu",
        "momentumInertia", "globalRelativeImbalance", "maximumSpeedRatio",
        "pressureRangeRatio", "wallForceX", "pressureForceX", "wallViscousForceX",
    )
    result = {key: raw.get(key) for key in keys}
    result["name"] = name
    result["prefix"] = str(path.with_suffix("" ).relative_to(ROOT))
    result["rawSha256"] = {
        suffix: sha256(RUN / f"{prefix}{suffix}")
        for suffix in (".json", ".cells.csv", ".faces.csv", ".residuals.csv")
    }
    return result


def cell_rows(prefix: str) -> list[dict[str, str]]:
    with (RUN / f"{prefix}.cells.csv").open(newline="") as stream:
        return list(csv.DictReader(stream))


def field_difference(a: str, b: str) -> dict:
    left, right = cell_rows(a), cell_rows(b)
    if len(left) != len(right):
        raise RuntimeError("field comparison requires identical final cells")
    result = {}
    for field in ("u", "v", "p"):
        result[f"maximumAbsolute{field.upper()}Difference"] = max(
            abs(float(x[field]) - float(y[field])) for x, y in zip(left, right)
        )
    return result


def topology(prefix: str) -> dict:
    path = RUN / f"{prefix}.topology-spectrum.csv"
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    converted = []
    for row in rows:
        converted.append({
            "cell": int(row["cell"]),
            "x": float(row["x"]),
            "y": float(row["y"]),
            "error": float(row["quadratic_consistency_error"]),
            "areaRatio": float(row["minimum_neighbour_area_ratio"]),
            "nonorthogonalRatio": float(row["maximum_nonorthogonal_ratio"]),
            "boundaryFaces": int(row["boundary_faces"]),
        })
    wall = sorted((row for row in converted if row["boundaryFaces"]),
                  key=lambda row: row["error"], reverse=True)
    interior = sorted((row for row in converted if not row["boundaryFaces"]),
                      key=lambda row: row["error"], reverse=True)
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256(path),
        "cells": len(rows),
        "maximumWallQuadraticConsistencyError": wall[0]["error"],
        "maximumInteriorQuadraticConsistencyError": interior[0]["error"],
        "eightWorstWallCells": wall[:8],
    }


def relative_change(a: float, b: float) -> float:
    return abs(a - b) / abs(b)


flat = summary("grid1-flat-initial", "far-20-fixed128-grid1-custom-flat")
stokes = summary("grid1-stokes", "far-20-fixed128-grid1-custom-stokes")
from_stokes = summary("grid1-from-stokes", "far-20-fixed128-grid1-custom-from-stokes")
grid1_stages = [
    summary("grid1-nu-1", "far-20-fixed128-grid1-custom-nu1"),
    summary("grid1-nu-0.5", "far-20-fixed128-grid1-custom-nu05"),
    summary("grid1-nu-0.2", "far-20-fixed128-grid1-custom-nu02"),
    summary("grid1-nu-0.1", "far-20-fixed128-grid1-custom-continuation"),
]
grid0_stages = [
    summary("grid0-nu-1", "far-20-fixed128-grid0-custom-nu1"),
    summary("grid0-nu-0.5", "far-20-fixed128-grid0-custom-nu05"),
    summary("grid0-nu-0.2", "far-20-fixed128-grid0-custom-nu02"),
    summary("grid0-nu-0.1", "far-20-fixed128-grid0-custom-continuation"),
]
grid0_direct = summary("grid0-direct", "far-20-fixed128-grid0-face-limited-linear")
grid2_direct = summary("grid2-direct", "far-20-face-limited-linear")
repaired = summary("grid1-volume-0.075-control",
                   "far-20-fixed128-grid1-volume-075-face-limited-linear")

external_prefix = "far-20-fixed128-grid1-face-limited-linear"
custom_prefix = "far-20-fixed128-grid1-custom-flat"
equivalent = {
    suffix: sha256(RUN / f"{external_prefix}{suffix}") == sha256(RUN / f"{custom_prefix}{suffix}")
    for suffix in (".cells.csv", ".faces.csv", ".residuals.csv")
}

grid1_final = grid1_stages[-1]
grid0_final = grid0_stages[-1]
evidence = {
    "schema": "cartmesh2d.native-laminar-branch-continuation.v1",
    "scope": (
        "Native fixed-geometry Re=20 cylinder branch selection and geometry-only "
        "quadratic consistency; no cells omitted and no product acceptance threshold added."
    ),
    "externalToCustomEquivalence": {
        "boundaryFile": str((RUN / "far-20-fixed128-grid1.external.boundaries").relative_to(ROOT)),
        "boundarySha256": sha256(RUN / "far-20-fixed128-grid1.external.boundaries"),
        "byteIdentical": equivalent,
        "allByteIdentical": all(equivalent.values()),
    },
    "grid1Branches": [flat, stokes, from_stokes, *grid1_stages],
    "grid1Continuation": {
        "viscosities": [stage["nu"] for stage in grid1_stages],
        "intermediateCoupledEvaluations": sum(
            stage["coupledEvaluations"] for stage in grid1_stages[:-1]
        ),
        "totalCoupledEvaluations": sum(stage["coupledEvaluations"] for stage in grid1_stages),
        "directFlatCoupledEvaluations": flat["coupledEvaluations"],
        "costRatioToDirectFlat": sum(
            stage["coupledEvaluations"] for stage in grid1_stages
        ) / flat["coupledEvaluations"],
        "finalWallForceBetweenAdjacentDirectMeshes": (
            min(grid0_direct["wallForceX"], grid2_direct["wallForceX"])
            <= grid1_final["wallForceX"]
            <= max(grid0_direct["wallForceX"], grid2_direct["wallForceX"])
        ),
        "finalWallForceRelativeDifferenceToFine": relative_change(
            grid1_final["wallForceX"], grid2_direct["wallForceX"]
        ),
    },
    "grid0NoBranchRegression": {
        "direct": grid0_direct,
        "continuationStages": grid0_stages,
        "finalFieldDifference": field_difference(
            "far-20-fixed128-grid0-face-limited-linear",
            "far-20-fixed128-grid0-custom-continuation",
        ),
        "wallForceAbsoluteDifference": abs(
            grid0_direct["wallForceX"] - grid0_final["wallForceX"]
        ),
    },
    "adjacentAndTopologyControls": {
        "grid2Direct": grid2_direct,
        "grid1VolumeRepair": repaired,
        "badGrid": topology("far-20-fixed128-grid1"),
        "volumeRepair": topology("far-20-fixed128-grid1-volume-075"),
        "finding": (
            "The geometry-only spectrum identifies all symmetric latent wall motifs, "
            "but its maximum remains large after the accepted field is normalized; "
            "it is not a sufficient standalone quality gate."
        ),
    },
    "artifacts": {
        "topologySource": "artifacts/current/native-laminar-topology-spectrum.cpp",
        "topologySourceSha256": sha256(ROOT / "artifacts/current/native-laminar-topology-spectrum.cpp"),
        "boundaryAdapterSource": "artifacts/current/native-laminar-external-boundary.cpp",
        "boundaryAdapterSourceSha256": sha256(ROOT / "artifacts/current/native-laminar-external-boundary.cpp"),
    },
    "qualification": {
        "productAlgorithmChanged": False,
        "defaultChanged": False,
        "qualityGateChanged": False,
        "relatedLinuxNativeTests": "3/3 passed: flow_boundary, solver_export, io",
        "fullNativeSuiteRerun": False,
        "frontendOrAppRerun": False,
        "conclusion": (
            "Viscosity continuation recovered a bounded branch on the known bad mesh "
            "and reproduced the existing bounded branch on the coarse mesh, but the "
            "2.21x evaluation cost and two-mesh evidence are insufficient for a default change."
        ),
    },
}

OUTPUT.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
print(OUTPUT)
