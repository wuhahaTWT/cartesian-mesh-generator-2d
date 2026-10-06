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


def normalized_field_distance(a: str, b: str, reference_speed: float = 1.0) -> dict:
    """Compare all cells while removing only the incompressible pressure gauge."""
    left, right = cell_rows(a), cell_rows(b)
    if len(left) != len(right):
        raise RuntimeError("field comparison requires identical final cells")
    total_area = sum(float(row["area"]) for row in left)
    pressure_difference = [float(x["p"]) - float(y["p"]) for x, y in zip(left, right)]
    pressure_gauge = sum(
        float(row["area"]) * value for row, value in zip(left, pressure_difference)
    ) / total_area
    velocity_square = 0.0
    pressure_square = 0.0
    maximum_velocity = 0.0
    maximum_pressure = 0.0
    for x, y, raw_pressure in zip(left, right, pressure_difference):
        area = float(x["area"])
        du = float(x["u"]) - float(y["u"])
        dv = float(x["v"]) - float(y["v"])
        dp = raw_pressure - pressure_gauge
        velocity_square += area * (du * du + dv * dv)
        pressure_square += area * dp * dp
        maximum_velocity = max(maximum_velocity, math.hypot(du, dv))
        maximum_pressure = max(maximum_pressure, abs(dp))
    return {
        "definition": "all-cell area-weighted; pressure difference has its area mean removed",
        "velocityRmsOverReferenceSpeed": math.sqrt(velocity_square / total_area) / reference_speed,
        "maximumVelocityVectorDifferenceOverReferenceSpeed": maximum_velocity / reference_speed,
        "pressureRmsOverReferenceSpeedSquared": math.sqrt(pressure_square / total_area) / reference_speed**2,
        "maximumGaugeAlignedPressureDifferenceOverReferenceSpeedSquared": maximum_pressure / reference_speed**2,
    }


def normalized_field_distance_paths(a: Path, b: Path, reference_speed: float = 1.0) -> dict:
    """Apply the same all-cell comparison to arbitrary exported native fields."""
    with a.open(newline="") as stream:
        left = list(csv.DictReader(stream))
    with b.open(newline="") as stream:
        right = list(csv.DictReader(stream))
    if len(left) != len(right):
        raise RuntimeError("field comparison requires identical final cells")
    total_area = sum(float(row["area"]) for row in left)
    pressure_difference = [float(x["p"]) - float(y["p"]) for x, y in zip(left, right)]
    pressure_gauge = sum(
        float(row["area"]) * value for row, value in zip(left, pressure_difference)
    ) / total_area
    velocity_square = pressure_square = maximum_velocity = maximum_pressure = 0.0
    for x, y, raw_pressure in zip(left, right, pressure_difference):
        area = float(x["area"])
        du = float(x["u"]) - float(y["u"])
        dv = float(x["v"]) - float(y["v"])
        dp = raw_pressure - pressure_gauge
        velocity_square += area * (du * du + dv * dv)
        pressure_square += area * dp * dp
        maximum_velocity = max(maximum_velocity, math.hypot(du, dv))
        maximum_pressure = max(maximum_pressure, abs(dp))
    return {
        "definition": "all-cell area-weighted; pressure difference has its area mean removed",
        "velocityRmsOverReferenceSpeed": math.sqrt(velocity_square / total_area) / reference_speed,
        "maximumVelocityVectorDifferenceOverReferenceSpeed": maximum_velocity / reference_speed,
        "pressureRmsOverReferenceSpeedSquared": math.sqrt(pressure_square / total_area) / reference_speed**2,
        "maximumGaugeAlignedPressureDifferenceOverReferenceSpeedSquared": maximum_pressure / reference_speed**2,
    }


def accuracy_summary(prefix: str) -> dict:
    path = ROOT / f"outputs/cloud-laminar/{prefix}.json"
    raw = json.loads(path.read_text())
    result = dict(raw)
    result["prefix"] = str(path.with_suffix("").relative_to(ROOT))
    result["rawSha256"] = {
        suffix: sha256(ROOT / f"outputs/cloud-laminar/{prefix}{suffix}")
        for suffix in (".json", ".cells.csv")
    }
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
grid1_short_nu02 = summary("grid1-flat-nu-0.2", "far-20-fixed128-grid1-short-nu02-flat")
grid1_short_from_nu02 = summary(
    "grid1-nu-0.1-from-flat-nu-0.2", "far-20-fixed128-grid1-short-nu01-from-flat02"
)
grid1_short_nu05 = summary("grid1-flat-nu-0.5", "far-20-fixed128-grid1-short-nu05-flat")
grid1_short_from_nu05 = summary(
    "grid1-nu-0.1-from-flat-nu-0.5", "far-20-fixed128-grid1-short-nu01-from-nu05"
)
grid1_short_from_nu1 = summary(
    "grid1-nu-0.1-from-flat-nu-1", "far-20-fixed128-grid1-short-nu01-from-nu1"
)
grid0_short_from_nu1 = summary(
    "grid0-nu-0.1-from-flat-nu-1", "far-20-fixed128-grid0-short-nu01-from-nu1"
)
grid0_direct = summary("grid0-direct", "far-20-fixed128-grid0-face-limited-linear")
grid2_direct = summary("grid2-direct", "far-20-face-limited-linear")
repaired = summary("grid1-volume-0.075-control",
                   "far-20-fixed128-grid1-volume-075-face-limited-linear")
grid1_anderson = summary("grid1-anderson-flat", "far-20-fixed128-grid1-anderson-flat")
grid1_anderson_from_nu1 = summary(
    "grid1-anderson-from-nu-1", "far-20-fixed128-grid1-anderson-from-nu1"
)
grid0_anderson = summary("grid0-anderson-flat", "far-20-fixed128-grid0-anderson-flat")
channel_default = accuracy_summary("channel-face-limited-linear-64")
channel_anderson = accuracy_summary("channel-face-limited-linear-64-anderson")
cavity_default = accuracy_summary("cavity-face-limited-linear-64")
cavity_anderson = accuracy_summary("cavity-face-limited-linear-64-anderson")

external_prefix = "far-20-fixed128-grid1-face-limited-linear"
custom_prefix = "far-20-fixed128-grid1-custom-flat"
equivalent = {
    suffix: sha256(RUN / f"{external_prefix}{suffix}") == sha256(RUN / f"{custom_prefix}{suffix}")
    for suffix in (".cells.csv", ".faces.csv", ".residuals.csv")
}

grid1_final = grid1_stages[-1]
grid0_final = grid0_stages[-1]
four_stage_cost = sum(stage["coupledEvaluations"] for stage in grid1_stages)
nu1_two_stage_cost = grid1_stages[0]["coupledEvaluations"] + grid1_short_from_nu1["coupledEvaluations"]
nu05_two_stage_cost = grid1_short_nu05["coupledEvaluations"] + grid1_short_from_nu05["coupledEvaluations"]
nu02_two_stage_cost = grid1_short_nu02["coupledEvaluations"] + grid1_short_from_nu02["coupledEvaluations"]
evidence = {
    "schema": "cartmesh2d.native-laminar-branch-continuation.v3",
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
    "grid1ShortContinuation": {
        "paths": [
            {
                "viscosities": [0.2, 0.1],
                "stages": [grid1_short_nu02, grid1_short_from_nu02],
                "totalCoupledEvaluations": nu02_two_stage_cost,
                "boundedTargetBranch": False,
                "finding": "Both stages pass every strict gate but retain a high-amplitude branch.",
            },
            {
                "viscosities": [0.5, 0.1],
                "stages": [grid1_short_nu05, grid1_short_from_nu05],
                "totalCoupledEvaluations": nu05_two_stage_cost,
                "boundedTargetBranch": True,
                "finalDistanceToFourStage": normalized_field_distance(
                    "far-20-fixed128-grid1-short-nu01-from-nu05",
                    "far-20-fixed128-grid1-custom-continuation",
                ),
            },
            {
                "viscosities": [1.0, 0.1],
                "stages": [grid1_stages[0], grid1_short_from_nu1],
                "totalCoupledEvaluations": nu1_two_stage_cost,
                "boundedTargetBranch": True,
                "finalDistanceToFourStage": normalized_field_distance(
                    "far-20-fixed128-grid1-short-nu01-from-nu1",
                    "far-20-fixed128-grid1-custom-continuation",
                ),
            },
        ],
        "lowestVerifiedBoundedCost": nu1_two_stage_cost,
        "costReductionRelativeToFourStage": 1.0 - nu1_two_stage_cost / four_stage_cost,
        "costRatioToDirectFlat": nu1_two_stage_cost / flat["coupledEvaluations"],
        "directVersusRecoveredBranchDistance": normalized_field_distance(
            "far-20-fixed128-grid1-custom-flat",
            "far-20-fixed128-grid1-short-nu01-from-nu1",
        ),
        "finding": (
            "A single nu=1 guide followed by the target equation reaches the same bounded branch "
            "as four-stage continuation with 35.8% fewer evaluations. Starting at nu=0.2 does "
            "not; this is a branch-disagreement certificate on the retained counterexample, not "
            "a universal viscosity threshold or an automatic branch-selection rule."
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
        "twoStage": {
            "viscosities": [1.0, 0.1],
            "stages": [grid0_stages[0], grid0_short_from_nu1],
            "totalCoupledEvaluations": (
                grid0_stages[0]["coupledEvaluations"] + grid0_short_from_nu1["coupledEvaluations"]
            ),
            "costRatioToDirect": (
                grid0_stages[0]["coupledEvaluations"] + grid0_short_from_nu1["coupledEvaluations"]
            ) / grid0_direct["coupledEvaluations"],
            "finalDistanceToDirect": normalized_field_distance(
                "far-20-fixed128-grid0-short-nu01-from-nu1",
                "far-20-fixed128-grid0-face-limited-linear",
            ),
        },
    },
    "independentAlgorithmCrosscheck": {
        "method": (
            "Existing explicit Anderson acceleration solves the identical target equation and "
            "uses the unchanged strict final gates; it is an independent nonlinear path, not a "
            "physical branch selector."
        ),
        "badGrid": {
            "defaultNewton": flat,
            "flatAnderson": grid1_anderson,
            "nu1InitializedAnderson": grid1_anderson_from_nu1,
            "newtonToFlatAndersonDistance": normalized_field_distance(
                "far-20-fixed128-grid1-custom-flat",
                "far-20-fixed128-grid1-anderson-flat",
            ),
            "recoveredNormalToFlatAndersonDistance": normalized_field_distance(
                "far-20-fixed128-grid1-short-nu01-from-nu1",
                "far-20-fixed128-grid1-anderson-flat",
            ),
            "dualFlatPathIterations": flat["iterations"] + grid1_anderson["iterations"],
            "continuationGuideAndTargetEvaluations": nu1_two_stage_cost,
            "finding": (
                "Flat-start Anderson converges faster than flat-start Newton but reaches an even "
                "larger-amplitude branch. Starting Anderson from the bounded nu=1 field reaches a "
                "third, load-cancelling target branch, so Anderson can expose non-uniqueness but "
                "cannot preserve or choose the bounded branch."
            ),
        },
        "normalGrid": {
            "defaultNewton": grid0_direct,
            "flatAnderson": grid0_anderson,
            "distance": normalized_field_distance(
                "far-20-fixed128-grid0-face-limited-linear",
                "far-20-fixed128-grid0-anderson-flat",
            ),
            "dualPathIterations": grid0_direct["iterations"] + grid0_anderson["iterations"],
        },
        "straightChannel64": {
            "defaultNewton": channel_default,
            "flatAnderson": channel_anderson,
            "distance": normalized_field_distance_paths(
                ROOT / "outputs/cloud-laminar/channel-face-limited-linear-64.cells.csv",
                ROOT / "outputs/cloud-laminar/channel-face-limited-linear-64-anderson.cells.csv",
            ),
        },
        "closedCavity64": {
            "defaultNewton": cavity_default,
            "flatAnderson": cavity_anderson,
            "distance": normalized_field_distance_paths(
                ROOT / "outputs/cloud-laminar/cavity-face-limited-linear-64.cells.csv",
                ROOT / "outputs/cloud-laminar/cavity-face-limited-linear-64-anderson.cells.csv",
            ),
        },
        "qualification": (
            "Agreement on three normal controls and disagreement on the retained counterexample "
            "supports an explicit branch-risk certificate. It does not identify which converged "
            "branch is physical, so no automatic acceptance, fallback, or default change is made."
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
            "Two-stage viscosity continuation recovered the same bounded branch on the known bad "
            "mesh with 35.8% fewer evaluations than the four-stage path and reproduced the existing "
            "bounded branch on the coarse mesh. Its 1.42x direct-solve cost, limited mesh evidence "
            "and unresolved branch-selection semantics remain insufficient for a default change."
        ),
    },
}

OUTPUT.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
print(OUTPUT)
