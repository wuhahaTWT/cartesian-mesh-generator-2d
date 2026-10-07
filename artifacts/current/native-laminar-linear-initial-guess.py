#!/usr/bin/env python3
"""Summarize native compatible linear-initial-guess A/B runs; no PDE reimplementation."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs/laminar-stability/compatible-linear-seed"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name):
    path = OUT / f"{name}.json"
    data = json.loads(path.read_text())
    assert data["completed"] and data["stop"] == 0
    return path, data


def pair(label, zero_name, current_name, compare_name):
    zero_path, zero = load(zero_name)
    current_path, current = load(current_name)
    compare_path = OUT / compare_name
    compare = json.loads(compare_path.read_text())
    assert zero["cells"] == current["cells"] == compare["cells"]
    assert len(zero["iterations"]) == len(current["iterations"])
    total = lambda data, key: sum(row[key] for row in data["iterations"])
    return {
        "label": label,
        "cells": zero["cells"],
        "faces": zero["faces"],
        "iterations": len(zero["iterations"]),
        "zero": {
            "seconds": zero["seconds"],
            "linearRestarts": total(zero, "restarts"),
            "matrixProducts": total(zero, "products"),
            "sha256": digest(zero_path),
        },
        "currentState": {
            "seconds": current["seconds"],
            "linearRestarts": total(current, "restarts"),
            "matrixProducts": total(current, "products"),
            "initialRelativeResiduals": [row["linearInitialRelativeResidual"] for row in current["iterations"]],
            "sha256": digest(current_path),
        },
        "relativeChange": {
            "seconds": current["seconds"] / zero["seconds"] - 1,
            "linearRestarts": total(current, "restarts") / total(zero, "restarts") - 1,
            "matrixProducts": total(current, "products") / total(zero, "products") - 1,
        },
        "fieldDifference": compare,
        "fieldDifferenceSha256": digest(compare_path),
        "perIteration": {
            "zero": [{key: row[key] for key in ("iteration", "linearInitialRelativeResidual", "restarts", "products", "linearRelativeResidual")} for row in zero["iterations"]],
            "currentState": [{key: row[key] for key in ("iteration", "linearInitialRelativeResidual", "restarts", "products", "linearRelativeResidual")} for row in current["iterations"]],
        },
    }


data = {
    "schemaVersion": 2,
    "scope": "Explicit algebraic current-state initial guess layered on the existing strict Krylov-cycle target, with FMA/TwoSum recomputation of the same assembled b-K*x residual. K, rhs, nonlinear equations, pseudo-time, preconditioner, tolerances, budgets and all original acceptance gates are unchanged.",
    "controls": {
        "equation": "NavierStokes",
        "viscosity": 0.1,
        "initialPseudoStep": 1.0,
        "maximumPseudoStep": 1e6,
        "linearTolerance": 1e-13,
        "krylovDirections": 60,
        "maximumLinearRestarts": 100,
        "pressureInverse": "diagonal-schur for pressure-reference outlet/cylinder; viscous-mass gauge inverse for closed flow",
        "quadratureOrder": 6,
    },
    "realCylinder4716": pair("fixed-128-segment 20D cylinder, original 4716-cell Linux mesh", "final-cylinder4716-zero", "final-cylinder4716-current", "final-cylinder4716-field-compare.json"),
    "openOutlet64": pair("64-cell analytic product-style open outlet", "final-outlet8-zero", "final-outlet8-current", "final-outlet8-field-compare.json"),
    "closedSheared64": pair("64-cell sheared stationary-wall manufactured flow", "final-sheared8-zero", "final-sheared8-current", "final-sheared8-field-compare.json"),
    "roundoffBoundary": {
        "ordinaryResidualRecompute": {
            "path": "outputs/laminar-stability/compatible-linear-seed/strict-cylinder4716-current.json",
            "completed": False,
            "firstLinearRestarts": 100,
            "firstMatrixProducts": 4164,
            "reportedRelativeResidual": 1.1263735857365017e-13,
            "sha256": digest(OUT / "strict-cylinder4716-current.json"),
        },
        "extendedAccumulationDiagnostic": {
            "path": "outputs/laminar-stability/compatible-linear-seed/strict-comp-cylinder4716-current.json",
            "acceptedStateSha256": digest(OUT / "strict-comp-cylinder4716-current.accepted.state"),
            "sha256": digest(OUT / "strict-comp-cylinder4716-current.json"),
        },
        "portableTwoSumFma": {
            "path": "outputs/laminar-stability/compatible-linear-seed/final-cylinder4716-current.json",
            "acceptedStateSha256": digest(OUT / "final-cylinder4716-current.accepted.state"),
            "sha256": digest(OUT / "final-cylinder4716-current.json"),
        },
        "extendedVsPortableFieldDifference": json.loads((OUT / "extended-vs-portable-cylinder4716-field-compare.json").read_text()),
        "extendedVsPortableFieldDifferenceSha256": digest(OUT / "extended-vs-portable-cylinder4716-field-compare.json"),
        "interpretation": "The ordinary separately rounded A*x path exhausted the unchanged 100-restart budget just above 1e-13. Extended accumulation and portable FMA/TwoSum independently crossed the unchanged target; their accepted fields agree to roundoff but are not byte-identical because FMA/TwoSum retains product error. This is an algebraic residual-evaluation repair, not a relaxed gate.",
    },
    "hashes": {
        "mesh4716": digest(OUT / "cylinder4716.solver.cm2d"),
        "binary": digest(ROOT / "build/native-laminar-compatible-solver-current"),
        **{str(path.relative_to(ROOT)): digest(path) for path in (
            ROOT / "include/cartmesh2d/fv/CompatibleIncompressible2D.hpp",
            ROOT / "include/cartmesh2d/fv/detail/CompatibleFlowLinear2D.hpp",
            ROOT / "src/fv/CompatibleIncompressible2D.cpp",
            ROOT / "artifacts/current/native-laminar-compatible-solver.cpp",
            ROOT / "tests/compatible_flow_test.cpp",
            ROOT / "tests/compatible_flow_checkpoint_test.cpp",
        )},
    },
    "linuxValidation": {
        "passed": 6,
        "scope": "compatible flow, compatible checkpoint, compatible boundary, legacy flow boundary, branch certificate, and real compatible CLI lifecycle",
        "outputsSha256": {
            name: digest(OUT / name) for name in (
                "final-compatible-flow-test-rerun.stdout",
                "final-compatible-checkpoint-test.stdout",
                "final-compatible-boundary-test.stdout",
                "final-flow-boundary-test.stdout",
                "final-flow-branch-test.stdout",
                "final-compatible-cli-test.stdout",
            )
        },
        "fullNativeSuiteRun": False,
        "frontendOrDesktopRun": False,
    },
    "limitations": [
        "The 4716-cell second nonlinear system regressed from 4291 to 4616 matrix products even though the current-state initial residual was 0.173 instead of 1; only complete same-source totals are used for the cost comparison.",
        "The real cylinder still required 17035 matrix products and 267.6 seconds with the current-state seed, so neither compensated residual evaluation nor this explicit seed by itself qualifies a scalable default backend.",
        "No 17260-cell rerun was made; no physical branch, curved-wall local pressure/traction, desktop or cross-platform qualification is claimed.",
    ],
}
print(json.dumps(data, indent=2, sort_keys=True))
