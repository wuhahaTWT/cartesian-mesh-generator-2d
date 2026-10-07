#!/usr/bin/env python3
"""Summarize native aggregation controls without implementing flow equations."""

import hashlib
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs/laminar-stability"
COMPARE = ROOT / "build/native-laminar-state-compare-crosscase"
SOLVER = ROOT / "build/native-laminar-compatible-solver-crosscase"
MERGED_SOLVER = ROOT / "build/native-laminar-compatible-solver-crosscase-merged"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def report(prefix):
    path = OUT / f"{prefix}.json"
    value = json.loads(path.read_text())
    iterations = value["iterations"]
    return {
        "completed": value["completed"],
        "stop": value["stop"],
        "reason": value["reason"],
        "cells": value["cells"],
        "faces": value["faces"],
        "seconds": value["seconds"],
        "acceptedIterations": sum(row["accepted"] for row in iterations),
        "totalRestarts": sum(row["restarts"] for row in iterations),
        "totalMatrixProducts": sum(row["products"] for row in iterations),
        "lastMetrics": {
            key: iterations[-1][key]
            for key in (
                "linearRelativeResidual", "cellMomentum", "faceMomentum",
                "divergence", "residualNorm", "stateChange"
            )
        },
        "iterations": iterations,
        "sha256": {
            "report": digest(path),
            "seed": digest(OUT / f"{prefix}.seed.state"),
            "accepted": digest(OUT / f"{prefix}.accepted.state"),
        },
    }


def compare(fixture, resolution, left, right):
    output = subprocess.check_output([
        str(COMPARE), "compare", fixture, str(resolution),
        str(OUT / f"{left}.accepted.state"),
        str(OUT / f"{right}.accepted.state"),
    ], text=True)
    return json.loads(output)


def pair(fixture, resolution, ic0_prefix, aggregation_prefix):
    ic0 = report(ic0_prefix)
    aggregation = report(aggregation_prefix)
    assert ic0["completed"] and aggregation["completed"]
    assert ic0["sha256"]["seed"] == aggregation["sha256"]["seed"]
    return {
        "fixture": fixture,
        "resolution": resolution,
        "cells": ic0["cells"],
        "faces": ic0["faces"],
        "controls": {
            "problem": "outlet-poiseuille",
            "viscosity": 0.1,
            "equation": "steady Navier-Stokes",
            "boundary": "product-like fixed static pressure via pseudo-traction",
            "initialState": "zero",
            "initialPseudoStep": 0.1,
            "maximumOuterIterations": 40,
            "maximumLinearRestartsPerSystem": 50,
            "linearTolerance": 1e-13,
        },
        "ic0": ic0,
        "aggregation": aggregation,
        "matrixProductReductionFraction": (
            1.0 - aggregation["totalMatrixProducts"] / ic0["totalMatrixProducts"]
        ),
        "elapsedReductionFraction": 1.0 - aggregation["seconds"] / ic0["seconds"],
        "acceptedStateComparison": compare(
            fixture, resolution, ic0_prefix, aggregation_prefix
        ),
    }


square16 = pair(
    "square", 16,
    "aggregation-crosscase-square16-schur-run",
    "aggregation-crosscase-square16-aggregation-run",
)
square32 = pair(
    "square", 32,
    "aggregation-crosscase-square32-schur-run",
    "aggregation-crosscase-square32-aggregation-run",
)
cut4 = pair(
    "cut", 4,
    "aggregation-crosscase-cut4-schur-run",
    "aggregation-crosscase-cut4-aggregation-run",
)
closed16 = pair(
    "square", 16,
    "aggregation-crosscase-closed16-schur-merged",
    "aggregation-crosscase-closed16-aggregation-merged",
)

failure_paths = {
    "sheared16PreMergeUnsupportedAdapter": OUT / "aggregation-crosscase-sheared16-adapter-failure.stderr",
    "cut16OriginalQualityGate": OUT / "aggregation-crosscase-cut16-quality-failure.stderr",
}
failures = {
    name: {
        "message": path.read_text().strip(),
        "sha256": digest(path),
        "bytes": path.stat().st_size,
    }
    for name, path in failure_paths.items()
}

result = {
    "schema": "cartmesh2d-compatible-aggregation-crosscase-v1",
    "scope": (
        "Matched native executable pairs, equations, zero seed, 1e-13 true linear "
        "target and original nonlinear acceptance gates. Only the explicit "
        "diagonal-Schur pressure inverse changes from IC0 to one fixed aggregation "
        "V-cycle. Open controls use the pre-merge executable; the closed control uses "
        "the ordinarily merged closed-Schur implementation. Python only summarizes "
        "native outputs and invokes the native state comparator."
    ),
    "commands": {
        "buildSolver": (
            "g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror "
            "-Wno-unused-parameter -I include artifacts/current/"
            "native-laminar-compatible-solver.cpp build/libcartmesh2d_fv.a "
            "build/libcartmesh2d.a -o build/native-laminar-compatible-solver-crosscase"
        ),
        "runTemplate": (
            "build/native-laminar-compatible-solver-crosscase FIXTURE N "
            "outlet-poiseuille .1 ns pressure FRESH_PREFIX .1 40 "
            "schur|schur-aggregation 50 zero"
        ),
        "buildMergedClosedSolver": (
            "g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror "
            "-Wno-unused-parameter -I include artifacts/current/"
            "native-laminar-compatible-solver.cpp src/fv/CompatibleIncompressible2D.cpp "
            "build/libcartmesh2d_fv.a build/libcartmesh2d.a -o "
            "build/native-laminar-compatible-solver-crosscase-merged"
        ),
    },
    "cases": {
        "square16": square16,
        "square32": square32,
        "trueLineCut14Cells": cut4,
        "closedSquare16Merged": closed16,
    },
    "failedControls": failures,
    "sha256": {
        "solverBinary": digest(SOLVER),
        "mergedClosedSolverBinary": digest(MERGED_SOLVER),
        "solverSource": digest(ROOT / "artifacts/current/native-laminar-compatible-solver.cpp"),
        "stateCompareBinary": digest(COMPARE),
        "stateCompareSource": digest(ROOT / "artifacts/current/native-laminar-state-compare.cpp"),
        "nativeLibrary": digest(ROOT / "build/libcartmesh2d_fv.a"),
        "mergedCompatibleSource": digest(ROOT / "src/fv/CompatibleIncompressible2D.cpp"),
        "mergedLinearHeader": digest(ROOT / "include/cartmesh2d/fv/detail/CompatibleFlowLinear2D.hpp"),
    },
    "interpretation": [
        "The 16x16 structured control reduces products by about 28 percent, but setup cost leaves wall time nearly unchanged.",
        "The 32x32 structured control reduces products by about 50 percent and native solve time by about 11 percent, showing scale-dependent benefit outside the cylinder cases.",
        "The accepted structured fields agree near roundoff, so the explicit pressure inverse does not select a different discrete state in these controls.",
        "The accepted 14-cell true line-cut control has essentially no product benefit and is slightly slower; aggregation is not universally cheaper on small Cut-cell systems.",
        "After the ordinary merge that enables closed-domain Schur inverses, the 16x16 no-slip control reduces products by about 39 percent but has nearly unchanged elapsed time; gauge-free accepted fields agree near roundoff.",
        "The finer line-cut fixture remains rejected by the pre-existing Solver quality gate. No cell was deleted and no gate was relaxed.",
        "The pre-merge sheared open-boundary attempt was rejected because this research adapter did not classify an oblique pressure outlet. Closed-domain Schur support does not add that missing outlet classification, so the run does not qualify or falsify aggregation on general oblique outlets.",
    ],
    "validation": {
        "nativeSolverProcessesAfterRuns": 0,
        "unrelatedSuitesRerun": False,
        "defaultChanged": False,
        "thresholdsChanged": False,
    },
}

target = ROOT / "artifacts/current/native-laminar-aggregation-crosscase.json"
target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(json.dumps({"written": str(target.relative_to(ROOT)), "sha256": digest(target)}))
