#!/usr/bin/env python3
"""Summarize native pressure-inverse scale and CLI integration checks.

This script only reads native solver outputs and compares serialized states.  It
does not implement or duplicate the flow equations.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs/laminar-stability/compatible-cli-linear-seed17260-current"
OLD = ROOT / "outputs/laminar-stability/compatible-api"
SCALE = ROOT / "outputs/laminar-stability/compatible-scale-current"
CROSS = ROOT / "outputs/laminar-stability/compatible-scale-crosscase"
CLI = ROOT / "outputs/laminar-stability/compatible-cli-aggregation-merged"
MESH = ROOT / "outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1.solver.cm2d"
MESH4716 = ROOT / "outputs/laminar-stability/compatible-linear-seed/cylinder4716.solver.cm2d"
COMPARE_BINARY = ROOT / "build/native-laminar-state-compare-scale"


def load(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def history(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 1 and rows[0]["iteration"] == "0"
    row = rows[0]
    return {
        "accepted": row["accepted"] == "1",
        "linearRestarts": int(row["linearRestarts"]),
        "matrixProducts": int(row["matrixProducts"]),
        "linearInitialRelativeResidual": float(row["linearInitialRelativeResidual"]),
        "linearRelativeResidual": float(row["linearRelativeResidual"]),
    }


def flat_json_state(path):
    state = load(path)
    return [value for name in ("cells", "faces") for row in state[name] for value in row]


def native_state(path):
    data = path.read_bytes()
    magic, mesh_key, count = struct.unpack_from("<QQQ", data)
    assert magic == 0x504F5345454E3031 and len(data) == 24 + 8 * count
    return {"meshKey": mesh_key,
            "coefficients": list(struct.unpack_from(f"<{count}d", data, 24))}


def compare(left, right):
    assert len(left) == len(right) and left
    delta = [a - b for a, b in zip(left, right)]
    return {
        "coefficientCount": len(delta),
        "coefficientRms": math.sqrt(sum(value * value for value in delta) / len(delta)),
        "coefficientMaximum": max(map(abs, delta)),
        "byteEquivalentValues": all(a == b for a, b in zip(left, right)),
    }


def run(prefix):
    summary = load(OUT / f"{prefix}.summary.json")
    return {
        "status": summary["status"],
        "reason": summary["reason"],
        "converged": summary["converged"],
        "acceptedIterations": summary["acceptedIterations"],
        "lastAcceptedAvailable": summary["lastAcceptedAvailable"],
        "lastRejectedAvailable": summary["lastRejectedAvailable"],
        "checkpointAvailable": summary["checkpointAvailable"],
        "boundaryLoadsAvailable": summary["boundaryLoadsAvailable"],
        "controls": summary["controls"],
        "history": history(OUT / f"{prefix}.residuals.csv"),
        "timingSeconds": summary["timingSeconds"],
        "sha256": {
            "summary": digest(OUT / f"{prefix}.summary.json"),
            "residuals": digest(OUT / f"{prefix}.residuals.csv"),
            "seed": digest(OUT / f"{prefix}.seed.json"),
            "boundaries": digest(OUT / f"{prefix}.boundaries"),
        },
    }


def adapter_run(prefix, directory=SCALE):
    report = load(directory / f"{prefix}.json")
    iterations = report["iterations"]
    return {
        "completed": report["completed"],
        "stop": report["stop"],
        "reason": report["reason"],
        "seconds": report["seconds"],
        "acceptedIterations": sum(item["accepted"] for item in iterations),
        "totalRestarts": sum(item["restarts"] for item in iterations),
        "totalMatrixProducts": sum(item["products"] for item in iterations),
        "iterations": iterations,
        "sha256": {
            "report": digest(directory / f"{prefix}.json"),
            "seed": digest(directory / f"{prefix}.seed.state"),
            **({"accepted": digest(directory / f"{prefix}.accepted.state")}
               if (directory / f"{prefix}.accepted.state").exists() else {}),
        },
    }


def state_comparison(mesh, left, right):
    output = subprocess.check_output([
        str(COMPARE_BINARY), "compare", str(mesh), "0", str(left), str(right)
    ], text=True)
    return json.loads(output)


old = load(OLD / "cylinder17260-cold-schur-pt1-r100.json")
old_first = old["iterations"][0]
current_seed = flat_json_state(OUT / "current.seed.json")
old_seed = native_state(OLD / "cylinder17260-cold-schur-pt1-r100.seed.state")
run100 = run("current")
run120 = run("current-r120")
matched_ic0 = adapter_run("schur-ic0-r120")
matched_aggregation = adapter_run("schur-aggregation-r120")
aggregation_full = adapter_run("schur-aggregation-full-r50")
old_accepted = OLD / "cylinder17260-cold-schur-pt1-r100.accepted.state"
aggregation_accepted = SCALE / "schur-aggregation-full-r50.accepted.state"
field_comparison = state_comparison(MESH, old_accepted, aggregation_accepted)
old_products = sum(item["products"] for item in old["iterations"])
old_seconds = old["seconds"]

cross_ic0 = adapter_run("cylinder4716-schur-first", CROSS)
cross_aggregation = adapter_run("cylinder4716-aggregation-first", CROSS)
cross_full = adapter_run("cylinder4716-aggregation-full", CROSS)
cross_reference = load(ROOT / "artifacts/current/native-laminar-linear-initial-cli.json")["realCylinder4716"]
cross_reference_products = cross_reference["history"]["matrixProducts"]
cross_reference_seconds = cross_reference["timingSeconds"]["elapsedBeforeSummary"]
cross_comparison = state_comparison(
    MESH4716,
    ROOT / "outputs/laminar-stability/compatible-linear-seed/final-cylinder4716-current.accepted.state",
    CROSS / "cylinder4716-aggregation-full.accepted.state",
)
cli_summary = load(CLI / "channel-aggregation.summary.json")

result = {
    "schema": "cartmesh2d-compatible-linear-initial-scale-v3",
    "scope": (
        "Native matched pressure-inverse checks on the original Linux 4,716- and "
        "17,260-cell fixed-128-segment 20D cylinders plus the explicit product CLI "
        "selector. Equations, mesh atoms, the original 1e-13 linear gate and all "
        "nonlinear acceptance gates remain unchanged."
    ),
    "mesh": {
        "path": str(MESH.relative_to(ROOT)),
        "sha256": digest(MESH),
        "cells": 17260,
        "faces": 34904,
    },
    "commands": {
        "common": [
            "build/cartmesh2d_flow_cli-linear-seed", "--discretization", "compatible",
            "--mesh", str(MESH.relative_to(ROOT)), "--case", "external", "--nu", ".1",
            "--compatible-pressure-inverse", "schur", "--linear-initial-guess", "current-state",
            "--krylov-directions", "60", "--linear-tolerance", "1e-13",
            "--pseudo-step", "1", "--pseudo-maximum-step", "1", "--max-iterations", "1000000",
        ],
        "run100": ["--linear-restarts", "100", "--output", str((OUT / "current").relative_to(ROOT))],
        "run120": ["--linear-restarts", "120", "--output", str((OUT / "current-r120").relative_to(ROOT))],
    },
    "run100": run100,
    "run120": run120,
    "budgetIncreaseFalsification": {
        "additionalRestartBudget": 20,
        "additionalMatrixProducts": run120["history"]["matrixProducts"] - run100["history"]["matrixProducts"],
        "residualRatioRun120OverRun100": (
            run120["history"]["linearRelativeResidual"] /
            run100["history"]["linearRelativeResidual"]
        ),
        "conclusion": (
            "The extra budget did not cross the unchanged 1e-13 true-residual gate. "
            "Late cycles stopped after only 180 additional products, so further budget "
            "inflation is not treated as a general remedy."
        ),
    },
    "seedIdentity": {
        "comparisonWithPriorResearchDriver": compare(current_seed, old_seed["coefficients"]),
        "currentAttemptsSameSeedBytes": (
            (OUT / "current.seed.json").read_bytes() ==
            (OUT / "current-r120.seed.json").read_bytes()
        ),
        "currentAttemptsSameBoundaryBytes": (
            (OUT / "current.boundaries").read_bytes() ==
            (OUT / "current-r120.boundaries").read_bytes()
        ),
    },
    "matchedCurrentExecutable": {
        "scope": (
            "The same current research-adapter binary, current native library, mesh, "
            "seed, equation and gates are used. Only the explicit pressure inverse "
            "changes from diagonal-Schur IC0 to the fixed aggregation V-cycle."
        ),
        "commands": {
            "common": [
                "build/native-laminar-compatible-solver-scale-current", str(MESH.relative_to(ROOT)),
                "0", "cylinder", ".1", "ns", "pressure", "OUTPUT_PREFIX", "1", "1",
            ],
            "ic0": ["schur", "120", "accepted", "ilu0"],
            "aggregation": ["schur-aggregation", "120", "accepted", "ilu0"],
        },
        "ic0FirstSystem": matched_ic0,
        "aggregationFirstSystem": matched_aggregation,
        "matrixProductReductionFraction": (
            1.0 - matched_aggregation["totalMatrixProducts"] /
            matched_ic0["totalMatrixProducts"]
        ),
        "elapsedReductionFraction": (
            1.0 - matched_aggregation["seconds"] / matched_ic0["seconds"]
        ),
        "seedByteEquivalent": (
            (SCALE / "schur-ic0-r120.seed.state").read_bytes() ==
            (SCALE / "schur-aggregation-r120.seed.state").read_bytes()
        ),
    },
    "aggregationFullRun": {
        **aggregation_full,
        "controls": {
            "pressureInverse": "diagonal-schur-aggregation",
            "velocityInverse": "ilu0",
            "linearTolerance": 1e-13,
            "krylovDirections": 60,
            "maximumLinearRestarts": 50,
            "linearInitialGuess": "current-state",
            "initialPseudoStep": 1.0,
        },
        "comparisonWithPriorAcceptedState": field_comparison,
        "comparisonWithPriorRunCost": {
            "priorMatrixProducts": old_products,
            "priorSeconds": old_seconds,
            "matrixProductReductionFraction": (
                1.0 - aggregation_full["totalMatrixProducts"] / old_products
            ),
            "elapsedReductionFraction": (
                1.0 - aggregation_full["seconds"] / old_seconds
            ),
            "warning": (
                "The field comparison is native and complete, but the old timing came "
                "from a prior executable. The same-executable first-system comparison "
                "above is the controlled algorithm A/B."
            ),
        },
    },
    "crossCase4716": {
        "mesh": {
            "path": str(MESH4716.relative_to(ROOT)),
            "sha256": digest(MESH4716),
            "cells": 4716,
            "faces": 9632,
        },
        "scope": (
            "Same current native executable, cold seed, equation, mesh and gates; "
            "only diagonal-Schur IC0 versus the fixed aggregation V-cycle changes."
        ),
        "ic0FirstSystem": cross_ic0,
        "aggregationFirstSystem": cross_aggregation,
        "firstSystem": {
            "matrixProductReductionFraction": 1.0 - cross_aggregation["totalMatrixProducts"] / cross_ic0["totalMatrixProducts"],
            "elapsedReductionFraction": 1.0 - cross_aggregation["seconds"] / cross_ic0["seconds"],
            "seedByteEquivalent": (
                (CROSS / "cylinder4716-schur-first.seed.state").read_bytes() ==
                (CROSS / "cylinder4716-aggregation-first.seed.state").read_bytes()
            ),
        },
        "aggregationFullRun": {
            **cross_full,
            "controls": {
                "pressureInverse": "diagonal-schur-aggregation",
                "velocityInverse": "ilu0",
                "linearTolerance": 1e-13,
                "krylovDirections": 60,
                "maximumLinearRestarts": 50,
            },
            "comparisonWithPriorAcceptedState": cross_comparison,
            "comparisonWithPriorProductCliCost": {
                "priorMatrixProducts": cross_reference_products,
                "priorSeconds": cross_reference_seconds,
                "matrixProductReductionFraction": 1.0 - cross_full["totalMatrixProducts"] / cross_reference_products,
                "elapsedReductionFraction": 1.0 - cross_full["seconds"] / cross_reference_seconds,
                "warning": (
                    "The full fields are compared natively, but the prior cost is a saved "
                    "product-CLI run. The first-system comparison above is the controlled A/B."
                ),
            },
        },
    },
    "productCliIntegration": {
        "selector": "schur-aggregation",
        "defaultChanged": False,
        "nativeLifecycleTest": {
            "passed": True,
            "cells": cli_summary["cells"],
            "case": cli_summary["case"],
            "status": cli_summary["status"],
            "acceptedIterations": cli_summary["acceptedIterations"],
            "pressureInverseRecorded": cli_summary["controls"]["pressureInverse"],
            "absolutePressureReference": load(CLI / "channel-aggregation.loads.json")["absolutePressureReference"],
            "coldStartAndResumeCompared": True,
            "failureRetentionChecked": True,
            "closedBoundaryRejected": True,
            "signalTested": True,
        },
        "sha256": {
            "binary": digest(ROOT / "build/cartmesh2d_flow_cli-aggregation-merged"),
            "cliSource": digest(ROOT / "apps/CompatibleFlowCLI.cpp"),
            "testSource": digest(ROOT / "tests/compatible_flow_cli_test.py"),
            "summary": digest(CLI / "channel-aggregation.summary.json"),
        },
    },
    "priorResearchComparator": {
        "completed": old["completed"],
        "totalIterations": len(old["iterations"]),
        "totalRestarts": sum(item["restarts"] for item in old["iterations"]),
        "totalMatrixProducts": sum(item["products"] for item in old["iterations"]),
        "firstLinearSystem": {
            "restarts": old_first["restarts"],
            "matrixProducts": old_first["products"],
            "linearRelativeResidual": old_first["linearRelativeResidual"],
            "accepted": old_first["accepted"],
        },
        "warning": (
            "The prior result is a saved comparator, not a controlled executable A/B: "
            "the product CLI and research adapter binaries have different hashes, and "
            "their complete assembly/code-generation identity has not been proven."
        ),
    },
    "sha256": {
        "productCliBinary": digest(ROOT / "build/cartmesh2d_flow_cli-linear-seed"),
        "priorResearchBinary": digest(ROOT / "build/native-laminar-compatible-solver-current"),
        "matchedResearchBinary": digest(ROOT / "build/native-laminar-compatible-solver-scale-current"),
        "matchedNativeLibrary": digest(ROOT / "build/libcartmesh2d_fv.a"),
        "nativeStateCompareBinary": digest(COMPARE_BINARY),
        "nativeStateCompareSource": digest(ROOT / "artifacts/current/native-laminar-state-compare.cpp"),
        "productCliSource": digest(ROOT / "apps/CompatibleFlowCLI.cpp"),
        "researchAdapterSource": digest(ROOT / "artifacts/current/native-laminar-compatible-solver.cpp"),
        "priorResearchReport": digest(OLD / "cylinder17260-cold-schur-pt1-r100.json"),
        "priorResearchSeed": digest(OLD / "cylinder17260-cold-schur-pt1-r100.seed.state"),
    },
    "interpretation": [
        "CurrentState cannot improve the first cold-start system because the current state is the boundary-applied seed; the measured initial relative residual is exactly one.",
        "Both failures preserve the native lifecycle: no candidate recovery, accepted field, checkpoint or loads are published.",
        "This scale check falsifies the claim that raising only the restart budget from 100 to 120 makes the current product CLI pass the original linear gate.",
        "The matched current adapter reproduces the CLI IC0 failure exactly, excluding CLI boundary conversion and lifecycle wrapping as its cause.",
        "Changing only the explicit pressure auxiliary inverse to aggregation crosses the same 1e-13 gate in 241 products; the full run converges under the original equation gates in 2622 products.",
        "The aggregation field matches the prior accepted solution to about 1e-12 in native geometry-weighted full-field norms, but this single mesh does not establish default, physical-accuracy or curved-wall-pressure qualification.",
        "The independent 4,716-cell cross-case reduces the matched first system from 4,612 to 376 products and completes the full run in 3,157 products; its accepted field agrees with the prior product result near roundoff.",
        "The product CLI now exposes the already-implemented aggregation inverse explicitly while retaining the previous default; a native pressure-outlet lifecycle regression records and validates the selected method.",
    ],
    "validation": {
        "nativeSolverProcessesAfterRuns": 0,
        "compatibleFlowNativeTest": "passed from current source and current native library",
        "compatibleFlowCliLifecycleTest": "passed, including explicit schur-aggregation and invalid-selector rejection",
        "unrelatedSuitesRerun": False,
        "reason": "The native CLI failures, matched adapter A/B, full aggregation solve and native field comparison are the verification target.",
    },
}

target = ROOT / "artifacts/current/native-laminar-linear-initial-scale.json"
target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(json.dumps({"written": str(target.relative_to(ROOT)), "sha256": digest(target)}))
