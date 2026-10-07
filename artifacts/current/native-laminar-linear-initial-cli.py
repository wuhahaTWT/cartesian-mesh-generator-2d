#!/usr/bin/env python3
"""Summarize native CLI initial-guess integration; no independent PDE audit."""
import csv
import hashlib
import json
import math
from pathlib import Path
import struct

ROOT = Path(__file__).resolve().parents[2]
SMALL = ROOT / "outputs/laminar-stability/compatible-cli-linear-seed-final"
REAL = ROOT / "outputs/laminar-stability/compatible-cli-linear-seed4716"
OLD = ROOT / "outputs/laminar-stability/compatible-linear-seed"


def load(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def history(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    return {
        "acceptedIterations": sum(row["accepted"] == "1" for row in rows),
        "linearRestarts": sum(int(row["linearRestarts"]) for row in rows),
        "matrixProducts": sum(int(row["matrixProducts"]) for row in rows),
        "linearInitialRelativeResiduals": [float(row["linearInitialRelativeResidual"]) for row in rows],
        "linearRelativeResiduals": [float(row["linearRelativeResidual"]) for row in rows],
    }


def flat(state):
    return [value for name in ("cells", "faces") for row in state[name] for value in row]


def compare(left, right):
    delta = [a - b for a, b in zip(left, right)]
    assert len(left) == len(right) and delta
    return {"coefficientRms": math.sqrt(sum(x*x for x in delta)/len(delta)),
            "coefficientMaximum": max(map(abs, delta))}


def native_state(path):
    data = path.read_bytes()
    magic, mesh_key, count = struct.unpack_from("<QQQ", data)
    assert magic == 0x504F5345454E3031 and len(data) == 24 + 8*count
    return {"meshKey": mesh_key, "coefficients": list(struct.unpack_from(f"<{count}d", data, 24))}


def load_summary(report, faces_path):
    with faces_path.open() as stream:
        faces = {int(row["face"]): row for row in csv.DictReader(stream)}
    wall = [entry for entry in report["boundaries"] if entry["embeddedBoundary"]]
    traction = [sum(entry["tractionMoments"][0][i] for entry in wall) for i in range(2)]
    pressure = [sum(entry["pressureMoments"][0][i] for entry in wall) for i in range(2)]
    weighted = []
    for entry in wall:
        face = faces[entry["face"]]
        length = math.hypot(float(face["Sx"]), float(face["Sy"]))
        weighted.append((length, float(face["x"]), float(face["y"])))
    perimeter = sum(row[0] for row in weighted)
    centre = [sum(row[0]*row[i+1] for row in weighted)/perimeter for i in range(2)]
    torque_origin = sum(entry["torqueOnFluid"] for entry in wall)
    torque_centre = torque_origin - centre[0]*traction[1] + centre[1]*traction[0]
    return {
        "absolutePressureReference": report["absolutePressureReference"],
        "boundaryFaces": len(report["boundaries"]),
        "embeddedBoundaryFaces": len(wall),
        "wallForceOnBody": [-value for value in traction],
        "wallPressureForceOnBody": [-value for value in pressure],
        "wallBoundaryLength": perimeter,
        "wallBoundaryLengthWeightedCentre": centre,
        "wallTorqueOnBodyAboutReportedOrigin": -torque_origin,
        "wallTorqueOnBodyAboutBoundaryCentre": -torque_centre,
        "globalBoundaryTractionOnFluid": report["boundaryTraction"],
        "globalBoundaryMomentumFlux": report["boundaryMomentumFlux"],
        "globalMomentumImbalance": report["momentumImbalance"],
        "globalVolumeFlux": report["boundaryVolumeFlux"],
    }


small_zero = load(SMALL / "cavity.accepted.json")
small_current = load(SMALL / "current-state.accepted.json")
small_zero_load = load(SMALL / "cavity.loads.json")
small_current_load = load(SMALL / "current-state.loads.json")
real_summary = load(REAL / "current.summary.json")
real_state = load(REAL / "current.accepted.json")
research_state = native_state(OLD / "final-cylinder4716-current.accepted.state")
test_log = REAL / "linux-tests.log"

result = {
    "schema": "cartmesh2d-compatible-linear-initial-cli-v1",
    "scope": "Explicit CLI exposure of the already validated algebraic current-state initial guess. The default remains zero; K, rhs, nonlinear equations, checkpoint context, conservative loads, tolerances, budgets and acceptance gates are unchanged.",
    "implementation": {
        "option": "--linear-initial-guess zero|current-state",
        "default": "zero",
        "historyAddition": "linearInitialRelativeResidual column",
        "checkpointBinding": "not bound: this is an algebraic strategy, and either choice may resume the same accepted physical checkpoint",
    },
    "smallNativeCli": {
        "cells": load(SMALL / "cavity.summary.json")["cells"],
        "defaultZero": history(SMALL / "cavity.residuals.csv"),
        "explicitZero": history(SMALL / "explicit-zero.residuals.csv"),
        "currentState": history(SMALL / "current-state.residuals.csv"),
        "defaultAndExplicitZeroCheckpointByteIdentical": (SMALL / "cavity.checkpoint").read_bytes() == (SMALL / "explicit-zero.checkpoint").read_bytes(),
        "currentStateDifference": compare(flat(small_zero), flat(small_current)),
        "boundaryTractionMaximumDifference": max(abs(a-b) for a, b in zip(small_zero_load["boundaryTraction"], small_current_load["boundaryTraction"])),
        "restartWithCurrentStateConverged": load(SMALL / "continued-current.summary.json")["converged"],
        "invalidChoiceRejected": (SMALL / "bad-linear-initial.stderr").read_text().strip(),
    },
    "realCylinder4716": {
        "mesh": "fixed 128-segment 20D Linux cylinder",
        "cells": real_summary["cells"],
        "faces": real_summary["faces"],
        "status": real_summary["status"],
        "timingSeconds": real_summary["timingSeconds"],
        "history": history(REAL / "current.residuals.csv"),
        "conservativeLoads": load_summary(load(REAL / "current.loads.json"), REAL / "current.faces.csv"),
        "differenceFromPriorResearchDriver": compare(flat(real_state), research_state["coefficients"]),
    },
    "linuxValidation": {
        "passed": 6,
        "total": 6,
        "tests": ["compatible loads", "compatible checkpoint", "compatible boundary", "compatible flow", "compatible element", "compatible CLI lifecycle"],
        "logSha256": digest(test_log),
    },
    "sha256": {
        "mesh4716": digest(OLD / "cylinder4716.solver.cm2d"),
        "binary": digest(ROOT / "build/cartmesh2d_flow_cli-linear-seed"),
        "cliSource": digest(ROOT / "apps/CompatibleFlowCLI.cpp"),
        "cliTest": digest(ROOT / "tests/compatible_flow_cli_test.py"),
        "checkpoint4716": digest(REAL / "current.checkpoint"),
        "loads4716": digest(REAL / "current.loads.json"),
        "accepted4716": digest(REAL / "current.accepted.json"),
    },
    "limitations": [
        "The current-state mode remains explicit because the real solve is still expensive and one 4716-cell nonlinear step previously required more products than zero initialization.",
        "This batch adds no spatial, curved-wall pressure or physical-branch acceptance criterion.",
        "No full native suite, frontend, Electron, macOS or packaged application validation was run.",
    ],
}
(ROOT / "artifacts/current/native-laminar-linear-initial-cli.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(json.dumps({"written": "artifacts/current/native-laminar-linear-initial-cli.json",
                  "sha256": digest(ROOT / "artifacts/current/native-laminar-linear-initial-cli.json")}))
