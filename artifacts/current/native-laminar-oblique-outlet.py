#!/usr/bin/env python3
"""Summarize native oblique-outlet runs; no PDE assembly or field recomputation."""
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "outputs/laminar-stability/aggregation-oblique"
OUT = ROOT / "artifacts/current/native-laminar-oblique-outlet.json"


def load(name):
    path = RAW / name
    return json.loads(path.read_text()), {
        "path": str(path.relative_to(ROOT)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def run(report, exact):
    data, report_id = load(report)
    error, exact_id = load(exact)
    return {
        "completed": data["completed"],
        "acceptedIterations": sum(i["accepted"] for i in data["iterations"]),
        "matrixProducts": sum(i["products"] for i in data["iterations"]),
        "linearRestarts": sum(i["restarts"] for i in data["iterations"]),
        "seconds": data["seconds"],
        "lastDivergencePerSecond": data["iterations"][-1]["divergence"],
        "nativeExactError": error,
        "files": {"solverReport": report_id, "nativeExactPostprocess": exact_id},
    }


pseudo = {
    "8": run("sheared8-pressure-aggregation.json", "sheared8-pressure-aggregation.exact.json"),
    "16": run("sheared16-aggregation.json", "sheared16-exact.json"),
    "32": run("sheared32-aggregation.json", "sheared32-exact.json"),
}
traction = {
    "8": run("sheared8-traction-aggregation.json", "sheared8-traction-aggregation.exact.json"),
    "16": run("sheared16-traction-aggregation.json", "sheared16-traction-exact.json"),
    "32": run("sheared32-traction-aggregation.json", "sheared32-traction-exact.json"),
}
for cases in (pseudo, traction):
    cases["orders"] = {}
    for field in ("velocityP1Rms", "faceVelocityRms"):
        cases["orders"][field] = [
            math.log(cases[a]["nativeExactError"][field] / cases[b]["nativeExactError"][field], 2)
            for a, b in (("8", "16"), ("16", "32"))
        ]

comparisons = {}
for n in (16, 32):
    ic0, ic0_id = load(f"sheared{n}-schur.json")
    agg, agg_id = load(f"sheared{n}-aggregation.json")
    field, field_id = load(f"sheared{n}-compare.json")
    pi = sum(i["products"] for i in ic0["iterations"])
    pa = sum(i["products"] for i in agg["iterations"])
    comparisons[str(n)] = {
        "schurIC0Products": pi,
        "aggregationProducts": pa,
        "productReductionFraction": 1 - pa / pi,
        "schurIC0Seconds": ic0["seconds"],
        "aggregationSeconds": agg["seconds"],
        "elapsedReductionFraction": 1 - agg["seconds"] / ic0["seconds"],
        "acceptedStateDifference": field,
        "files": {"schurIC0": ic0_id, "aggregation": agg_id, "comparison": field_id},
    }

binary = ROOT / "build/native-laminar-compatible-solver-oblique"
post = ROOT / "build/native-laminar-state-compare-oblique"
solver_source = ROOT / "artifacts/current/native-laminar-compatible-solver.cpp"
post_source = ROOT / "artifacts/current/native-laminar-state-compare.cpp"
result = {
    "scope": "Native sheared parallelogram manufactured outlet control; product defaults unchanged",
    "geometry": {"shear": 0.7, "cells": [64, 256, 1024], "realPolygonMesh": True},
    "equation": {"problem": "outlet-poiseuille", "viscosity": 0.1, "equation": "steady Navier-Stokes"},
    "boundaryIdentity": {
        "pseudoTraction": "(nu G-pI)n=-p_D n",
        "exactTraction": "(nu(G+G^T)-pI)n=t_exact",
        "obstruction": "For u=(4y(1-y),0), the slanted outlet has n_y!=0 and G n=((4-8y)n_y,0), so p_D=p does not satisfy pseudo-traction except at y=0.5.",
    },
    "pseudoTraction": pseudo,
    "exactTraction": traction,
    "pressureInverseAB": comparisons,
    "build": {
        "solverSha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "postprocessorSha256": hashlib.sha256(post.read_bytes()).hexdigest(),
        "solverSourceSha256": hashlib.sha256(solver_source.read_bytes()).hexdigest(),
        "postprocessorSourceSha256": hashlib.sha256(post_source.read_bytes()).hexdigest(),
        "solverCommand": "g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -I include artifacts/current/native-laminar-compatible-solver.cpp src/fv/CompatibleIncompressible2D.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -o build/native-laminar-compatible-solver-oblique",
        "postprocessorCommand": "g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -I include artifacts/current/native-laminar-state-compare.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -o build/native-laminar-state-compare-oblique",
        "postprocessorRole": "Native polynomial field integration only; no PDE assembly or acceptance gate",
    },
    "runControls": {
        "solverTemplate": "native-laminar-compatible-solver-oblique sheared N outlet-poiseuille .1 ns pressure|traction PREFIX .1 40 schur|schur-aggregation 50 zero",
        "initialState": "zero",
        "initialPseudoStep": 0.1,
        "maximumOuterIterations": 40,
        "maximumLinearRestarts": 50,
        "productQualityGate": "unchanged",
    },
    "conclusion": "The compatible discretization and aggregation inverse solve the oblique system, but product-like static-pressure pseudo-traction is not a coordinate-invariant exact outlet model. Full traction restores second-order velocity convergence on this control.",
    "limitations": [
        "Three manufactured parallelogram meshes are not a curved-wall or backflow qualification.",
        "Pressure errors under exact traction are roundoff-scale and do not define an observed pressure order.",
        "No product default, quality gate, threshold, or boundary input format changed.",
    ],
}
OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(json.dumps({"output": str(OUT.relative_to(ROOT)), "sha256": hashlib.sha256(OUT.read_bytes()).hexdigest()}))
