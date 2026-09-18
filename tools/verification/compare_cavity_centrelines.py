#!/usr/bin/env python3
"""Compare existing cavity centre-lines with independent FD references.

This is a diagnostic comparison only.  It does not launch a solver and does
not add or alter any cavity acceptance gate.
"""

from __future__ import annotations

import argparse
import csv
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


SIZES = (33, 65, 129)
SAMPLES = np.linspace(0.0, 1.0, 101)
LEVELS = {196: 4, 900: 5, 3844: 6}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fail(message: str) -> None:
    raise ValueError(message)


def finite(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not numeric") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} is not finite")
    return result


def load_fd_reference(root: Path) -> dict[str, Any]:
    summary_path = root / "summary.json"
    if not summary_path.is_file():
        fail(f"missing FD reference summary: {summary_path}")
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid FD reference JSON: {summary_path}: {exc}") from exc
    cases = summary.get("cases")
    if not isinstance(cases, dict):
        fail("FD reference summary cases must be an object")
    for n in SIZES:
        case = cases.get(f"N{n}")
        if not isinstance(case, dict) or case.get("converged") is not True:
            fail(f"FD reference N{n} is not converged")
    poisson = summary.get("poissonValidation", {})
    poisson_ok = (summary.get("poissonCheck") is True
                  or isinstance(poisson, dict) and poisson.get("passedRoundoff") is True)
    if not poisson_ok:
        fail("FD reference Poisson check is not true")

    profiles: dict[int, dict[str, list[float]]] = {}
    metadata: dict[str, Any] = {}
    for n in SIZES:
        arrays: dict[str, np.ndarray] = {}
        for component in ("u", "v"):
            path = root / f"N{n}_{component}.npy"
            if not path.is_file():
                fail(f"missing FD reference array: {path}")
            array = np.asarray(np.load(path, allow_pickle=False), dtype=float)
            if array.shape != (n, n):
                fail(f"N{n}_{component} has shape {array.shape}, expected {(n, n)}")
            if not np.isfinite(array).all():
                fail(f"N{n}_{component} contains non-finite values")
            declared = cases[f"N{n}"].get("arrays", {}).get(component, {})
            if isinstance(declared, dict) and declared.get("sha256"):
                actual_hash = sha256_file(path)
                if actual_hash != declared["sha256"]:
                    fail(f"N{n}_{component} hash differs from reference summary")
            arrays[component] = array
        axis = np.linspace(0.0, 1.0, n)
        # The FD arrays include their true wall values.  The vertical u line
        # uses the x=.5 column; horizontal v uses the y=.5 row.
        u = np.interp(SAMPLES, axis, arrays["u"][:, n // 2])
        v = np.interp(SAMPLES, axis, arrays["v"][n // 2, :])
        profiles[n] = {"u": u.tolist(), "v": v.tolist()}
        metadata[str(n)] = {
            "converged": True,
            "arraySha256": {component: sha256_file(root / f"N{n}_{component}.npy")
                             for component in ("u", "v")},
            "declaredArraySha256": {
                component: cases[f"N{n}"].get("arrays", {}).get(component, {}).get("sha256")
                for component in ("u", "v")
            },
            "shape": [n, n],
            "sourceCase": copy.deepcopy(cases[f"N{n}"]),
        }

    differences = {}
    for coarse, fine in ((33, 65), (65, 129)):
        cu, cv = np.asarray(profiles[coarse]["u"]), np.asarray(profiles[coarse]["v"])
        fu, fv = np.asarray(profiles[fine]["u"]), np.asarray(profiles[fine]["v"])
        differences[f"{coarse}vs{fine}"] = {
            "uRms": float(np.sqrt(np.mean((cu - fu) ** 2))),
            "vRms": float(np.sqrt(np.mean((cv - fv) ** 2))),
            "combinedRms": float(np.sqrt(np.mean(np.r_[cu - fu, cv - fv] ** 2))),
        }
    return {
        "root": str(root.resolve()),
        "summarySha256": sha256_file(summary_path),
        "summary": {
            "path": str(summary_path.resolve()),
            "format": summary.get("format"),
            "diagnosticOnly": summary.get("diagnosticOnly"),
            "poissonCheck": True,
            "assessment": copy.deepcopy(summary.get("assessment")),
        },
        "metadata": metadata,
        "profiles": profiles,
        "fdDifferenceRms": differences,
    }


def path_from_summary(raw: str, summary_path: Path) -> Path:
    path = Path(raw)
    return path if path.is_absolute() else (summary_path.parent / path).resolve()


def read_cell_profile(case: dict[str, Any], summary_path: Path) -> dict[str, Any]:
    if case.get("case") != "cavity" or case.get("valid") is not True:
        fail(f"internal error: selected case is not valid cavity: {case.get('label')}")
    native = case.get("native")
    measurement = case.get("meshMeasurement")
    if not isinstance(native, dict) or native.get("case") != "cavity":
        fail(f"{case.get('label')}: missing cavity native metadata")
    if native.get("converged") is not True or native.get("status") != "converged":
        fail(f"{case.get('label')}: native case is not converged")
    if not math.isclose(finite(native.get("speed"), "native speed"), 1.0, abs_tol=1e-12):
        fail(f"{case.get('label')}: only speed=1 cavity cases are supported")
    if not math.isclose(finite(native.get("nu"), "native nu"), 0.01, abs_tol=1e-14):
        fail(f"{case.get('label')}: only nu=0.01 cavity cases are supported")
    if not isinstance(measurement, dict):
        fail(f"{case.get('label')}: missing mesh measurement")
    bounds = measurement.get("bounds")
    if not isinstance(bounds, list) or len(bounds) != 4 or any(
            not math.isclose(finite(value, "mesh bound"), target, abs_tol=1e-10)
            for value, target in zip(bounds, (0.0, 0.0, 1.0, 1.0))):
        fail(f"{case.get('label')}: only unit-square cavity meshes are supported")

    prefix = path_from_summary(str(case.get("prefix")), summary_path)
    cells_path = Path(str(prefix) + ".cells.csv")
    if not cells_path.is_file():
        fail(f"{case.get('label')}: missing cells CSV {cells_path}")
    with cells_path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    expected_cells = int(case.get("counts", {}).get("cells", len(rows)))
    if len(rows) != expected_cells:
        fail(f"{case.get('label')}: CSV has {len(rows)} rows, expected {expected_cells}")
    if expected_cells not in LEVELS:
        fail(f"{case.get('label')}: unsupported cavity cell count {expected_cells}")

    xs = sorted({round(finite(row.get("x"), "cell x"), 12) for row in rows})
    ys = sorted({round(finite(row.get("y"), "cell y"), 12) for row in rows})
    if len(xs) * len(ys) != len(rows):
        fail(f"{case.get('label')}: rounded coordinates do not form a complete rectangle")
    u = np.full((len(ys), len(xs)), np.nan)
    v = np.full((len(ys), len(xs)), np.nan)
    x_index, y_index = {x: i for i, x in enumerate(xs)}, {y: i for i, y in enumerate(ys)}
    seen: set[tuple[float, float]] = set()
    for row in rows:
        x, y = round(finite(row.get("x"), "cell x"), 12), round(finite(row.get("y"), "cell y"), 12)
        key = (x, y)
        if key in seen:
            fail(f"{case.get('label')}: duplicate rounded cell coordinate {key}")
        seen.add(key)
        if not (0.0 < x < 1.0 and 0.0 < y < 1.0):
            fail(f"{case.get('label')}: cell centre is outside the unit square: {key}")
        u_value, v_value = finite(row.get("u"), "cell u"), finite(row.get("v"), "cell v")
        if not math.isclose(finite(row.get("speed"), "cell speed"),
                            math.hypot(u_value, v_value), abs_tol=1e-10, rel_tol=1e-9):
            fail(f"{case.get('label')}: cell speed does not match velocity magnitude")
        u[y_index[y], x_index[x]] = u_value
        v[y_index[y], x_index[x]] = v_value
    if len(seen) != len(rows) or not np.isfinite(u).all() or not np.isfinite(v).all():
        fail(f"{case.get('label')}: rectangular cell arrays are incomplete or non-finite")
    actual_cells_hash = sha256_file(cells_path)
    declared_cells_hash = None
    for raw_path, declared_hash in (case.get("artifactSha256", {}) or {}).items():
        declared_path = Path(raw_path)
        if not declared_path.is_absolute():
            declared_path = summary_path.parent / declared_path
        if declared_path.resolve() == cells_path.resolve():
            declared_cells_hash = declared_hash
            break
    if declared_cells_hash is not None and actual_cells_hash != declared_cells_hash:
        fail(f"{case.get('label')}: cells CSV hash differs from verifier artifact metadata")

    x_axis, y_axis = np.asarray(xs), np.asarray(ys)
    # First interpolate cell centres to x=.5/y=.5, then add the true wall
    # values before the final 101-point centre-line interpolation.
    u_at_centre = np.asarray([np.interp(0.5, x_axis, row) for row in u])
    v_at_centre = np.asarray([np.interp(0.5, y_axis, column) for column in v.T])
    u_profile = np.interp(SAMPLES, np.r_[0.0, y_axis, 1.0],
                          np.r_[0.0, u_at_centre, 1.0])
    v_profile = np.interp(SAMPLES, np.r_[0.0, x_axis, 1.0],
                          np.r_[0.0, v_at_centre, 0.0])
    if not np.isfinite(u_profile).all() or not np.isfinite(v_profile).all():
        fail(f"{case.get('label')}: interpolated profile is non-finite")
    return {
        "label": case.get("label"),
        "level": LEVELS[expected_cells],
        "cells": len(rows),
        "prefix": str(prefix),
        "meshSha256": case.get("meshSha256"),
        "cellsCsvSha256": actual_cells_hash,
        "declaredCellsCsvSha256": declared_cells_hash,
        "integrityScope": ("cells CSV hash matched verifier artifact metadata"
                            if declared_cells_hash is not None
                            else "cells CSV hash computed; verifier artifact metadata did not declare it"),
        "artifactSha256": copy.deepcopy(case.get("artifactSha256", {})),
        "cellGrid": {"nx": len(xs), "ny": len(ys), "coordinateRoundingDigits": 12},
        "profile": {"u": u_profile.tolist(), "v": v_profile.tolist()},
        "metadata": {
            "native": copy.deepcopy(native),
            "stage": copy.deepcopy(case.get("stage")),
            "command": copy.deepcopy(case.get("stage", {}).get("command")),
        },
    }


def compare_case(case: dict[str, Any], summary_path: Path,
                 fd_profiles: dict[int, dict[str, list[float]]]) -> dict[str, Any]:
    profile = read_cell_profile(case, summary_path)
    fvm_u = np.asarray(profile["profile"]["u"])
    fvm_v = np.asarray(profile["profile"]["v"])
    errors: dict[str, Any] = {}
    for n, reference in fd_profiles.items():
        ref_u, ref_v = np.asarray(reference["u"]), np.asarray(reference["v"])
        errors[f"FD{n}"] = {
            "uRms": float(np.sqrt(np.mean((fvm_u - ref_u) ** 2))),
            "vRms": float(np.sqrt(np.mean((fvm_v - ref_v) ** 2))),
            "combinedRms": float(np.sqrt(np.mean(np.r_[fvm_u - ref_u, fvm_v - ref_v] ** 2))),
        }
    profile["errorsVsFD"] = errors
    return profile


def infer_scheme(summary_path: Path, summary: dict[str, Any], cases: list[dict[str, Any]]) -> str:
    schemes = {case.get("native", {}).get("convection") for case in cases}
    schemes.discard(None)
    if len(schemes) != 1 or next(iter(schemes)) not in {"upwind", "limited-linear"}:
        fail(f"{summary_path}: cavity cases do not identify one supported scheme")
    return str(next(iter(schemes)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", type=Path, required=True,
                        help="FD directory containing summary.json and N33/N65/N129 u/v arrays")
    parser.add_argument("--flow-summary", type=Path, action="append", required=True,
                        help="formal independent verifier summary; repeat per scheme")
    parser.add_argument("--output", type=Path, required=True, help="diagnostic JSON output path")
    args = parser.parse_args()

    reference_root = args.reference_root.resolve()
    reference = load_fd_reference(reference_root)
    source_summaries: list[dict[str, Any]] = []
    comparison_cases: list[dict[str, Any]] = []
    for raw_path in args.flow_summary:
        summary_path = raw_path.resolve()
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"cannot read flow summary {summary_path}: {exc}") from exc
        cases = summary.get("cases")
        if not isinstance(cases, list):
            fail(f"{summary_path}: flow summary cases must be a list")
        cavity_cases = [case for case in cases
                        if isinstance(case, dict) and case.get("case") == "cavity"
                        and case.get("valid") is True]
        if not cavity_cases:
            fail(f"{summary_path}: no valid cavity cases to compare")
        scheme = infer_scheme(summary_path, summary, cavity_cases)
        source = {
            "path": str(summary_path),
            "sha256": sha256_file(summary_path),
            "format": summary.get("format"),
            "valid": summary.get("valid"),
            "issues": copy.deepcopy(summary.get("issues", [])),
            # Preserve this exactly, including a failing sequence gate.
            "originalSequenceChecks": copy.deepcopy(summary.get("sequenceChecks")),
            "executables": copy.deepcopy(summary.get("executables")),
            "scheme": scheme,
            "selectedLabels": [case.get("label") for case in cavity_cases],
        }
        source_summaries.append(source)
        for case in cavity_cases:
            selected = compare_case(case, summary_path, reference["profiles"])
            selected["scheme"] = scheme
            selected["sourceSummary"] = str(summary_path)
            selected["originalSequenceChecks"] = copy.deepcopy(summary.get("sequenceChecks"))
            comparison_cases.append(selected)

    comparison_cases.sort(key=lambda case: (case["scheme"], case["level"]))
    output = {
        "format": "cartmesh2d-cavity-centreline-diagnostic-v1",
        "diagnosticOnly": True,
        "scriptSha256": sha256_file(Path(__file__).resolve()),
        "scope": "Existing cavity cells CSV compared with independent finite-grid references; no physical acceptance gate",
        "sampleCoordinates": SAMPLES.tolist(),
        "limitations": ["unit-square cavity only", "speed=1 only", "101 equally spaced centre-line samples"],
        "reference": reference,
        "sourceSummaries": source_summaries,
        "cases": comparison_cases,
    }
    output_path = args.output.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({"output": str(output_path), "cases": len(comparison_cases),
                      "fdDifferenceRms": reference["fdDifferenceRms"],
                      "errors": [(case["scheme"], case["level"], case["errorsVsFD"]["FD129"]["combinedRms"])
                                 for case in comparison_cases]},
                     indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
