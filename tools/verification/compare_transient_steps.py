#!/usr/bin/env python3
"""Compare completed fixed-mesh transient runs at strictly halved time steps.

This is a consistency report.  A finite self-convergence order is
reported when available, but no order or accuracy threshold is an acceptance
criterion.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_transient_flow as audit  # noqa: E402


def fail(message: str) -> None:
    raise audit.native.VerificationError(message)


def _command_value(command: list[str], option: str) -> str:
    if command.count(option) != 1:
        fail(f"run command must contain exactly one {option}")
    try:
        return command[command.index(option) + 1]
    except (ValueError, IndexError) as exc:
        fail(f"run command is missing {option}")
        raise AssertionError from exc


def _prefix(command: list[str]) -> Path:
    return Path(_command_value(command, "--output"))


def _optional_command_value(command: list[str], option: str) -> str | None:
    try:
        return command[command.index(option) + 1]
    except (ValueError, IndexError):
        return None


def _read_cells(prefix: Path) -> dict[int, tuple[float, float, float]]:
    path = Path(str(prefix) + ".cells.csv")
    if not path.is_file():
        fail(f"missing transient cell field: {path}")
    try:
        with path.open(newline="", encoding="utf-8-sig") as stream:
            rows = list(csv.DictReader(stream))
    except (OSError, csv.Error) as exc:
        fail(f"cannot read {path}: {exc}")
    required = {"cell", "area", "u", "v"}
    if not rows or not required.issubset(rows[0]):
        fail(f"{path}: missing cell/area/u/v columns")
    values: dict[int, tuple[float, float, float]] = {}
    for line, row in enumerate(rows, 2):
        try:
            cell = audit.native.integer(row["cell"], f"{path}:{line} cell")
            area = audit.native.finite(row["area"], f"{path}:{line} area")
            u = audit.native.finite(row["u"], f"{path}:{line} u")
            v = audit.native.finite(row["v"], f"{path}:{line} v")
        except (KeyError, audit.native.VerificationError):
            raise
        if cell in values:
            fail(f"{path}:{line}: duplicate cell id {cell}")
        if area <= 0:
            fail(f"{path}:{line}: cell area must be positive")
        values[cell] = (area, u, v)
    return values


def _same(a: Any, b: Any, label: str) -> None:
    if isinstance(a, (int, float)) or isinstance(b, (int, float)):
        try:
            if not audit.native.close(float(a), float(b), 1e-14, 1e-12):
                fail(f"runset {label} differs: {a!r} vs {b!r}")
            return
        except (TypeError, ValueError):
            pass
    if a != b:
        fail(f"runset {label} differs: {a!r} vs {b!r}")


def compare(runset: Path, output: Path) -> dict[str, Any]:
    try:
        data = json.loads(runset.read_text(encoding="utf-8-sig"), parse_constant=lambda x: fail(f"non-finite JSON token {x}"))
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"cannot read runset {runset}: {exc}")
    if not isinstance(data, dict) or not isinstance(data.get("runs"), list):
        fail("runset must be an object containing runs")
    entries = data["runs"]
    if len(entries) < 3:
        fail("runset requires at least three completed runs")
    mesh = Path(str(data.get("mesh", ""))).expanduser()
    if not mesh.is_file():
        fail(f"runset mesh is missing: {mesh}")
    mesh_hash = audit.native.sha256_file(mesh)
    if data.get("meshSha256") != mesh_hash:
        fail("runset mesh hash does not match current source mesh")
    expected_binary_hash = data.get("binarySha256")
    if not isinstance(expected_binary_hash, str):
        fail("runset is missing binarySha256")

    records: list[dict[str, Any]] = []
    controls: dict[str, Any] | None = None
    expected_end: float | None = None
    mesh_ids: set[int] | None = None
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict) or entry.get("returnCode") != 0 or entry.get("valid") is not True:
            fail(f"run {index} is not a completed valid run")
        command = entry.get("command")
        if not isinstance(command, list) or not command or not all(isinstance(x, str) for x in command):
            fail(f"run {index} has no usable command")
        if "--restart" in command:
            fail(f"run {index} uses --restart; comparison requires fresh runs")
        binary = Path(command[0]).expanduser()
        if not binary.is_file():
            fail(f"run {index} binary is missing: {binary}")
        binary_hash = audit.native.sha256_file(binary)
        if binary_hash != expected_binary_hash:
            fail(f"run {index} binary hash differs from runset")
        command_mesh = Path(_command_value(command, "--mesh")).expanduser()
        if command_mesh.resolve() != mesh.resolve():
            fail(f"run {index} command mesh differs from runset mesh")
        dt = audit.native.finite(entry.get("dt"), f"run {index} dt")
        requested_end = audit.native.finite(entry.get("endTime"), f"run {index} endTime")
        if dt <= 0 or requested_end <= 0:
            fail(f"run {index} dt/endTime must be positive")
        command_dt = audit.native.finite(_command_value(command, "--time-step"), f"run {index} command dt")
        command_steps = audit.native.integer(_command_value(command, "--steps"), f"run {index} command steps")
        if command_steps < 1 or not audit.native.close(command_dt, dt, 1e-14, 1e-12):
            fail(f"run {index} entry dt differs from command --time-step")
        if not audit.native.close(command_steps * command_dt, requested_end, 1e-11, 1e-10):
            fail(f"run {index} command steps*dt does not equal endTime")
        prefix = _prefix(command)
        audit_path = output.parent / (output.stem + f".run{index}.audit.json")
        result = audit.verify(mesh, prefix, audit_path)
        if result.get("valid") is not True:
            fail(f"run {index} independent verification did not produce valid=true")
        history = result.get('history')
        if not isinstance(history, list) or len(history) != command_steps or not history:
            fail(f"run {index} accepted history length differs from command steps")
        if not audit.native.close(audit.native.finite(history[0].get('time'), 'first time'), dt, 1e-12, 1e-10):
            fail(f"run {index} history does not start at physical time zero")
        run_controls = dict(result.get("controls") or {})
        run_controls["case"] = result.get("case")
        required_controls = ("case", "nu", "speed", "convection", "viscousStress",
                             "temporalFaceInterpolation", "velocityRelaxation", "tolerance")
        if any(key not in run_controls for key in required_controls):
            fail(f"run {index} independent verifier omitted required controls")
        command_controls = {
            "case": _command_value(command, "--case"),
            "nu": audit.native.finite(_command_value(command, "--nu"), f"run {index} command nu"),
            "speed": audit.native.finite(_command_value(command, "--speed"), f"run {index} command speed"),
            "convection": _command_value(command, "--convection"),
            "tolerance": audit.native.finite(_command_value(command, "--tolerance"), f"run {index} command tolerance"),
        }
        for key, value in command_controls.items():
            _same(value, run_controls.get(key), f"run {index} command {key}")
        if _optional_command_value(command, "--viscous-stress") is not None:
            _same(_optional_command_value(command, "--viscous-stress"), run_controls.get("viscousStress"),
                  f"run {index} command viscousStress")
        if not audit.native.close(dt, audit.native.finite(result.get("dt"), f"run {index} verified dt"), 1e-14, 1e-12):
            fail(f"run {index} entry dt differs from independently verified dt")
        if controls is None:
            controls = dict(run_controls)
            expected_end = requested_end
        else:
            for key in required_controls:
                _same(controls.get(key), run_controls.get(key), key)
            _same(expected_end, requested_end, "endTime")
        accepted_time = audit.native.finite(result.get("time"), f"run {index} accepted end time")
        if not audit.native.close(accepted_time, requested_end, 1e-11, 1e-9):
            fail(f"run {index} accepted end time differs from requested endTime")
        fields = _read_cells(prefix)
        ids = set(fields)
        if mesh_ids is None:
            mesh_ids = ids
        elif ids != mesh_ids:
            fail(f"run {index} cell id coverage differs from other runs")
        records.append({"index": index, "dt": dt, "prefix": prefix, "acceptedTime": accepted_time,
                        "binarySha256": binary_hash, "audit": str(audit_path),
                        "auditSha256": audit.native.sha256_file(audit_path), "cells": fields})
    records.sort(key=lambda item: item["dt"], reverse=True)
    if controls is None or expected_end is None or mesh_ids is None:
        fail("runset contains no comparable runs")
    for coarse, fine in zip(records, records[1:]):
        if not audit.native.close(fine["dt"], coarse["dt"] / 2.0, 1e-14, 1e-12):
            fail("time steps are not strictly halved")
    distances = []
    for coarse, fine in zip(records, records[1:]):
        area_sum = 0.0
        du2 = dv2 = 0.0
        for cell in mesh_ids:
            area, u0, v0 = coarse["cells"][cell]
            area1, u1, v1 = fine["cells"][cell]
            if not audit.native.close(area, area1, 1e-14, 1e-12):
                fail(f"cell {cell} area differs between runs")
            area_sum += area
            du2 += area * (u1 - u0) ** 2
            dv2 += area * (v1 - v0) ** 2
        if not math.isfinite(area_sum) or area_sum <= 0:
            fail("non-positive total comparison area")
        if not math.isfinite(du2) or not math.isfinite(dv2) or du2 < 0 or dv2 < 0:
            fail("non-finite or invalid area-weighted velocity distance")
        u_l2, v_l2 = math.sqrt(du2 / area_sum), math.sqrt(dv2 / area_sum)
        vector_l2 = math.hypot(u_l2, v_l2)
        if not math.isfinite(vector_l2):
            fail("velocity distance overflowed")
        distances.append({"coarseDt": coarse["dt"], "fineDt": fine["dt"], "uL2": u_l2,
                          "vL2": v_l2, "vectorL2": vector_l2})
    for i, item in enumerate(distances[1:]):
        previous = distances[i]["vectorL2"]
        current = item["vectorL2"]
        order = math.log(previous / current, 2.0) if previous > 0 and current > 0 and math.isfinite(previous) and math.isfinite(current) else None
        item["observedOrder"] = order if order is None or math.isfinite(order) else None
    result = {"valid": True, "scope": "fixed-mesh transient artifact consistency; self-convergence is diagnostic, not an accuracy qualification",
              "runsetSha256": audit.native.sha256_file(runset),
              "mesh": str(mesh.resolve()), "meshSha256": mesh_hash, "binarySha256": expected_binary_hash,
              "controls": controls, "acceptedEndTime": expected_end, "cellIds": len(mesh_ids),
              "runs": [{"dt": r["dt"], "acceptedTime": r["acceptedTime"], "binarySha256": r["binarySha256"], "audit": r["audit"], "auditSha256": r["auditSha256"]} for r in records],
              "adjacentDistances": distances}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = compare(args.runs, args.output)
    except (audit.native.VerificationError, OSError, ValueError, json.JSONDecodeError) as exc:
        failure = {"valid": False, "issues": [str(exc)]}
        audit.native.write_json(args.output, failure)
        print(json.dumps(failure, indent=2))
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
