#!/usr/bin/env python3
"""Fail-closed acceptance gate for the V1c OpenFOAM verification products."""

import argparse
import json
import math
import re
from pathlib import Path


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def last_residual(path: Path) -> float:
    matches = re.findall(r"Solving for T, Initial residual = [^,]+, Final residual = ([^,]+)",
                         path.read_text(encoding="utf-8"))
    if not matches:
        return math.inf
    return float(matches[-1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mms", type=Path, action="append", required=True)
    parser.add_argument("--constant", type=Path, required=True)
    parser.add_argument("--linear", type=Path, required=True)
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--check", type=Path, action="append", required=True)
    parser.add_argument("--solver-log", type=Path, action="append", required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    issues = []
    mms = [load(path) for path in args.mms]
    if len(mms) < 3:
        issues.append("at least three manufactured-solution grids are required")
    for report in mms:
        if not report.get("valid") or report.get("linear_system_final_residual", math.inf) > 1e-12:
            issues.append("manufactured solution or linear-system residual failed")
    orders = {}
    for norm in ("l1", "l2", "linf"):
        values = [report[norm] for report in mms]
        orders[norm] = [math.log(values[i] / values[i + 1], 2.0)
                        for i in range(len(values) - 1)]
        if any(values[i + 1] >= values[i] for i in range(len(values) - 1)):
            issues.append(f"{norm} is not strictly decreasing")
        if any(order < 0.9 for order in orders[norm]):
            issues.append(f"{norm} observed order falls below 0.9")
    constant = load(args.constant)
    linear = load(args.linear)
    if not constant.get("valid") or constant.get("l2", math.inf) > 1e-12:
        issues.append("constant-field free-stream preservation exceeds 1e-12 L2")
    if not linear.get("valid") or linear.get("l2", math.inf) > 5e-4:
        issues.append("linear-field preservation exceeds 5e-4 L2")
    geometry = load(args.geometry)
    if (not geometry.get("valid") or
            abs(geometry.get("global_constant_flux_balance", math.inf)) > 1e-12 or
            geometry.get("max_closure_residual", math.inf) > 1e-12):
        issues.append("independent geometric conservation gate failed")
    checks = []
    for path in args.check:
        text = path.read_text(encoding="utf-8")
        valid = "Mesh OK." in text and "Failed " not in text and "FOAM FATAL" not in text
        checks.append({"path": str(path), "valid": valid})
        if not valid:
            issues.append(f"{path}: checkMesh lacks unqualified Mesh OK")
    residuals = [last_residual(path) for path in args.solver_log]
    if any(value > 1e-12 for value in residuals):
        issues.append("a solver log final residual exceeds 1e-12")
    report = {"valid": not issues, "mms": mms, "observed_orders": orders,
              "constant": constant, "linear": linear, "geometry": geometry,
              "check_mesh": checks, "solver_final_residuals": residuals,
              "issues": issues}
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    raise SystemExit(0 if report["valid"] else 1)


if __name__ == "__main__":
    main()
