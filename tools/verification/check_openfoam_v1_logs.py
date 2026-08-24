#!/usr/bin/env python3
"""Fail closed on OpenFOAM checkMesh and simpleFoam acceptance logs."""

import argparse
import json
import re
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", type=Path, action="append", required=True)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    issues = []
    checks = []
    for path in args.check:
        text = path.read_text(encoding="utf-8")
        valid = "Mesh OK." in text and "Failed " not in text and "FOAM FATAL" not in text
        if not valid:
            issues.append(f"{path.name}: checkMesh did not report an unqualified Mesh OK")
        checks.append({"path": str(path), "valid": valid})
    solver_text = args.solver.read_text(encoding="utf-8")
    converged = re.search(r"SIMPLE solution converged in (\d+) iterations", solver_text)
    continuity = re.findall(
        r"time step continuity errors\s*:\s*sum local = ([^,]+), global = ([^,]+), cumulative = ([^\n]+)",
        solver_text,
    )
    if not converged or not continuity or "FOAM FATAL" in solver_text or "[stack trace]" in solver_text:
        issues.append("simpleFoam log lacks clean convergence evidence")
    local, global_error, cumulative = map(float, continuity[-1]) if continuity else (float("inf"),) * 3
    if abs(local) > 1.0e-8 or abs(global_error) > 1.0e-8:
        issues.append("final local/global continuity errors exceed 1e-8")
    report = {
        "valid": not issues,
        "checks": checks,
        "solver": {
            "path": str(args.solver),
            "converged": bool(converged),
            "iterations": int(converged.group(1)) if converged else None,
            "final_continuity_local": local,
            "final_continuity_global": global_error,
            "final_continuity_cumulative": cumulative,
        },
        "issues": issues,
    }
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    raise SystemExit(0 if report["valid"] else 1)


if __name__ == "__main__":
    main()
