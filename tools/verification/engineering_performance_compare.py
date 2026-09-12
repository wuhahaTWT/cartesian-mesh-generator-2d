#!/usr/bin/env python3
"""Reproducible bounded before/after performance pair runner.

The output directory must be empty or absent.  Native timing comes from the
CLI's timing_total_seconds; /usr/bin/time's real duration is reported separately.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import signal
import subprocess
from pathlib import Path

from engineering_scale_survey import parse_native_seconds, parse_rss, measure, check_openfoam

REPO = Path(__file__).resolve().parents[2]
CASES = (("circle-r02", "examples/acceptance/circle.xy", "2"),
         ("naca-r02", "examples/complex/naca2412_dense.xy", "1"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--phase", choices=("before", "after"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=180.0)
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("timeout must be finite and positive")
    args.cli = args.cli.resolve(strict=True)
    args.output_dir = args.output_dir.resolve()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise SystemExit(f"refusing to overwrite non-empty output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for name, input_path, reference in CASES:
        destination = args.output_dir / name
        destination.mkdir(parents=True, exist_ok=True)
        prefix, case = destination / "mesh", destination / "openfoam"
        command = ["/usr/bin/time", "-l", str(args.cli), str(REPO / input_path), str(prefix),
                   "11", "0.25", "0.1", "exterior", str(case), "0", "--size-field",
                   "--reference-length", reference, "--wall-relative-size", "0.005",
                   "--background-relative-size", "0.0125", "--far-field-spans", "0.5",
                   "--cells-per-level", "3", "--max-safe-wall-level", "11"]
        (destination / "command.json").write_text(json.dumps(command, indent=2) + "\n")
        process = subprocess.Popen(command, cwd=REPO, text=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            output, _ = process.communicate(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            output, _ = process.communicate()
            process.returncode = 124
        (destination / "command.log").write_text(output or "", encoding="utf-8")
        elapsed = re.search(r"([0-9.]+)\s+real\s+", output or "")
        result = {"case": name, "phase": args.phase, "command": command,
                  "wrapper_returncode": process.returncode,
                  "native_seconds": parse_native_seconds(output or ""),
                  "wrapper_elapsed_seconds": float(elapsed.group(1)) if elapsed else None,
                  "peak_rss_bytes": parse_rss(output or ""),
                  "native_pass": "cartmesh2d end-to-end PASS" in (output or ""),
                  "timeout_seconds": args.timeout, "issues": []}
        for label, pattern in (("solver", r"^solver_quality=(\w+)"),):
            match = re.search(pattern, output or "", re.M)
            result[label] = match.group(1) if match else "not_reported"
        result["external_checkMesh"] = "not_run"
        cm2d, report = Path(str(prefix) + ".solver.cm2d"), Path(str(prefix) + ".resolution.json")
        if result["native_pass"]:
            try:
                result["cm2d_sha256"] = hashlib.sha256(cm2d.read_bytes()).hexdigest()
                result["resolution_independent"] = measure(cm2d, report)
                result["openfoam_independent"] = check_openfoam(case)
                result["cell_count"] = result["openfoam_independent"]["cell_count"]
            except (OSError, ValueError, KeyError, TypeError) as exc:
                result["issues"].append(str(exc))
        result["measurement_complete"] = (process.returncode == 0 and result["native_pass"]
            and all(result[k] is not None for k in
                    ("native_seconds", "wrapper_elapsed_seconds", "peak_rss_bytes"))
            and result.get("resolution_independent", {}).get("valid", False)
            and result.get("openfoam_independent", {}).get("valid", False)
            and not result["issues"])
        (destination / "measurement.json").write_text(json.dumps(result, indent=2) + "\n")
        results.append(result)
        print(f"{args.phase}/{name}: measurement_complete={result['measurement_complete']}")
    summary = {"cli": str(args.cli), "cli_sha256": hashlib.sha256(args.cli.read_bytes()).hexdigest(),
               "phase": args.phase, "cases": results,
               "notes": ["Measurement completeness is separate from external solver acceptance.",
                         "Compare the two phase summaries by case and CM2D SHA256; do not mix native and wrapper times."]}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0 if all(r["measurement_complete"] for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
