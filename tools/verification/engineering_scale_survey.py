#!/usr/bin/env python3
"""Bounded engineering survey of fixed dimensionless size requests.

This is a measurement harness, not a success preset.  Every request is run
exactly once with a 120 second limit and failures are retained verbatim.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import signal
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from tools.verification.check_mesh_resolution import measure  # noqa: E402
from tools.verification.check_openfoam2d import check as check_openfoam  # noqa: E402


REQUESTS = (("r01", 0.01, 0.025, 3),
            ("r02", 0.005, 0.0125, 3),
            ("r03", 0.0025, 0.0125, 40))
CASES = (
    ("circle", REPO / "examples/acceptance/circle.xy", "exterior", 2.0),
    ("nozzle", REPO / "examples/complex/nozzle_profile.xy", "interior", None),
    ("naca", REPO / "examples/complex/naca2412_dense.xy", "exterior", 1.0),
)


def bbox_span(path: Path) -> float:
    points = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 2:
            points.append((float(fields[0]), float(fields[1])))
    if not points:
        raise ValueError(f"{path}: no points")
    return max(max(x for x, _ in points) - min(x for x, _ in points),
               max(y for _, y in points) - min(y for _, y in points))


def parse_rss(text: str) -> int | None:
    match = re.search(r"(?:maximum resident set size:\s*(\d+)|(\d+)\s+maximum resident set size)",
                      text, re.I)
    if match:
        return int(match.group(1) or match.group(2))
    return None


def parse_native_seconds(text: str) -> float | None:
    match = re.search(r"timing_total_seconds=([0-9.eE+-]+)", text)
    return float(match.group(1)) if match else None


def run_one(root: Path, cli: Path, case_name: str, boundary: Path,
            fluid_region: str, reference: float, request: tuple[str, float, float, int]) -> dict:
    request_name, wall, background, band = request
    output = root / case_name / request_name
    output.mkdir(parents=True, exist_ok=True)
    prefix = output / "mesh"
    foam_case = output / "openfoam"
    command = ["/usr/bin/time", "-l", str(cli), str(boundary), str(prefix), "11", "0.25", "0.1",
               fluid_region, str(foam_case), "0", "--size-field", "--reference-length",
               f"{reference:.17g}", "--wall-relative-size", f"{wall:.17g}",
               "--background-relative-size", f"{background:.17g}", "--far-field-spans", "0.5",
               "--cells-per-level", str(band), "--max-safe-wall-level", "11"]
    process = subprocess.Popen(command, cwd=REPO, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, start_new_session=True)
    try:
        output_text, _ = process.communicate(timeout=120)
        completed = subprocess.CompletedProcess(command, process.returncode, output_text or "")
        status = "success" if completed.returncode == 0 else "failed"
        failure = None if status == "success" else f"wrapper_returncode={completed.returncode}"
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        output_text, _ = process.communicate()
        completed = subprocess.CompletedProcess(command, 124, output_text or "")
        status, failure = "timeout", "timeout after 120 seconds"
    log_path = output / "command.log"
    log_path.write_text(completed.stdout or "", encoding="utf-8")
    native_completed = "cartmesh2d end-to-end PASS" in (completed.stdout or "")
    result = {"case": case_name, "request": {"name": request_name, "wall_h_over_lref": wall,
              "background_h_over_lref": background, "cells_per_level": band},
              "reference_length": reference, "fluid_region": fluid_region, "status": status,
              "failure": failure, "command": command, "log": str(log_path),
              "wrapper_returncode": completed.returncode,
              "native_completed_marker": native_completed,
              "peak_rss_bytes": parse_rss(completed.stdout or ""),
              "native_timing_seconds": parse_native_seconds(completed.stdout or "")}
    result["missing_metrics"] = []
    if result["peak_rss_bytes"] is None:
        result["missing_metrics"].append("peak_rss_bytes")
    if result["native_timing_seconds"] is None:
        result["missing_metrics"].append("native_timing_seconds")
    if status == "success" and result["missing_metrics"]:
        result["status"] = "measurement_incomplete"
        result["failure"] = "wrapper completed but required measurement metrics are missing"
    if status == "failed" and native_completed:
        result["status"] = "measurement_incomplete"
        result["failure"] = "native completion marker present but /usr/bin/time wrapper returned nonzero"
    report_path = Path(str(prefix) + ".resolution.json")
    cm2d_path = Path(str(prefix) + ".solver.cm2d")
    if result["status"] in ("success", "measurement_incomplete"):
        if not report_path.exists() or not cm2d_path.exists() or not foam_case.exists():
            result["status"] = "measurement_incomplete"
            result["failure"] = "success exit but required solver/report/OpenFOAM output is missing"
        else:
            try:
                resolution = measure(cm2d_path, report_path)
                openfoam = check_openfoam(foam_case)
                result["resolution_independent"] = resolution
                result["openfoam_independent"] = openfoam
                result["final_cell_count"] = resolution["measured"]["cell_count"]
                result["requested_wall_actual"] = resolution["measured"]["actual"].get(
                    "wall_owner_tangential_extent_over_reference")
                if not resolution["valid"] or not openfoam["valid"]:
                    result["status"] = "measurement_incomplete"
                    result["failure"] = "independent CM2D or OpenFOAM reader failed"
            except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
                result["status"] = "measurement_incomplete"
                result["failure"] = f"independent reader exception: {exc}"
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", type=Path, default=REPO / "build/cartmesh2d_cli")
    parser.add_argument("--output-dir", type=Path, default=REPO / "outputs/engineering-scale")
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise SystemExit(f"refusing to overwrite non-empty output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for case_name, boundary, fluid_region, configured_reference in CASES:
        reference = configured_reference if configured_reference is not None else bbox_span(boundary)
        for request in REQUESTS:
            results.append(run_one(args.output_dir, args.cli, case_name, boundary,
                                   fluid_region, reference, request))
    summary = {"format": "cartmesh2d-engineering-scale-v1", "timeout_seconds": 120,
               "cases": results,
               "success_count": sum(result["status"] == "success" for result in results),
               "failure_count": sum(result["status"] not in ("success", "measurement_incomplete")
                                     for result in results),
               "measurement_incomplete_count": sum(result["status"] == "measurement_incomplete"
                                                    for result in results),
               "missing_metrics": sorted({metric for result in results
                                           for metric in result.get("missing_metrics", [])}),
               "first_failure": next((result for result in results if result["status"] != "success"), None)}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                                                    encoding="utf-8")
    print(json.dumps({"success_count": summary["success_count"],
                      "failure_count": summary["failure_count"],
                      "first_failure": summary["first_failure"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
