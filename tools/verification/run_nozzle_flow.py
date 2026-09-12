#!/usr/bin/env python3
"""Generate and solve the fixed three-grid nozzle fixture in a fresh directory.

Flow-pipeline completion and expanded checkMesh are separate results.
The driver never overwrites an earlier run or leaves a timed-out container.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import signal
import subprocess
import sys
import time
import uuid
from pathlib import Path
from types import SimpleNamespace

import openfoam_nozzle_flow as nozzle
from check_openfoam2d import check

ROOT = Path(__file__).resolve().parents[2]
GRIDS = (("r01", ".01", ".025", "3"),
         ("r02", ".005", ".0125", "3"),
         ("r03", ".0025", ".0125", "40"))
IMAGE = "opencfd/openfoam-run:2606"


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def run_stage(directory, command, timeout, container=None):
    directory.mkdir(parents=True, exist_ok=False)
    write_json(directory / "command.json", command)
    started = time.monotonic()
    result = {"command": command, "timeout_seconds": timeout,
              "status": "not_started", "returncode": None}
    with (directory / "stdout.log").open("w") as log:
        try:
            process = subprocess.Popen(command, cwd=ROOT, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            try:
                process.wait(timeout=timeout)
                result.update(status="passed" if process.returncode == 0 else "failed",
                              returncode=process.returncode)
            except subprocess.TimeoutExpired:
                result.update(status="timeout", returncode=124)
                if container:
                    try:
                        cleanup = subprocess.run(["docker", "rm", "-f", container],
                            capture_output=True, text=True, timeout=15)
                        result["container_cleanup"] = {"returncode": cleanup.returncode,
                            "stdout": cleanup.stdout, "stderr": cleanup.stderr}
                    except (OSError, subprocess.TimeoutExpired) as exc:
                        result["container_cleanup"] = {"error": str(exc)}
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
        except OSError as exc:
            result.update(status="failed", error=str(exc))
    result["elapsed_seconds"] = time.monotonic() - started
    write_json(directory / "stage.json", result)
    return result


def docker_stage(directory, case, program, extra=()):
    name = "cartmesh-nozzle-" + uuid.uuid4().hex
    command = ["docker", "run", "--name", name, "--rm", "--network", "none", "-v",
        f"{case.resolve()}:/home/openfoam/workingDir/case", IMAGE, program,
        "-case", "/home/openfoam/workingDir/case", *extra]
    result = run_stage(directory, command, 600 if program == "simpleFoam" else 120, name)
    log = (directory / "stdout.log").read_text(errors="replace")
    if program == "checkMesh":
        result["mesh_ok"] = (result["status"] == "passed"
                             and "Mesh OK." in log and "Failed " not in log)
        match = re.search(r"Writing (\d+) concave cells", log)
        result["concave_cells"] = int(match.group(1)) if match else None
    else:
        result.update(nozzle.solver_convergence(log))
    write_json(directory / "stage.json", result)
    return result



def run_case(root, cli, grid):
    name, wall, background, band = grid
    generated = root / "generated" / name
    generated.mkdir(parents=True)
    prefix, source = generated / "nozzle", generated / "openfoam"
    stages_root = root / "stages" / name
    command = [str(cli), str(nozzle.INPUT), str(prefix), "11", ".25", ".1", "interior",
        str(source), "0", "--size-field", "--reference-length", "6",
        "--wall-relative-size", wall, "--background-relative-size", background,
        "--far-field-spans", ".5", "--cells-per-level", band, "--max-safe-wall-level", "11"]
    item = {"name": name, "request": {"wall": wall, "background": background, "band": band},
            "stages": {}, "flow_pipeline_valid": False}
    generated_stage = run_stage(stages_root / "generate", command, 180)
    item["stages"]["generate"] = generated_stage
    if generated_stage["status"] != "passed":
        return item
    log = (stages_root / "generate/stdout.log").read_text()
    for key, expression in (("solver_quality", r"^solver_quality=(\w+)"),):
        match = re.search(expression, log, re.M)
        item[key] = match.group(1) if match else "not_reported"
    config_report = stages_root / "configuration.json"
    nozzle.configure(SimpleNamespace(output_root=root / "cases", case=[(name, source)], report=config_report))
    case = root / "cases" / name
    item["case"] = str(case)
    independent = check(case)
    write_json(stages_root / "independent.json", independent)
    item["independent_valid"] = independent["valid"]
    if not independent["valid"]:
        return item
    standard = docker_stage(stages_root / "checkMesh-standard", case, "checkMesh")
    expanded = docker_stage(stages_root / "checkMesh-expanded", case, "checkMesh", ("-allGeometry", "-allTopology"))
    item["stages"].update(checkMesh_standard=standard, checkMesh_expanded=expanded)
    if not standard["mesh_ok"]:
        return item
    solve = docker_stage(stages_root / "solve", case, "simpleFoam")
    item["stages"]["solve"] = solve
    if solve["status"] == "passed":
        item["evaluation"] = nozzle.evaluate_case(case)
        item["flow_pipeline_valid"] = (item["evaluation"]["valid"]
                                          and solve["converged_at_requested_residual"])
    return item


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--cli", type=Path, default=ROOT / "build/cartmesh2d_cli")
    parser.add_argument("--cases", nargs="+", choices=[g[0] for g in GRIDS],
                        default=[g[0] for g in GRIDS], help="explicit subset for a small pipeline check")
    args = parser.parse_args()
    root, cli = args.output_root.resolve(), args.cli.resolve(strict=True)
    if root.exists() and any(root.iterdir()):
        parser.error(f"refusing non-empty run directory: {root}")
    nozzle.source_points()  # Verify the fixed source before launching any process.
    root.mkdir(parents=True, exist_ok=True)
    cases = []
    for grid in GRIDS:
        if grid[0] not in args.cases:
            continue
        try:
            item = run_case(root, cli, grid)
        except (OSError, ValueError, KeyError, IndexError, TypeError) as exc:
            item = {"name": grid[0], "flow_pipeline_valid": False, "error": str(exc)}
        cases.append(item)
        write_json(root / "progress.json", cases)
    result = {"cases": cases, "grid_scope": [c["name"] for c in cases],
        "flow_pipeline_valid": all(c["flow_pipeline_valid"] for c in cases),
        "expanded_checks_passed": all(c.get("stages", {}).get("checkMesh_expanded", {}).get("mesh_ok", False) for c in cases),
        "notes": ["Pipeline completion is not general CFD accuracy acceptance.",
                  "Pressure is kinematic (m2/s2), phi is volumetric (m3/s)."]}
    if all(c.get("evaluation", {}).get("valid", False) for c in cases):
        nozzle.evaluate(SimpleNamespace(case=[Path(c["case"]) for c in cases], report=root / "comparison.json",
            solve_logs=[root / "stages" / c["name"] / "solve/stdout.log" for c in cases]))
        result["comparison"] = str(root / "comparison.json")
    write_json(root / "run-summary.json", result)
    print(json.dumps({k:v for k,v in result.items() if k != "cases"}, indent=2))
    return 0 if result["flow_pipeline_valid"] else 1


if __name__ == "__main__":
    sys.exit(main())
