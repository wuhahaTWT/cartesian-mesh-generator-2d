#!/usr/bin/env python3
"""Reproduce the fixed three-grid OpenFOAM manufactured-solution run.

This driver deliberately has no parameter search.  It records every command,
return code, timeout, and output in outputs/engineering-cfd/mms/<case>/ and
fails closed: a failed stage makes all dependent stages ``not_run``.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[2]
INPUT = REPO / "examples/acceptance/circle.xy"
OUTPUT_ROOT = REPO / "outputs/engineering-cfd/mms"
IMAGE = "opencfd/openfoam-run:2606"
TIMEOUT_SECONDS = 120


@dataclass(frozen=True)
class Grid:
    name: str
    wall: str
    background: str


GRIDS = (
    Grid("h04", "0.04", "0.1"),
    Grid("h02", "0.02", "0.05"),
    Grid("h01", "0.01", "0.025"),
)


def command_text(command: list[str]) -> str:
    return shlex.join(command)


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def output_text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def run_stage(case_dir: Path, name: str, command: list[str], timeout: int) -> dict[str, Any]:
    stage_dir = case_dir / name
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "command.txt").write_text(command_text(command) + "\n", encoding="utf-8")
    result: dict[str, Any] = {
        "stage": name,
        "command": command,
        "command_text": command_text(command),
        "timeout_seconds": timeout,
        "status": "running",
        "returncode": None,
        "timed_out": False,
        "stdout": str(stage_dir / "stdout.log"),
        "stderr": str(stage_dir / "stderr.log"),
    }
    try:
        completed = subprocess.run(
            command, cwd=REPO, text=True, capture_output=True, timeout=timeout,
            check=False,
        )
        (stage_dir / "stdout.log").write_text(completed.stdout, encoding="utf-8")
        (stage_dir / "stderr.log").write_text(completed.stderr, encoding="utf-8")
        result["returncode"] = completed.returncode
        result["status"] = "passed" if completed.returncode == 0 else "failed"
    except subprocess.TimeoutExpired as exc:
        (stage_dir / "stdout.log").write_text(output_text(exc.stdout), encoding="utf-8")
        (stage_dir / "stderr.log").write_text(output_text(exc.stderr), encoding="utf-8")
        result["status"] = "timeout"
        result["timed_out"] = True
        result["returncode"] = 124
    except OSError as exc:
        (stage_dir / "stdout.log").write_text("", encoding="utf-8")
        (stage_dir / "stderr.log").write_text(str(exc) + "\n", encoding="utf-8")
        result["status"] = "failed"
        result["returncode"] = 127
    write_json(stage_dir / "stage.json", result)
    return result


def docker(case: Path, *args: str) -> list[str]:
    return [
        "docker", "run", "--rm", "--network", "none", "-v",
        f"{case}:/home/openfoam/workingDir/case", IMAGE, *args,
    ]


def stage_not_run(case_dir: Path, name: str, command: list[str], reason: str) -> dict[str, Any]:
    stage_dir = case_dir / name
    stage_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "stage": name, "command": command, "command_text": command_text(command),
        "timeout_seconds": TIMEOUT_SECONDS, "status": "not_run", "returncode": None,
        "timed_out": False, "reason": reason,
        "stdout": str(stage_dir / "stdout.log"),
        "stderr": str(stage_dir / "stderr.log"),
    }
    (stage_dir / "command.txt").write_text(command_text(command) + "\n", encoding="utf-8")
    (stage_dir / "stdout.log").write_text("", encoding="utf-8")
    (stage_dir / "stderr.log").write_text("", encoding="utf-8")
    write_json(stage_dir / "stage.json", result)
    return result


def fixed_commands(grid: Grid, case_dir: Path) -> dict[str, list[str]]:
    case = case_dir / "case"
    prefix = case_dir / "circle"
    cli = REPO / "build/cartmesh2d_cli"
    independent = REPO / "tools/verification/check_openfoam2d.py"
    mms = REPO / "tools/verification/openfoam_harmonic_mms.py"
    return {
        "generate": [
            str(cli), str(INPUT), str(prefix), "8", "0.25", "0.1", "exterior",
            str(case), "0", "0", "--size-field", "--reference-length", "2",
            "--wall-relative-size", grid.wall, "--background-relative-size", grid.background,
            "--far-field-spans", "0.5", "--cells-per-level", "3",
        ],
        "independent": [
            sys.executable, str(independent), str(case), "--report",
            str(case_dir / "independent.json"),
        ],
        "check_standard": docker(case, "checkMesh", "-case", "/home/openfoam/workingDir/case"),
        "check_all": docker(case, "checkMesh", "-case", "/home/openfoam/workingDir/case",
                             "-allGeometry", "-allTopology"),
        "configure": [
            sys.executable, str(mms), "configure", str(case), "--field", "quadratic",
        ],
        "solve": docker(case, "laplacianFoam", "-case", "/home/openfoam/workingDir/case"),
        "evaluate": [
            sys.executable, str(mms), "evaluate", str(case), "--field", "quadratic",
            "--time", "30", "--log", str(case_dir / "solve" / "stdout.log"),
            "--report", str(case_dir / "mms.json"),
        ],
    }


def checkmesh_gate(stage: dict[str, Any], strict: bool) -> dict[str, Any]:
    """Apply the textual checkMesh gate, keeping extended diagnostics separate."""
    if stage["status"] not in ("passed", "failed"):
        return stage
    stage_dir = Path(stage["stdout"]).parent
    text = "\n".join(
        path.read_text(encoding="utf-8", errors="replace")
        for path in (stage_dir / "stdout.log", stage_dir / "stderr.log")
        if path.exists()
    )
    has_mesh_ok = "Mesh OK." in text
    failed_match = re.search(r"Failed\s+([1-9]\d*)\s+mesh checks", text)
    stage["mesh_ok"] = has_mesh_ok
    stage["failed_mesh_checks"] = int(failed_match.group(1)) if failed_match else 0
    if has_mesh_ok and not failed_match and stage["returncode"] == 0:
        stage["status"] = "passed"
        stage["gate"] = "mesh_ok"
    elif not strict and failed_match and stage["returncode"] != 124:
        stage["status"] = "diagnostic_failed"
        stage["gate"] = "extended_diagnostic"
    else:
        stage["status"] = "failed"
        stage["gate"] = "mesh_check_failed"
    write_json(stage_dir / "stage.json", stage)
    return stage


def stage_allowed_to_continue(name: str, status: str) -> bool:
    if name == "check_all":
        return status in ("passed", "diagnostic_failed")
    return status == "passed"


def run_grid(grid: Grid, dry_run: bool) -> dict[str, Any]:
    case_dir = OUTPUT_ROOT / grid.name
    commands = fixed_commands(grid, case_dir)
    if dry_run:
        stages = [{
            "stage": name,
            "command": command,
            "command_text": command_text(command),
            "timeout_seconds": TIMEOUT_SECONDS,
            "status": "dry_run",
            "returncode": None,
            "timed_out": False,
        } for name, command in commands.items()]
        return {
            "grid": grid.name,
            "parameters": {
                "input": str(INPUT), "reference_length": 2.0,
                "wall_relative_size": float(grid.wall),
                "background_relative_size": float(grid.background),
                "far_field_spans": 0.5, "cells_per_level": 3,
                "small_cell_alpha": 0.1, "fluid_region": "exterior",
            },
            "dry_run": True,
            "mms_pipeline_valid": False,
            "all_mesh_checks_passed": False,
            "valid": False,
            "stages": stages,
        }
    case_dir.mkdir(parents=True, exist_ok=True)
    stages: list[dict[str, Any]] = []
    failed = False
    for name in ("generate", "independent", "check_standard", "check_all",
                 "configure", "solve", "evaluate"):
        command = commands[name]
        if failed:
            stages.append(stage_not_run(case_dir, name, command, "upstream stage failed"))
            continue
        stage = run_stage(case_dir, name, command, TIMEOUT_SECONDS)
        if name == "check_standard":
            stage = checkmesh_gate(stage, strict=True)
        elif name == "check_all":
            stage = checkmesh_gate(stage, strict=False)
        stages.append(stage)
        if not stage_allowed_to_continue(name, stage["status"]):
            failed = True
    stage_status = {stage["stage"]: stage["status"] for stage in stages}
    mms_pipeline_valid = all(
        stage_status.get(name) == "passed"
        for name in ("generate", "independent", "configure", "solve", "evaluate")
    ) and stage_status.get("check_standard") == "passed"
    all_mesh_checks_passed = all(
        stage_status.get(name) == "passed" for name in ("check_standard", "check_all")
    )
    summary = {
        "grid": grid.name,
        "parameters": {
            "input": str(INPUT), "reference_length": 2.0,
            "wall_relative_size": float(grid.wall),
            "background_relative_size": float(grid.background),
            "far_field_spans": 0.5, "cells_per_level": 3,
            "small_cell_alpha": 0.1, "fluid_region": "exterior",
        },
        "dry_run": dry_run,
        "mms_pipeline_valid": mms_pipeline_valid,
        "all_mesh_checks_passed": all_mesh_checks_passed,
        "valid": mms_pipeline_valid and all_mesh_checks_passed,
        "stages": stages,
    }
    write_json(case_dir / "summary.json", summary)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=[grid.name for grid in GRIDS] + ["all"], default="all")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    selected = GRIDS if args.case == "all" else tuple(grid for grid in GRIDS if grid.name == args.case)
    summaries = [run_grid(grid, args.dry_run) for grid in selected]
    overall = {
        "image": IMAGE, "timeout_seconds": TIMEOUT_SECONDS,
        "dry_run": args.dry_run, "grids": summaries,
        "mms_pipeline_valid": all(summary["mms_pipeline_valid"] for summary in summaries),
        "all_mesh_checks_passed": all(summary["all_mesh_checks_passed"] for summary in summaries),
        "valid": all(summary["valid"] for summary in summaries),
    }
    if not args.dry_run:
        write_json(OUTPUT_ROOT / "summary.json", overall)
    print(json.dumps(overall, indent=2, sort_keys=True))
    return 0 if args.dry_run or overall["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
