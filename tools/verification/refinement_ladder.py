#!/usr/bin/env python3
"""Measure and gate refinement-ladder robustness for the native 2D generator.

CI exercises every acceptance case at a single refinement level, so nothing
observed what happens when the same geometry is refined further.  This tool runs
a level ladder for one or more cases and records, per level, the deterministic
mesh facts plus the wall-clock attribution the CLIs now report.

The gate is deliberately about *monotonicity*, not about absolute numbers: no
level in a ladder may fail, and the hard solver-quality limits already declared
by ``SolverQualityPolicy2D`` may not be crossed.  It never relaxes a threshold.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import platform
import subprocess
import time
from typing import Any


FORMAT_VERSION = "cartmesh2d-refinement-ladder-v1"

# Hard limits mirrored from include/cartmesh2d/quality/SolverQuality2D.hpp
# (SolverQualityPolicy2D) which in turn match OpenFOAM's
# etc/caseDicts/meshQualityDict.  Mirrored, never redefined: the C++ gate stays
# authoritative and this tool only refuses to accept a ladder that crossed it.
HARD_MIN_FACE_WEIGHT = 0.05
HARD_MIN_VOLUME_RATIO = 0.01
HARD_MAX_NON_ORTHOGONALITY_DEG = 70.0

CASES: dict[str, str] = {
    "circle": "examples/acceptance/circle.xy",
    "airfoil_like": "examples/acceptance/airfoil_like.xy",
    "superellipse": "examples/complex/superellipse_24.xy",
    "concave_l": "examples/h4_3/concave_l.xy",
    "narrow_gap": "examples/h4_3/narrow_gap.xy",
    "sharp_trailing_edge": "examples/h4_3/sharp_trailing_edge.xy",
}

# Per-case hybrid boundary-layer parameters, taken from the existing acceptance
# commands in .github/workflows so a ladder rung at the CI level reproduces CI.
HYBRID_LAYERS: dict[str, tuple[int, float, float, float]] = {
    "circle": (4, 0.02, 1.2, 1.0),
    "airfoil_like": (4, 0.012, 1.15, 1.0),
    "superellipse": (3, 0.015, 1.15, 1.0),
    "concave_l": (4, 0.012, 1.15, 1.0),
    "narrow_gap": (4, 0.012, 1.15, 1.0),
    "sharp_trailing_edge": (4, 0.012, 1.15, 1.0),
}

# Pure Cut-cell ladders need a global floor to be a controlled sequence; a
# boundary-only ladder leaves the far field coarse.  STAGE2DV1C used
# (minimum, boundary) = (level - 2, level).
CUTCELL_MINIMUM_LEVEL_OFFSET = 2
CUTCELL_PADDING_FRACTION = 0.25
CUTCELL_SMALL_ALPHA = 0.10


def parse_key_values(text: str) -> dict[str, str]:
    """Parse the ``key=value`` summaries both CLIs print on stdout."""
    values: dict[str, str] = {}
    for token in text.split():
        key, separator, value = token.partition("=")
        if separator and key:
            values[key] = value
    return values


def number(text: str | None) -> float | None:
    if text is None:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    return parsed


def hybrid_command(executable: pathlib.Path, boundary: pathlib.Path,
                   prefix: pathlib.Path, case: str, level: int,
                   minimum_level: int) -> list[str]:
    layers, thickness, growth, padding = HYBRID_LAYERS[case]
    return [str(executable), str(boundary), str(prefix),
            str(level), str(minimum_level), str(level),
            str(layers), repr(thickness), repr(growth), repr(padding),
            f"{prefix}-case", "0.01"]


def cutcell_command(executable: pathlib.Path, boundary: pathlib.Path,
                    prefix: pathlib.Path, level: int,
                    minimum_level: int) -> list[str]:
    return [str(executable), str(boundary), str(prefix),
            str(level), repr(CUTCELL_PADDING_FRACTION), repr(CUTCELL_SMALL_ALPHA),
            "exterior", f"{prefix}-case", str(minimum_level)]


def hard_issue_counts(path: pathlib.Path) -> dict[str, int] | None:
    """Per-metric Q1 hard counts, from the full runtime contract report."""
    if not path.exists():
        return None
    report = json.loads(path.read_text(encoding="utf-8"))
    counts: dict[str, int] = {}
    for issue in report.get("issues", []):
        if issue.get("level") != "hard":
            continue
        metric = issue.get("metric", "unknown")
        counts[metric] = counts.get(metric, 0) + 1
    return dict(sorted(counts.items()))


def split_hybrid_rejection(stderr: str) -> tuple[dict[str, str], dict[str, str]]:
    """Separate the rejected hybrid candidate from the pure Cut-cell fallback.

    Both halves of the rejection line use the same key names
    (``min_face_weight``, ``max_nonorthogonality``, ...), so a single pass would
    silently record the fallback's numbers as the rung's and hide the hybrid's
    real violation. The fallback usually looks *better* than what was rejected,
    which is exactly the wrong direction for a gate to be wrong in.
    """
    marker = "fallback_failure="
    index = stderr.find(marker)
    if index < 0:
        return parse_key_values(stderr), {}
    return (parse_key_values(stderr[:index]),
            parse_key_values(stderr[index + len(marker):]))


def collect_hybrid(prefix: pathlib.Path, stdout: str, stderr: str) -> dict[str, Any]:
    rung: dict[str, Any] = {}
    # The rejection summary goes to stderr and the timing summary to stdout, so a
    # failing rung is only fully described by both streams.
    rejected, fallback = split_hybrid_rejection(stderr)
    keys = parse_key_values(stdout)
    keys.update(rejected)
    report_path = pathlib.Path(f"{prefix}.hybrid.json")
    if report_path.exists():
        report = json.loads(report_path.read_text(encoding="utf-8"))
        for key in ("hybrid_status", "failure_reason", "quadtree_leaf_count",
                    "cell_count", "solver_cell_count", "boundary_layer_cell_count",
                    "zero_layer_column_count", "termination_cell_count",
                    "remainder_cut_cell_count", "remainder_cartesian_cell_count",
                    "transition_polygon_count", "topology_valid",
                    "mesh_quality_valid", "solver_quality_valid",
                    "quality_contract_status", "solver_quality_issue_count",
                    "area_error", "interface_length_error",
                    "max_non_orthogonality_deg", "min_face_weight",
                    "min_volume_ratio"):
            if key in report:
                rung[key] = report[key]
    rung["hard_issue_counts"] = hard_issue_counts(
        pathlib.Path(f"{prefix}.hybrid.quality-contract.json"))
    if rung["hard_issue_counts"] is not None:
        rung["hard_issue_total"] = sum(rung["hard_issue_counts"].values())
    # The failure path writes no report; the stderr keys are the only record, and
    # they must describe the *rejected* candidate, not the fallback.
    if "hybrid_status" not in rung:
        rung["hybrid_status"] = keys.get("h4_status", "unknown")
        rung["failure_reason"] = keys.get("hybrid_failure", "unknown")
        rung["max_non_orthogonality_deg"] = number(keys.get("max_nonorthogonality"))
        rung["min_face_weight"] = number(keys.get("min_face_weight"))
        rung["min_volume_ratio"] = number(keys.get("min_volume_ratio"))
        rung["rejected_issue_count"] = number(keys.get("issues"))
        if fallback:
            rung["fallback_min_face_weight"] = number(
                fallback.get("min_face_weight"))
            rung["fallback_min_volume_ratio"] = number(
                fallback.get("min_volume_ratio"))
            rung["fallback_max_non_orthogonality_deg"] = number(
                fallback.get("max_nonorthogonality"))
            rung["fallback_issue_count"] = number(fallback.get("issues"))
    rung["fallback_stage"] = keys.get("fallback_stage")
    rung["mesh_mode"] = keys.get("mesh_mode", "hybrid")
    return rung


def collect_cutcell(prefix: pathlib.Path, stdout: str, stderr: str) -> dict[str, Any]:
    rung: dict[str, Any] = {}
    keys = parse_key_values(stdout)
    for key, target in (("leaf_count", "quadtree_leaf_count"),
                        ("cut_cells", "cut_cell_count"),
                        ("source_cells", "source_cell_count"),
                        ("stabilized_cells", "cell_count"),
                        ("small_cells", "small_cell_count"),
                        ("solver_max_nonorthogonality_deg",
                         "max_non_orthogonality_deg"),
                        ("solver_min_face_weight", "min_face_weight"),
                        ("solver_min_volume_ratio", "min_volume_ratio"),
                        ("solver_min_face_length", "min_face_length"),
                        ("solver_max_internal_skewness", "max_internal_skewness")):
        if key in keys:
            rung[target] = number(keys[key]) if "." in keys[key] or "e" in keys[key] \
                else int(keys[key])
    quality_path = pathlib.Path(f"{prefix}-case") / "solver_quality.json"
    if quality_path.exists():
        quality = json.loads(quality_path.read_text(encoding="utf-8"))
        rung["solver_quality_valid"] = bool(quality.get("valid"))
        rung["solver_quality_issue_count"] = int(quality.get("issue_count", 0))
        for metric, value in (quality.get("metrics") or {}).items():
            rung.setdefault(metric, value)
    source = number(keys.get("source_fluid_area"))
    expected = number(keys.get("expected_fluid_area"))
    if source is not None and expected is not None:
        rung["area_error"] = source - expected
    rung["hybrid_status"] = ("success" if "cartmesh2d end-to-end PASS" in stdout
                             else "failed")
    if rung["hybrid_status"] != "success":
        rung["failure_reason"] = cutcell_failure_reason(stderr)
    return rung


# Solver-quality issue codes, in SolverQualityIssueCode2D declaration order
# (include/cartmesh2d/quality/SolverQuality2D.hpp). Mirrored for reporting only.
SOLVER_QUALITY_ISSUE_CODES = (
    "invalid_topology", "invalid_cell", "short_face",
    "excessive_non_orthogonality", "excessive_skewness",
    "excessive_boundary_skewness", "excessive_concavity", "excessive_aspect",
    "small_interior_angle", "low_face_weight", "low_volume_ratio",
)


def cutcell_failure_reason(stderr: str) -> str:
    """Name the first rejection the plain CLI reported, not just 'unknown'."""
    for line in stderr.splitlines():
        if "solver_quality_issue[" not in line:
            continue
        keys = parse_key_values(line)
        raw = keys.get("code")
        if raw is not None and raw.isdigit():
            index = int(raw)
            if index < len(SOLVER_QUALITY_ISSUE_CODES):
                return f"solver_quality:{SOLVER_QUALITY_ISSUE_CODES[index]}"
        return "solver_quality:unclassified"
    for line in reversed(stderr.splitlines()):
        if line.strip():
            return line.strip()[:200]
    return "unknown"


TIMING_KEYS = (
    "timing_input_seconds", "timing_build_seconds", "timing_export_seconds",
    "timing_total_seconds", "timing_refinement_seconds", "timing_balance_seconds",
    "timing_cut_cell_seconds", "timing_source_topology_seconds",
    "timing_agglomeration_seconds", "timing_solver_topology_seconds",
    "timing_serialization_export_seconds",
    "h4_total_seconds", "h4_requested_layer_seconds",
    "h4_requested_hybrid_seconds", "h4_local_layer_seconds",
    "h4_local_hybrid_seconds", "h4_pure_cutcell_fallback_seconds",
    "h4_unattributed_seconds",
)

DETERMINISTIC_COUNTER_KEYS = (
    "h4_requested_hybrid_attempts", "h4_local_hybrid_attempts",
    "h4_pure_cutcell_fallback_attempts", "h4_conformal_hybrid_build_calls",
)


def collect_attribution(stdout: str) -> dict[str, Any]:
    """Wall-clock attribution plus the deterministic call counters.

    Kept in its own section because seconds are not reproducible and must never
    take part in a determinism or gate comparison.
    """
    keys = parse_key_values(stdout)
    seconds = {key: number(keys[key]) for key in TIMING_KEYS if key in keys}
    counters = {key: int(keys[key]) for key in DETERMINISTIC_COUNTER_KEYS
                if key in keys and keys[key].isdigit()}
    return {"measurement_class": "wall_time", "reproducible": False,
            "seconds": seconds, "deterministic_counters": counters}


def run_rung(command: list[str], repo: pathlib.Path, timeout: float,
             env: dict[str, str] | None) -> dict[str, Any]:
    started = time.monotonic()
    try:
        completed = subprocess.run(command, cwd=repo, env=env, timeout=timeout,
                                   capture_output=True, text=True, check=False)
    except subprocess.TimeoutExpired as expired:
        return {"exit_code": None, "timed_out": True,
                "wall_seconds": time.monotonic() - started,
                "stdout": expired.stdout or "", "stderr": expired.stderr or ""}
    wall = time.monotonic() - started
    if "Library not loaded" in completed.stderr:
        raise RuntimeError(
            "the CLI could not start: export DYLD_LIBRARY_PATH to the toolchain "
            "library directory (see --dyld-library-path) before running the ladder")
    return {"exit_code": completed.returncode, "timed_out": False,
            "wall_seconds": wall, "stdout": completed.stdout,
            "stderr": completed.stderr}


def gate(ladders: list[dict[str, Any]]) -> list[str]:
    """Monotonicity gate. Returns the list of violations, empty when clean."""
    violations: list[str] = []
    for ladder in ladders:
        label = f"{ladder['case']}:{ladder['mode']}"
        previous_hard: int | None = None
        for rung in ladder["rungs"]:
            level = rung["level"]
            where = f"{label}:level {level}"
            if rung["timed_out"]:
                violations.append(f"{where} timed out")
                continue
            if rung["exit_code"] != 0:
                violations.append(
                    f"{where} exited {rung['exit_code']} "
                    f"({rung['mesh'].get('failure_reason', 'unknown')})")
                continue
            mesh = rung["mesh"]
            if mesh.get("hybrid_status") != "success":
                violations.append(f"{where} status {mesh.get('hybrid_status')}")
            for key in ("topology_valid", "mesh_quality_valid",
                        "solver_quality_valid"):
                if key in mesh and not mesh[key]:
                    violations.append(f"{where} {key} is false")
            weight = mesh.get("min_face_weight")
            if weight is not None and weight < HARD_MIN_FACE_WEIGHT:
                violations.append(
                    f"{where} min_face_weight {weight} < {HARD_MIN_FACE_WEIGHT}")
            ratio = mesh.get("min_volume_ratio")
            if ratio is not None and ratio < HARD_MIN_VOLUME_RATIO:
                violations.append(
                    f"{where} min_volume_ratio {ratio} < {HARD_MIN_VOLUME_RATIO}")
            angle = mesh.get("max_non_orthogonality_deg")
            if angle is not None and angle > HARD_MAX_NON_ORTHOGONALITY_DEG:
                violations.append(
                    f"{where} max_non_orthogonality_deg {angle} > "
                    f"{HARD_MAX_NON_ORTHOGONALITY_DEG}")
            hard = mesh.get("hard_issue_total")
            if hard is not None:
                if previous_hard is not None and hard > previous_hard:
                    violations.append(
                        f"{where} Q1 hard issues grew {previous_hard} -> {hard}")
                previous_hard = hard
    return violations


def parse_ladder(specification: str) -> tuple[str, str, list[int]]:
    parts = specification.split(":")
    if len(parts) != 3:
        raise ValueError(
            f"expected CASE:MODE:LEVELS, got {specification!r}")
    case, mode, levels = parts
    if case not in CASES:
        raise ValueError(f"unknown case {case!r}; known: {sorted(CASES)}")
    if mode not in ("hybrid", "cutcell"):
        raise ValueError(f"unknown mode {mode!r}; expected hybrid or cutcell")
    parsed = [int(value) for value in levels.split(",") if value]
    if not parsed:
        raise ValueError(f"ladder {specification!r} has no levels")
    if parsed != sorted(parsed):
        raise ValueError(f"ladder {specification!r} levels must ascend")
    return case, mode, parsed


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run and gate refinement ladders for the native 2D generator.")
    parser.add_argument("--repo", type=pathlib.Path, default=pathlib.Path("."))
    parser.add_argument("--build-dir", type=pathlib.Path,
                        default=pathlib.Path("build"))
    parser.add_argument("--evidence-dir", type=pathlib.Path,
                        default=pathlib.Path("build/r2_ladder"))
    parser.add_argument("--output-dir", type=pathlib.Path,
                        default=pathlib.Path("artifacts/r2"))
    parser.add_argument("--manifest-name", default="refinement-ladder-manifest.json")
    parser.add_argument("--ladder", action="append", required=True,
                        metavar="CASE:MODE:LEVELS",
                        help="e.g. circle:cutcell:6,7,8,9 or narrow_gap:hybrid:8,9")
    parser.add_argument("--timeout-seconds", type=float, default=1800.0)
    parser.add_argument("--dyld-library-path",
                        help="prepended to DYLD_LIBRARY_PATH; macOS toolchain builds "
                             "need it or every CLI dies at launch")
    parser.add_argument("--gate", action="store_true",
                        help="fail when the ladder is not monotone")
    parser.add_argument("--stop-on-failure", action="store_true",
                        help="stop a ladder at its first failing level")
    args = parser.parse_args()

    repo = args.repo.resolve()
    build_dir = (repo / args.build_dir).resolve()
    evidence_dir = (repo / args.evidence_dir).resolve()
    output_dir = (repo / args.output_dir).resolve()
    evidence_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    env: dict[str, str] | None = None
    if args.dyld_library_path:
        env = dict(os.environ)
        existing = env.get("DYLD_LIBRARY_PATH")
        env["DYLD_LIBRARY_PATH"] = (f"{args.dyld_library_path}:{existing}"
                                    if existing else args.dyld_library_path)

    hybrid_exe = build_dir / "cartmesh2d_hybrid_cli"
    cutcell_exe = build_dir / "cartmesh2d_cli"
    ladders: list[dict[str, Any]] = []
    for specification in args.ladder:
        case, mode, levels = parse_ladder(specification)
        boundary = repo / CASES[case]
        rungs: list[dict[str, Any]] = []
        for level in levels:
            minimum_level = (level - CUTCELL_MINIMUM_LEVEL_OFFSET if mode == "cutcell"
                             else 3)
            if minimum_level < 0:
                raise ValueError(f"{specification}: level {level} is too coarse")
            prefix = evidence_dir / f"{case}-{mode}-L{level}"
            if mode == "hybrid":
                command = hybrid_command(hybrid_exe, boundary, prefix, case,
                                         level, minimum_level)
            else:
                command = cutcell_command(cutcell_exe, boundary, prefix, level,
                                          minimum_level)
            outcome = run_rung(command, repo, args.timeout_seconds, env)
            collect = collect_hybrid if mode == "hybrid" else collect_cutcell
            rung = {
                "level": level,
                "minimum_level": minimum_level,
                "command": [pathlib.Path(part).name if part == command[0] else part
                            for part in command],
                "exit_code": outcome["exit_code"],
                "timed_out": outcome["timed_out"],
                "mesh": collect(prefix, outcome["stdout"], outcome["stderr"]),
                "attribution": collect_attribution(outcome["stdout"]),
            }
            rung["attribution"]["seconds"]["measured_wall_seconds"] = \
                outcome["wall_seconds"]
            if outcome["exit_code"] != 0 or outcome["timed_out"]:
                rung["stderr_tail"] = outcome["stderr"][-2000:]
            rungs.append(rung)
            if args.stop_on_failure and (outcome["timed_out"] or
                                         outcome["exit_code"] != 0):
                break
        ladders.append({"case": case, "mode": mode, "input": CASES[case],
                        "levels": levels, "rungs": rungs})

    violations = gate(ladders)
    manifest = {
        "format_version": FORMAT_VERSION,
        "platform": platform.platform(),
        "hard_limits": {
            "min_face_weight": HARD_MIN_FACE_WEIGHT,
            "min_volume_ratio": HARD_MIN_VOLUME_RATIO,
            "max_non_orthogonality_deg": HARD_MAX_NON_ORTHOGONALITY_DEG,
            "source": "include/cartmesh2d/quality/SolverQuality2D.hpp",
        },
        "ladders": ladders,
        "gate_violations": violations,
    }
    manifest_path = output_dir / args.manifest_name
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")

    for ladder in ladders:
        for rung in ladder["rungs"]:
            mesh = rung["mesh"]
            seconds = rung["attribution"]["seconds"].get("measured_wall_seconds", 0.0)
            print(f"{ladder['case']}:{ladder['mode']}:L{rung['level']} "
                  f"exit={rung['exit_code']} {seconds:8.2f}s "
                  f"cells={mesh.get('solver_cell_count', mesh.get('cell_count', '-'))} "
                  f"status={mesh.get('hybrid_status', '-')} "
                  f"minFW={mesh.get('min_face_weight', '-')} "
                  f"minVR={mesh.get('min_volume_ratio', '-')} "
                  f"hard={mesh.get('hard_issue_total', '-')}")
    print(f"manifest={manifest_path}")
    if violations:
        for violation in violations:
            print(f"ladder violation: {violation}")
        if args.gate:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
