#!/usr/bin/env python3
"""Measure the size-field configurations against the pre-existing default sizing.

The refinement ladder in `refinement_ladder.py` pins `minimum_level = level - 2`,
so its numbers describe a near-uniform stress mesh rather than a product mesh.
This survey is the product-side counterpart: it runs each acceptance geometry at
the CLI defaults and at each size-field configuration and records the four solver
metrics side by side, so a sizing change can be judged instead of assumed.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import subprocess
import sys

CASES = {
    "circle": "examples/acceptance/circle.xy",
    "airfoil_like": "examples/acceptance/airfoil_like.xy",
    "superellipse": "examples/complex/superellipse_24.xy",
    "naca2412_dense": "examples/complex/naca2412_dense.xy",
    "concave_l": "examples/h4_3/concave_l.xy",
    "narrow_gap": "examples/h4_3/narrow_gap.xy",
    "sharp_trailing_edge": "examples/h4_3/sharp_trailing_edge.xy",
}

# Every configuration keeps the same requested wall cell size, so a change in the
# table is a change in how the field is graded rather than in how fine it is.
CONFIGURATIONS = {
    "legacy-default": [],
    "boundary-only": ["--size-field", "--cells-per-level", "0"],
    "graded-1": ["--size-field", "--cells-per-level", "1"],
    "graded-3": ["--size-field", "--cells-per-level", "3"],
    "graded-3-curv-gap": ["--size-field", "--cells-per-level", "3",
                          "--curvature-cells-per-radius", "8", "--gap-cells", "4"],
    "graded-3-wake": ["--size-field", "--cells-per-level", "3", "--wake", "0", "6", "0.6", "4"],
}

METRICS = ("leaf_count", "stabilized_cells", "solver_min_volume_ratio",
           "solver_min_face_weight", "solver_max_nonorthogonality_deg",
           "solver_max_internal_skewness", "solver_max_boundary_skewness",
           "timing_total_seconds", "size_field_max_level", "size_field_wall_cell_size")


def parse_keys(stdout: str) -> dict[str, str]:
    return dict(re.findall(r"^([a-z0-9_]+)=(.+)$", stdout, re.M))


def run_case(cli: pathlib.Path, repo: pathlib.Path, boundary: str,
             output: pathlib.Path, wall_cells: float, legacy_level: int,
             extra: list[str], env: dict[str, str]) -> dict[str, object]:
    command = [str(cli), str(repo / boundary), str(output)]
    if extra:
        # max-level and padding-fraction are ignored while the field is active, but
        # they are positional so they still have to be supplied.
        command += [str(legacy_level), "0.25", "0.15", "exterior", f"{output}-case", "0", "0"]
        command += extra + ["--wall-cells-per-span", repr(wall_cells)]
    else:
        command += [str(legacy_level), "0.25", "0.15", "exterior", f"{output}-case", "0"]
    completed = subprocess.run(command, capture_output=True, text=True, env=env)
    record: dict[str, object] = {"command": command, "exit_code": completed.returncode}
    if completed.returncode != 0:
        record["failure"] = (completed.stderr.strip().splitlines() or ["unknown"])[0]
        return record
    keys = parse_keys(completed.stdout)
    for metric in METRICS:
        if metric in keys:
            record[metric] = keys[metric]
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--build-dir", type=pathlib.Path, default=pathlib.Path("build"))
    parser.add_argument("--output-dir", type=pathlib.Path,
                        default=pathlib.Path("build/size-field-survey"))
    parser.add_argument("--manifest", type=pathlib.Path,
                        default=pathlib.Path("artifacts/wb/size-field-survey.json"))
    parser.add_argument("--wall-cells-per-span", type=float, default=128.0)
    parser.add_argument("--legacy-level", type=int, default=8)
    parser.add_argument("--dyld-library-path", default="/Applications/mesasdk/lib")
    parser.add_argument("--case", action="append", default=[])
    arguments = parser.parse_args()

    repo = arguments.repo.resolve()
    cli = (repo / arguments.build_dir / "cartmesh2d_cli").resolve()
    if not cli.exists():
        print(f"missing {cli}", file=sys.stderr)
        return 1
    output_root = (repo / arguments.output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    import os
    env = dict(os.environ)
    if arguments.dyld_library_path:
        existing = env.get("DYLD_LIBRARY_PATH", "")
        env["DYLD_LIBRARY_PATH"] = (f"{arguments.dyld_library_path}:{existing}"
                                    if existing else arguments.dyld_library_path)

    selected = arguments.case or list(CASES)
    results: dict[str, dict[str, object]] = {}
    for case in selected:
        if case not in CASES:
            print(f"unknown case {case}", file=sys.stderr)
            return 1
        results[case] = {}
        for name, extra in CONFIGURATIONS.items():
            output = output_root / f"{case}-{name}"
            results[case][name] = run_case(cli, repo, CASES[case], output,
                                           arguments.wall_cells_per_span,
                                           arguments.legacy_level, extra, env)
            entry = results[case][name]
            status = entry.get("failure", "ok")
            print(f"{case:22s} {name:18s} exit={entry['exit_code']} "
                  f"leaves={entry.get('leaf_count','-'):>7} "
                  f"cells={entry.get('stabilized_cells','-'):>7} "
                  f"volR={entry.get('solver_min_volume_ratio','-'):>12} "
                  f"faceW={entry.get('solver_min_face_weight','-'):>12} "
                  f"{'' if status == 'ok' else status[:60]}")

    manifest = {
        "format": "cartmesh2d-size-field-survey-v1",
        "wall_cells_per_span": arguments.wall_cells_per_span,
        "legacy_level": arguments.legacy_level,
        "configurations": CONFIGURATIONS,
        "results": results,
    }
    destination = repo / arguments.manifest
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"\nwrote {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
