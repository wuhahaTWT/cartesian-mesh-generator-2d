#!/usr/bin/env python3
"""Measure how sensitive the solver-quality gate is to grid alignment alone.

Every other measurement in this repository fixes the geometry and varies the
sizing, which cannot separate "this sizing is better" from "this sizing happened
to move the boundary vertices off the grid lines".  This tool holds the sizing
request fixed and perturbs the body scale by a fraction of a percent, so the mesh
resolution is unchanged and the only thing that moves is where the wall vertices
land relative to the dyadic grid.

Report the pass rate, not just the metrics: a gate that flips on a sub-percent
geometry change is sampling alignment, not quality.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import re
import statistics
import subprocess
import sys

METRICS = ("solver_min_volume_ratio", "solver_min_face_weight",
           "solver_max_nonorthogonality_deg", "solver_max_internal_skewness",
           "solver_max_boundary_skewness", "leaf_count", "stabilized_cells")


def write_polygon(path: pathlib.Path, vertices: int, scale: float) -> None:
    with path.open("w") as handle:
        handle.write(f"# {vertices}-gon scale {scale!r}\n")
        for index in range(vertices):
            angle = 2.0 * math.pi * index / vertices
            handle.write(f"{scale * math.cos(angle):.12f} {scale * math.sin(angle):.12f}\n")


def scale_source(path: pathlib.Path, destination: pathlib.Path, scale: float) -> None:
    rows = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            rows.append(line)
            continue
        parts = stripped.split()
        rows.append(f"{float(parts[0]) * scale:.12f} {float(parts[1]) * scale:.12f}")
    destination.write_text("\n".join(rows) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--build-dir", type=pathlib.Path, default=pathlib.Path("build"))
    parser.add_argument("--work-dir", type=pathlib.Path,
                        default=pathlib.Path("build/alignment-sensitivity"))
    parser.add_argument("--manifest", type=pathlib.Path,
                        default=pathlib.Path("artifacts/wb/alignment-sensitivity.json"))
    parser.add_argument("--boundary", default="", help="source .xy; default is a generated 32-gon")
    parser.add_argument("--polygon-vertices", type=int, default=32)
    parser.add_argument("--samples", type=int, default=20)
    parser.add_argument("--scale-step", type=float, default=0.001)
    parser.add_argument("--wall-cells-per-span", type=float, default=0.0,
                        help="0 leaves it unset so the wall lands on the safe ceiling")
    parser.add_argument("--far-field-spans", type=float, default=10.0)
    parser.add_argument("--cells-per-level", type=int, default=1)
    parser.add_argument("--max-safe-wall-level", type=int, default=0,
                        help="0 keeps the built-in ceiling")
    parser.add_argument("--allow-unsafe-wall-level", action="store_true")
    parser.add_argument("--small-alpha", type=float, default=0.15)
    parser.add_argument("--dyld-library-path", default="/Applications/mesasdk/lib")
    arguments = parser.parse_args()

    repo = arguments.repo.resolve()
    cli = (repo / arguments.build_dir / "cartmesh2d_cli").resolve()
    if not cli.exists():
        print(f"missing {cli}", file=sys.stderr)
        return 1
    work = (repo / arguments.work_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)

    env = dict(os.environ)
    if arguments.dyld_library_path:
        existing = env.get("DYLD_LIBRARY_PATH", "")
        env["DYLD_LIBRARY_PATH"] = (f"{arguments.dyld_library_path}:{existing}"
                                    if existing else arguments.dyld_library_path)

    samples: list[dict[str, object]] = []
    for index in range(arguments.samples):
        scale = 1.0 + arguments.scale_step * index
        boundary = work / f"sample{index:03d}.xy"
        if arguments.boundary:
            scale_source(repo / arguments.boundary, boundary, scale)
        else:
            write_polygon(boundary, arguments.polygon_vertices, scale)
        prefix = work / f"sample{index:03d}"
        command = [str(cli), str(boundary), str(prefix), "8", "0.25",
                   repr(arguments.small_alpha), "exterior", f"{prefix}-case", "0", "0",
                   "--size-field", "--cells-per-level", str(arguments.cells_per_level),
                   "--far-field-spans", repr(arguments.far_field_spans)]
        if arguments.wall_cells_per_span > 0.0:
            command += ["--wall-cells-per-span", repr(arguments.wall_cells_per_span)]
        if arguments.max_safe_wall_level > 0:
            command += ["--max-safe-wall-level", str(arguments.max_safe_wall_level)]
        if arguments.allow_unsafe_wall_level:
            command += ["--allow-unsafe-wall-level"]
        completed = subprocess.run(command, capture_output=True, text=True, env=env)
        record: dict[str, object] = {"scale": scale, "exit_code": completed.returncode}
        if completed.returncode == 0:
            keys = dict(re.findall(r"^([a-z0-9_]+)=(.+)$", completed.stdout, re.M))
            for metric in METRICS:
                if metric in keys:
                    record[metric] = float(keys[metric])
        else:
            failure = (completed.stderr.strip().splitlines() or ["unknown"])[0]
            record["failure"] = failure
            # A refusal by the safe-level guard is not an alignment failure; counting
            # the two together would hide exactly what this tool exists to measure.
            record["refused_by_guard"] = "maxSafeWallLevel" in failure
        samples.append(record)
        if completed.returncode == 0:
            status = "PASS"
        elif record.get("refused_by_guard"):
            status = "REFUSED"
        else:
            status = "GATE FAIL"
        print(f"scale={scale:<8.4f} {status:9s} "
              f"volR={record.get('solver_min_volume_ratio','-')!s:>14} "
              f"faceW={record.get('solver_min_face_weight','-')!s:>14}")

    passed = [sample for sample in samples if sample["exit_code"] == 0]
    refused = [sample for sample in samples if sample.get("refused_by_guard")]
    graded = [sample for sample in samples if not sample.get("refused_by_guard")]
    summary: dict[str, object] = {
        "samples": len(samples),
        "passed": len(passed),
        "refused_by_guard": len(refused),
        "alignment_pass_rate": len(passed) / len(graded) if graded else 0.0,
    }
    for metric in ("solver_min_volume_ratio", "solver_min_face_weight"):
        values = [float(sample[metric]) for sample in passed if metric in sample]
        if values:
            summary[metric] = {"min": min(values), "max": max(values),
                               "median": statistics.median(values),
                               "distinct": len(set(values))}
    print(f"\npass {summary['passed']}/{summary['samples']} "
          f"(refused by guard {summary['refused_by_guard']}, "
          f"alignment pass rate {summary['alignment_pass_rate']:.2f})")
    for metric in ("solver_min_volume_ratio", "solver_min_face_weight"):
        if metric in summary:
            print(f"{metric}: {summary[metric]}")

    destination = repo / arguments.manifest
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(
        {"format": "cartmesh2d-alignment-sensitivity-v1",
         "request": {"wall_cells_per_span": arguments.wall_cells_per_span or "unset",
                     "far_field_spans": arguments.far_field_spans,
                     "cells_per_level": arguments.cells_per_level,
                     "max_safe_wall_level": arguments.max_safe_wall_level or "default",
                     "allow_unsafe_wall_level": arguments.allow_unsafe_wall_level,
                     "small_alpha": arguments.small_alpha,
                     "boundary": arguments.boundary or f"generated-{arguments.polygon_vertices}-gon"},
         "summary": summary, "samples": samples}, indent=2, sort_keys=True) + "\n")
    print(f"wrote {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
