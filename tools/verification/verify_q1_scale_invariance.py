#!/usr/bin/env python3
"""Replay H4 cases at a uniform scale and compare Q1 dimensionless conclusions."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import subprocess

from generate_q1_baselines import CASES


def scale_boundary(source: pathlib.Path, destination: pathlib.Path, factor: float) -> None:
    scaled: list[str] = []
    for line in source.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            scaled.append(line)
            continue
        fields = stripped.split()
        if len(fields) != 2:
            raise ValueError(f"expected x y boundary line in {source}: {line}")
        scaled.append(f"{float(fields[0]) * factor:.17g} "
                      f"{float(fields[1]) * factor:.17g}")
    destination.write_text("\n".join(scaled) + "\n", encoding="utf-8")


def close(first: float, second: float, metric: str) -> bool:
    absolute = 2.0e-5 if metric.endswith("_deg") else 1.0e-9
    return math.isclose(first, second, rel_tol=1.0e-8, abs_tol=absolute)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=pathlib.Path, default=pathlib.Path("."))
    parser.add_argument("--build-dir", type=pathlib.Path, default=pathlib.Path("build"))
    parser.add_argument("--reference-dir", type=pathlib.Path,
                        default=pathlib.Path("build/q1_evidence"))
    parser.add_argument("--output-dir", type=pathlib.Path,
                        default=pathlib.Path("build/q1_scale"))
    parser.add_argument("--scale", type=float, default=1000.0)
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    args = parser.parse_args()
    if not math.isfinite(args.scale) or args.scale <= 0.0:
        raise ValueError("scale must be finite and positive")

    repo = args.repo.resolve()
    executable = (repo / args.build_dir / "cartmesh2d_hybrid_cli").resolve()
    reference_dir = (repo / args.reference_dir).resolve()
    output_dir = (repo / args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    # Circle and superellipse cover a regular and the known micro-face case.
    # The H4-3 local-termination generator itself is not scale invariant at
    # extreme factors under its current absolute geometry tolerance; changing
    # that algorithm is explicitly outside Q1.
    selected = args.case or ["circle", "superellipse"]
    result = {"scale": args.scale, "cases": {}}

    for name in selected:
        input_name, raw_parameters = CASES[name]
        scaled_input = output_dir / f"{name}.scaled.xy"
        scale_boundary(repo / input_name, scaled_input, args.scale)
        parameters = list(raw_parameters)
        parameters[4] = float(parameters[4]) * args.scale
        parameters[6] = float(parameters[6]) * args.scale
        prefix = output_dir / name
        subprocess.run([str(executable), str(scaled_input), str(prefix),
                        *[str(value) for value in parameters]], cwd=repo, check=True)
        reference = json.loads((reference_dir /
            f"{name}.hybrid.quality-contract.json").read_text(encoding="utf-8"))
        scaled = json.loads(pathlib.Path(
            f"{prefix}.hybrid.quality-contract.json").read_text(encoding="utf-8"))
        if scaled["status"] != reference["status"]:
            raise AssertionError(f"{name}: overall status changed")
        reference_types = {key: value["status"]
                           for key, value in reference["by_cell_type"].items()}
        scaled_types = {key: value["status"]
                        for key, value in scaled["by_cell_type"].items()}
        if scaled_types != reference_types:
            raise AssertionError(f"{name}: typed status changed")
        maximum_delta = 0.0
        for section in ("ordinary_metrics", "boundary_layer_metrics"):
            if scaled[section].keys() != reference[section].keys():
                raise AssertionError(f"{name}: metric set changed")
            for metric, summary in reference[section].items():
                other = scaled[section][metric]
                for statistic in ("p50", "p95", "p99", "worst"):
                    if not close(float(summary[statistic]),
                                 float(other[statistic]), metric):
                        raise AssertionError(
                            f"{name}: {metric}.{statistic} changed from "
                            f"{summary[statistic]} to {other[statistic]}")
                    maximum_delta = max(maximum_delta,
                        abs(float(summary[statistic]) - float(other[statistic])))
        result["cases"][name] = {
            "status": scaled["status"],
            "by_cell_type": scaled_types,
            "maximum_absolute_metric_delta": maximum_delta,
            "legacy_min_face_reference":
                reference["legacy_hard_safety"]["min_face_length_absolute"],
            "legacy_min_face_scaled":
                scaled["legacy_hard_safety"]["min_face_length_absolute"],
        }

    (output_dir / "scale-invariance.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
