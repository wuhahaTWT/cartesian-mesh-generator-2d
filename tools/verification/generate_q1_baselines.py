#!/usr/bin/env python3
"""Generate compact Q1 baselines from full runtime quality-contract evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import subprocess
from typing import Any


FORMAT_VERSION = "cartmesh2d-q1-quality-contract-v1"
CASES = {
    "circle": ("examples/acceptance/circle.xy", [6, 3, 6, 4, 0.02, 1.2, 1.0]),
    "superellipse": ("examples/complex/superellipse_24.xy", [6, 3, 6, 3, 0.015, 1.15, 1.0]),
    "concave_l": ("examples/h4_3/concave_l.xy", [8, 3, 8, 4, 0.012, 1.15, 1.0]),
    "narrow_gap": ("examples/h4_3/narrow_gap.xy", [8, 3, 8, 4, 0.012, 1.15, 1.0]),
    "sharp_trailing_edge": ("examples/h4_3/sharp_trailing_edge.xy", [8, 3, 8, 4, 0.012, 1.15, 1.0]),
}


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def severity(issue: dict[str, Any]) -> float:
    measured = abs(float(issue["measured"]))
    limit = abs(float(issue["limit"]))
    lower_bound = issue["metric"] in {
        "face_weight", "volume_ratio", "minimum_interior_angle_deg",
        "face_length_over_local_background_h",
        "face_length_over_sqrt_owner_area",
        "face_length_over_sqrt_neighbour_area",
    }
    if lower_bound:
        return float("inf") if measured == 0.0 else limit / measured
    return float("inf") if limit == 0.0 else measured / limit


def compact(report: dict[str, Any], case: str, command: list[str],
            source_commit: str) -> dict[str, Any]:
    hard = [issue for issue in report["issues"] if issue["level"] == "hard"]
    counts = {metric: sum(issue["metric"] == metric for issue in hard)
              for metric in sorted({issue["metric"] for issue in hard})}
    primary = max(hard, key=severity) if hard else None
    return {
        "format_version": FORMAT_VERSION,
        "case": case,
        "source_commit": source_commit,
        "generation_command": command,
        "status": report["status"],
        "dimensionless": report["dimensionless"],
        "contract": report["contract"],
        "by_cell_type": report["by_cell_type"],
        "ordinary_metrics": report["ordinary_metrics"],
        "boundary_layer_metrics": report["boundary_layer_metrics"],
        "legacy_hard_safety": report["legacy_hard_safety"],
        "hard_issue_counts": counts,
        "primary_hard_failure": primary,
        "input_issues": report["input_issues"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=pathlib.Path, default=pathlib.Path("."))
    parser.add_argument("--build-dir", type=pathlib.Path, default=pathlib.Path("build"))
    parser.add_argument("--evidence-dir", type=pathlib.Path,
                        default=pathlib.Path("build/q1_evidence"))
    parser.add_argument("--output-dir", type=pathlib.Path,
                        default=pathlib.Path("artifacts/q1"))
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--collect-only", action="store_true")
    parser.add_argument("--expect-superellipse-short-faces", choices=("present", "absent"),
                        default="present", help="present reproduces historical Q1; absent validates the Q2 superellipse fix")
    args = parser.parse_args()

    repo = args.repo.resolve()
    build_dir = (repo / args.build_dir).resolve()
    evidence_dir = (repo / args.evidence_dir).resolve()
    output_dir = (repo / args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    executable = build_dir / "cartmesh2d_hybrid_cli"
    evidence_dir.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, Any] = {
        "format_version": FORMAT_VERSION,
        "generator_commit": args.source_commit,
        "cases": {},
    }
    for name, (input_name, parameters) in CASES.items():
        input_path = repo / input_name
        prefix = evidence_dir / name
        command = [str(executable), str(input_path), str(prefix),
                   *[str(value) for value in parameters]]
        recorded_command = [str(args.build_dir / "cartmesh2d_hybrid_cli"),
                            input_name, str(args.evidence_dir / name),
                            *[str(value) for value in parameters]]
        if not args.collect_only:
            subprocess.run(command, cwd=repo, check=True)
        full_path = pathlib.Path(f"{prefix}.hybrid.quality-contract.json")
        report = json.loads(full_path.read_text(encoding="utf-8"))
        if report.get("quality_class") != "solver_quality_contract":
            raise ValueError(f"not a Q1 contract report: {full_path}")
        baseline = compact(report, name, recorded_command, args.source_commit)
        baseline_path = output_dir / f"{name}.quality-contract-baseline.json"
        baseline_path.write_text(json.dumps(baseline, indent=2, sort_keys=True) + "\n",
                                 encoding="utf-8")
        manifest["cases"][name] = {
            "input": input_name,
            "input_sha256": sha256(input_path),
            "status": baseline["status"],
            "baseline_sha256": sha256(baseline_path),
            "primary_hard_failure": baseline["primary_hard_failure"],
        }

    superellipse = json.loads((output_dir /
        "superellipse.quality-contract-baseline.json").read_text(encoding="utf-8"))
    micro = superellipse["ordinary_metrics"]["face_length_over_local_background_h"]
    if args.expect_superellipse_short_faces == "present":
        if not micro["worst"] < 0.01 or not any(
                key.startswith("face_length_over_")
                for key in superellipse["hard_issue_counts"]):
            raise RuntimeError("superellipse micro internal face is not a Q1 hard failure")
    elif micro["worst"] < 0.01 or any(
            key.startswith("face_length_over_") for key in superellipse["hard_issue_counts"]):
        raise RuntimeError("Q2 superellipse still contains a Q1 hard short-face violation")
    manifest_path = output_dir / "baseline-manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
