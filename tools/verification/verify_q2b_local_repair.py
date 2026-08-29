#!/usr/bin/env python3
"""Verify the two Q2-B constrained local-repair cases without running full acceptance."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

from generate_q0_baselines import read_cm2d


CASES = ("narrow_gap", "sharp_trailing_edge")
SUFFIXES = (".solver.cm2d", ".quality-contract.json", ".json",
            ".construction.json", ".intersections.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def boundary_patches(mesh) -> Counter[int]:
    return Counter(edge.patch for edge in mesh.edges if edge.neighbour is None)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, default=Path("build-q2a/evidence"))
    parser.add_argument("--after", type=Path, default=Path("build-q2b/final-targeted"))
    parser.add_argument("--output", type=Path,
                        default=Path("build-q2b/q2b-targeted-verification.json"))
    args = parser.parse_args()
    result = {"format": "cartmesh2d-q2b-targeted-v1", "cases": {}}
    for case in CASES:
        before = args.baseline / f"shared-{case}.hybrid"
        after = args.after / f"{case}.hybrid"
        repeat = args.after / f"repeat-{case}.hybrid"
        before_quality = json.loads(Path(str(before)+".quality-contract.json").read_text())
        after_quality = json.loads(Path(str(after)+".quality-contract.json").read_text())
        geometry = json.loads(Path(str(after)+".json").read_text())
        before_mesh = read_cm2d(Path(str(before)+".solver.cm2d"))
        after_mesh = read_cm2d(Path(str(after)+".solver.cm2d"))
        regressions = {}
        for section in ("ordinary_metrics", "boundary_layer_metrics"):
            for metric, summary in before_quality[section].items():
                old = summary["worst"]
                new = after_quality[section][metric]["worst"]
                worse = new < old if summary["worst_direction"] == "min" else new > old
                if worse and abs(new-old) > 1.e-8*max(1.0, abs(old)):
                    regressions[metric] = [old, new]
        hard_faces = [issue for issue in after_quality["issues"]
                      if issue["level"] == "hard" and
                      issue["metric"].startswith("face_length_over_")]
        deterministic = all(sha256(Path(str(after)+suffix)) ==
                            sha256(Path(str(repeat)+suffix)) for suffix in SUFFIXES)
        entry = {
            "min_face_over_local_h_before":
                before_quality["ordinary_metrics"]
                ["face_length_over_local_background_h"]["worst"],
            "min_face_over_local_h_after":
                after_quality["ordinary_metrics"]
                ["face_length_over_local_background_h"]["worst"],
            "min_absolute_face_before":
                before_quality["legacy_hard_safety"]["min_face_length_absolute"],
            "min_absolute_face_after":
                after_quality["legacy_hard_safety"]["min_face_length_absolute"],
            "hard_dimensionless_face_issues_after": len(hard_faces),
            "solver_cells_before": len(before_mesh.cells),
            "solver_cells_after": len(after_mesh.cells),
            "boundary_patch_counts_before": dict(boundary_patches(before_mesh)),
            "boundary_patch_counts_after": dict(boundary_patches(after_mesh)),
            "q2b_candidate_count": geometry["q2b_short_face_candidate_count"],
            "q2b_accepted_transactions": geometry["q2b_accepted_transaction_count"],
            "q2b_solver_cell_reduction": geometry["q2b_solver_cell_reduction"],
            "area_error": geometry["area_error"],
            "worst_metric_regressions": regressions,
            "legacy_solver_quality_valid": after_quality["legacy_hard_safety"]["valid"],
            "topology_valid": geometry["topology_valid"],
            "deterministic": deterministic,
        }
        assert entry["min_face_over_local_h_after"] >= 0.01
        assert entry["hard_dimensionless_face_issues_after"] == 0
        assert entry["legacy_solver_quality_valid"] and entry["topology_valid"]
        assert not regressions and deterministic
        assert entry["boundary_patch_counts_before"] == entry["boundary_patch_counts_after"]
        assert entry["q2b_accepted_transactions"] > 0
        assert entry["solver_cells_before"]-entry["solver_cells_after"] == \
               entry["q2b_solver_cell_reduction"]
        assert abs(entry["area_error"]) < 1.e-10*geometry["expected_fluid_area"]
        result["cases"][case] = entry
    result["status"] = "PASS"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print("Q2-B targeted PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
