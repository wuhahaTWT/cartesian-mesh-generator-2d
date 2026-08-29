#!/usr/bin/env python3
"""Compare shared construction against the retained Q2 legacy implementation."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
from generate_q1_baselines import CASES
from generate_q0_baselines import read_cm2d


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(command, log):
    start = time.perf_counter()
    result = subprocess.run(list(map(str, command)), text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=600)
    log.write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"failed: {command}; see {log}: {result.stdout[-2000:]}")
    return time.perf_counter()-start, result.stdout


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--executable", type=Path, default=Path("build-q2a/cartmesh2d_hybrid_cli"))
    p.add_argument("--evidence", type=Path, default=Path("build-q2a/evidence"))
    p.add_argument("--baseline", type=Path, default=Path("build-q2/after"))
    p.add_argument("--output", type=Path, default=Path("artifacts/q2a/comparison.json"))
    p.add_argument("--collect-only", action="store_true")
    a = p.parse_args()
    a.evidence.mkdir(parents=True, exist_ok=True)
    result = {"format": "cartmesh2d-q2a-construction-comparison-v1",
              "base_commit": "dde1c26eceede55364dca53d3c9ecac2d700b825",
              "executable_sha256": sha(a.executable), "cases": {}}
    sources = [Path("CMakeLists.txt")]
    for root in ("src", "include", "apps", "tests"):
        sources.extend(p for p in Path(root).rglob("*") if p.is_file())
    result["source_files_sha256"] = {str(p): sha(p) for p in sorted(sources)}
    for name, (boundary, params) in CASES.items():
        timings = {}
        commands = {}
        for mode in ("legacy", "shared", "repeat"):
            prefix = a.evidence/f"{mode}-{name}"
            command = [a.executable, boundary, prefix, *params]
            if mode == "legacy":
                command += [a.evidence/f"legacy-{name}-case", ".01", "--legacy-construction"]
            elif mode == "shared":
                command += [a.evidence/f"{name}-case", ".01"]
            commands[mode] = list(map(str, command))
            if not a.collect_only:
                timings[mode], stdout = run(command, a.evidence/f"{mode}-{name}.generation.log")
                assert "hybrid_status=success" in stdout, f"{name}: fallback is not acceptance"
                print(name, mode, round(timings[mode], 3), flush=True)
        timing_file = a.evidence/f"{name}.timings.json"
        if a.collect_only:
            timings = json.loads(timing_file.read_text())
        else:
            timing_file.write_text(json.dumps(timings, indent=2)+"\n")
        def path(mode, suffix):
            return a.evidence/f"{mode}-{name}.hybrid{suffix}"
        before = json.loads(path("legacy", ".quality-contract.json").read_text())
        after = json.loads(path("shared", ".quality-contract.json").read_text())
        geom = json.loads(path("shared", ".json").read_text())
        old_geom = json.loads(path("legacy", ".json").read_text())
        registry = json.loads(path("shared", ".construction.json").read_text())
        solver = read_cm2d(path("shared", ".solver.cm2d"))
        old_solver = read_cm2d(path("legacy", ".solver.cm2d"))
        deterministic = all(sha(path("shared", s)) == sha(path("repeat", s)) for s in
                            (".cm2d", ".solver.cm2d", ".quality-contract.json", ".intersections.json", ".construction.json"))
        baseline_identical = sha(path("legacy", ".solver.cm2d")) == sha(a.baseline/f"{name}.hybrid.solver.cm2d")
        hard = lambda report: Counter(i["metric"] for i in report["issues"] if i["level"] == "hard")
        old_hard, new_hard = hard(before), hard(after)
        regressions = {}
        for section in ("ordinary_metrics", "boundary_layer_metrics"):
            for metric, summary in before[section].items():
                old, new = summary["worst"], after[section][metric]["worst"]
                worse = new < old if summary["worst_direction"] == "min" else new > old
                if worse and abs(new-old) > 1.e-8*max(1., abs(old)):
                    regressions[metric] = [old, new]
        assert len(registry["solver_vertex_handles"]) == len(solver.vertices)
        assert registry["intersection_evaluations"] > 0 and registry["intersection_cache_hits"] > 0
        assert registry["solver_partition_cache_hits"] > 0
        handle_points = dict(zip(registry["solver_vertex_handles"], solver.vertices))
        for e in registry["events"]:
            assert e["local_h"] > 0 and e["feature_classification"] != "none"
            delta = sum((x-y)**2 for x, y in zip(e["original_position"], e["canonical_position"]))**.5
            assert abs(delta-e["displacement"]) <= e["local_h"]*1.e-14
            if e["canonical_handle"] in handle_points:
                assert tuple(e["canonical_position"]) == handle_points[e["canonical_handle"]]
        case = a.evidence/f"{name}-case"
        run([sys.executable, "tools/verification/check_openfoam2d.py", case, "--report",
             a.evidence/f"{name}.independent-foam.json"], a.evidence/f"{name}.independent-foam.log")
        if name in ("circle", "superellipse"):
            run([sys.executable, "tools/verification/check_hybrid_mesh2d.py", path("shared", ".vtk"),
                 path("shared", ".json"), "--output", a.evidence/f"{name}.independent-hybrid.json"],
                a.evidence/f"{name}.independent-hybrid.log")
        foam_log = a.evidence/f"{name}.checkMesh.log"
        foam = "NOT_RUN"
        if foam_log.exists():
            text = foam_log.read_text()
            assert "Mesh OK" in text and "Failed " not in text
            foam = "PASS"
        short = "face_length_over_local_background_h"
        entry = {"input": boundary, "input_sha256": sha(Path(boundary)), "commands": commands,
                 "wall_seconds": timings, "deterministic": deterministic, "legacy_matches_base_solver": baseline_identical,
                 "min_face_over_local_h_before": before["ordinary_metrics"][short]["worst"],
                 "min_face_over_local_h_after": after["ordinary_metrics"][short]["worst"],
                 "min_absolute_face_after": after["legacy_hard_safety"]["min_face_length_absolute"],
                 "hard_counts_before": dict(old_hard), "hard_counts_after": dict(new_hard),
                 "increased_hard_counts": dict(new_hard-old_hard), "worst_metric_regressions": regressions,
                 "area_error_before": old_geom["area_error"], "area_error_after": geom["area_error"],
                 "solver_cells_before": len(old_solver.cells), "solver_cells_after": len(solver.cells),
                 "patch_counts_before": dict(Counter(e.patch for e in old_solver.edges)),
                 "patch_counts_after": dict(Counter(e.patch for e in solver.edges)),
                 "intersection_evaluations": registry["intersection_evaluations"],
                 "intersection_cache_hits": registry["intersection_cache_hits"],
                 "solver_partition_count": registry["solver_partition_count"],
                 "solver_partition_cache_hits": registry["solver_partition_cache_hits"],
                 "independent_reader": "PASS", "checkMesh": foam,
                 "quality_contract_status": after["status"],
                 "solver_sha256": sha(path("shared", ".solver.cm2d")),
                 "construction_sha256": sha(path("shared", ".construction.json"))}
        assert deterministic and baseline_identical and after["legacy_hard_safety"]["valid"]
        assert geom["topology_valid"] and abs(geom["area_error"]) < 1.e-10*geom["expected_fluid_area"]
        assert not entry["increased_hard_counts"] and not regressions, f"{name}: quality regression: {entry}"
        assert entry["patch_counts_before"] == entry["patch_counts_after"]
        assert entry["solver_cells_before"] == entry["solver_cells_after"]
        result["cases"][name] = entry
    result["q2a_status"] = "PASS" if all(c["checkMesh"] == "PASS" for c in result["cases"].values()) else "CHECKMESH_PENDING"
    result["q2_full_status"] = "PARTIAL_NOT_ACCEPTED"
    assert sha(a.executable) == result["executable_sha256"], "executable changed during verification"
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(result["q2a_status"], flush=True)


if __name__ == "__main__":
    main()
