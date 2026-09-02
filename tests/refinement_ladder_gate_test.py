#!/usr/bin/env python3
"""Unit tests for the refinement-ladder gate.

The gate decides whether a refinement ladder regressed, so it must be tested on
synthetic rungs rather than only by running the real CLIs: a gate that silently
accepts a failing rung would let exactly the class of defect R2 exists to catch
through. No mesh is generated here.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] /
                       "tools" / "verification"))

import refinement_ladder as ladder  # noqa: E402


FAILURES: list[str] = []


def check(condition: bool, message: str) -> None:
    if not condition:
        FAILURES.append(message)


def rung(level: int, **mesh: object) -> dict[str, object]:
    return {"level": level, "exit_code": 0, "timed_out": False, "mesh": dict(mesh)}


def clean_mesh(**overrides: object) -> dict[str, object]:
    mesh: dict[str, object] = {
        "hybrid_status": "success",
        "topology_valid": True,
        "mesh_quality_valid": True,
        "solver_quality_valid": True,
        "min_face_weight": 0.20,
        "min_volume_ratio": 0.20,
        "max_non_orthogonality_deg": 40.0,
        "hard_issue_total": 10,
    }
    mesh.update(overrides)
    return mesh


def ladder_of(*rungs: dict[str, object]) -> list[dict[str, object]]:
    return [{"case": "circle", "mode": "hybrid", "rungs": list(rungs)}]


def test_monotone_ladder_passes() -> None:
    violations = ladder.gate(ladder_of(
        rung(6, **clean_mesh(hard_issue_total=20)),
        rung(7, **clean_mesh(hard_issue_total=20)),
        rung(8, **clean_mesh(hard_issue_total=12)),
    ))
    check(violations == [], f"a monotone ladder must pass, got {violations}")


def test_growing_hard_issue_count_is_a_violation() -> None:
    violations = ladder.gate(ladder_of(
        rung(6, **clean_mesh(hard_issue_total=80)),
        rung(7, **clean_mesh(hard_issue_total=132)),
    ))
    check(len(violations) == 1, f"expected one violation, got {violations}")
    check(any("hard issues grew 80 -> 132" in item for item in violations),
          f"the growth must be named explicitly, got {violations}")


def test_hard_limits_are_enforced_at_the_documented_values() -> None:
    # These are the SolverQualityPolicy2D limits, which also match OpenFOAM's
    # meshQualityDict. Crossing them is a failure, not a warning.
    check(ladder.HARD_MIN_FACE_WEIGHT == 0.05, "min face weight limit is 0.05")
    check(ladder.HARD_MIN_VOLUME_RATIO == 0.01, "min volume ratio limit is 0.01")
    check(ladder.HARD_MAX_NON_ORTHOGONALITY_DEG == 70.0,
          "max non-orthogonality limit is 70 degrees")

    violations = ladder.gate(ladder_of(
        rung(9, **clean_mesh(min_face_weight=0.0244, min_volume_ratio=0.0053,
                             max_non_orthogonality_deg=78.6)),
    ))
    check(any("min_face_weight" in item for item in violations),
          f"a low face weight must be reported, got {violations}")
    check(any("min_volume_ratio" in item for item in violations),
          f"a low volume ratio must be reported, got {violations}")
    check(any("max_non_orthogonality_deg" in item for item in violations),
          f"an excessive non-orthogonality must be reported, got {violations}")


def test_failed_and_timed_out_rungs_are_violations() -> None:
    failed = {"level": 9, "exit_code": 1, "timed_out": False,
              "mesh": {"failure_reason": "solver_quality_failed"}}
    timed_out = {"level": 10, "exit_code": None, "timed_out": True, "mesh": {}}
    violations = ladder.gate(ladder_of(failed, timed_out))
    check(any("exited 1" in item and "solver_quality_failed" in item
              for item in violations),
          f"a failing rung must name its reason, got {violations}")
    check(any("timed out" in item for item in violations),
          f"a timed-out rung must be a violation, got {violations}")


def test_invalid_flags_are_violations_even_when_the_exit_code_is_zero() -> None:
    violations = ladder.gate(ladder_of(
        rung(8, **clean_mesh(solver_quality_valid=False)),
    ))
    check(any("solver_quality_valid is false" in item for item in violations),
          f"an invalid flag must be a violation, got {violations}")


def test_ladder_specifications_are_validated() -> None:
    case, mode, levels = ladder.parse_ladder("circle:cutcell:6,7,8")
    check((case, mode, levels) == ("circle", "cutcell", [6, 7, 8]),
          "a well-formed specification parses")
    for bad in ("circle:cutcell", "nope:cutcell:6", "circle:wrong:6",
                "circle:cutcell:", "circle:cutcell:8,7"):
        try:
            ladder.parse_ladder(bad)
        except ValueError:
            continue
        FAILURES.append(f"{bad!r} must be rejected")


def test_key_value_parsing_matches_the_cli_summaries() -> None:
    parsed = ladder.parse_key_values(
        "h4_status=failed min_face_weight=0.0244 issue[0]=(code=9,cell=1) noise")
    check(parsed["h4_status"] == "failed", "plain keys parse")
    check(parsed["min_face_weight"] == "0.0244", "numeric keys parse")
    check("noise" not in parsed, "bare tokens are not keys")
    check(ladder.number("1e-6") == 1e-06, "exponent notation parses")
    check(ladder.number("nope") is None, "non-numbers become None")
    check(ladder.number(None) is None, "missing values stay None")


# Real stderr shape from circle:hybrid level 9. The rejected candidate and the
# pure Cut-cell fallback reuse the same key names, and the fallback looks better,
# so a single-pass parse would record 0.110669 as the rung's face weight and let
# the actual 0.0244111 violation through the gate.
HYBRID_LEVEL9_STDERR = (
    "h4_status=failed mesh_mode=failed fallback_stage=hybrid_candidate "
    "requested_layer_failure=none local_layer_failure=none "
    "hybrid_failure=solver_quality_failed "
    "hybrid_detail=solver quality remains invalid after constrained repair: "
    "issues=233 max_nonorthogonality=78.6117 min_face_weight=0.0244111 "
    "min_volume_ratio=0.00228732 issue[0]=(code=9,cell=1723,edge=486) "
    "fallback_failure=pure Cut-cell fallback solver quality failed: issues=16 "
    "max_nonorthogonality=45.859 min_face_weight=0.110669 "
    "min_volume_ratio=0.0787994\n")


def test_rejected_candidate_metrics_are_not_shadowed_by_the_fallback() -> None:
    rejected, fallback = ladder.split_hybrid_rejection(HYBRID_LEVEL9_STDERR)
    check(rejected["min_face_weight"] == "0.0244111",
          f"the rejected face weight must survive, got {rejected.get('min_face_weight')}")
    check(fallback["min_face_weight"] == "0.110669",
          f"the fallback face weight is recorded separately, got {fallback}")

    mesh = ladder.collect_hybrid(pathlib.Path("/nonexistent-prefix"), "",
                                 HYBRID_LEVEL9_STDERR)
    check(mesh["min_face_weight"] == 0.0244111,
          f"the rung reports the rejected candidate, got {mesh['min_face_weight']}")
    check(mesh["min_volume_ratio"] == 0.00228732, "the rejected volume ratio is kept")
    check(mesh["max_non_orthogonality_deg"] == 78.6117,
          "the rejected non-orthogonality is kept")
    check(mesh["fallback_min_face_weight"] == 0.110669,
          "the fallback metrics stay under their own keys")
    check(mesh["failure_reason"] == "solver_quality_failed",
          "the failure reason comes from the rejected half")

    # And the gate must then actually fire on those numbers.
    violations = ladder.gate([{"case": "circle", "mode": "hybrid", "rungs": [
        {"level": 9, "exit_code": 0, "timed_out": False, "mesh": mesh}]}])
    check(any("min_face_weight" in item for item in violations),
          f"the rejected face weight must trip the gate, got {violations}")


def test_stderr_without_a_fallback_half_still_parses() -> None:
    rejected, fallback = ladder.split_hybrid_rejection(
        "h4_status=failed min_face_weight=0.01\n")
    check(rejected["min_face_weight"] == "0.01", "a single half parses")
    check(fallback == {}, "there is no fallback half to report")


def main() -> int:
    for name, function in sorted(globals().items()):
        if name.startswith("test_") and callable(function):
            function()
    if FAILURES:
        for failure in FAILURES:
            print(f"FAIL: {failure}")
        print(f"{len(FAILURES)} refinement-ladder gate checks failed")
        return 1
    print("refinement ladder gate tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
