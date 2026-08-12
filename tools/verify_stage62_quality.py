#!/usr/bin/env python3
"""Cross-check native Stage 6.2 diagnostics against independent readers/checkMesh."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def load(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def close(actual: float, expected: float, absolute: float = 1.0e-9) -> bool:
    return abs(actual - expected) <= absolute * max(1.0, abs(expected))


def case_result(name: str, quality_path: Path, repeat_quality_path: Path,
                reader_path: Path, checkmesh_path: Path) -> dict:
    quality = load(quality_path)
    reader = load(reader_path)
    checkmesh = load(checkmesh_path)
    summary = quality["summary"]
    failures: list[str] = []
    if quality.get("schema") != "cartmesh-solver-mesh-quality-v1":
        failures.append("native_quality_schema")
    if not quality.get("topologyPass") or not quality.get("qualityPass"):
        failures.append("native_quality_gate")
    if summary.get("issueCount") != 0:
        failures.append("native_issue_count")
    if reader.get("status") != "pass" or reader.get("failures"):
        failures.append("independent_reader")
    if checkmesh.get("status") != "pass" or not checkmesh.get("meshOkMarker"):
        failures.append("openfoam_checkmesh")
    for key, native_key, reader_key, check_key in (
        ("cell", "cellCount", "cellCount", "cells"),
        ("face", "faceCount", "faceCount", "faces"),
        ("internal_face", "internalFaceCount", "internalFaceCount", "internalFaces"),
    ):
        values = (summary[native_key], reader[reader_key], checkmesh[check_key])
        if len(set(values)) != 1:
            failures.append(f"{key}_count")
    if not close(summary["minimumCellVolume"], reader["minimumCellVolume"]):
        failures.append("native_reader_minimum_volume")
    if not close(summary["minimumCellVolume"], checkmesh["minimumVolume"]):
        failures.append("native_checkmesh_minimum_volume")
    if abs(summary["maximumNonOrthogonalityDegrees"] -
           checkmesh["maximumNonOrthogonality"]) > 1.0e-3:
        failures.append("non_orthogonality_correspondence")
    quality_sha = sha256(quality_path)
    repeat_sha = sha256(repeat_quality_path)
    if quality_sha != repeat_sha:
        failures.append("quality_report_determinism")
    return {
        "name": name,
        "status": "pass" if not failures else "fail",
        "failures": failures,
        "nativeIssueCount": summary["issueCount"],
        "cellCount": summary["cellCount"],
        "faceCount": summary["faceCount"],
        "minimumCellVolume": summary["minimumCellVolume"],
        "nativeMaximumNonOrthogonalityDegrees":
            summary["maximumNonOrthogonalityDegrees"],
        "openFoamMaximumNonOrthogonalityDegrees":
            checkmesh["maximumNonOrthogonality"],
        "nativeMaximumSkewness": summary["maximumSkewness"],
        "openFoamMaximumSkewness": checkmesh["maximumSkewness"],
        "qualityReportSha256": quality_sha,
        "repeatQualityReportSha256": repeat_sha,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    for name in ("cube", "l"):
        parser.add_argument(f"--{name}-quality", required=True, type=Path)
        parser.add_argument(f"--{name}-repeat-quality", required=True, type=Path)
        parser.add_argument(f"--{name}-reader", required=True, type=Path)
        parser.add_argument(f"--{name}-checkmesh", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    cases = [
        case_result("adaptive_cube", args.cube_quality,
                    args.cube_repeat_quality, args.cube_reader,
                    args.cube_checkmesh),
        case_result("adaptive_l_prism", args.l_quality,
                    args.l_repeat_quality, args.l_reader,
                    args.l_checkmesh),
    ]
    native_skew = [case["nativeMaximumSkewness"] for case in cases]
    openfoam_skew = [case["openFoamMaximumSkewness"] for case in cases]
    same_skewness_ranking = ((native_skew[1] > native_skew[0]) ==
                             (openfoam_skew[1] > openfoam_skew[0]))
    failures = [failure for case in cases for failure in case["failures"]]
    if not same_skewness_ranking:
        failures.append("skewness_ranking")
    result = {
        "schema": "cartmesh-stage62-quality-acceptance-v1",
        "status": "pass" if not failures else "fail",
        "nativeAndOpenFoamUseDifferentSkewnessNormalizations": True,
        "skewnessTrendSameRanking": same_skewness_ranking,
        "cases": cases,
        "failures": failures,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
