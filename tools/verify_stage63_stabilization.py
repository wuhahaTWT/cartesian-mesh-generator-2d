#!/usr/bin/env python3
"""Cross-check Stage 6.3 stabilization, independent reader, and checkMesh."""

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


def close(first: float, second: float, tolerance: float = 1.0e-12) -> bool:
    return abs(first - second) <= tolerance * max(1.0, abs(first), abs(second))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stabilization", required=True, type=Path)
    parser.add_argument("--repeat-stabilization", required=True, type=Path)
    parser.add_argument("--quality", required=True, type=Path)
    parser.add_argument("--repeat-quality", required=True, type=Path)
    parser.add_argument("--reader", required=True, type=Path)
    parser.add_argument("--repeat-reader", required=True, type=Path)
    parser.add_argument("--checkmesh", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--repeat-mapping", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    stabilization = load(args.stabilization)
    quality = load(args.quality)
    reader = load(args.reader)
    repeat_reader = load(args.repeat_reader)
    checkmesh = load(args.checkmesh)
    mapping = load(args.mapping)
    failures: list[str] = []

    if stabilization.get("schema") != "cartmesh-solver-mesh-stabilization-v1":
        failures.append("stabilization_schema")
    if (not stabilization.get("pass") or
            not stabilization.get("conservationPass") or
            stabilization.get("unresolvedStableIds")):
        failures.append("stabilization_gate")
    if stabilization.get("agglomerationCount", 0) <= 0:
        failures.append("agglomeration_not_exercised")
    if stabilization.get("finalCellCount") != (
            stabilization.get("initialCellCount", 0) -
            stabilization.get("agglomerationCount", 0)):
        failures.append("cell_count_reduction")
    if not close(stabilization["initialVolume"],
                 stabilization["finalVolume"]):
        failures.append("volume_conservation")
    if any(not close(first, second) for first, second in zip(
            stabilization["initialFirstMoment"],
            stabilization["finalFirstMoment"])):
        failures.append("first_moment_conservation")

    if quality.get("schema") != "cartmesh-solver-mesh-quality-v1" or not (
            quality.get("topologyPass") and quality.get("qualityPass")):
        failures.append("native_quality_gate")
    summary = quality.get("summary", {})
    if summary.get("issueCount") != 0:
        failures.append("native_quality_issues")
    minimum_fraction = summary.get("minimumVolumeFraction", 0.0)
    required_fraction = quality.get("thresholds", {}).get(
        "minimumVolumeFraction", 0.0)
    if minimum_fraction < required_fraction:
        failures.append("minimum_volume_fraction")

    if reader.get("status") != "pass" or reader.get("failures"):
        failures.append("independent_reader")
    if checkmesh.get("status") != "pass" or not checkmesh.get("meshOkMarker"):
        failures.append("openfoam_checkmesh")
    for native_key, reader_key, check_key in (
            ("cellCount", "cellCount", "cells"),
            ("faceCount", "faceCount", "faces"),
            ("internalFaceCount", "internalFaceCount", "internalFaces")):
        if len({summary.get(native_key), reader.get(reader_key),
                checkmesh.get(check_key)}) != 1:
            failures.append(f"count_{native_key}")
    if not close(summary.get("minimumCellVolume", 0.0),
                 reader.get("minimumCellVolume", 0.0), 1.0e-9):
        failures.append("reader_minimum_volume")
    if not close(summary.get("minimumCellVolume", 0.0),
                 checkmesh.get("minimumVolume", 0.0), 1.0e-9):
        failures.append("checkmesh_minimum_volume")

    if mapping.get("schema") != "cartmesh-openfoam-cell-mapping-v2":
        failures.append("mapping_schema")
    mapped_cells = mapping.get("cells", [])
    aggregated_cells = sum(
        len(cell.get("sourceMembers", [])) > 1 for cell in mapped_cells)
    if len(mapped_cells) != summary.get("cellCount") or aggregated_cells <= 0:
        failures.append("source_provenance")

    deterministic = {
        "stabilizationReport": sha256(args.stabilization) ==
            sha256(args.repeat_stabilization),
        "qualityReport": sha256(args.quality) == sha256(args.repeat_quality),
        "polyMesh": reader.get("polyMeshSha256") ==
            repeat_reader.get("polyMeshSha256"),
        "mapping": sha256(args.mapping) == sha256(args.repeat_mapping),
    }
    failures.extend(f"determinism_{name}" for name, passed in
                    deterministic.items() if not passed)

    result = {
        "schema": "cartmesh-stage63-stabilization-acceptance-v1",
        "status": "pass" if not failures else "fail",
        "agglomerationCount": stabilization.get("agglomerationCount"),
        "initialCellCount": stabilization.get("initialCellCount"),
        "finalCellCount": stabilization.get("finalCellCount"),
        "minimumVolumeFraction": minimum_fraction,
        "requiredMinimumVolumeFraction": required_fraction,
        "aggregatedCellCountWithMultipleSourceMembers": aggregated_cells,
        "nativeIssueCount": summary.get("issueCount"),
        "readerStatus": reader.get("status"),
        "openFoamMeshOk": checkmesh.get("meshOkMarker"),
        "deterministic": deterministic,
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
