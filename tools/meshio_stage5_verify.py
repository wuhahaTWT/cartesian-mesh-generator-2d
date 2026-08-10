#!/usr/bin/env python3
"""使用 meshio/NumPy 独立读取并验证阶段 5 增量输出。"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import meshio
import numpy as np


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as source:
        return json.load(source)


def cell_data(mesh: meshio.Mesh, name: str) -> np.ndarray:
    if name not in mesh.cell_data:
        raise RuntimeError(f"missing cell data field: {name}")
    arrays = mesh.cell_data[name]
    if len(arrays) != len(mesh.cells):
        raise RuntimeError(f"cell data block mismatch: {name}")
    return np.concatenate([np.asarray(array, dtype=np.float64) for array in arrays])


def verify(args: argparse.Namespace) -> dict[str, Any]:
    old_mesh = meshio.read(args.old_mesh)
    new_mesh = meshio.read(args.new_mesh)
    full_mesh = meshio.read(args.full_mesh)
    report = load_json(args.report)
    mapping = load_json(args.mapping)

    old_count = sum(len(block.data) for block in old_mesh.cells)
    new_count = sum(len(block.data) for block in new_mesh.cells)
    full_count = sum(len(block.data) for block in full_mesh.cells)
    failures: list[str] = []
    if old_count != report["oldLeafCount"]:
        failures.append("old_leaf_count")
    if new_count != report["newLeafCount"]:
        failures.append("new_leaf_count")
    if full_count != report["newLeafCount"]:
        failures.append("full_leaf_count")

    old_fraction = cell_data(old_mesh, "fluid_volume_fraction")
    new_fraction = cell_data(new_mesh, "fluid_volume_fraction")
    rebuilt = cell_data(new_mesh, "stage5_rebuilt")
    reused = cell_data(new_mesh, "stage5_reused")
    low = np.rint(cell_data(new_mesh, "octree_node_code_low32")).astype(np.uint64)
    high = np.rint(cell_data(new_mesh, "octree_node_code_high32")).astype(np.uint64)
    codes = low | (high << np.uint64(32))
    full_fraction = cell_data(full_mesh, "fluid_volume_fraction")
    full_cut = cell_data(full_mesh, "cut_cell")
    full_low = np.rint(cell_data(full_mesh, "octree_node_code_low32")).astype(np.uint64)
    full_high = np.rint(cell_data(full_mesh, "octree_node_code_high32")).astype(np.uint64)
    full_codes = full_low | (full_high << np.uint64(32))
    if not np.array_equal(codes, full_codes):
        failures.append("external_stable_code_equivalence")
    if not np.array_equal(new_fraction, full_fraction):
        failures.append("external_fluid_fraction_equivalence")
    if not np.array_equal(cell_data(new_mesh, "cut_cell"), full_cut):
        failures.append("external_cut_mask_equivalence")

    if np.any(old_fraction < -1.0e-12) or np.any(old_fraction > 1.0 + 1.0e-12):
        failures.append("old_fluid_fraction_range")
    if np.any(new_fraction < -1.0e-12) or np.any(new_fraction > 1.0 + 1.0e-12):
        failures.append("new_fluid_fraction_range")
    if not np.all(np.isin(rebuilt, [0.0, 1.0])):
        failures.append("rebuilt_mask_not_binary")
    if not np.all(np.isin(reused, [0.0, 1.0])):
        failures.append("reused_mask_not_binary")
    if not np.allclose(rebuilt + reused, 1.0, atol=0.0, rtol=0.0):
        failures.append("rebuilt_reused_not_complementary")
    if int(np.count_nonzero(reused)) != report["geometryReusedLeafCount"]:
        failures.append("reused_count")
    if int(np.count_nonzero(rebuilt)) != report["geometryRebuiltLeafCount"]:
        failures.append("rebuilt_count")
    if len(np.unique(codes)) != new_count:
        failures.append("stable_node_code_uniqueness")

    entries = mapping.get("entries", [])
    if mapping.get("schema") != "cartmesh-stage5-mapping-v1":
        failures.append("mapping_schema")
    if len(entries) != report["mappingEntryCount"]:
        failures.append("mapping_entry_count")
    background_sum = sum(float(entry["backgroundOverlapVolume"]) for entry in entries)
    if not np.isclose(background_sum, report["mappingBackgroundVolume"],
                      rtol=1.0e-12, atol=1.0e-12):
        failures.append("mapping_background_volume")
    for entry in entries:
        if not entry.get("exactFluidOverlap", False):
            failures.append("mapping_not_exact")
            break
        old_weight = float(entry["oldVolumePreservedFraction"])
        new_weight = float(entry["newVolumeFromOldFraction"])
        if not (-1.0e-12 <= old_weight <= 1.0 + 1.0e-10 and
                -1.0e-12 <= new_weight <= 1.0 + 1.0e-10):
            failures.append("mapping_weight_range")
            break

    required_true = (
        "topologyEqualsFullRebuild",
        "geometryEqualsFullRebuild",
        "partitionValid",
        "faceBalanceValid",
    )
    for field in required_true:
        if report.get(field) is not True:
            failures.append(field)
    required_zero = (
        "nonclosedCellCount",
        "negativeVolumeCellCount",
        "sharedFaceMismatchCount",
        "classificationConflictCount",
    )
    for field in required_zero:
        if report.get(field) != 0:
            failures.append(field)

    return {
        "schema": "cartmesh-stage5-meshio-v1",
        "status": "pass" if not failures else "fail",
        "stage5CaseExternalPass": not failures,
        "externalIndependentReaderAccepted": not failures,
        "reader": f"meshio {meshio.__version__}",
        "oldLeafCount": old_count,
        "newLeafCount": new_count,
        "fullRebuildLeafCount": full_count,
        "externalBackgroundFieldsEqualFullRebuild": not any(
            failure.startswith("external_") for failure in failures
        ),
        "uniqueStableNodeCodeCount": int(len(np.unique(codes))),
        "rebuiltLeafCount": int(np.count_nonzero(rebuilt)),
        "reusedLeafCount": int(np.count_nonzero(reused)),
        "minimumFluidVolumeFraction": float(np.min(new_fraction)),
        "maximumFluidVolumeFraction": float(np.max(new_fraction)),
        "mappingEntryCount": len(entries),
        "mappingBackgroundVolume": background_sum,
        "failures": failures,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old-mesh", type=Path, required=True)
    parser.add_argument("--new-mesh", type=Path, required=True)
    parser.add_argument("--full-mesh", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        json.dump(result, output, ensure_ascii=False, indent=2, sort_keys=True)
        output.write("\n")
    raise SystemExit(0 if result["status"] == "pass" else 2)


if __name__ == "__main__":
    main()
