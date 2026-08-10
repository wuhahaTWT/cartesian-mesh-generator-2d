#!/usr/bin/env python3
"""Use meshio/NumPy to independently verify a Stage-2 adaptive-octree VTU."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
from collections import defaultdict

import meshio
import numpy as np


def decode_morton(code: int) -> tuple[int, int, int]:
    coordinates = [0, 0, 0]
    for bit in range(21):
        for axis in range(3):
            coordinates[axis] |= ((code >> (3 * bit + axis)) & 1) << bit
    return tuple(coordinates)


def decode_node(code: int) -> tuple[int, int, int, int]:
    if code <= 0:
        raise ValueError(f"invalid zero/negative node code {code}")
    highest = code.bit_length() - 1
    if highest % 3:
        raise ValueError(f"node code {code} has no level sentinel")
    level = highest // 3
    morton = code ^ (1 << highest)
    x, y, z = decode_morton(morton)
    limit = 1 << level
    if level > 21 or x >= limit or y >= limit or z >= limit:
        raise ValueError(f"node code {code} decodes outside level {level}")
    return level, x, y, z


def exact_cube_classification(minimum: np.ndarray, maximum: np.ndarray) -> int:
    overlaps = np.all(maximum >= 0.0) and np.all(minimum <= 1.0)
    surface = False
    if overlaps:
        for axis in range(3):
            others = [candidate for candidate in range(3) if candidate != axis]
            other_overlap = all(
                maximum[other] >= 0.0 and minimum[other] <= 1.0 for other in others
            )
            if other_overlap and (
                minimum[axis] <= 0.0 <= maximum[axis]
                or minimum[axis] <= 1.0 <= maximum[axis]
            ):
                surface = True
                break
    if surface:
        return 2
    center = 0.5 * (minimum + maximum)
    return 1 if np.all(center > 0.0) and np.all(center < 1.0) else 0


def face_balance(cells: list[tuple[int, int, int, int]], maximum_level: int) -> dict[str, int | bool]:
    faces: list[dict[int, tuple[list[tuple[int, int, int, int, int]], list[tuple[int, int, int, int, int]]]]] = [
        defaultdict(lambda: ([], [])) for _ in range(3)
    ]
    for level, x, y, z in cells:
        scale = 1 << (maximum_level - level)
        starts = [x * scale, y * scale, z * scale]
        ends = [(x + 1) * scale, (y + 1) * scale, (z + 1) * scale]
        for axis in range(3):
            transverse = [candidate for candidate in range(3) if candidate != axis]
            rectangle = (
                starts[transverse[0]],
                ends[transverse[0]],
                starts[transverse[1]],
                ends[transverse[1]],
                level,
            )
            faces[axis][starts[axis]][1].append(rectangle)
            faces[axis][ends[axis]][0].append(rectangle)
    pairs = 0
    violations = 0
    maximum_difference = 0
    for axis_faces in faces:
        for positive, negative in axis_faces.values():
            for first in positive:
                for second in negative:
                    overlap_u = min(first[1], second[1]) > max(first[0], second[0])
                    overlap_v = min(first[3], second[3]) > max(first[2], second[2])
                    if overlap_u and overlap_v:
                        difference = abs(first[4] - second[4])
                        pairs += 1
                        maximum_difference = max(maximum_difference, difference)
                        violations += difference > 1
    return {
        "faceNeighborPairs": pairs,
        "maximumLevelDifference": maximum_difference,
        "violatingFacePairCount": violations,
        "balanced": violations == 0,
    }


def parse_expected_point(text: str) -> tuple[np.ndarray, int]:
    try:
        coordinates, level = text.rsplit(":", 1)
        point = np.asarray([float(value) for value in coordinates.split(",")], dtype=float)
        if point.shape != (3,):
            raise ValueError
        return point, int(level)
    except ValueError as error:
        raise ValueError(f"invalid --expect-point {text!r}; expected x,y,z:level") from error


def verify(
    mesh_path: pathlib.Path,
    report_path: pathlib.Path,
    shape: str,
    expected_points: list[str],
    expected_classifications: list[str],
) -> dict[str, object]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    mesh = meshio.read(mesh_path)
    hexa_blocks = [block for block in mesh.cells if block.type == "hexahedron"]
    if len(hexa_blocks) != 1 or len(mesh.cells) != 1:
        raise ValueError("VTU must contain exactly one hexahedron block")
    hexahedra = np.asarray(hexa_blocks[0].data, dtype=np.int64)
    cell_count = len(hexahedra)
    if cell_count != report["leafCount"]:
        raise ValueError("VTU cell count differs from report leafCount")
    required = {
        "octree_level",
        "octree_node_code_low32",
        "octree_node_code_high32",
        "stl_cell_classification",
        "inside_stl_center_sample",
    }
    if not required.issubset(mesh.cell_data):
        raise ValueError(f"missing cell fields: {sorted(required - set(mesh.cell_data))}")

    def field(name: str) -> np.ndarray:
        blocks = mesh.cell_data[name]
        if len(blocks) != 1:
            raise ValueError(f"field {name} must have one block")
        return np.asarray(blocks[0]).reshape(-1)

    levels = field("octree_level").astype(np.int64)
    lows = field("octree_node_code_low32").astype(np.uint64)
    highs = field("octree_node_code_high32").astype(np.uint64)
    classifications = field("stl_cell_classification").astype(np.int64)
    codes = lows | (highs << np.uint64(32))
    domain_minimum = np.asarray(report["domain"]["minimum"], dtype=float)
    domain_maximum = np.asarray(report["domain"]["maximum"], dtype=float)
    domain_extent = domain_maximum - domain_minimum
    maximum_level = int(report["maximumLevel"])
    decoded: list[tuple[int, int, int, int]] = []
    anchors: list[int] = []
    geometry_mismatches: list[int] = []
    analytic_mismatches: list[int] = []
    volume_sum = 0.0
    inside_volume = 0.0
    inside_intersected_volume = 0.0
    physical_bounds: list[tuple[np.ndarray, np.ndarray]] = []
    surface_target_mismatches: list[int] = []
    surface_target_level = int(report["refinementRules"]["surfaceTargetLevel"])
    for cell_id, (code_value, connectivity) in enumerate(zip(codes.tolist(), hexahedra)):
        level, x, y, z = decode_node(int(code_value))
        decoded.append((level, x, y, z))
        if level != int(levels[cell_id]):
            geometry_mismatches.append(cell_id)
        shift = maximum_level - level
        morton = int(code_value) ^ (1 << (3 * level))
        anchors.append(morton << (3 * shift))
        expected_minimum = domain_minimum + domain_extent * np.array([x, y, z]) / (1 << level)
        expected_maximum = domain_minimum + domain_extent * np.array([x + 1, y + 1, z + 1]) / (1 << level)
        physical_bounds.append((expected_minimum, expected_maximum))
        coordinates = mesh.points[np.asarray(connectivity, dtype=np.int64)]
        if len(set(int(value) for value in connectivity)) != 8:
            geometry_mismatches.append(cell_id)
        actual_minimum = coordinates.min(axis=0)
        actual_maximum = coordinates.max(axis=0)
        if not (
            np.allclose(actual_minimum, expected_minimum, rtol=0.0, atol=2e-14)
            and np.allclose(actual_maximum, expected_maximum, rtol=0.0, atol=2e-14)
        ):
            geometry_mismatches.append(cell_id)
        expected_corners = np.asarray(
            [
                [expected_minimum[0], expected_minimum[1], expected_minimum[2]],
                [expected_maximum[0], expected_minimum[1], expected_minimum[2]],
                [expected_maximum[0], expected_maximum[1], expected_minimum[2]],
                [expected_minimum[0], expected_maximum[1], expected_minimum[2]],
                [expected_minimum[0], expected_minimum[1], expected_maximum[2]],
                [expected_maximum[0], expected_minimum[1], expected_maximum[2]],
                [expected_maximum[0], expected_maximum[1], expected_maximum[2]],
                [expected_minimum[0], expected_maximum[1], expected_maximum[2]],
            ]
        )
        if not np.allclose(coordinates, expected_corners, rtol=0.0, atol=2e-14):
            geometry_mismatches.append(cell_id)
        jacobian = float(
            np.dot(
                coordinates[1] - coordinates[0],
                np.cross(coordinates[3] - coordinates[0], coordinates[4] - coordinates[0]),
            )
        )
        if not jacobian > 0.0:
            geometry_mismatches.append(cell_id)
        volume = float(np.prod(expected_maximum - expected_minimum))
        volume_sum += volume
        if classifications[cell_id] == 1:
            inside_volume += volume
            inside_intersected_volume += volume
        elif classifications[cell_id] == 2:
            inside_intersected_volume += volume
        if shape == "cube":
            expected_classification = exact_cube_classification(expected_minimum, expected_maximum)
            if int(classifications[cell_id]) != expected_classification:
                analytic_mismatches.append(cell_id)
            if expected_classification == 2 and level < surface_target_level:
                surface_target_mismatches.append(cell_id)
    if len(set(map(int, codes))) != cell_count:
        raise ValueError("Morton node codes are not unique")
    if anchors != sorted(anchors):
        raise ValueError("leaf cells are not in deterministic Morton-anchor order")
    if not np.all(np.isin(classifications, [0, 1, 2, 3])):
        raise ValueError("classification contains values outside 0..3")
    expected_domain_volume = float(np.prod(domain_extent))
    if not math.isclose(volume_sum, expected_domain_volume, rel_tol=2e-12, abs_tol=2e-12):
        raise ValueError("adaptive leaf volumes do not partition the root domain")
    counts = {
        name: int(np.count_nonzero(classifications == value))
        for name, value in (("outside", 0), ("inside", 1), ("intersected", 2), ("conflict", 3))
    }
    if counts != report["classificationCounts"]:
        raise ValueError("classification counts differ from report")
    balance = face_balance(decoded, maximum_level)
    if report["faceBalance"]["enforced"] and not balance["balanced"]:
        raise ValueError("independent face 2:1 check failed")
    exact_volume = float(report["exactPolyhedralVolume"])
    volume_bracket_pass = inside_volume <= exact_volume <= inside_intersected_volume
    expected_point_results: list[dict[str, object]] = []
    for text in expected_points:
        point, expected_level = parse_expected_point(text)
        matches = [
            cell_id
            for cell_id, (minimum, maximum) in enumerate(physical_bounds)
            if np.all(point >= minimum) and np.all(point < maximum)
        ]
        actual_level = int(levels[matches[0]]) if len(matches) == 1 else None
        morton_code = int(codes[matches[0]]) if len(matches) == 1 else None
        expected_point_results.append(
            {
                "point": point.tolist(),
                "expectedLevel": expected_level,
                "matchingLeafCount": len(matches),
                "actualLevel": actual_level,
                "mortonCode": morton_code,
                "pass": len(matches) == 1 and actual_level == expected_level,
            }
        )
    expected_classification_results: list[dict[str, object]] = []
    for text in expected_classifications:
        point, expected_value = parse_expected_point(text)
        matches = [
            cell_id
            for cell_id, (minimum, maximum) in enumerate(physical_bounds)
            if np.all(point >= minimum) and np.all(point < maximum)
        ]
        actual_value = int(classifications[matches[0]]) if len(matches) == 1 else None
        morton_code = int(codes[matches[0]]) if len(matches) == 1 else None
        expected_classification_results.append(
            {
                "point": point.tolist(),
                "expectedClassification": expected_value,
                "matchingLeafCount": len(matches),
                "actualClassification": actual_value,
                "mortonCode": morton_code,
                "pass": len(matches) == 1 and actual_value == expected_value,
            }
        )
    rules = report["refinementRules"]
    boxes = rules.get("boxes")
    spheres = rules.get("spheres")
    cylinders = rules.get("cylinders")
    if not isinstance(boxes, list) or not isinstance(spheres, list) or not isinstance(cylinders, list):
        raise ValueError("refinementRules 必须完整保存 boxes/spheres/cylinders 数组")
    if (
        int(rules.get("boxCount", -1)) != len(boxes)
        or int(rules.get("sphereCount", -1)) != len(spheres)
        or int(rules.get("cylinderCount", -1)) != len(cylinders)
    ):
        raise ValueError("用户细化区数量与配置数组不一致")

    reported_region_checks: list[dict[str, object]] = []
    region_points: list[tuple[str, int, np.ndarray, int]] = []
    for index, region in enumerate(boxes):
        point = 0.5 * (
            np.asarray(region["minimum"], dtype=float)
            + np.asarray(region["maximum"], dtype=float)
        )
        region_points.append(("box", index, point, int(region["targetLevel"])))
    for index, region in enumerate(spheres):
        region_points.append(
            ("sphere", index, np.asarray(region["center"], dtype=float), int(region["targetLevel"]))
        )
    for index, region in enumerate(cylinders):
        point = 0.5 * (
            np.asarray(region["firstAxisPoint"], dtype=float)
            + np.asarray(region["secondAxisPoint"], dtype=float)
        )
        region_points.append(("cylinder", index, point, int(region["targetLevel"])))
    for kind, index, point, target_level in region_points:
        matches = [
            cell_id
            for cell_id, (minimum, maximum) in enumerate(physical_bounds)
            if np.all(point >= minimum) and np.all(point < maximum)
        ]
        actual_level = int(levels[matches[0]]) if len(matches) == 1 else None
        reported_region_checks.append(
            {
                "kind": kind,
                "index": index,
                "point": point.tolist(),
                "targetLevel": target_level,
                "matchingLeafCount": len(matches),
                "actualLevel": actual_level,
                "pass": len(matches) == 1 and actual_level is not None and actual_level >= target_level,
            }
        )
    passed = (
        not geometry_mismatches
        and not analytic_mismatches
        and not surface_target_mismatches
        and counts["conflict"] == 0
        and volume_bracket_pass
        and balance["balanced"]
        and all(result["pass"] for result in expected_point_results)
        and all(result["pass"] for result in expected_classification_results)
        and all(result["pass"] for result in reported_region_checks)
    )
    return {
        "schemaVersion": 1,
        "checker": "meshio_numpy_independent_stage2",
        "meshioVersion": meshio.__version__,
        "mesh": str(mesh_path),
        "report": str(report_path),
        "shape": shape,
        "pointCount": len(mesh.points),
        "leafCount": cell_count,
        "classificationCounts": counts,
        "uniqueMortonCodeCount": len(set(map(int, codes))),
        "mortonOrderStable": anchors == sorted(anchors),
        "domainVolume": expected_domain_volume,
        "leafVolumeSum": volume_sum,
        "insideVolumeLowerBound": inside_volume,
        "insideIntersectedVolumeUpperBound": inside_intersected_volume,
        "exactVolumeWithinBounds": volume_bracket_pass,
        "geometryMismatchCount": len(set(geometry_mismatches)),
        "geometryMismatchExamples": sorted(set(geometry_mismatches))[:16],
        "analyticClassificationMismatchCount": len(analytic_mismatches),
        "analyticClassificationMismatchExamples": analytic_mismatches[:16],
        "surfaceTargetMismatchCount": len(surface_target_mismatches),
        "surfaceTargetMismatchExamples": surface_target_mismatches[:16],
        "faceBalance": balance,
        "expectedPointChecks": expected_point_results,
        "expectedClassificationChecks": expected_classification_results,
        "reportedUserRegionChecks": reported_region_checks,
        "status": "pass" if passed else "fail",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", type=pathlib.Path)
    parser.add_argument("report", type=pathlib.Path)
    parser.add_argument("--shape", choices=("cube", "none"), default="none")
    parser.add_argument("--expect-point", action="append", default=[])
    parser.add_argument("--expect-classification", action="append", default=[])
    parser.add_argument("--output", type=pathlib.Path)
    arguments = parser.parse_args()
    try:
        result = verify(
            arguments.mesh,
            arguments.report,
            arguments.shape,
            arguments.expect_point,
            arguments.expect_classification,
        )
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        result = {"status": "fail", "error": str(error)}
    text = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0 if result.get("status") == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
