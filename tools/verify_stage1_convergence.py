#!/usr/bin/env python3
"""独立复核阶段 1 解析 STL 的中心采样误差与分类体积包络收敛。"""

from __future__ import annotations

import argparse
import json
import pathlib
from typing import Any


def load_report(path: pathlib.Path) -> dict[str, Any]:
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("projectStage") != 1 or report.get("status") != "pass":
        raise ValueError(f"不是通过的阶段 1 报告：{path}")
    counts = report["classificationCounts"]
    if sum(int(counts[name]) for name in ("outside", "inside", "intersected", "conflict")) != int(
        report["cellCount"]
    ):
        raise ValueError(f"四类计数未完整分割网格：{path}")
    if int(counts["conflict"]) != 0:
        raise ValueError(f"解析案例存在分类冲突：{path}")
    if report.get("centerSampleVolumeValid") is not True:
        raise ValueError(f"中心采样体积未标记为有效：{path}")
    return report


def metrics(path: pathlib.Path, report: dict[str, Any]) -> dict[str, Any]:
    exact = float(report["exactPolyhedralVolume"])
    estimate = float(report["estimatedCenterSampleVolume"])
    lower = float(report["definitelyInsideVolumeLowerBound"])
    upper = float(report["insidePlusIntersectedVolumeUpperBound"])
    if not lower <= exact <= upper:
        raise ValueError(f"多面体体积没有落在分类包络中：{path}")
    if report.get("polyhedralVolumeInsideClassificationBounds") is not True:
        raise ValueError(f"报告自身未声明体积包络通过：{path}")
    return {
        "report": str(path),
        "dimensions": report["dimensions"],
        "cellCount": report["cellCount"],
        "classificationCounts": report["classificationCounts"],
        "centerPointCounts": report["centerPointCounts"],
        "exactPolyhedralVolume": exact,
        "estimatedCenterSampleVolume": estimate,
        "centerSampleAbsoluteError": abs(estimate - exact),
        "definitelyInsideVolumeLowerBound": lower,
        "insidePlusIntersectedVolumeUpperBound": upper,
        "classificationVolumeBracketWidth": upper - lower,
        "resultHashFnv1a64": report["resultHashFnv1a64"],
    }


def verify_case(name: str, coarse_path: pathlib.Path, fine_path: pathlib.Path) -> dict[str, Any]:
    coarse_report = load_report(coarse_path)
    fine_report = load_report(fine_path)
    if coarse_report["input"] != fine_report["input"]:
        raise ValueError(f"{name} 的粗细网格不是同一 STL 输入")
    coarse = metrics(coarse_path, coarse_report)
    fine = metrics(fine_path, fine_report)
    exact_difference = abs(
        coarse["exactPolyhedralVolume"] - fine["exactPolyhedralVolume"]
    )
    exact_scale = max(1.0, abs(coarse["exactPolyhedralVolume"]))
    if exact_difference > 1.0e-13 * exact_scale:
        raise ValueError(f"{name} 的粗细网格多面体体积不一致")
    if fine["centerSampleAbsoluteError"] >= coarse["centerSampleAbsoluteError"]:
        raise ValueError(f"{name} 的中心采样误差没有随加密下降")
    if (
        fine["classificationVolumeBracketWidth"]
        >= coarse["classificationVolumeBracketWidth"]
    ):
        raise ValueError(f"{name} 的 inside/intersected 体积包络没有随加密收缩")
    return {
        "name": name,
        "input": coarse_report["input"],
        "coarse": coarse,
        "fine": fine,
        "centerSampleErrorDecreased": True,
        "classificationVolumeBracketShrank": True,
        "exactVolumeInsideBothBrackets": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        nargs=3,
        required=True,
        metavar=("NAME", "COARSE_REPORT", "FINE_REPORT"),
        help="可重复指定名称、粗网格 JSON 和细网格 JSON",
    )
    parser.add_argument("--output", required=True, type=pathlib.Path)
    arguments = parser.parse_args()
    try:
        cases = [
            verify_case(name, pathlib.Path(coarse), pathlib.Path(fine))
            for name, coarse, fine in arguments.case
        ]
        result = {
            "schemaVersion": 1,
            "projectStage": 1,
            "checker": "independent JSON arithmetic and convergence invariants",
            "cases": cases,
            "solverReadyCutCellMesh": False,
            "status": "pass",
        }
        serialized = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True)
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(serialized + "\n", encoding="utf-8")
        print(serialized)
        return 0
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(json.dumps({"status": "fail", "error": str(error)}, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
