#!/usr/bin/env python3
"""汇总阶段 5 重复性能测量，不隐藏慢运行或联合峰值 RSS 范围。"""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path
from typing import Any


def load(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as source:
        return json.load(source)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reports", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    reports = [load(path) for path in args.reports]
    failures: list[str] = []
    if len(reports) < 3:
        failures.append("fewer_than_three_runs")

    reference_fields = (
        "oldLeafCount", "newLeafCount", "geometryReusedLeafCount",
        "geometryRebuiltLeafCount", "incrementalResultHashFnv1a64",
        "fullRebuildResultHashFnv1a64", "hardwareModel",
        "hardwareMemoryBytes", "hardwareLogicalCpu", "buildType",
        "compiler", "runtimeThreads",
    )
    for field in reference_fields:
        values = {json.dumps(report.get(field), sort_keys=True) for report in reports}
        if len(values) != 1:
            failures.append(f"inconsistent_{field}")
    for index, report in enumerate(reports):
        if report.get("status") != "internal_pass_external_pending":
            failures.append(f"run_{index}_status")
        if not report.get("topologyEqualsFullRebuild"):
            failures.append(f"run_{index}_topology")
        if not report.get("geometryEqualsFullRebuild"):
            failures.append(f"run_{index}_geometry")
        for field in ("nonclosedCellCount", "negativeVolumeCellCount",
                      "sharedFaceMismatchCount", "classificationConflictCount"):
            if report.get(field) != 0:
                failures.append(f"run_{index}_{field}")

    incremental = [float(report["incrementalTotalSeconds"]) for report in reports]
    full = [float(report["fullRebuildSeconds"]) for report in reports]
    speedup = [float(report["measuredSpeedup"]) for report in reports]
    peak_rss = [int(report["peakRssBytes"]) for report in reports]
    median_speedup = statistics.median(speedup)
    non_micro = reports[0]["newLeafCount"] >= 5000
    performance_pass = not failures and non_micro and median_speedup > 1.0

    result = {
        "schema": "cartmesh-stage5-benchmark-v1",
        "status": "pass" if performance_pass else "fail",
        "stage5PerformancePass": performance_pass,
        "runCount": len(reports),
        "oldLeafCount": reports[0]["oldLeafCount"],
        "newLeafCount": reports[0]["newLeafCount"],
        "geometryReusedLeafCount": reports[0]["geometryReusedLeafCount"],
        "geometryRebuiltLeafCount": reports[0]["geometryRebuiltLeafCount"],
        "geometryReuseFraction": reports[0]["geometryReuseFraction"],
        "incrementalSeconds": {
            "minimum": min(incremental),
            "median": statistics.median(incremental),
            "maximum": max(incremental),
        },
        "fullRebuildSeconds": {
            "minimum": min(full),
            "median": statistics.median(full),
            "maximum": max(full),
        },
        "speedup": {
            "minimum": min(speedup),
            "median": median_speedup,
            "maximum": max(speedup),
        },
        "maximumObservedPeakRssBytes": max(peak_rss),
        "peakRssScope": reports[0]["peakRssScope"],
        "runtimeThreads": reports[0]["runtimeThreads"],
        "hardwareModel": reports[0]["hardwareModel"],
        "hardwareMemoryBytes": reports[0]["hardwareMemoryBytes"],
        "hardwareLogicalCpu": reports[0]["hardwareLogicalCpu"],
        "buildType": reports[0]["buildType"],
        "compiler": reports[0]["compiler"],
        "resultHashFnv1a64": reports[0]["incrementalResultHashFnv1a64"],
        "failures": failures,
        "notes": [
            "Each process built the old baseline, then separately timed incremental and full-new paths.",
            "Peak RSS is the conservative maximum for old baseline plus incremental plus full validation plus export.",
            "Export time is excluded from both compared compute timings and reported separately by each run.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        json.dump(result, output, ensure_ascii=False, indent=2, sort_keys=True)
        output.write("\n")
    raise SystemExit(0 if performance_pass else 2)


if __name__ == "__main__":
    main()
