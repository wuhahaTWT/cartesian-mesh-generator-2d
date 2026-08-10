#!/usr/bin/env python3
"""汇总三类增量案例、独立读取和性能记录，生成阶段 5 终态契约。"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as source:
        return json.load(source)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", type=Path, required=True)
    parser.add_argument("--external", action="append", type=Path, required=True)
    parser.add_argument("--performance", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if len(args.case) != 3 or len(args.external) != 3:
        raise SystemExit("exactly three --case and three --external inputs are required")

    cases = [load(path) for path in args.case]
    external = [load(path) for path in args.external]
    performance = load(args.performance)
    failures: list[str] = []
    case_summaries = []
    for index, (case, check) in enumerate(zip(cases, external, strict=True)):
        if case.get("status") != "internal_pass_external_pending":
            failures.append(f"case_{index}_internal")
        if not case.get("topologyEqualsFullRebuild"):
            failures.append(f"case_{index}_topology")
        if not case.get("geometryEqualsFullRebuild"):
            failures.append(f"case_{index}_geometry")
        if case.get("geometryReusedLeafCount", 0) <= 0:
            failures.append(f"case_{index}_reuse")
        for field in ("nonclosedCellCount", "negativeVolumeCellCount",
                      "sharedFaceMismatchCount", "classificationConflictCount"):
            if case.get(field) != 0:
                failures.append(f"case_{index}_{field}")
        if not check.get("stage5CaseExternalPass"):
            failures.append(f"case_{index}_external")
        if not check.get("externalBackgroundFieldsEqualFullRebuild"):
            failures.append(f"case_{index}_external_equivalence")
        case_summaries.append({
            "oldInput": case["oldInput"],
            "newInput": case["newInput"],
            "oldInputHashFnv1a64": case["oldInputHashFnv1a64"],
            "newInputHashFnv1a64": case["newInputHashFnv1a64"],
            "oldLeafCount": case["oldLeafCount"],
            "newLeafCount": case["newLeafCount"],
            "geometryReuseFraction": case["geometryReuseFraction"],
            "resultHashFnv1a64": case["incrementalResultHashFnv1a64"],
            "measuredSpeedup": case["measuredSpeedup"],
            "externalReader": check["reader"],
        })
    if not performance.get("stage5PerformancePass"):
        failures.append("performance")

    passed = not failures
    result = {
        "schema": "cartmesh-stage5-acceptance-v1",
        "status": "pass" if passed else "fail",
        "stage5Complete": passed,
        "externalIndependentReaderAccepted": passed,
        "incrementalResultsEquivalentToFullRebuild": passed,
        "stableCellIdsVerified": passed,
        "exactFluidOverlapMappingVerified": passed,
        "performanceEvidenceComplete": passed,
        "caseCount": len(case_summaries),
        "cases": case_summaries,
        "performance": {
            "runCount": performance["runCount"],
            "newLeafCount": performance["newLeafCount"],
            "geometryReuseFraction": performance["geometryReuseFraction"],
            "incrementalSeconds": performance["incrementalSeconds"],
            "fullRebuildSeconds": performance["fullRebuildSeconds"],
            "speedup": performance["speedup"],
            "maximumObservedPeakRssBytes": performance["maximumObservedPeakRssBytes"],
            "peakRssScope": performance["peakRssScope"],
            "runtimeThreads": performance["runtimeThreads"],
            "hardwareModel": performance["hardwareModel"],
            "hardwareMemoryBytes": performance["hardwareMemoryBytes"],
            "hardwareLogicalCpu": performance["hardwareLogicalCpu"],
            "buildType": performance["buildType"],
            "compiler": performance["compiler"],
        },
        "acceptanceBlockers": failures,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        json.dump(result, output, ensure_ascii=False, indent=2, sort_keys=True)
        output.write("\n")
    raise SystemExit(0 if passed else 2)


if __name__ == "__main__":
    main()
