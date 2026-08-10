#!/usr/bin/env python3
"""复核重复运行及 ASCII/二进制 STL 的阶段 1 VTU 与分类哈希一致性。"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run",
        action="append",
        nargs=3,
        required=True,
        metavar=("LABEL", "VTU", "REPORT"),
    )
    parser.add_argument("--output", required=True, type=pathlib.Path)
    arguments = parser.parse_args()
    try:
        runs = []
        for label, vtu_value, report_value in arguments.run:
            vtu = pathlib.Path(vtu_value)
            report_path = pathlib.Path(report_value)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            if report.get("projectStage") != 1 or report.get("status") != "pass":
                raise ValueError(f"{label} 不是通过的阶段 1 报告")
            runs.append(
                {
                    "label": label,
                    "input": report["input"],
                    "inputFormat": report["inputFormat"],
                    "vtu": str(vtu),
                    "vtuBytes": vtu.stat().st_size,
                    "vtuSha256": sha256(vtu),
                    "report": str(report_path),
                    "resultHashFnv1a64": report["resultHashFnv1a64"],
                    "classificationCounts": report["classificationCounts"],
                    "centerPointCounts": report["centerPointCounts"],
                }
            )
        if len(runs) < 2:
            raise ValueError("至少需要两次运行才能验证确定性")
        for field in (
            "vtuBytes",
            "vtuSha256",
            "resultHashFnv1a64",
            "classificationCounts",
            "centerPointCounts",
        ):
            if any(run[field] != runs[0][field] for run in runs[1:]):
                raise ValueError(f"不同运行的 {field} 不一致")
        result = {
            "schemaVersion": 1,
            "projectStage": 1,
            "runs": runs,
            "vtuByteIdentical": True,
            "semanticClassificationHashIdentical": True,
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
