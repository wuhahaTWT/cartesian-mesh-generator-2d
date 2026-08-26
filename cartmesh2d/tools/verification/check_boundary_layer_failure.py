#!/usr/bin/env python3
"""Independently validate a structured H4-1 fail-closed JSON report."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path)
    parser.add_argument("--reason", required=True)
    args = parser.parse_args()
    try:
        data = json.loads(args.report.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"boundary-layer failure check: FAIL: {exc}")
        return 1
    issues: list[str] = []
    if data.get("layer_status") != "failed":
        issues.append("layer_status is not failed")
    if data.get("failure_reason") != args.reason:
        issues.append("failure_reason does not match expected value")
    if not isinstance(data.get("message"), str) or not data["message"]:
        issues.append("message is missing")
    if not isinstance(data.get("chain_id"), int):
        issues.append("chain_id is missing")
    for key in ("requested_thickness", "growth_ratio"):
        value = data.get(key)
        if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            issues.append(f"{key} is not finite")
    if not isinstance(data.get("n_layers"), int):
        issues.append("n_layers is missing")
    if issues:
        print("boundary-layer failure check: FAIL: " + "; ".join(issues))
        return 1
    print(
        "boundary-layer failure check: PASS "
        f"reason={data['failure_reason']} chain={data['chain_id']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
