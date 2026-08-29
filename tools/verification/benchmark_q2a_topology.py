#!/usr/bin/env python3
"""Isolated common-partition benchmark, not an end-to-end mesher speed claim."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import statistics
import subprocess

p = argparse.ArgumentParser()
p.add_argument("--executable", type=Path, default=Path("build-q2a/cartmesh2d_shared_construction_tests"))
p.add_argument("--output", type=Path, default=Path("artifacts/q2a/topology-benchmark.json"))
a = p.parse_args()
result = {"scope": "isolated topology; uniform rectangles; shared includes exact interning and index construction",
          "order": "legacy then shared; three process runs per size", "platform": platform.platform(),
          "executable_sha256": hashlib.sha256(a.executable.read_bytes()).hexdigest(), "cases": []}
for n in (100, 316):
    runs = [json.loads(subprocess.check_output([str(a.executable), str(n)], text=True)) for _ in range(3)]
    assert all(r["identical"] for r in runs)
    entry = {"cells": n*n, "runs": runs,
             "legacy_median_seconds": statistics.median(r["legacy_seconds"] for r in runs),
             "shared_median_seconds": statistics.median(r["shared_seconds"] for r in runs)}
    result["cases"].append(entry)
    print(json.dumps(entry), flush=True)
a.output.parent.mkdir(parents=True, exist_ok=True)
a.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
