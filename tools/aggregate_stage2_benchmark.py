#!/usr/bin/env python3
"""Aggregate three /usr/bin/time Stage-2 ten-million-leaf benchmark runs."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import pathlib
import platform
import re
import statistics
from typing import Any


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_manifest(root: pathlib.Path) -> dict[str, Any]:
    patterns = (
        "AGENTS.md",
        "CARTESIAN_MESH_GENERATOR_PROJECT_BRIEF_CN.md",
        "CMakeLists.txt",
        "CMakePresets.json",
        "apps/**/*.cpp",
        "include/**/*.hpp",
        "src/**/*.cpp",
        "tests/**/*.cpp",
        "tests/data/**/*.stl",
        "tools/**/*.py",
    )
    paths: set[pathlib.Path] = set()
    for pattern in patterns:
        paths.update(path for path in root.glob(pattern) if path.is_file())
    digest = hashlib.sha256()
    entries: list[dict[str, Any]] = []
    for path in sorted(paths):
        relative = path.relative_to(root).as_posix()
        contents = path.read_bytes()
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(contents)
        digest.update(b"\0")
        entries.append(
            {
                "path": relative,
                "bytes": len(contents),
                "sha256": hashlib.sha256(contents).hexdigest(),
            }
        )
    return {"sha256": digest.hexdigest(), "fileCount": len(entries), "files": entries}


def parse_time(path: pathlib.Path) -> dict[str, float | int]:
    text = path.read_text(encoding="utf-8")
    timing = re.search(
        r"\s*([0-9.]+) real\s+([0-9.]+) user\s+([0-9.]+) sys", text
    )
    rss = re.search(r"\s*(\d+)\s+maximum resident set size", text)
    footprint = re.search(r"\s*(\d+)\s+peak memory footprint", text)
    if timing is None or rss is None:
        raise ValueError(f"cannot parse macOS time -l output: {path}")
    return {
        "wallSeconds": float(timing.group(1)),
        "userSeconds": float(timing.group(2)),
        "systemSeconds": float(timing.group(3)),
        "maximumResidentSetBytes": int(rss.group(1)),
        "peakMemoryFootprintBytes": int(footprint.group(1)) if footprint else 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path("."))
    parser.add_argument("--executable", type=pathlib.Path, required=True)
    parser.add_argument("--run-json", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--time-log", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--hardware-model", required=True)
    parser.add_argument("--chip", required=True)
    parser.add_argument("--physical-cores", type=int, required=True)
    parser.add_argument("--logical-cores", type=int, required=True)
    parser.add_argument("--memory-bytes", type=int, required=True)
    parser.add_argument("--os-version", required=True)
    parser.add_argument("--compiler", required=True)
    parser.add_argument("--cmake-version", required=True)
    arguments = parser.parse_args()
    if len(arguments.run_json) != len(arguments.time_log) or len(arguments.run_json) < 3:
        raise ValueError("provide matching run/time paths for at least three measured runs")
    root = arguments.root.resolve()
    runs: list[dict[str, Any]] = []
    for index, (json_path, time_path) in enumerate(
        zip(arguments.run_json, arguments.time_log), start=1
    ):
        internal = json.loads(json_path.read_text(encoding="utf-8"))
        external = parse_time(time_path)
        runs.append(
            {
                "run": index,
                "internal": internal,
                "external": external,
                "rawReport": json_path.as_posix(),
                "rawTimeLog": time_path.as_posix(),
            }
        )
    hashes = {run["internal"]["resultHashFnv1a64Decimal"] for run in runs}
    leaf_counts = {run["internal"]["finalLeafCount"] for run in runs}
    if len(hashes) != 1 or len(leaf_counts) != 1:
        raise ValueError("measured Stage-2 runs are not deterministic")
    wall = [run["external"]["wallSeconds"] for run in runs]
    rss = [run["external"]["maximumResidentSetBytes"] for run in runs]
    manifest = source_manifest(root)
    result = {
        "schemaVersion": 1,
        "projectStage": 2,
        "benchmark": "compact_linear_octree_ten_million_leaf_construction",
        "date": dt.datetime.now().astimezone().isoformat(),
        "solverReadyCutCellMesh": False,
        "meshKind": "adaptive_cartesian_background_leaf_codes_only",
        "sourceRevision": {"gitRepository": False, "sourceManifest": manifest},
        "environment": {
            "operatingSystem": {
                "name": "macOS",
                "version": arguments.os_version,
                "architecture": platform.machine(),
            },
            "hardware": {
                "model": arguments.hardware_model,
                "chip": arguments.chip,
                "physicalCores": arguments.physical_cores,
                "logicalCores": arguments.logical_cores,
                "memoryBytes": arguments.memory_bytes,
            },
            "build": {
                "type": "Release",
                "compiler": arguments.compiler,
                "cmake": arguments.cmake_version,
                "sanitizers": False,
                "runtimeThreads": 1,
            },
        },
        "executable": {
            "path": arguments.executable.as_posix(),
            "sha256": sha256_file(arguments.executable),
        },
        "measuredRuns": runs,
        "summary": {
            "runCount": len(runs),
            "requestedLeafCount": runs[0]["internal"]["requestedLeafCount"],
            "finalLeafCount": runs[0]["internal"]["finalLeafCount"],
            "compactLeafStorageBytes": runs[0]["internal"]["compactLeafStorageBytes"],
            "resultHashFnv1a64Decimal": runs[0]["internal"]["resultHashFnv1a64Decimal"],
            "wallSecondsMedian": statistics.median(wall),
            "wallSecondsMinimum": min(wall),
            "wallSecondsMaximum": max(wall),
            "maximumResidentSetBytesMaximum": max(rss),
            "partitionValidEveryRun": all(
                run["internal"]["partitionValid"] for run in runs
            ),
            "faceBalanced2To1EveryRun": all(
                run["internal"]["faceBalanced2To1"] for run in runs
            ),
        },
        "notes": [
            "The benchmark stores one uint64 node code per leaf in the compact leaf array.",
            "The measured tree has only two adjacent levels, so global level span proves face 2:1 before the generic neighbor tests remain independently covered by unit and meshio checks.",
            "This is a background octree memory/construction benchmark, not Cut-cell generation or solver-ready mesh output.",
        ],
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result["summary"], ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
