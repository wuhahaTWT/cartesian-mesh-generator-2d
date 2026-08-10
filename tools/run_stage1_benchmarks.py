#!/usr/bin/env python3
"""重复运行阶段 1 精确单元分类，并记录独立墙钟、子进程 RSS 和解析真值。"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import pathlib
import platform
import statistics
import subprocess
import sys
import time
from typing import Any
from zoneinfo import ZoneInfo


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
        file_digest = hashlib.sha256(contents).hexdigest()
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(contents)
        digest.update(b"\0")
        entries.append({"path": relative, "bytes": len(contents), "sha256": file_digest})
    return {"sha256": digest.hexdigest(), "fileCount": len(entries), "files": entries}


def axis_counts(
    minimum: float, maximum: float, count: int, solid_minimum: float, solid_maximum: float
) -> tuple[int, int, int]:
    spacing = (maximum - minimum) / count
    overlap = 0
    strict_interior = 0
    center_inside = 0
    for index in range(count):
        lower = minimum + spacing * index
        upper = minimum + spacing * (index + 1)
        overlaps_solid = lower <= solid_maximum and upper >= solid_minimum
        touches_boundary = (
            lower <= solid_minimum <= upper or lower <= solid_maximum <= upper
        )
        overlap += int(overlaps_solid)
        strict_interior += int(overlaps_solid and not touches_boundary)
        center = 0.5 * (lower + upper)
        center_inside += int(solid_minimum <= center <= solid_maximum)
    return overlap, strict_interior, center_inside


def analytic_cube_counts(report: dict[str, Any]) -> dict[str, int]:
    dimensions = report["dimensions"]
    domain = report["domain"]
    minima = domain["minimum"]
    maxima = domain["maximum"]
    axis = [
        axis_counts(minima[i], maxima[i], dimensions[name], 0.0, 1.0)
        for i, name in enumerate(("nx", "ny", "nz"))
    ]
    overlap = axis[0][0] * axis[1][0] * axis[2][0]
    inside = axis[0][1] * axis[1][1] * axis[2][1]
    intersected = overlap - inside
    total = dimensions["nx"] * dimensions["ny"] * dimensions["nz"]
    center_inside = axis[0][2] * axis[1][2] * axis[2][2]
    return {
        "outside": total - inside - intersected,
        "inside": inside,
        "intersected": intersected,
        "conflict": 0,
        "centerInside": center_inside,
    }


def wait4_run(command: list[str], stdout_path: pathlib.Path, stderr_path: pathlib.Path) -> dict[str, Any]:
    if not hasattr(os, "wait4") or not hasattr(os, "posix_spawn"):
        raise RuntimeError("该外部基准需要 POSIX posix_spawn/wait4 以取得子进程资源用量")
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("wb") as stdout_file, stderr_path.open("wb") as stderr_file:
        actions = [
            (os.POSIX_SPAWN_DUP2, stdout_file.fileno(), 1),
            (os.POSIX_SPAWN_DUP2, stderr_file.fileno(), 2),
        ]
        start = time.perf_counter()
        pid = os.posix_spawn(command[0], command, os.environ.copy(), file_actions=actions)
        waited_pid, status, usage = os.wait4(pid, 0)
        wall = time.perf_counter() - start
    if waited_pid != pid:
        raise RuntimeError("wait4 返回了非预期子进程")
    exit_code = os.waitstatus_to_exitcode(status)
    rss = int(usage.ru_maxrss)
    if sys.platform != "darwin":
        rss *= 1024
    return {
        "exitCode": exit_code,
        "wallSeconds": wall,
        "userSeconds": float(usage.ru_utime),
        "systemSeconds": float(usage.ru_stime),
        "maximumResidentSetBytes": rss,
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
    }


def cmake_version() -> str:
    completed = subprocess.run(
        ["cmake", "--version"], check=True, capture_output=True, text=True
    )
    return completed.stdout.splitlines()[0].removeprefix("cmake version ")


def execute_case(
    root: pathlib.Path,
    executable: pathlib.Path,
    surface: pathlib.Path,
    output_directory: pathlib.Path,
    resolution: int,
    label: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    stem = f"stage1_cube_{resolution}_exact_{label}"
    report_path = output_directory / f"{stem}.json"
    command = [
        str(executable),
        "--stl",
        str(surface),
        "--resolution",
        str(resolution),
        "--no-vtk",
        "--report",
        str(report_path),
    ]
    external = wait4_run(
        command,
        output_directory / f"{stem}.stdout.txt",
        output_directory / f"{stem}.stderr.txt",
    )
    if external["exitCode"] != 0:
        raise RuntimeError(
            f"基准子进程失败：resolution={resolution}, label={label}, "
            f"exit={external['exitCode']}"
        )
    project = json.loads(report_path.read_text(encoding="utf-8"))
    expected = analytic_cube_counts(project)
    if project["classificationCounts"] != {
        name: expected[name] for name in ("outside", "inside", "intersected", "conflict")
    }:
        raise RuntimeError(
            f"{resolution}^3 分类与独立解析计数不一致："
            f"{project['classificationCounts']} != {expected}"
        )
    if project["centerPointCounts"]["inside"] != expected["centerInside"]:
        raise RuntimeError(f"{resolution}^3 中心采样计数与独立解析计数不一致")
    if project["status"] != "pass" or project["solverReadyCutCellMesh"] is not False:
        raise RuntimeError("阶段 1 报告状态或非 Cut-cell 边界声明不正确")
    external["projectReport"] = str(report_path)
    external["internalTotalSeconds"] = project["timingsSeconds"]["total"]
    external["internalClassificationSeconds"] = project["timingsSeconds"]["classification"]
    external["internalPeakRssBytes"] = project["peakRssBytes"]
    external["resultHashFnv1a64"] = project["resultHashFnv1a64"]
    return external, project


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--executable", type=pathlib.Path, default=pathlib.Path("build/release/cartmesh_cli"))
    parser.add_argument("--surface", type=pathlib.Path, default=pathlib.Path("tests/data/closed_unit_cube_ascii.stl"))
    parser.add_argument("--output-dir", type=pathlib.Path, default=pathlib.Path("artifacts/stage1_benchmarks"))
    parser.add_argument("--summary", type=pathlib.Path)
    parser.add_argument("--resolutions", type=int, nargs="+", default=(100, 216))
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--hardware-model", default="unknown")
    parser.add_argument("--hardware-chip", default="unknown")
    parser.add_argument("--physical-cores", type=int, default=0)
    parser.add_argument("--logical-cores", type=int, default=0)
    parser.add_argument("--memory-bytes", type=int, default=0)
    arguments = parser.parse_args()

    timezone_name = "Asia/Shanghai"
    measurement_date = datetime.datetime.now(ZoneInfo(timezone_name)).date().isoformat()
    root = arguments.root.resolve()
    executable = (root / arguments.executable).resolve() if not arguments.executable.is_absolute() else arguments.executable.resolve()
    surface = (root / arguments.surface).resolve() if not arguments.surface.is_absolute() else arguments.surface.resolve()
    output_directory = (root / arguments.output_dir).resolve() if not arguments.output_dir.is_absolute() else arguments.output_dir.resolve()
    requested_summary = arguments.summary or pathlib.Path(
        f"benchmarks/baselines/stage1_uniform_exact_m1_{measurement_date}.json"
    )
    summary_path = (
        (root / requested_summary).resolve()
        if not requested_summary.is_absolute()
        else requested_summary.resolve()
    )
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise SystemExit(f"Release 可执行文件不存在或不可执行：{executable}")
    if not surface.is_file():
        raise SystemExit(f"STL 输入不存在：{surface}")
    if arguments.warmups < 0 or arguments.repeats < 1:
        raise SystemExit("warmups 必须非负且 repeats 至少为 1")
    if any(value < 1 for value in arguments.resolutions):
        raise SystemExit("分辨率必须为正整数")

    output_directory.mkdir(parents=True, exist_ok=True)
    cases: list[dict[str, Any]] = []
    compiler = None
    for resolution in arguments.resolutions:
        for warmup in range(arguments.warmups):
            execute_case(
                root,
                executable,
                surface,
                output_directory,
                resolution,
                f"warmup{warmup + 1}",
            )
        runs: list[dict[str, Any]] = []
        representative: dict[str, Any] | None = None
        for repeat in range(arguments.repeats):
            external, project = execute_case(
                root,
                executable,
                surface,
                output_directory,
                resolution,
                f"run{repeat + 1}",
            )
            runs.append(external)
            representative = project
            compiler = project["compiler"]
        assert representative is not None
        hashes = {run["resultHashFnv1a64"] for run in runs}
        if len(hashes) != 1:
            raise RuntimeError(f"{resolution}^3 正式重复运行的结果哈希不一致")
        expected = analytic_cube_counts(representative)
        cases.append(
            {
                "resolution": resolution,
                "dimensions": representative["dimensions"],
                "domain": representative["domain"],
                "cellSpacing": representative["cellSpacing"],
                "cellCount": representative["cellCount"],
                "classificationCounts": representative["classificationCounts"],
                "centerPointCounts": representative["centerPointCounts"],
                "analyticExpected": expected,
                "analyticCountsMatch": True,
                "bvh": representative["bvh"],
                "runs": runs,
                "summary": {
                    "externalWallMedianSeconds": statistics.median(
                        run["wallSeconds"] for run in runs
                    ),
                    "externalWallMinimumSeconds": min(run["wallSeconds"] for run in runs),
                    "externalWallMaximumSeconds": max(run["wallSeconds"] for run in runs),
                    "externalMaximumRssBytes": max(
                        run["maximumResidentSetBytes"] for run in runs
                    ),
                    "resultHashFnv1a64": next(iter(hashes)),
                },
            }
        )

    manifest = source_manifest(root)
    mac_version = platform.mac_ver()[0]
    payload = {
        "schemaVersion": 3,
        "projectStage": 1,
        "date": measurement_date,
        "timezone": timezone_name,
        "benchmark": "uniform Cartesian exact triangle-AABB surface intersection plus center parity",
        "sourceRevision": {"gitRepository": False, "sourceManifest": manifest},
        "hardware": {
            "model": arguments.hardware_model,
            "chip": arguments.hardware_chip,
            "physicalCores": arguments.physical_cores,
            "logicalCores": arguments.logical_cores,
            "memoryBytes": arguments.memory_bytes,
        },
        "operatingSystem": {
            "name": "macOS" if mac_version else platform.system(),
            "version": mac_version or platform.release(),
            "architecture": platform.machine(),
        },
        "build": {
            "type": "Release",
            "compiler": compiler,
            "cmake": cmake_version(),
            "sanitizers": False,
            "runtimeThreads": 1,
        },
        "input": {
            "path": surface.relative_to(root).as_posix() if surface.is_relative_to(root) else str(surface),
            "format": "ascii_stl",
            "sha256": sha256_file(surface),
            "triangleCount": 12,
            "bounds": {"minimum": [0.0, 0.0, 0.0], "maximum": [1.0, 1.0, 1.0]},
        },
        "warmupRunsPerCase": arguments.warmups,
        "measuredRunsPerCase": arguments.repeats,
        "cases": cases,
        "outputEnabled": False,
        "solverReadyCutCellMesh": False,
        "notes": [
            "相交标签来自 BVH 粗筛后的精确 triangle-AABB SAT；inside/outside 来自非相交单元中心奇偶射线。",
            "这是均匀背景单元分类，不是 Cut-cell 或求解器可用网格。",
            "外部资源数据由独立 Python 父进程的 POSIX wait4 采集；macOS ru_maxrss 单位为字节。",
            "关闭 VTU 只测分类核心；小规模 VTU 另由 meshio 和 ParaView 验证。",
        ],
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
