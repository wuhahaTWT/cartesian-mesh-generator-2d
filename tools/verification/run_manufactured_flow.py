#!/usr/bin/env python3
"""Reproducibly run the native manufactured-flow refinement ladder.

The default output directory is deliberately separate from historical
``outputs/native-flow/manufactured`` artifacts.  A normal invocation creates
fresh meshes with cartmesh2d_cli, optionally applies the boundary-preserving
skew map, runs the flow CLI, and then calls verify_native_flow.verify_case.
``--reuse`` never launches either native executable; it verifies artifacts
already present in the selected output directory and labels their provenance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shlex
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any


REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = REPO / "outputs/native-flow/manufactured-repro"
LEVELS = {
    4: (14, 196),
    5: (30, 900),
    6: (62, 3844),
}
KINDS = ("cartesian", "warped")
SCHEMES = ("upwind", "limited-linear")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_text(command: list[str]) -> str:
    return shlex.join(command)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
                    encoding="utf-8")


def square_geometry(root: Path) -> Path:
    path = root / "geometry" / "unit-square.xy"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("0 0\n1 0\n1 1\n0 1\n", encoding="utf-8")
    return path


def mesh_command(mesh_cli: Path, geometry: Path, prefix: Path, level: int,
                 padding: float, openfoam_case: Path) -> list[str]:
    return [str(mesh_cli), str(geometry), str(prefix), str(level),
            f"{padding:.17g}", "0.1", "interior", str(openfoam_case),
            str(level), "0"]


def expected_mesh_path(prefix: Path) -> Path:
    return Path(str(prefix) + ".solver.cm2d")


def warp_mesh(source: Path, destination: Path) -> dict[str, Any]:
    """Apply the fixed boundary-preserving map used by the manufactured study."""
    sys.path.insert(0, str(REPO / "tools" / "verification"))
    from verify_native_flow import read_cm2d  # pylint: disable=import-outside-toplevel

    mesh = read_cm2d(source)
    vertices: list[tuple[float, float]] = []
    for x, y in mesh.vertices:
        offset = 0.06 * math.sin(math.pi * x) * math.sin(math.pi * y)
        vertices.append((x + offset, y + 0.6 * offset))

    lines = ["CM2D 1", f"VERTICES {len(vertices)}"]
    lines.extend(f"{index} {x:.17g} {y:.17g}"
                 for index, (x, y) in enumerate(vertices))
    lines.append(f"EDGES {len(mesh.edges)}")
    lines.extend(f"{edge.id} {edge.v0} {edge.v1} {edge.owner} "
                 f"{edge.neighbour} {edge.patch}" for edge in mesh.edges)
    lines.append(f"CELLS {len(mesh.cells)}")
    for cell in mesh.cells:
        points = [vertices[index] for index in cell.vertices]
        area = math.fsum(a[0] * b[1] - b[0] * a[1]
                         for a, b in zip(points, points[1:] + points[:1])) / 2.0
        if not math.isfinite(area) or area <= 0.0:
            raise ValueError(f"warped cell {cell.id} has non-positive area {area}")
        values = [cell.id, 0, 0, format(area, ".17g"), len(cell.vertices),
                  *cell.vertices, len(cell.edges), *cell.edges]
        lines.append(" ".join(map(str, values)))
    lines.extend(("AUDIT 0 0 0 0 0 0 0", "END"))
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {"source": str(source.resolve()), "sourceSha256": sha256_file(source),
            "map": "x += 0.06*sin(pi*x)*sin(pi*y); y += 0.6*offset",
            "destination": str(destination.resolve()),
            "destinationSha256": sha256_file(destination)}


def run_command(command: list[str], log: Path, timeout: int) -> dict[str, Any]:
    log.parent.mkdir(parents=True, exist_ok=True)
    stdout_path = log.with_suffix(".stdout.log")
    stderr_path = log.with_suffix(".stderr.log")
    record: dict[str, Any] = {
        "command": command,
        "commandText": command_text(command),
        "returncode": None,
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
        "timedOut": False,
    }
    start = time.monotonic()
    try:
        with stdout_path.open("w", encoding="utf-8") as stdout, \
                stderr_path.open("w", encoding="utf-8") as stderr:
            completed = subprocess.run(command, stdout=stdout, stderr=stderr,
                                       timeout=timeout, check=False)
        record["returncode"] = completed.returncode
        record["status"] = "passed" if completed.returncode == 0 else "failed"
    except subprocess.TimeoutExpired:
        # The redirected files contain the partial native logs.  Preserve
        # them; replacing them here would erase the most useful failure data.
        record.update(status="failed", timedOut=True, error=f"timed out after {timeout}s")
    except OSError as exc:
        stdout_path.write_text("", encoding="utf-8")
        stderr_path.write_text(str(exc) + "\n", encoding="utf-8")
        record.update(status="failed", error=str(exc))
    record["wallSeconds"] = time.monotonic() - start
    return record


def verifier_namespace(args: argparse.Namespace) -> SimpleNamespace:
    # These are the only verifier controls used by verify_case for the
    # manufactured case.  Keeping them explicit prevents this runner from
    # silently inheriting acceptance gates for other physical cases.
    return SimpleNamespace(
        geometry_absolute_tolerance=1e-11,
        geometry_relative_tolerance=1e-9,
        continuity_absolute_tolerance=1e-10,
        continuity_relative_tolerance=1e-7,
        max_reported_continuity=1e-8,
        max_iterations=args.max_iterations,
        manufactured_pressure_slope=args.manufactured_pressure_slope,
    )


def validate_mesh(path: Path, expected_cells: int) -> dict[str, Any]:
    sys.path.insert(0, str(REPO / "tools" / "verification"))
    from verify_native_flow import measure, read_cm2d  # pylint: disable=import-outside-toplevel

    mesh = read_cm2d(path)
    measured = measure(mesh, 1e-11, 1e-9)
    issues = list(measured.issues)
    if len(mesh.cells) != expected_cells:
        issues.append(f"expected {expected_cells} cells, got {len(mesh.cells)}")
    if not all(abs(value - target) <= 1e-10 for value, target in
               zip(measured.bounds, (0.0, 0.0, 1.0, 1.0))):
        issues.append(f"mesh bounds are not [0,1]^2: {measured.bounds}")
    return {"valid": not issues, "issues": issues, "counts": {"cells": len(mesh.cells),
            "faces": len(mesh.edges)}, "measurement": {"area": measured.total_area,
            "characteristicH": measured.characteristic_h, "bounds": measured.bounds},
            "sha256": sha256_file(path)}


def verify_flow(mesh: Path, prefix: Path, args: argparse.Namespace) -> dict[str, Any]:
    sys.path.insert(0, str(REPO / "tools" / "verification"))
    from verify_native_flow import verify_case  # pylint: disable=import-outside-toplevel

    return verify_case(mesh, prefix, "manufactured", args.nu, args.speed,
                       verifier_namespace(args))


def refinement_groups(cases: list[dict[str, Any]]) -> dict[str, Any]:
    sys.path.insert(0, str(REPO / "tools" / "verification"))
    from verify_native_flow import sequence_checks  # pylint: disable=import-outside-toplevel

    groups: dict[str, list[dict[str, Any]]] = {}
    for item in cases:
        key = f"{item['meshKind']}/{item['scheme']}"
        benchmark = item.get("verification", {}).get("benchmark", {})
        measurement = item.get("verification", {}).get("meshMeasurement", {})
        if "velocityL2Relative" in benchmark and "characteristicH" in measurement:
            groups.setdefault(key, []).append({
                "label": item["label"], "h": measurement["characteristicH"],
                "error": benchmark["velocityL2Relative"],
            })
    result: dict[str, Any] = {"valid": True, "issues": [], "series": {}, "groups": {}}
    for key, values in groups.items():
        selected = [{"case": "manufactured", "label": value["label"],
                     "benchmark": {"velocityL2Relative": value["error"]},
                     "meshMeasurement": {"characteristicH": value["h"]}}
                    for value in values]
        check = sequence_checks(selected)
        result["groups"][key] = check
        result["series"][key] = check["series"]["manufactured"]
        if not check["valid"]:
            result["valid"] = False
            result["issues"].extend(f"{key}: {issue}" for issue in check["issues"])
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--mesh-cli", type=Path, default=REPO / "build/cartmesh2d_cli")
    parser.add_argument("--flow-cli", type=Path, default=REPO / "build/cartmesh2d_flow_cli")
    parser.add_argument("--levels", nargs="+", type=int, choices=tuple(LEVELS), default=list(LEVELS))
    parser.add_argument("--kind", choices=KINDS, action="append", dest="kinds",
                        help="mesh kind; repeat to select both (default: both)")
    parser.add_argument("--scheme", choices=SCHEMES, action="append", dest="schemes",
                        help="convection scheme; repeat to select both (default: both)")
    parser.add_argument("--limit", type=int, help="run at most this many selected flow cases")
    parser.add_argument("--reuse", action="store_true",
                        help="verify existing artifacts; do not launch either native CLI")
    parser.add_argument("--mesh-only", action="store_true", help="generate/validate meshes without flow runs")
    parser.add_argument("--dry-run", action="store_true", help="write planned commands without launching them")
    parser.add_argument("--timeout", type=int, default=180, help="per native command timeout in seconds")
    parser.add_argument("--max-iterations", type=int, default=7000)
    parser.add_argument("--tolerance", type=float, default=1e-8)
    parser.add_argument("--nu", type=float, default=0.1)
    parser.add_argument("--speed", type=float, default=1.0)
    parser.add_argument("--manufactured-pressure-slope", type=float, default=0.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.kinds is None:
        args.kinds = list(KINDS)
    if args.schemes is None:
        args.schemes = list(SCHEMES)
    if args.timeout <= 0 or args.max_iterations <= 0 or args.limit is not None and args.limit <= 0:
        raise SystemExit("timeout, max-iterations and limit must be positive")
    if not math.isfinite(args.tolerance) or args.tolerance <= 0.0 or not math.isfinite(args.nu) or args.nu <= 0.0:
        raise SystemExit("tolerance and nu must be positive finite values")
    if not math.isfinite(args.speed) or args.speed <= 0.0:
        raise SystemExit("speed must be positive and finite")
    if not math.isfinite(args.manufactured_pressure_slope):
        raise SystemExit("manufactured pressure slope must be finite")
    if args.reuse and args.dry_run:
        raise SystemExit("--reuse and --dry-run are mutually exclusive")

    output_root = args.output_root.resolve()
    if args.reuse:
        if not output_root.is_dir():
            raise SystemExit(f"--reuse requires an existing output directory: {output_root}")
        geometry = output_root / "geometry" / "unit-square.xy"
    else:
        if output_root.exists() and any(output_root.iterdir()):
            raise SystemExit(f"fresh run requires an empty output directory: {output_root}")
        output_root.mkdir(parents=True, exist_ok=True)
        geometry = square_geometry(output_root)
    mesh_cli = args.mesh_cli.resolve()
    flow_cli = args.flow_cli.resolve()
    summary: dict[str, Any] = {
        "format": "cartmesh2d-manufactured-flow-runner-v1", "valid": False,
        "provenance": "reused-existing-artifacts" if args.reuse else "fresh-native-processes",
        "outputRoot": str(output_root), "geometry": str(geometry),
        "geometrySha256": sha256_file(geometry) if geometry.is_file() else None,
        "parameters": {"levels": args.levels, "kinds": args.kinds, "schemes": args.schemes,
                       "timeout": args.timeout, "maxIterations": args.max_iterations,
                       "tolerance": args.tolerance, "nu": args.nu, "speed": args.speed,
                       "manufacturedPressureSlope": args.manufactured_pressure_slope},
        "meshGeneration": [], "runs": [], "cases": [],
        "refinementGroups": None, "issues": [],
    }
    try:
        if not args.reuse:
            for path, name in ((mesh_cli, "mesh-cli"), (flow_cli, "flow-cli")):
                if not path.is_file():
                    raise FileNotFoundError(f"{name} not found: {path}")
        summary["verificationEnvironmentExecutables"] = {
            "meshCli": {"path": str(mesh_cli), "sha256": sha256_file(mesh_cli)} if mesh_cli.is_file() else None,
            "flowCli": {"path": str(flow_cli), "sha256": sha256_file(flow_cli)} if flow_cli.is_file() else None,
        }
        mesh_paths: dict[tuple[str, int], Path] = {}
        mesh_ok = True
        for level in args.levels:
            divisor, expected_cells = LEVELS[level]
            cart_dir = output_root / "meshes" / f"cartesian-l{level}"
            cart_prefix = cart_dir / f"cartesian-l{level}"
            cart_mesh = expected_mesh_path(cart_prefix)
            stage: dict[str, Any] = {"meshKind": "cartesian", "level": level,
                                     "expectedCells": expected_cells,
                                     "path": str(cart_mesh.resolve())}
            if not args.reuse and mesh_cli.is_file():
                stage["binarySha256"] = sha256_file(mesh_cli)
            if args.reuse:
                stage.update(status="reused-existing-artifacts", execution="reused-existing-artifacts")
            elif args.dry_run:
                command = mesh_command(mesh_cli, geometry, cart_prefix, level, 1.0 / divisor,
                                       output_root / "openfoam" / f"cartesian-l{level}")
                stage.update(status="planned", execution="dry-run", expectedCommand=command,
                             expectedCommandText=command_text(command))
            else:
                command = mesh_command(mesh_cli, geometry, cart_prefix, level, 1.0 / divisor,
                                       output_root / "openfoam" / f"cartesian-l{level}")
                stage.update(run_command(command, output_root / "logs" / f"mesh-cartesian-l{level}", args.timeout))
            command_succeeded = args.reuse or stage.get("status") == "passed"
            if not args.dry_run and command_succeeded and cart_mesh.is_file():
                stage["validation"] = validate_mesh(cart_mesh, expected_cells)
                mesh_ok = mesh_ok and stage["validation"]["valid"]
            elif not args.dry_run:
                issue = (f"mesh command failed with status {stage.get('status')}"
                         if not command_succeeded else f"missing mesh: {cart_mesh}")
                stage["validation"] = {"valid": False, "issues": [issue]}
                mesh_ok = False
            summary["meshGeneration"].append(stage)
            mesh_paths[("cartesian", level)] = cart_mesh

            warped_mesh = output_root / "meshes" / f"warped-l{level}" / f"warped-l{level}.solver.cm2d"
            warp_stage: dict[str, Any] = {"meshKind": "warped", "level": level,
                                          "expectedCells": expected_cells,
                                          "path": str(warped_mesh.resolve())}
            if args.reuse:
                warp_stage.update(status="reused-existing-artifacts", execution="reused-existing-artifacts")
            elif args.dry_run:
                warp_stage.update(status="planned", execution="dry-run", source=str(cart_mesh.resolve()),
                                  expectedTransform="x += 0.06*sin(pi*x)*sin(pi*y); y += 0.6*offset")
            elif cart_mesh.is_file() and stage.get("validation", {}).get("valid"):
                warp_stage.update(status="passed", execution="local-deterministic-transform",
                                   transform=warp_mesh(cart_mesh, warped_mesh))
            else:
                warp_stage.update(status="failed", execution="local-deterministic-transform",
                                  issues=["cartesian mesh unavailable or invalid"])
            if not args.dry_run and warped_mesh.is_file():
                warp_stage["validation"] = validate_mesh(warped_mesh, expected_cells)
                mesh_ok = mesh_ok and warp_stage["validation"]["valid"]
            elif not args.dry_run:
                warp_stage["validation"] = {"valid": False, "issues": [f"missing mesh: {warped_mesh}"]}
                mesh_ok = False
            summary["meshGeneration"].append(warp_stage)
            mesh_paths[("warped", level)] = warped_mesh

        if not args.mesh_only and not args.dry_run:
            selected = [(kind, scheme, level) for kind in args.kinds for scheme in args.schemes
                        for level in args.levels]
            if args.limit is not None:
                selected = selected[:args.limit]
            for kind, scheme, level in selected:
                mesh = mesh_paths[(kind, level)]
                label = f"{kind}-{scheme}-l{level}"
                prefix = output_root / "runs" / scheme / f"{kind}-l{level}" / "manufactured"
                command = [str(flow_cli), "--mesh", str(mesh), "--output", str(prefix),
                           "--case", "manufactured", "--nu", f"{args.nu:.17g}",
                           "--speed", f"{args.speed:.17g}", "--max-iterations", str(args.max_iterations),
                           "--tolerance", f"{args.tolerance:.17g}", "--convection", scheme, "--profile",
                           "--manufactured-pressure-slope", f"{args.manufactured_pressure_slope:.17g}"]
                stage: dict[str, Any] = {"label": label, "meshKind": kind, "scheme": scheme,
                                         "level": level, "mesh": str(mesh.resolve()),
                                         "meshSha256": sha256_file(mesh) if mesh.is_file() else None}
                if not args.reuse and flow_cli.is_file():
                    stage["binarySha256"] = sha256_file(flow_cli)
                if args.reuse:
                    stage.update(status="reused-existing-artifacts", execution="reused-existing-artifacts",
                                 returncode=None, expectedCommand=command,
                                 expectedCommandText=command_text(command))
                else:
                    stage["command"] = command
                    stage["commandText"] = command_text(command)
                    stage.update(run_command(command, output_root / "logs" / f"flow-{label}", args.timeout))
                item: dict[str, Any] = {"label": label, "meshKind": kind, "scheme": scheme,
                                        "level": level, "stage": stage}
                if stage.get("status") == "passed" or args.reuse:
                    try:
                        item["verification"] = verify_flow(mesh, prefix, args)
                    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
                        item["verification"] = {"valid": False, "issues": [str(exc)]}
                else:
                    item["verification"] = {"valid": False, "issues": ["native flow command failed"]}
                summary["runs"].append(stage)
                summary["cases"].append(item)
        summary["refinementGroups"] = refinement_groups(summary["cases"])
        summary["valid"] = None if args.dry_run else (mesh_ok and (args.mesh_only or
                                          bool(summary["cases"]) and
                                          all(item.get("verification", {}).get("valid") is True
                                              for item in summary["cases"]) and
                                          summary["refinementGroups"]["valid"]))
        summary["status"] = ("planned" if args.dry_run else
                              "failed" if not summary["valid"] else
                              "mesh-only-validated" if args.mesh_only else
                              "verified-reused" if args.reuse and summary["valid"] else
                              "verified" if summary["valid"] else "failed")
        if not summary["valid"] and not args.dry_run:
            summary["issues"].append("one or more mesh or manufactured-flow checks failed")
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        summary["issues"].append(str(exc))
    # A read-back must not destroy the original commands, return codes and
    # producer hashes retained by a fresh run.
    summary_path = output_root / ("runner-reverification.json" if args.reuse else "runner-summary.json")
    write_json(summary_path, summary)
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0 if summary["valid"] or (args.dry_run and not summary["issues"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
