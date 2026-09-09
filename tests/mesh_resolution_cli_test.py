#!/usr/bin/env python3
"""Small scale-invariance regression for the dimensionless CLI controls."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import shutil
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.verification.check_mesh_resolution import measure
from tools.verification.check_layer_resolution import measure as measure_layers
from tools.verification.check_directional_connectivity import (
    measure as measure_directional, verify_native_report)


def run(command: list[str], log, allow_failure: bool = False,
        timeout_seconds: float = 20.0) -> subprocess.CompletedProcess[str]:
    try:
        completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   timeout=timeout_seconds)
    except subprocess.TimeoutExpired as exc:
        output = (exc.output or "")
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        log.write(f"$ {' '.join(command)}\nreturncode=124 timeout={timeout_seconds}s\n{output}\n")
        if not allow_failure:
            raise AssertionError(f"command timed out after {timeout_seconds}s: {' '.join(command)}")
        return subprocess.CompletedProcess(command, 124, output)
    log.write(f"$ {' '.join(command)}\nreturncode={completed.returncode}\n{completed.stdout}\n")
    if completed.returncode and not allow_failure:
        raise AssertionError(f"command failed ({completed.returncode}): {' '.join(command)}\n"
                             f"{completed.stdout}")
    return completed


def parse_kv(text: str) -> dict[str, str]:
    values = {}
    for line in text.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def rectangle(path: Path, scale: float) -> None:
    path.write_text("\n".join(f"{scale * x:.17g} {scale * y:.17g}"
                              for x, y in ((0, 0), (1, 0), (1, 0.5), (0, 0.5))) + "\n",
                    encoding="utf-8")


def assert_close(a: float, b: float, message: str) -> None:
    if not math.isclose(a, b, rel_tol=3e-10, abs_tol=3e-12):
        raise AssertionError(f"{message}: {a} != {b}")


def sizing_args(reference: float, hybrid: bool = False) -> list[str]:
    values = ["--size-field", "--reference-length", f"{reference:.17g}",
              "--wall-relative-size", "0.25", "--background-relative-size", "0.5",
              "--far-field-spans", "0.5"]
    if hybrid:
        values.extend(["--first-layer-relative-size", "0.05"])
    return values


def pure_size_only(cli: Path, boundary: Path, prefix: Path, reference: float, log) -> dict[str, str]:
    command = [str(cli), str(boundary), str(prefix), "5", "0.25", "0.1", "exterior", "-", "0",
               *sizing_args(reference), "--size-field-only"]
    completed = run(command, log)
    return parse_kv(completed.stdout)


def hybrid_size_only(cli: Path, boundary: Path, prefix: Path, reference: float, log) -> dict[str, str]:
    command = [str(cli), str(boundary), str(prefix), "5", "2", "5", "2", "0.05", "1.1", "0.5",
               "-", f"{reference * 0.1:.17g}", *sizing_args(reference, hybrid=True), "--size-field-only"]
    completed = run(command, log)
    report = json.loads(Path(str(prefix) + ".size-field.json").read_text(encoding="utf-8"))
    return {"size_field_domain_span": str(report["domain_span"]),
            "size_field_wall_cell_size": str(report["wall_cell_size"]),
            "size_field_wall_level": str(report["wall_level"]),
            "size_field_max_level": str(report["max_level"]),
            "reference_source": report["reference_source"]}


def actual_pure(cli: Path, boundary: Path, prefix: Path, case: Path, reference: float, log) -> dict:
    command = [str(cli), str(boundary), str(prefix), "5", "0.25", "0.1", "exterior", str(case), "0",
               *sizing_args(reference)]
    completed = run(command, log)
    return {"mode": "pure", "status": "success", "stdout": completed.stdout,
            "prefix": str(prefix), "resolution": str(prefix) + ".resolution.json"}


def actual_hybrid(cli: Path, boundary: Path, prefix: Path, case: Path, reference: float, log) -> dict:
    command = [str(cli), str(boundary), str(prefix), "5", "2", "5", "2", "0.05", "1.1", "0.5",
               str(case), f"{reference * 0.1:.17g}", *sizing_args(reference, hybrid=True)]
    completed = run(command, log, allow_failure=True)
    output = {"mode": "hybrid", "status": "success" if completed.returncode == 0 else "failed",
              "returncode": completed.returncode, "stdout": completed.stdout,
              "prefix": str(prefix), "resolution": str(prefix) + ".resolution.json"}
    if completed.returncode == 0:
        output["hybrid_status"] = "hybrid_status=success" in completed.stdout
        if not output["hybrid_status"]:
            output["status"] = "fallback_or_unlabelled"
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--hybrid-cli", type=Path, required=True)
    parser.add_argument("--log", type=Path,
                        default=Path("outputs/engineering-resolution/cli-regression.log"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("outputs/engineering-resolution/scale-regression"))
    args = parser.parse_args()
    args.log.parent.mkdir(parents=True, exist_ok=True)
    results = {"valid": False, "scales": [], "issues": [], "log": str(args.log)}
    try:
        root = args.output_dir
        if root.exists():
            shutil.rmtree(root)
        root.mkdir(parents=True, exist_ok=True)
        args.log.write_text("", encoding="utf-8")
        invalid_boundary = root / "invalid-control.xy"
        rectangle(invalid_boundary, 1.0)
        with args.log.open("a", encoding="utf-8") as log:
            invalid_pure = run(
                [str(args.cli), str(invalid_boundary), str(root / "invalid-pure"), "5", "0.25",
                 "0.1", "exterior", "-", "0", "--size-field", "--reference-length", "1",
                 "--wall-relative-size", "0", "--size-field-only"], log, allow_failure=True)
            invalid_hybrid = run(
                [str(args.hybrid_cli), str(invalid_boundary), str(root / "invalid-hybrid"), "5", "2",
                 "5", "2", "0.05", "1.1", "0.5", "-", "0.1", "--size-field",
                 "--reference-length", "1", "--wall-relative-size", "0", "--size-field-only"],
                log, allow_failure=True)
        if invalid_pure.returncode == 0 or invalid_hybrid.returncode == 0:
            raise AssertionError("non-positive relative wall size was not rejected")
        results["invalid_parameter_rejection"] = {
            "pure_returncode": invalid_pure.returncode,
            "hybrid_returncode": invalid_hybrid.returncode,
        }
        for scale in (0.001, 1.0, 1000.0):
            boundary = root / f"rectangle-{scale:g}.xy"
            rectangle(boundary, scale)
            with args.log.open("a", encoding="utf-8") as log:
                pure = pure_size_only(args.cli, boundary, root / f"pure-{scale:g}", scale, log)
                hybrid = hybrid_size_only(args.hybrid_cli, boundary, root / f"hybrid-{scale:g}", scale, log)
            expected_domain = scale * 2.0  # body span = scale, padding = 0.5 Lref
            for mode, values in (("pure", pure), ("hybrid", hybrid)):
                assert_close(float(values["size_field_domain_span"]), expected_domain,
                             f"{mode} domain span at scale {scale}")
                assert_close(float(values["size_field_wall_cell_size"]), scale * 0.25,
                             f"{mode} wall size at scale {scale}")
                if values.get("size_field_wall_level") != values.get("size_field_max_level"):
                    raise AssertionError(f"{mode} resolved levels diverge at scale {scale}")
            if pure["size_field_domain_span"] != hybrid["size_field_domain_span"]:
                raise AssertionError(f"pure/hybrid domain mismatch at scale {scale}")
            with args.log.open("a", encoding="utf-8") as log:
                actual = actual_pure(args.cli, boundary, root / f"pure-mesh-{scale:g}",
                                     root / f"pure-case-{scale:g}", scale, log)
                hybrid_actual = actual_hybrid(args.hybrid_cli, boundary,
                                              root / f"hybrid-mesh-{scale:g}",
                                              root / f"hybrid-case-{scale:g}", scale, log)
            actual["size_only"] = pure
            hybrid_actual["size_only"] = hybrid
            for product in (actual, hybrid_actual):
                report_path = Path(product["resolution"])
                if product["status"] == "success" and report_path.exists():
                    report = json.loads(report_path.read_text(encoding="utf-8"))
                    product["native_report"] = report
                    if report.get("reference_source") != "explicit":
                        raise AssertionError(f"{product['mode']} lost explicit reference length")
                    cm2d_path = (Path(product["prefix"] + ".solver.cm2d")
                                 if product["mode"] == "pure"
                                 else Path(product["prefix"] + ".hybrid.solver.cm2d"))
                    independent = measure(cm2d_path, report_path)
                    product["independent_resolution"] = independent
                    directional = measure_directional(cm2d_path)
                    verify_native_report(directional, report)
                    product["independent_directional"] = directional
                    if not independent["valid"]:
                        raise AssertionError(f"independent {product['mode']} resolution check failed: "
                                             f"{independent['issues']}")
                    if product["mode"] == "hybrid":
                        product["independent_layers"] = measure_layers(
                            Path(product["prefix"]+".hybrid.vtk"),cm2d_path,report_path)
                    if product["mode"] == "pure" and scale == 1.0:
                        forged_path = root / "forged-resolution.json"
                        forged = json.loads(json.dumps(report))
                        forged["actual"]["sqrt_area_over_reference"]["p50"] += 0.125
                        forged_path.write_text(json.dumps(forged), encoding="utf-8")
                        forged_check = measure(cm2d_path, forged_path)
                        results["forged_report_rejected"] = not forged_check["valid"]
                        if forged_check["valid"]:
                            raise AssertionError("independent reader accepted forged native statistics")
                        legacy_path = root / "legacy-no-wall-resolution.json"
                        legacy = json.loads(json.dumps(report))
                        legacy["requested"]["wall_h_over_reference"] = None
                        legacy["wall_owner_tangential_exceedance_length_fraction"] = None
                        legacy_path.write_text(json.dumps(legacy), encoding="utf-8")
                        legacy_check = measure(cm2d_path, legacy_path)
                        results["legacy_no_wall_request"] = legacy_check["valid"]
                        if not legacy_check["valid"]:
                            raise AssertionError("legacy no-wall request was not accepted")
                elif product["mode"] == "hybrid":
                    product["fallback_recorded"] = True
            results["scales"].append({"scale": scale, "pure": actual,
                                      "hybrid": hybrid_actual})
        for row in results["scales"]:
            for mode in ("pure", "hybrid"):
                product = row[mode]
                if product.get("status") != "success":
                    results["issues"].append(
                        f"{mode} scale {row['scale']} did not complete successfully")
                if not product.get("independent_resolution", {}).get("valid", False):
                    results["issues"].append(
                        f"{mode} scale {row['scale']} lacks an independent resolution pass")
        # A sparse input polyline must not fix the layer spacing when an
        # explicit wall size was requested. The original 32 corners are kept.
        circle = Path(__file__).resolve().parents[1] / "examples/acceptance/circle.xy"
        prefix = root / "tangential-circle"
        with args.log.open("a", encoding="utf-8") as log:
            completed = run([str(args.hybrid_cli),str(circle),str(prefix),"6","3","6","4",
                "0.02","1.2","1","-","0.02","--size-field","--wall-relative-size","0.03125",
                "--background-relative-size","0.25","--far-field-spans","0.5",
                "--cells-per-level","0","--first-layer-relative-size","0.01"],log)
        if "hybrid_status=success" not in completed.stdout:
            raise AssertionError("tangential circle fell back instead of producing layers")
        report_path=Path(str(prefix)+".resolution.json")
        report=json.loads(report_path.read_text())
        if report["wall_owner_tangential_exceedance_length_fraction"]>1.e-9:
            raise AssertionError("hybrid wall size failed to refine the sparse circle polyline")
        results["tangential_circle"]=measure_layers(Path(str(prefix)+".hybrid.vtk"),
            Path(str(prefix)+".hybrid.solver.cm2d"),report_path)
        results["valid"] = not results["issues"]
    except (AssertionError, KeyError, ValueError, OSError, json.JSONDecodeError) as exc:
        results["issues"].append(str(exc))
    summary = {
        "valid": results["valid"],
        "issues": results["issues"],
        "scales": [
            {"scale": row["scale"],
             "pure_status": row["pure"].get("status"),
             "pure_independent": row["pure"].get("independent_resolution", {}).get("valid", False),
             "hybrid_status": row["hybrid"].get("status"),
             "hybrid_independent": row["hybrid"].get("independent_resolution", {}).get("valid", False)}
            for row in results["scales"]
        ],
        "invalid_parameter_rejection": results.get("invalid_parameter_rejection"),
        "forged_report_rejected": results.get("forged_report_rejected", False),
        "legacy_no_wall_request": results.get("legacy_no_wall_request", False),
        "log": results["log"],
        "output_dir": str(args.output_dir),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    encoded = json.dumps(results, indent=2, sort_keys=True) + "\n"
    print(encoded, end="")
    return 0 if results["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
