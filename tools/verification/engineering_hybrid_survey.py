#!/usr/bin/env python3
"""Bounded, measurement-only survey of fixed hybrid relative-size requests."""
from __future__ import annotations
import argparse, hashlib, json, math, os, re, signal, subprocess
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools/verification"))
from engineering_scale_survey import parse_native_seconds, parse_rss
from check_mesh_resolution import measure
from check_openfoam2d import check as check_openfoam
from check_layer_resolution import measure as measure_layers

# Freeze requests before running the survey. Exterior requests match the
# measured pure-path survey; the smaller interior region starts one step finer.
# Band growth resolves a larger near-wall region without padding the far field.
REQUESTS = {
    name: {
        "r01": dict(wall=.01, background=.025, first=.0025, layers=4, band=3),
        "r02": dict(wall=.005, background=.0125, first=.00125, layers=4, band=3),
        "r03": dict(wall=.0025, background=.0125, first=.000625, layers=4, band=40),
    } for name in ("circle", "naca2412_dense")
}
REQUESTS["nozzle"] = {
    "r01": dict(wall=.005, background=.0125, first=.00125, layers=4, band=3),
    "r02": dict(wall=.0025, background=.0125, first=.000625, layers=4, band=40),
    "r03": dict(wall=.00125, background=.0125, first=.0003125, layers=4, band=40),
}
GRID_NAMES = ("r01", "r02", "r03")
CASES = {
    "circle": ("examples/acceptance/circle.xy", "exterior", 2.0),
    "nozzle": ("examples/complex/nozzle_profile.xy", "interior", 6.0),
    "naca2412_dense": ("examples/complex/naca2412_dense.xy", "exterior", 1.0),
}

def run_one(output, cli, case_name, request_name, timeout):
    boundary, region, lref = CASES[case_name]; q = REQUESTS[case_name][request_name]
    output.mkdir(parents=True, exist_ok=True); prefix, foam = output / "mesh", output / "openfoam"
    command = ["/usr/bin/time", "-l", str(cli), str(REPO / boundary), str(prefix), "11", "0", "11",
               str(q["layers"]), "0.01", "1.2", "0.5", str(foam), "0.02", "--small-alpha=0.1",
               "--fluid-region=" + region, "--reference-length", str(lref),
               "--wall-relative-size", str(q["wall"]), "--background-relative-size", str(q["background"]),
               "--far-field-spans", "0.5", "--first-layer-relative-size", str(q["first"]),
               "--cells-per-level", str(q["band"]), "--max-safe-wall-level", "11"]
    (output / "command.json").write_text(json.dumps(command, indent=2, allow_nan=False) + "\n")
    proc = subprocess.Popen(command, cwd=REPO, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, start_new_session=True)
    timed = False
    try:
        text, _ = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed = True
        try: os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError: pass
        text, _ = proc.communicate()
    text = text or ""; (output / "command.log").write_text(text, encoding="utf-8")
    elapsed = re.search(r"([0-9.]+)\s+real\s+", text)
    hybrid_marker = bool(re.search(r"^hybrid_status=success\b", text, re.M))
    fallback_marker = bool(re.search(r"^h4_status=success\s+mesh_mode=pure_cutcell_fallback\b", text, re.M))
    result = {"case": case_name, "request": request_name, "reference_length": lref,
              "fluid_region": region, "parameters": q, "command": command,
              "wrapper_returncode": 124 if timed else proc.returncode, "timeout_seconds": timeout,
              "timed_out": timed, "native_completed_marker": hybrid_marker or fallback_marker,
              "native_mode": "hybrid" if hybrid_marker else ("pure_cutcell_fallback" if fallback_marker else "unknown"),
              "native_elapsed_seconds": parse_native_seconds(text),
              "wrapper_elapsed_seconds": float(elapsed.group(1)) if elapsed else None,
              "peak_rss_bytes": parse_rss(text), "solver": "not_reported",
              "external_checkMesh": "not_run", "issues": []}
    for key, pattern in (("solver", r"solver_quality=(\w+)"),):
        match = re.search(pattern, text, re.M)
        if match: result[key] = match.group(1)
    if timed: result["status"] = "timeout"
    elif proc.returncode != 0: result["status"] = "failed"
    else: result["status"] = "success"
    suffix = ".fallback.solver" if fallback_marker else ".hybrid.solver"
    cm2d, report, vtk = Path(str(prefix) + suffix + ".cm2d"), Path(str(prefix) + ".resolution.json"), Path(str(prefix) + (".fallback.vtk" if fallback_marker else ".hybrid.vtk"))
    if result["status"] == "success":
        try:
            result["resolution_independent"] = measure(cm2d, report)
            result["openfoam_independent"] = check_openfoam(foam) if foam.exists() else {"valid": False, "issues": ["missing OpenFOAM case"]}
            result["cell_count"] = result["resolution_independent"]["measured"]["cell_count"]
            native_report = json.loads(report.read_text(encoding="utf-8"))
            result["boundary_layer_coverage_status"] = native_report.get("boundary_layer_coverage_status", "missing")
            bl = native_report.get("boundary_layers") if isinstance(native_report.get("boundary_layers"), dict) else None
            result["boundary_layer_summary"] = ({k: bl.get(k) for k in ("status", "requested_cells", "constructed_cells", "retained_cells", "first_layer_wall_length_fraction", "full_requested_layers_wall_length_fraction", "first_layer_height_exceedance_wall_length_fraction")} if bl else None)
            result["pure_fallback"] = fallback_marker
            if fallback_marker and result["boundary_layer_coverage_status"] != "no_layers_retained":
                result["issues"].append("fallback marker disagrees with final layer report")
            if hybrid_marker and vtk.exists() and result["resolution_independent"].get("valid"):
                layer_out = output / "layer-independent.json"
                result["layer_independent"] = measure_layers(vtk, Path(str(prefix) + ".hybrid.solver.cm2d"), report)
                layer_out.write_text(json.dumps(result["layer_independent"], indent=2, allow_nan=False) + "\n")
            else: result["layer_independent"] = {"valid": fallback_marker, "status": "not_applicable_fallback" if fallback_marker else "not_available"}
        except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
            result["issues"].append(str(exc))
    finite_metrics = all(isinstance(result.get(k), (int, float)) and math.isfinite(result[k]) for k in ("native_elapsed_seconds", "wrapper_elapsed_seconds", "peak_rss_bytes"))
    result["measurement_complete"] = (result["status"] == "success" and result["wrapper_returncode"] == 0
        and result["native_completed_marker"] and finite_metrics and not result["issues"]
        and result.get("resolution_independent", {}).get("valid", False)
        and result.get("openfoam_independent", {}).get("valid", False)
        and result.get("layer_independent", {}).get("valid", False))
    (output / "measurement.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--cli", type=Path, default=REPO / "build/cartmesh2d_hybrid_cli")
    ap.add_argument("--output-dir", type=Path, default=REPO / "outputs/engineering-hybrid-survey"); ap.add_argument("--timeout", type=float, default=180.0)
    ap.add_argument("--cases", nargs="+", choices=tuple(CASES), default=list(CASES)); ap.add_argument("--grids", nargs="+", choices=GRID_NAMES, default=list(GRID_NAMES))
    args=ap.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0: ap.error("timeout must be finite and positive")
    args.cli=args.cli.resolve(strict=True); args.output_dir=args.output_dir.resolve()
    if args.output_dir.exists() and any(args.output_dir.iterdir()): raise SystemExit(f"refusing to overwrite non-empty output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True); results=[]
    for case_name in args.cases:
        for request_name in args.grids:
            try: results.append(run_one(args.output_dir / case_name / ("request-" + request_name), args.cli, case_name, request_name, args.timeout))
            except Exception as exc:
                results.append({"case":case_name,"request":request_name,"status":"failed","measurement_complete":False,"issues":[str(exc)]})
    summary={"format":"cartmesh2d-engineering-hybrid-survey-v1","scope":{"cases":list(CASES),"requests":REQUESTS,"selected_cases":args.cases,"selected_requests":args.grids},"cli":str(args.cli),"cli_sha256":hashlib.sha256(args.cli.read_bytes()).hexdigest(),"size_field_auto_enabled":True,"external_checkMesh":"not_run","cases":results,"process_success_count":sum(r["status"]=="success" for r in results),"hybrid_certified_count":sum(r.get("native_mode")=="hybrid" and r["measurement_complete"] for r in results),"fallback_count":sum(r.get("native_mode")=="pure_cutcell_fallback" for r in results),"measurement_complete_count":sum(r["measurement_complete"] for r in results),"notes":["Requested quantities are measurements, not claims of 1e4/5e4/1e5 cells.","Solver, external checkMesh and layer authentication are separate reports.","Relative options implicitly create the size-field policy in cartmesh2d_hybrid_cli; --size-field is not required."]}
    (args.output_dir/"summary.json").write_text(json.dumps(summary,indent=2, allow_nan=False)+"\n"); print(json.dumps({k:summary[k] for k in ('process_success_count','hybrid_certified_count','fallback_count','measurement_complete_count','external_checkMesh')},indent=2, allow_nan=False)); return 0 if all(r["measurement_complete"] for r in results) else 1
if __name__ == "__main__": raise SystemExit(main())
