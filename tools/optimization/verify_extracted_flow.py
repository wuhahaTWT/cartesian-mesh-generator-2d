#!/usr/bin/env python3
"""Send an extracted topology boundary through the actual native Cut-cell chain.

Preserves every mesh/solver failure. Uses explicit interior semantics, original
quality gates, an independent CM2D reader, and an optional low-Re native flow.
This checks the SHARP extracted shape, not the finite-resistance analysis mesh.
Native outlets use p=0 rather than the prescribed Stokes outflow profile, so
this is a separate physical configuration and does not validate equal objectives.
"""
import argparse
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"tools"/"verification"))
import verify_native_flow as native
from brinkman import port_average


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def execute(command, prefix, timeout):
    command = list(map(str, command))
    start = time.monotonic()
    environment = dict(os.environ, VECLIB_MAXIMUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    try:
        result = subprocess.run(command, capture_output=True, timeout=timeout, env=environment)
        code, stdout, stderr, timed_out = result.returncode, result.stdout, result.stderr, False
    except subprocess.TimeoutExpired as exc:
        code, stdout, stderr, timed_out = None, exc.stdout or b"", exc.stderr or b"", True
    Path(str(prefix)+".stdout.log").write_bytes(stdout)
    Path(str(prefix)+".stderr.log").write_bytes(stderr)
    return dict(command=command, returncode=code, timedOut=timed_out,
                elapsedSeconds=time.monotonic()-start,
                stdout=str(prefix)+".stdout.log", stderr=str(prefix)+".stderr.log")


def parabolic_boundaries(source, target, problem, speed):
    centres = [problem["height"]/4]
    if problem["case"] == "double-pipe":
        centres.append(3*problem["height"]/4)
    right_centres = centres if problem["case"] == "double-pipe" else [3*problem["height"]/4]
    coordinate_tolerance = 1e-11+1e-9*max(problem["width"], problem["height"])
    counts = dict(inletFaces=0, outletFaces=0, closedArtificialOpenings=0)
    rows = []
    for line in source.read_text().splitlines():
        if not line.startswith("BOUNDARY "):
            rows.append(line)
            continue
        fields = shlex.split(line)
        if len(fields) != 12:
            raise ValueError("invalid native boundary export")
        if fields[7] in ("velocity-inlet", "pressure-outlet"):
            y, length = float(fields[4]), abs(float(fields[5]))
            if length <= 0 or abs(float(fields[6])) > 1e-12*length:
                raise ValueError("topology prototype only supports vertical port faces")
            x = float(fields[3])
            inlet = fields[7] == "velocity-inlet"
            port_centres = centres if inlet else right_centres
            actual_side = abs(x-(0 if inlet else problem["width"])) <= coordinate_tolerance
            profile = sum(port_average(y-length/2, y+length/2, c, problem["port_width"]) for c in port_centres)
            if actual_side and profile > 0 and inlet:
                velocity = speed*profile
                fields[9] = format(velocity, ".17g")
                fields[8] = "inlet_"+str(min(range(len(centres)), key=lambda i:abs(y-centres[i])))
                counts["inletFaces"] += 1
            elif actual_side and profile > 0:
                fields[8] = "outlet_"+str(min(range(len(right_centres)), key=lambda i:abs(y-right_centres[i])))
                counts["outletFaces"] += 1
            else:
                fields[7:12] = ["wall", "closed_end", "0", "0", "0"]
                counts["closedArtificialOpenings"] += 1
        fields[8] = json.dumps(fields[8])
        rows.append(" ".join(fields))
    target.write_text("\n".join(rows)+"\n")
    return counts


def pressure_metrics(flow, boundary, problem, speed, viscosity):
    with Path(str(flow)+".faces.csv").open() as stream:
        faces = {int(row["face"]):row for row in csv.DictReader(stream)}
    inlet_power, outlet_power, inlet_flux, outlet_flux = 0.0, 0.0, 0.0, 0.0
    for line in boundary.read_text().splitlines():
        if not line.startswith("BOUNDARY "):
            continue
        fields = shlex.split(line)
        face = faces[int(fields[1])]
        flux, pressure = float(face["flux"]), float(face["pressure"])
        if fields[7] == "velocity-inlet":
            inlet_flux -= flux; inlet_power -= pressure*flux
        elif fields[7] == "pressure-outlet":
            outlet_flux += flux; outlet_power += pressure*flux
    if inlet_flux <= 0:
        raise ValueError("pressure-drop comparison requires a positive inlet throughput")
    return dict(inletFlux=inlet_flux, outletFlux=outlet_flux,
                pressurePower=inlet_power-outlet_power,
                fluxWeightedPressureDrop=(inlet_power-outlet_power)/inlet_flux,
                speed=speed, viscosity=viscosity,
                nominalReynolds=speed*problem["port_width"]/viscosity)


def run(args):
    directory = args.directory.resolve(strict=True)
    root = args.output.resolve()
    if root.exists():
        raise ValueError("native output exists; preserve it and choose a fresh directory")
    extraction = json.loads((directory/"extraction.json").read_text())
    research = json.loads((directory/"summary.json").read_text())
    if args.component is None and extraction["components"] > 1:
        # Every region is retained and solved. Each disconnected fluid region
        # needs its own pressure reference; the native solver rejects a combined
        # disconnected solve. This is decomposition, never deletion of cells.
        root.mkdir(parents=True)
        components = []
        for index in range(extraction["components"]):
            child = copy.copy(args)
            child.component, child.output = index, root/f"component-{index}"
            components.append(run(child))
        report = dict(schema="cartmesh2d-topology-native-v1", components=components,
                      scope="All disconnected fluid regions solved separately with independent pressure references.",
                      physicalAccuracyQualified=False,
                      originalBoundarySha256=sha(directory/"fluid.xy"),
                      expectedArea=extraction["fluidArea"],
                      externalCheckMesh=[r["externalCheckMesh"] for r in components])
        for key in ("meshAccepted", "independentTopologyPassed", "solverQualityPassed",
                    "nativeFlowConverged", "independentFlowPassed"):
            report[key] = all(r[key] for r in components)
        report["totalAcceptedArea"] = sum(r.get("acceptedMesh", {}).get("area", 0) for r in components)
        report["areaError"] = report["totalAcceptedArea"]-report["expectedArea"]
        if abs(report["areaError"]) > 1e-11+1e-9*abs(report["expectedArea"]):
            report["independentTopologyPassed"] = False
        (root/"summary.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        return report
    if args.component is not None:
        if not 0 <= args.component < len(extraction.get("regions", [])):
            raise ValueError("invalid or missing extracted component")
        selected = extraction["regions"][args.component]
        source = directory/selected["file"]
        extraction = dict(extraction, fluidArea=selected["fluidArea"], boundarySha256=selected["boundarySha256"])
    else:
        source = directory/"fluid.xy"
    if sha(source) != extraction["boundarySha256"]:
        raise ValueError("extracted boundary changed after provenance was recorded")
    mesh_cli, flow_cli = args.mesh_cli.resolve(strict=True), args.flow_cli.resolve(strict=True)
    root.mkdir(parents=True)
    report = dict(schema="cartmesh2d-topology-native-v1", source=str(source), sourceSha256=sha(source),
                  executables={str(p):sha(p) for p in (mesh_cli, flow_cli)},
                  smallCellAggregationFraction=getattr(args, "small_alpha", .1),
                  meshAccepted=False, independentTopologyPassed=False, solverQualityPassed=False,
                  nativeFlowConverged=False, independentFlowPassed=False, physicalAccuracyQualified=False,
                  externalCheckMesh="not-run", cases=[], issues=[],
                  scope="Extracted sharp-wall interior mesh; separate low-Re Navier-Stokes flow with pressure outlets.")
    def save():
        (root/"summary.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    save()
    for level in args.levels:
        case = root/f"level-{level}"
        case.mkdir()
        prefix = case/"mesh"
        item = dict(level=level, runs=[], status="meshing")
        report["cases"].append(item)
        save()
        command = [mesh_cli, source, prefix, level, 1/30, getattr(args, "small_alpha", .1),
                   "interior", case/"openfoam", level, 0]
        meshing = execute(command, case/"mesh-run", args.timeout)
        item["runs"].append(meshing)
        mesh_path = Path(str(prefix)+".solver.cm2d")
        if meshing["returncode"] != 0 or not mesh_path.exists():
            item["status"] = "mesh-rejected"
            save()
            continue
        mesh = native.read_cm2d(mesh_path)
        measured = native.measure(mesh, 1e-11, 1e-9)
        area_error = measured.total_area-extraction["fluidArea"]
        # Same existing reader geometry tolerances; no new quality threshold.
        area_tolerance = 1e-11+1e-9*abs(extraction["fluidArea"])
        item["mesh"] = dict(path=str(mesh_path), sha256=sha(mesh_path), cells=len(mesh.cells),
                            faces=len(mesh.edges), area=measured.total_area,
                            expectedArea=extraction["fluidArea"], areaError=area_error,
                            independentIssues=measured.issues)
        if measured.issues or abs(area_error) > area_tolerance:
            item["status"] = "independent-topology-rejected"
            save()
            continue
        report["meshAccepted"] = report["independentTopologyPassed"] = report["solverQualityPassed"] = True
        report["acceptedMesh"] = item["mesh"]
        item["status"] = "mesh-accepted"
        check_mesh = shutil.which("checkMesh")
        if check_mesh:
            check = execute([check_mesh, "-case", case/"openfoam"], case/"checkMesh", args.timeout)
            item["runs"].append(check)
            text = Path(check["stdout"]).read_text(errors="replace")
            report["externalCheckMesh"] = "passed-standard" if check["returncode"] == 0 and "Mesh OK" in text else "failed"
        if args.mesh_only:
            save()
            break
        template, boundary = case/"template.boundaries", case/"flow.boundaries"
        boundary_run = execute([flow_cli, "--mesh", mesh_path, "--case", "duct", "--speed", args.speed,
                                "--export-boundaries", template], case/"boundary-run", args.timeout)
        item["runs"].append(boundary_run)
        if boundary_run["returncode"] != 0:
            item["status"] = "boundary-rejected"; save(); continue
        ports = parabolic_boundaries(template, boundary, research["problem"], args.speed)
        item["ports"] = ports
        if not ports["inletFaces"] or not ports["outletFaces"]:
            item["status"] = "ports-rejected"
            item["issue"] = "A component without both global inlet and outlet ports is retained but not qualified by this through-flow bridge."
            save()
            continue
        flow = case/"flow"
        command = [flow_cli, "--mesh", mesh_path, "--case", "custom", "--boundary", boundary,
                   "--output", flow, "--nu", args.nu, "--speed", args.speed,
                   "--max-iterations", args.iterations, "--tolerance", args.tolerance,
                   "--velocity-relaxation", .2, "--pressure-corrections", 1,
                   "--steady-acceleration", "anderson", "--linear-policy", "adaptive",
                   "--pressure-preconditioner", "aggregation"]
        item["status"] = "solving"
        save()
        solve = execute(command, case/"flow-run", args.timeout)
        item["runs"].append(solve)
        summary_path = Path(str(flow)+".json")
        summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
        item["flow"] = summary
        report["nativeFlowConverged"] = solve["returncode"] == 0 and summary.get("converged") is True
        if not report["nativeFlowConverged"]:
            item["status"] = "flow-not-converged"; save(); continue
        options = native.argument_parser().parse_args(["--max-iterations", str(args.iterations)])
        audit = native.verify_case(mesh_path, flow, "custom", args.nu, args.speed, options)
        (case/"independent-flow.json").write_text(json.dumps(audit, indent=2, allow_nan=False)+"\n")
        report["independentFlowPassed"] = bool(audit["valid"])
        item["status"] = "flow-audited" if audit["valid"] else "flow-audit-failed"
        item["independentFlowIssues"] = audit["issues"]
        item["metrics"] = pressure_metrics(flow, boundary, research["problem"], args.speed, args.nu)
        save()
        break
    if not report["meshAccepted"]:
        report["issues"].append("All requested native mesh attempts were explicitly rejected.")
    elif not args.mesh_only and not report["independentFlowPassed"]:
        report["issues"].append("No requested mesh completed the native flow and independent equation audit.")
    save()
    print(json.dumps({key:report[key] for key in ("meshAccepted", "independentTopologyPassed", "solverQualityPassed",
                                                  "nativeFlowConverged", "independentFlowPassed", "externalCheckMesh")}, indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mesh-cli", type=Path, default=ROOT/"build/cartmesh2d_cli")
    parser.add_argument("--flow-cli", type=Path, default=ROOT/"build/cartmesh2d_flow_cli")
    parser.add_argument("--levels", nargs="+", type=int, default=[5, 6])
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--iterations", type=int, default=1500)
    parser.add_argument("--speed", type=float, default=.02)
    parser.add_argument("--nu", type=float, default=1.)
    parser.add_argument("--tolerance", type=float, default=1e-6,
                        help="existing native normalised stopping tolerance, separate from physical accuracy")
    parser.add_argument("--mesh-only", action="store_true")
    parser.add_argument("--component", type=int, help="one extracted region; default processes ALL regions")
    parser.add_argument("--small-alpha", type=float, default=.1,
                        help="existing conservative small-cell aggregation fraction; does not change quality gates")
    args = parser.parse_args()
    if (not args.levels or not all(3 <= level <= 9 for level in args.levels) or
            not 1 <= args.iterations <= 20000 or not 0 < args.timeout <= 600 or not 0 < args.small_alpha <= .5 or
            not all(math.isfinite(v) and v > 0 for v in (args.speed, args.nu, args.tolerance))):
        parser.error("invalid bounded mesh or flow controls")
    try:
        result = run(args)
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        parser.exit(1, f"native topology verification failed: {exc}\n")
    checks = result["externalCheckMesh"]
    checks = checks if isinstance(checks, list) else [checks]
    good = (result["meshAccepted"] and result["independentTopologyPassed"] and
            "failed" not in checks and (args.mesh_only or result["independentFlowPassed"]))
    raise SystemExit(0 if good else 1)
