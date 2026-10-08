#!/usr/bin/env python3
"""Generate and solve extracted sharp-wall designs with the native tools.

The native mesher and solver decide geometry quality and convergence.
Matched inlet/outlet profiles use the porous model's Reynolds number.
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
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"tools"/"flow"))
import native_mesh as native
from brinkman import port_average
from topology_cases import ports_for


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


def matched_controls(problem,reynolds,speed):
    if not math.isfinite(reynolds) or reynolds<0 or not math.isfinite(speed) or speed<=0:
        raise ValueError("matched controls require nonnegative Re and positive peak speed")
    return dict(speed=speed,nu=(2*speed*problem["port_width"]/3/reynolds if reynolds else 1.),
                momentumInertia=1 if reynolds else 0,viscousStress="laplacian",
                outletMode="prescribed",reynolds=reynolds,
                definition="Re=mean inlet velocity * inlet width / nu; Re=0 denotes exact Stokes equations")


def prescribed_boundaries(source,target,problem,speed):
    ports=ports_for(problem)
    tolerance=1e-11+1e-9*max(problem["width"],problem["height"])
    counts=dict(inletFaces=0,outletFaces=0,closedArtificialOpenings=0)
    rows=[]
    for line in source.read_text().splitlines():
        if not line.startswith("BOUNDARY "):
            rows.append(line);continue
        f=shlex.split(line)
        x,y,sx,sy=map(float,f[3:7]);length=math.hypot(sx,sy)
        selected=None
        for port in ports:
            vertical=port.side in ("left","right")
            coordinate=x if vertical else y
            side=0 if port.side in ("left","bottom") else problem["width"] if port.side=="right" else problem["height"]
            if abs(coordinate-side)>tolerance or abs(sy if vertical else sx)>1e-12*length:
                continue
            tangent=y if vertical else x
            profile=port.peak*port_average(tangent-length/2,tangent+length/2,port.centre,port.width)
            if profile>0:
                if selected is not None:
                    raise ValueError("overlapping prescribed ports")
                selected=(port,profile)
        if selected is None:
            counts["closedArtificialOpenings"]+=int(f[7] in ("velocity-inlet","pressure-outlet"))
            f[7:12]=["wall","wall","0","0","0"]
        else:
            port,profile=selected
            sign=-1 if port.role=="inlet" else 1
            f[7:12]=["velocity-"+port.role,port.name,format(sign*speed*profile*sx/length,".17g"),
                     format(sign*speed*profile*sy/length,".17g"),"0"]
            counts[port.role+"Faces"]+=1
        f[8]=json.dumps(f[8]);rows.append(" ".join(f))
    target.write_text("\n".join(rows)+"\n")
    return counts


def neutral_boundary_template(mesh,measured,target,geometries=None):
    """Explicit all-wall template from actual CM2D geometry, not a duct preset.

    A valid elbow/U-return need not have global left AND right openings. The
    native custom-boundary reader still validates all IDs, owners and geometry
    before it solves the subsequently assigned physical port conditions.
    """
    if geometries is None:
        geometries=native.face_geometry(mesh,measured)
    edges=[e for e in mesh.edges if e.neighbour<0]
    rows=["CARTMESH2D_FLOW_BOUNDARIES 1",f"COUNTS {len(mesh.cells)} {len(mesh.edges)} {len(edges)}"]
    for edge in edges:
        g=geometries[edge.id]
        values=" ".join(format(v,".17g") for v in (*g.centre,*g.area_vector))
        rows.append(f'BOUNDARY {edge.id} {edge.owner} {values} wall "wall" 0 0 0')
    target.write_text("\n".join(rows)+"\nEND\n")


def parabolic_boundaries(source, target, problem, speed, prescribed=False):
    if prescribed:
        return prescribed_boundaries(source,target,problem,speed)
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


def pressure_metrics(flow, boundary, problem, speed, viscosity, momentum_inertia=1):
    with Path(str(flow)+".faces.csv").open() as stream:
        faces = {int(row["face"]):row for row in csv.DictReader(stream)}
    inlet_power, outlet_power, inlet_flux, outlet_flux = 0.0, 0.0, 0.0, 0.0
    kinetic_power=0.
    prescribed_outlet=True
    ports={p.name:p for p in ports_for(problem)}
    for line in boundary.read_text().splitlines():
        if not line.startswith("BOUNDARY "):
            continue
        fields = shlex.split(line)
        face = faces[int(fields[1])]
        flux, pressure = float(face["flux"]), float(face["pressure"])
        if fields[7] == "velocity-inlet":
            inlet_flux -= flux; inlet_power -= pressure*flux
        elif fields[7] in ("pressure-outlet","velocity-outlet"):
            outlet_flux += flux; outlet_power += pressure*flux
            prescribed_outlet &= fields[7]=="velocity-outlet"
        if fields[7] in ("velocity-inlet","velocity-outlet") and fields[8] in ports:
            port=ports[fields[8]]
            tangent=float(fields[4] if port.side in ("left","right") else fields[3])
            length=math.hypot(float(fields[5]),float(fields[6]))
            a=max(tangent-length/2-port.centre,-port.width/2)
            b=min(tangent+length/2-port.centre,port.width/2)
            w=port.width
            primitive=lambda s:s-4*s**3/w**2+48*s**5/(5*w**4)-64*s**7/(7*w**6)
            kinetic_power+=(1 if port.role=="inlet" else -1)*.5*(speed*port.peak)**3*(primitive(b)-primitive(a))
    if inlet_flux <= 0:
        raise ValueError("pressure-drop comparison requires a positive inlet throughput")
    return dict(inletFlux=inlet_flux, outletFlux=outlet_flux,
                pressurePower=inlet_power-outlet_power,
                fluxWeightedPressureDrop=(inlet_power-outlet_power)/inlet_flux,
                speed=speed, viscosity=viscosity,
                nominalReynolds=momentum_inertia*(2/3)*speed*problem["port_width"]/viscosity,
                totalPressurePower=(inlet_power-outlet_power+momentum_inertia*kinetic_power) if prescribed_outlet else None,
                kineticPower=momentum_inertia*kinetic_power if prescribed_outlet else None,
                normalization=viscosity*speed**2/problem.get("viscosity",1),
                dimensionlessTotalPower=(inlet_power-outlet_power+momentum_inertia*kinetic_power)/
                                       (viscosity*speed**2/problem.get("viscosity",1)) if prescribed_outlet else None)


def run(args):
    directory = args.directory.resolve(strict=True)
    root = args.output.resolve()
    extraction = json.loads((directory/"extraction.json").read_text())
    research = json.loads((directory/"summary.json").read_text())
    matching = None
    if getattr(args, "reynolds", None) is not None:
        matching = matched_controls(research["problem"], args.reynolds, args.speed)
        args = copy.copy(args)
        args.nu = matching["nu"]
    root.mkdir(parents=True)
    if args.component is None and extraction["components"] > 1:
        components = []
        for index in range(extraction["components"]):
            child = copy.copy(args)
            child.component, child.output = index, root/f"component-{index}"
            components.append(run(child))
        report = dict(schema="cartmesh2d-topology-native-v2", components=components,
            meshPaddingFraction=getattr(args, "padding_fraction", 1/30),
            scope="Native solves for each disconnected fluid region.",
            physicalAccuracyQualified=False, expectedArea=extraction["fluidArea"],
            totalAcceptedArea=sum(c.get("acceptedMesh", {}).get("area", 0) for c in components))
        for key in ("meshAccepted", "solverQualityPassed", "nativeFlowConverged"):
            report[key] = all(c[key] for c in components)
        (root/"summary.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        return report
    source = directory/"fluid.xy"
    if args.component is not None:
        selected = extraction["regions"][args.component]
        source = directory/selected["file"]
        extraction = dict(extraction, fluidArea=selected["fluidArea"])
    mesh_cli, flow_cli = args.mesh_cli.resolve(strict=True), args.flow_cli.resolve(strict=True)
    report = dict(schema="cartmesh2d-topology-native-v2", source=str(source), sourceSha256=sha(source),
        executables={str(p):sha(p) for p in (mesh_cli, flow_cli)},
        smallCellAggregationFraction=getattr(args, "small_alpha", .1),
        meshPaddingFraction=getattr(args, "padding_fraction", 1/30),
        meshAccepted=False, solverQualityPassed=False, nativeFlowConverged=False,
        physicalAccuracyQualified=False, cases=[], issues=[],
        scope="Native sharp-wall solve; acceptance uses native mesh and solver status.")
    if matching:
        report["matching"] = matching
    def save():
        (root/"summary.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    for level in args.levels:
        case = root/f"level-{level}"
        case.mkdir()
        prefix = case/"mesh"
        item = dict(level=level, runs=[], status="meshing")
        report["cases"].append(item)
        save()
        command = [mesh_cli, source, prefix, level, report["meshPaddingFraction"],
                   report["smallCellAggregationFraction"], "interior", case/"openfoam", level, 0]
        meshing = execute(command, case/"mesh-run", args.timeout)
        item["runs"].append(meshing)
        if meshing["returncode"] != 0:
            item["status"] = "mesh-rejected"
            save()
            continue
        mesh_path = Path(str(prefix)+".solver.cm2d")
        mesh = native.read_cm2d(mesh_path)
        measured = native.measure(mesh)
        item["mesh"] = dict(path=str(mesh_path), sha256=sha(mesh_path), cells=len(mesh.cells),
                            faces=len(mesh.edges), area=measured.total_area)
        report["meshAccepted"] = report["solverQualityPassed"] = True
        report["acceptedMesh"] = item["mesh"]
        item["status"] = "mesh-accepted"
        if args.mesh_only:
            break
        template, boundary = case/"template.boundaries", case/"flow.boundaries"
        if matching:
            neutral_boundary_template(mesh, measured, template)
        else:
            boundary_run = execute([flow_cli, "--mesh", mesh_path, "--case", "duct", "--speed", args.speed,
                                    "--export-boundaries", template], case/"boundary-run", args.timeout)
            item["runs"].append(boundary_run)
            if boundary_run["returncode"] != 0:
                item["status"] = "boundary-rejected"
                save()
                continue
        item["ports"] = parabolic_boundaries(template, boundary, research["problem"], args.speed,
                                               prescribed=matching is not None)
        if not item["ports"]["inletFaces"] or not item["ports"]["outletFaces"]:
            item["status"] = "ports-rejected"
            save()
            continue
        flow = case/"flow"
        command = [flow_cli, "--mesh", mesh_path, "--case", "custom", "--boundary", boundary,
            "--output", flow, "--nu", args.nu, "--speed", args.speed,
            "--max-iterations", args.iterations, "--tolerance", args.tolerance,
            "--velocity-relaxation", (.2 if getattr(args, "velocity_relaxation", None) is None else args.velocity_relaxation),
            "--pressure-corrections", 1,
            "--steady-acceleration", getattr(args, "steady_acceleration", None) or "anderson",
            "--linear-policy", "adaptive", "--pressure-preconditioner", "aggregation"]
        if matching:
            command += ["--momentum-inertia", matching["momentumInertia"],
                        "--viscous-stress", "laplacian", "--convection", "limited-linear"]
        if getattr(args, "convergence", None) is not None:
            command += ["--convergence", args.convergence]
        item["status"] = "solving"
        save()
        solve = execute(command, case/"flow-run", args.timeout)
        item["runs"].append(solve)
        summary_path = Path(str(flow)+".json")
        item["flow"] = json.loads(summary_path.read_text()) if summary_path.exists() else {}
        report["nativeFlowConverged"] = solve["returncode"] == 0 and item["flow"].get("converged") is True
        item["status"] = "flow-converged" if report["nativeFlowConverged"] else "flow-not-converged"
        if report["nativeFlowConverged"]:
            item["metrics"] = pressure_metrics(flow, boundary, research["problem"], args.speed, args.nu,
                                                matching["momentumInertia"] if matching else 1)
            break
        save()
    if not report["meshAccepted"]:
        report["issues"].append("Native mesh generation failed at all requested levels.")
    elif not args.mesh_only and not report["nativeFlowConverged"]:
        report["issues"].append("Native flow did not converge at any requested level.")
    save()
    print(json.dumps({key:report[key] for key in ("meshAccepted", "nativeFlowConverged")}, indent=2))
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
    parser.add_argument("--reynolds",type=float,help="match mean-port Re and prescribed outlet; overrides nu, exact Stokes at 0")
    parser.add_argument("--convergence",choices=["strict","engineering"])
    parser.add_argument("--velocity-relaxation",type=float)
    parser.add_argument("--steady-acceleration",choices=["none","anderson"])
    parser.add_argument("--tolerance", type=float, default=1e-6,
                        help="existing native normalised stopping tolerance, separate from physical accuracy")
    parser.add_argument("--mesh-only", action="store_true")
    parser.add_argument("--component", type=int, help="one extracted region; default processes ALL regions")
    parser.add_argument("--small-alpha", type=float, default=.1,
                        help="existing conservative small-cell aggregation fraction; does not change quality gates")
    parser.add_argument("--padding-fraction",type=float,default=1/30,
                        help="existing mesh CLI background-domain padding; changes grid phase, preserves imported geometry")
    args = parser.parse_args()
    result = run(args)
    raise SystemExit(0 if result["meshAccepted"] and (args.mesh_only or result["nativeFlowConverged"]) else 1)
