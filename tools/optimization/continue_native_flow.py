#!/usr/bin/env python3
"""Additional bounded steady iterations on an identical verified topology mesh.

An iteration-limit field is only an initial guess, never a physical checkpoint.
The original failure is preserved. The new result must converge and pass a fresh
independent audit; all original equation, boundary and stopping controls remain.
"""
import argparse
import copy
import csv
import json
from pathlib import Path
import shutil

from optimize_flow import write_json
import verify_extracted_flow as bridge
import extract_steady_iterate as steady


def continuation_command(command, source, target, iterations, recorded_continuation=False):
    result = list(map(str, command))
    if result[result.index("--case")+1] != "custom":
        raise ValueError("continuation only supports this bridge's steady custom case")
    if any(flag in result for flag in ("--restart", "--time-step")):
        raise ValueError("source is not a single cold steady run: physical-time restart/time stepping is unsupported")
    for flag in ("--initial-guess", "--initial-flux"):
        if flag in result:
            if not recorded_continuation or result.count(flag) != 1:
                raise ValueError("source is not a single cold steady run or a recorded continuation")
            index = result.index(flag)
            del result[index:index+2]
    for flag, expected in (("--mesh", source/"mesh.solver.cm2d"),
                           ("--boundary", source/"flow.boundaries"), ("--output", source/"flow")):
        index = result.index(flag)+1
        if Path(result[index]).resolve() != expected.resolve():
            raise ValueError(f"source {flag} does not match the recorded mesh/field")
        result[index] = str(target/expected.name)
    result[result.index("--max-iterations")+1] = str(iterations)
    result += ["--initial-guess", str(target/"initial.csv"), "--initial-flux", str(target/"initial-flux.csv")]
    return result


def initial_files(mesh, source, target):
    """Use the existing verified same-mesh extractor, including face centres."""
    cells, faces = steady.extract(mesh, steady.rows(source/"flow.cells.csv"), steady.rows(source/"flow.faces.csv"))
    for name, header, rows in (("initial.csv", ["cell", "x", "y", "u", "v", "p"], cells),
                              ("initial-flux.csv", ["face", "owner", "neighbour", "x", "y", "flux"], faces)):
        with (target/name).open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)


def run(args):
    directory, output = args.directory.resolve(strict=True), args.output.resolve()
    if output.exists():
        raise ValueError("preserve earlier results; continuation output must be new")
    report = json.loads((directory/"summary.json").read_text())
    if "components" in report or len(report.get("cases", [])) != 1:
        raise ValueError("choose one recorded native component and one grid")
    if not all(report.get(key) for key in ("meshAccepted", "independentTopologyPassed", "solverQualityPassed")):
        raise ValueError("source mesh must have passed topology, area and solver quality")
    original = report["cases"][0]
    if original.get("flow", {}).get("status") != "iteration_limit" or original["flow"].get("case") != "custom":
        raise ValueError("only a finite steady iteration-limit field may initialise this continuation")
    mesh = Path(original["mesh"]["path"])
    source, level = mesh.parent, original["level"]
    if mesh.name != "mesh.solver.cm2d" or bridge.sha(mesh) != original["mesh"]["sha256"]:
        raise ValueError("source mesh changed")
    flow_cli = bridge.ROOT/"build/cartmesh2d_flow_cli"
    if report["executables"].get(str(flow_cli)) != bridge.sha(flow_cli):
        raise ValueError("source solver binary changed")
    command = original["runs"][-1]["command"]
    if Path(command[0]).resolve() != flow_cli:
        raise ValueError("source command is not the recorded native solver")
    nu, speed = original["flow"]["nu"], original["flow"]["speed"]
    options = bridge.native.argument_parser().parse_args(["--max-iterations", str(original["flow"]["iterations"])])
    audit = bridge.native.verify_case(mesh, source/"flow", "custom", nu, speed, options)
    expected = {"native JSON status is not converged", "native solver did not report convergence"}
    if set(audit["issues"])-expected:
        raise ValueError("source failed independent reconstruction beyond its iteration-limit status")
    prior = report.get("steadyContinuation")
    previous_iterations = prior["totalIterations"] if prior else original["flow"]["iterations"]
    if not isinstance(previous_iterations, int) or previous_iterations < original["flow"]["iterations"]:
        raise ValueError("source continuation lacks a complete iteration count")
    original_directory = (prior.get("originalNativeDirectory", prior["source"]) if prior else str(directory))
    if prior:
        for path, fingerprint in prior["initialFiles"].items():
            if bridge.sha(path) != fingerprint:
                raise ValueError("recorded source initialization changed")
    target = output/f"level-{level}"
    command = continuation_command(command, source, target, args.iterations, recorded_continuation=bool(prior))
    target.mkdir(parents=True)
    inputs = [directory/"summary.json", mesh, source/"flow.cells.csv", source/"flow.faces.csv",
              source/"flow.json", source/"flow.boundaries", Path(__file__), Path(steady.__file__)]
    for name in ("mesh.solver.cm2d", "flow.boundaries", "template.boundaries", "openfoam/solver_quality.json"):
        destination = target/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source/name, destination)
    initial_files(bridge.native.read_cm2d(mesh), source, target)
    report = copy.deepcopy(report)
    case = report["cases"][0]
    case["mesh"]["path"] = str(target/"mesh.solver.cm2d")
    report["acceptedMesh"] = copy.deepcopy(case["mesh"])
    report["nativeFlowConverged"] = report["independentFlowPassed"] = False
    report["issues"] = []
    report["steadyContinuation"] = dict(source=str(directory), inputs={str(p):bridge.sha(p) for p in inputs},
                originalNativeDirectory=original_directory,
                previousIterations=previous_iterations, additionalIterationBudget=args.iterations,
                initialFiles={str(target/name):bridge.sha(target/name) for name in ("initial.csv", "initial-flux.csv")},
                meshRegenerated=False, physicalCheckpoint=False, sourceAuditIssues=audit["issues"],
                note=__doc__)
    case.update(status="solving", flow={})
    case.pop("metrics", None)
    write_json(output/"summary.json", report)
    trial = bridge.execute(command, target/"flow-continuation-run", args.timeout)
    case["runs"].append(trial)
    path = target/"flow.json"
    case["flow"] = json.loads(path.read_text()) if path.exists() else {}
    additional = case["flow"].get("iterations")
    report["steadyContinuation"]["totalIterations"] = (previous_iterations+additional
                                                       if additional is not None else None)
    report["nativeFlowConverged"] = trial["returncode"] == 0 and case["flow"].get("converged") is True
    case["status"] = "flow-not-converged"
    if report["nativeFlowConverged"]:
        options = bridge.native.argument_parser().parse_args(["--max-iterations", str(args.iterations)])
        audit = bridge.native.verify_case(target/"mesh.solver.cm2d", target/"flow", "custom", nu, speed, options)
        write_json(target/"independent-flow.json", audit)
        report["independentFlowPassed"] = audit["valid"]
        case["independentFlowIssues"] = audit["issues"]
        case["status"] = "flow-audited" if audit["valid"] else "flow-audit-failed"
        problem = json.loads((Path(report["source"]).parent/"summary.json").read_text())["problem"]
        case["metrics"] = bridge.pressure_metrics(target/"flow", target/"flow.boundaries", problem, speed, nu)
    if not report["independentFlowPassed"]:
        report["issues"] = ["Continuation did not complete all native convergence and independent audit gates."]
    write_json(output/"summary.json", report)
    print(json.dumps(dict(nativeFlowConverged=report["nativeFlowConverged"],
                         independentFlowPassed=report["independentFlowPassed"],
                         totalIterations=report["steadyContinuation"]["totalIterations"]), indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=2000)
    parser.add_argument("--timeout", type=float, default=300)
    args = parser.parse_args()
    if not 1 <= args.iterations <= 4000 or not 0 < args.timeout <= 600:
        parser.error("invalid bounded continuation budget")
    result = run(args)
    raise SystemExit(0 if result["independentFlowPassed"] else 1)
