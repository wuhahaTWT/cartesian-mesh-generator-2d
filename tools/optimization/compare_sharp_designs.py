#!/usr/bin/env python3
"""Equal-area, equal-port comparison of a density design and its geometric seed.

The scalar interior contour offset enforces the SHARP fluid area. It does not
smooth imported XY, discard islands or change passive port collars. Both designs
are recalculated with identical native equations, boundary profiles and controls.
Every requested resolution is retained, including rejected meshes and solves.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np

from brinkman import Problem, StokesBrinkman
from optimize_flow import write_json
from topology_artifacts import contours, extract, signed_area
import verify_extracted_flow as native_bridge


def area(groups):
    return sum(abs(signed_area(group[0]))-sum(abs(signed_area(hole)) for hole in group[1:])
               for group in groups)


def match_sharp_area(rho, model, target):
    """Monotone interior level-set shift with fixed passive collars and walls.

    Geometry tolerance is the existing native-reader area tolerance, in squared
    design length units. It is unrelated to any CFD or optimization precision.
    """
    rho = np.array(rho, copy=True, dtype=float)
    if (rho.shape != (model.ny, model.nx) or not np.isfinite(rho).all() or
            np.any(rho < 0) or np.any(rho > 1)):
        raise ValueError("invalid density field")
    if not np.isfinite(target) or target <= 0:
        raise ValueError("target sharp area must be positive and finite")
    mask = model.design.reshape(rho.shape)
    fixed = model.fixed_design.reshape(rho.shape)
    if not np.array_equal(rho[~mask], fixed[~mask]):
        raise ValueError("input density does not preserve the problem's passive cells")
    def sample(shift):
        field = rho.copy()
        field[mask] = np.clip(rho[mask]+shift, 0, 1)
        groups = contours(field, model.problem.width, model.problem.height)
        return field, area(groups)
    lo, hi = -1., 1.
    if not sample(lo)[1] < target < sample(hi)[1]:
        raise ValueError("sharp area is infeasible with the fixed collars")
    tolerance = 1e-11+1e-9*abs(target)
    before = sample(0)[1]
    for _ in range(55):
        shift = (lo+hi)/2
        field, measured = sample(shift)
        if abs(measured-target) <= tolerance:
            return field, dict(targetArea=target, originalArea=before, matchedArea=measured,
                               interiorDensityOffset=shift, areaTolerance=tolerance,
                               passiveCellsUnchanged=bool(np.array_equal(field[~mask], rho[~mask])))
        if measured < target:
            lo = shift
        else:
            hi = shift
    raise ArithmeticError("sharp-area contour projection did not reach its geometry tolerance")


def prepare_design(directory, source_field, model, target, label):
    with np.load(source_field) as archive:
        rho = archive["rho"].copy()
        if float(archive["width"]) != model.problem.width or float(archive["height"]) != model.problem.height:
            raise ValueError("field domain differs from the recorded problem")
    field, matching = match_sharp_area(rho, model, target)
    directory.mkdir(parents=True)
    # This is an extraction field, not a new porous solution: do not copy u/p.
    np.savez_compressed(directory/"final.npz", rho=field, width=model.problem.width, height=model.problem.height)
    extracted = extract(directory)
    metadata = dict(schema="cartmesh2d-sharp-design-v1", label=label, problem=asdict(model.problem),
                    sourceField=str(source_field.resolve()),
                    sourceFieldSha256=hashlib.sha256(source_field.read_bytes()).hexdigest(),
                    extractionOnly=True, matching=matching, extraction=extracted)
    write_json(directory/"summary.json", metadata)
    return metadata


def flow_metrics(report):
    if not (report.get("independentTopologyPassed") and report.get("nativeFlowConverged") and
            report.get("independentFlowPassed")):
        return None
    leaves = report.get("components", [report])
    check = [leaf["externalCheckMesh"] for leaf in leaves]
    if "failed" in check:
        return None
    cases = [next(case for case in leaf["cases"] if case["status"] == "flow-audited") for leaf in leaves]
    flux = sum(case["metrics"]["inletFlux"] for case in cases)
    power = sum(case["metrics"]["pressurePower"] for case in cases)
    return dict(pressurePower=power, inletFlux=flux, fluxWeightedPressureDrop=power/flux,
                cells=sum(case["mesh"]["cells"] for case in cases),
                area=sum(case["mesh"]["area"] for case in cases),
                iterations=[case["flow"]["iterations"] for case in cases],
                externalCheckMesh=check)


def assess(rows):
    """Observed ranking stability, NOT a rigorous discretization error bound."""
    paired = []
    for row in rows:
        baseline, candidate = row.get("baseline"), row.get("candidate")
        if baseline is None or candidate is None:
            continue
        values = [baseline["inletFlux"], candidate["inletFlux"], baseline["fluxWeightedPressureDrop"],
                  candidate["fluxWeightedPressureDrop"]]
        if not all(np.isfinite(value) and value > 0 for value in values):
            raise ValueError("comparison requires positive finite flow and pressure drop")
        # Equal inlet flow is part of the experiment, not an objective benefit.
        if abs(baseline["inletFlux"]-candidate["inletFlux"]) > 1e-10*baseline["inletFlux"]:
            raise ValueError("unequal inlet flux invalidates the design comparison")
        b, c = baseline["fluxWeightedPressureDrop"], candidate["fluxWeightedPressureDrop"]
        paired.append(dict(level=row["level"], baselineDrop=b, candidateDrop=c,
                           relativeReduction=1-c/b, difference=b-c))
    unpaired = [row["level"] for row in rows if row.get("baseline") is None or row.get("candidate") is None]
    result = dict(paired=paired, unpairedLevels=unpaired, meshRobustImprovementObserved=False,
                  physicalAccuracyQualified=False, status="incomplete-native-comparison")
    if not paired:
        return result
    result["status"] = "single-grid-comparison"
    if len(paired) >= 2:
        previous, final = paired[-2:]
        change_b = abs(final["baselineDrop"]-previous["baselineDrop"])
        change_c = abs(final["candidateDrop"]-previous["candidateDrop"])
        margin = final["difference"]-change_b-change_c
        result.update(observedBaselineGridChange=change_b, observedCandidateGridChange=change_c,
                      improvementBeyondObservedGridChanges=margin,
                      meshRobustImprovementObserved=bool(margin > 0 and previous["difference"] > 0),
                      status="improvement-observed" if margin > 0 and previous["difference"] > 0 else "no-robust-improvement",
                      interpretation="Two-grid changes are observed sensitivity, not an error estimator or grid-independence proof.")
    if unpaired:
        result["status"] = "partial-native-comparison"
        result["meshRobustImprovementObserved"] = False
    return result


def run(args):
    source, root = args.directory.resolve(strict=True), args.output.resolve()
    if root.exists():
        raise ValueError("comparison output exists; choose a fresh directory")
    research = json.loads((source/"summary.json").read_text())
    model = StokesBrinkman(Problem(**research["problem"]))
    budget = model.problem.width*model.problem.height*model.problem.volume_fraction
    target = budget
    if args.area_target == "candidate":
        with np.load(source/args.candidate) as archive:
            target = area(contours(archive["rho"], model.problem.width, model.problem.height))
        if target > budget+1e-11+1e-9*budget:
            raise ValueError("candidate sharp area exceeds the physical budget; use budget projection explicitly")
    root.mkdir(parents=True)
    report = dict(schema="cartmesh2d-sharp-comparison-v1", status="running", problem=research["problem"],
                  targetArea=target, areaBudget=budget, areaTarget=args.area_target,
                  controls=dict(levels=args.levels, speed=args.speed, nu=args.nu,
                                tolerance=args.tolerance, iterations=args.iterations, smallAlpha=args.small_alpha),
                  metric="inlet/outlet flux-weighted static kinematic pressure drop; pressure power per density and depth",
                  sourceHashes={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in
                                [Path(__file__), Path(native_bridge.__file__)]},
                  physicalAccuracyQualified=False, designs={}, rows=[], issues=[])
    for label, filename in (("baseline", args.baseline), ("candidate", args.candidate)):
        report["designs"][label] = prepare_design(root/label, source/filename, model, target, label)
    write_json(root/"summary.json", report)
    for level in args.levels:
        row = dict(level=level)
        report["rows"].append(row)
        for label in ("baseline", "candidate"):
            output = root/f"{label}-level-{level}"
            controls = argparse.Namespace(directory=root/label, output=output, mesh_cli=args.mesh_cli,
                        flow_cli=args.flow_cli, levels=[level], timeout=args.timeout, iterations=args.iterations,
                        speed=args.speed, nu=args.nu, tolerance=args.tolerance, small_alpha=args.small_alpha,
                        mesh_only=False, component=None)
            try:
                native = native_bridge.run(controls)
                row[label] = flow_metrics(native)
            except (ValueError, RuntimeError, OSError, ArithmeticError) as exc:
                row[label] = None
                report["issues"].append(dict(design=label, level=level, issue=str(exc)))
            row[label+"Output"] = str(output)
            write_json(root/"summary.json", report)
        report["assessment"] = assess(report["rows"])
        write_json(root/"summary.json", report)
    report["status"] = report["assessment"]["status"]
    write_json(root/"summary.json", report)
    print(json.dumps(report["assessment"], indent=2))
    return report


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--baseline", default="reference-geometricSeed.npz")
    p.add_argument("--candidate", default="final.npz")
    p.add_argument("--area-target", choices=["candidate", "budget"], default="candidate",
                   help="match the baseline to the unchanged candidate area, or project both to the area budget")
    p.add_argument("--levels", type=int, nargs="+", default=[6, 7])
    p.add_argument("--mesh-cli", type=Path, default=native_bridge.ROOT/"build/cartmesh2d_cli")
    p.add_argument("--flow-cli", type=Path, default=native_bridge.ROOT/"build/cartmesh2d_flow_cli")
    p.add_argument("--timeout", type=float, default=180)
    p.add_argument("--iterations", type=int, default=2500)
    p.add_argument("--speed", type=float, default=.02)
    p.add_argument("--nu", type=float, default=1)
    p.add_argument("--tolerance", type=float, default=1e-8)
    p.add_argument("--small-alpha", type=float, default=.1)
    return p


if __name__ == "__main__":
    p = parser(); args = p.parse_args()
    if (args.levels != sorted(set(args.levels)) or not all(3 <= level <= 9 for level in args.levels) or
            not 1 <= args.iterations <= 20000 or not 0 < args.timeout <= 600 or not 0 < args.small_alpha <= .5 or
            not all(np.isfinite(v) and v > 0 for v in (args.nu, args.speed, args.tolerance))):
        p.error("invalid bounded comparison controls")
    try:
        run(args)
    except (ValueError, RuntimeError, OSError, ArithmeticError) as exc:
        p.exit(1, f"sharp-wall comparison failed: {exc}\n")
