#!/usr/bin/env python3
"""Research prototype: constrained 2D fluid topology optimisation.

Example (optional dependencies in requirements.txt):
  python tools/optimization/optimize_flow.py --output outputs/topology/double-pipe

This implements a density filter, smooth projection, discrete adjoint, bounded
optimality-criteria updates and objective backtracking. Continuation stages use
different problems: objective values must only be compared WITHIN a stage or
after reevaluating both designs with the same final parameters.
"""
import argparse
import csv
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
import scipy

from brinkman import Problem, StokesBrinkman


def write_json(path, data):
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
    temporary.replace(path)


def stationarity(model, x, e):
    """Normalised projected KKT residual for the current volume tangent plane.

    Unit step after scaling the objective gradient by its infinity norm; the
    reported norm is in [0,1] design-variable units. This is not a CFD tolerance
    or a global-optimum certificate.
    """
    g = e.gradient/max(np.max(np.abs(e.gradient[model.design])), 1e-30)
    v = e.volume_gradient
    def projected(lam):
        return model.enforce_passive(np.clip(x-g-lam*v, 0, 1))
    candidate = projected(0)
    if e.volume >= model.problem.volume_fraction-1e-10 and v@(candidate-x) > 0:
        lo, hi = 0.0, 1.0
        for _ in range(80):
            if v@(projected(hi)-x) <= 0:
                break
            hi *= 2
        for _ in range(50):
            mid = (lo+hi)/2
            if v@(projected(mid)-x) > 0:
                lo = mid
            else:
                hi = mid
        candidate = projected(hi)
    return float(np.max(np.abs(candidate-x)))


def snapshot(path, model, x, e, q, beta):
    u, v = model.cell_velocity(e.velocity)
    temporary = path.with_suffix(path.suffix+".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, design=x.reshape(model.ny, model.nx),
                            rho=e.rho.reshape(model.ny, model.nx),
                            u=u.reshape(model.ny, model.nx), v=v.reshape(model.ny, model.nx),
                            p=e.pressure.reshape(model.ny, model.nx), face_velocity=e.velocity,
                            gradient=e.gradient.reshape(model.ny, model.nx), q=q, beta=beta,
                            objective=e.objective, volume=e.volume,
                            width=model.problem.width, height=model.problem.height)
    temporary.replace(path)


def vtk(path, model, e):
    """Analysis output only; explicitly not a product fluid mesh / CM2D input."""
    uc, vc = model.cell_velocity(e.velocity)
    with path.open("w") as f:
        f.write("# vtk DataFile Version 3.0\nBrinkman analysis: finite resistance solid, not solver mesh\nASCII\n")
        f.write(f"DATASET STRUCTURED_POINTS\nDIMENSIONS {model.nx+1} {model.ny+1} 1\n")
        f.write(f"ORIGIN 0 0 0\nSPACING {model.dx:.17g} {model.dy:.17g} 1\nCELL_DATA {model.cells}\n")
        for name, values in (("fluid_fraction", e.rho), ("kinematic_pressure", e.pressure)):
            f.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n")
            f.writelines(f"{a:.17g}\n" for a in values)
        f.write("VECTORS velocity double\n")
        f.writelines(f"{u:.17g} {v:.17g} 0\n" for u, v in zip(uc, vc))


def metrics(model, x, e):
    return dict(objective=e.objective, volumeFraction=e.volume, grayness=e.grayness,
                dissipation=e.dissipation, pressurePower=e.pressure_power,
                fluxWeightedPressureDrop=e.pressure_power/model.inflow,
                maxSolidToGlobalSpeed=e.solid_speed_fraction,
                linearResidual=e.linear_residual, continuity=e.continuity,
                adjointResidual=e.adjoint_residual, projectedKkt=stationarity(model, x, e))


def reference_design(model, beta):
    """Reproducible geometric seed, reevaluated at the SAME final model."""
    p = model.problem
    x = np.zeros(model.cells)
    half_width = p.volume_fraction*p.height/(4 if p.case == "double-pipe" else 2)
    for j in range(model.ny):
        y = (j+0.5)*model.dy
        for i in range(model.nx):
            t = (i+0.5)/model.nx
            centres = model.left_ports if p.case == "double-pipe" else [p.height*(.25+.5*t)]
            x[model.c(i, j)] = float(any(abs(y-c) <= half_width for c in centres))
    return model.feasible_design(x, beta)


def run(args):
    root = args.output.resolve()
    if root.exists():
        raise ValueError("output already exists; choose a fresh directory to preserve evidence")
    problem = Problem(nx=args.nx, ny=args.ny, width=args.width, height=1,
                      volume_fraction=args.volume, filter_radius=args.filter_radius,
                      alpha_max=args.alpha_max, case=args.case, port_width=args.port_width)
    model = StokesBrinkman(problem)
    # Construct/validate the initial design before creating any output directory.
    x = model.uniform_design()
    root.mkdir(parents=True)
    (root/"snapshots").mkdir()
    os.environ.setdefault("MPLCONFIGDIR", str(root/"plot-cache"))
    stages = [(0.01, 0.0), (0.1, 2.0), (0.1, 6.0)]
    budgets = args.iterations if len(args.iterations) == 3 else args.iterations*3
    start = time.monotonic()
    report = dict(schema="cartmesh2d-fluid-topology-v1", experimental=True,
                  model="nondimensional 2D Stokes-Brinkman, staggered finite volume",
                  objective=args.objective, problem=asdict(problem),
                  optimizer="filtered/projected OC with feasible-volume bisection and backtracking",
                  controls=dict(stages=stages, stageIterations=budgets, move=args.move,
                                stationarityTolerance=args.stationarity_tolerance,
                                maxSeconds=args.max_seconds),
                  analysisIsProductFluidMesh=False, physicalAccuracyQualified=False,
                  optimizationConverged=False, status="running", stages=[], history=[],
                  environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                  sourceHashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                                [Path(__file__), Path(__file__).with_name("brinkman.py")]})
    global_iteration = 0
    current = None
    interrupted = False
    failure = None
    try:
        for stage, ((q, beta), budget) in enumerate(zip(stages, budgets)):
            stage_x = model.feasible_design(x, beta)
            stage_evaluation = model.evaluate(stage_x, q, beta, args.objective)
            x, current = stage_x, stage_evaluation
            record = dict(stage=stage, q=q, beta=beta, initial=metrics(model, x, current),
                          status="iteration-limit", accepted=0, rejectedTrials=0)
            report["stages"].append(record)
            snapshot(root/"accepted-design.npz", model, x, current, q, beta)
            snapshot(root/"snapshots"/f"stage-{stage}-initial.npz", model, x, current, q, beta)
            for iteration in range(budget):
                row = dict(iteration=global_iteration, stage=stage, stageIteration=iteration,
                           elapsedSeconds=time.monotonic()-start, **metrics(model, x, current))
                report["history"].append(row)
                write_json(root/"progress.json", report)
                if row["projectedKkt"] <= args.stationarity_tolerance:
                    record["status"] = "stationary"
                    break
                if time.monotonic()-start >= args.max_seconds:
                    record["status"] = "time-limit"
                    break
                move = args.move
                accepted = False
                for trial in range(10):
                    proposal = model.oc_candidate(x, current, beta, move)
                    candidate = model.evaluate(proposal, q, beta, args.objective)
                    # Small allowance is floating-point objective repeatability,
                    # not a relaxation of any physical or solver acceptance gate.
                    if (candidate.volume <= problem.volume_fraction+1e-10 and
                            candidate.objective <= current.objective+1e-11*max(current.objective, 1.0)):
                        accepted = True
                        break
                    move *= .5
                    record["rejectedTrials"] += 1
                if not accepted:
                    record["status"] = "line-search-stalled"
                    break
                row["acceptedChange"] = float(np.max(np.abs(proposal-x)))
                x, current = proposal, candidate
                snapshot(root/"accepted-design.npz", model, x, current, q, beta)
                record["accepted"] += 1
                global_iteration += 1
                if global_iteration % args.snapshot_every == 0:
                    snapshot(root/"snapshots"/f"accepted-{global_iteration:04d}.npz", model, x, current, q, beta)
                print(f"stage={stage} step={iteration+1} J={current.objective:.8g} "
                      f"volume={current.volume:.7f} gray={current.grayness:.4f}", flush=True)
            record["final"] = metrics(model, x, current)
            if record["final"]["projectedKkt"] <= args.stationarity_tolerance:
                record["status"] = "stationary"
            snapshot(root/"snapshots"/f"stage-{stage}-final.npz", model, x, current, q, beta)
            if record["status"] == "time-limit":
                break
    except KeyboardInterrupt:
        interrupted = True
    except (ArithmeticError, RuntimeError) as exc:
        failure = str(exc)
    if current is None:
        report.update(status="initial-analysis-failed", issues=[failure or "interrupted before first accepted state"])
        write_json(root/"summary.json", report)
        raise RuntimeError("no accepted analysis state available")
    stage = report["stages"][-1]
    q, beta = stage["q"], stage["beta"]
    report["status"] = "analysis-failed" if failure else ("interrupted" if interrupted else stage["status"])
    report["optimizationConverged"] = (not interrupted and not failure and len(report["stages"]) == len(stages)
                                        and report["status"] == "stationary")
    report["final"] = metrics(model, x, current)
    stage["final"] = report["final"]
    report["finalParameters"] = dict(q=q, beta=beta)
    snapshot(root/"final.npz", model, x, current, q, beta)
    vtk(root/"final.analysis.vtk", model, current)
    if failure or interrupted:
        report["issues"] = [failure] if failure else ["User interrupted; last accepted design preserved."]
        report["elapsedSeconds"] = time.monotonic()-start
        write_json(root/"summary.json", report)
        write_json(root/"progress.json", report)
        return report
    # Equal-volume, equal-parameter comparisons; continuation improvements are
    # not reported as if they had used one unchanged objective.
    comparisons = {}
    for name, baseline in (("uniformPorous", model.uniform_design(beta)),
                           ("geometricSeed", reference_design(model, beta))):
        evaluation = model.evaluate(baseline, q, beta, args.objective)
        comparisons[name] = metrics(model, baseline, evaluation)
        comparisons[name]["relativeObjectiveReduction"] = 1-current.objective/evaluation.objective
        snapshot(root/f"reference-{name}.npz", model, baseline, evaluation, q, beta)
    report["sameModelReferences"] = comparisons
    report["elapsedSeconds"] = time.monotonic()-start
    report["limits"] = ["Finite-resistance porous analysis; solid leakage is reported, not suppressed.",
                        "Stokes creeping flow; no convective inertia, turbulence or compressibility.",
                        "Filter radius is not a certified manufacturing minimum width.",
                        "Local stationarity is not a global optimum or physical accuracy certificate.",
                        "Threshold extraction changes the model and requires separate sharp-wall CFD."]
    with (root/"history.csv").open("w", newline="") as f:
        keys = sorted(set().union(*(row.keys() for row in report["history"])))
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader(); writer.writerows(report["history"])
    write_json(root/"summary.json", report)
    write_json(root/"progress.json", report)
    if not args.no_plot:
        from topology_artifacts import render, extract
        extracted = extract(root, args.threshold)
        report["extraction"] = extracted
        write_json(root/"summary.json", report)
        render(root)
    print(json.dumps({k:report[k] for k in ("status", "optimizationConverged", "elapsedSeconds", "final")}, indent=2))
    return report


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--case", choices=["double-pipe", "bend"], default="double-pipe")
    p.add_argument("--nx", type=int, default=48)
    p.add_argument("--ny", type=int, default=32)
    p.add_argument("--width", type=float, default=1.5)
    p.add_argument("--port-width", type=float, default=1/6)
    p.add_argument("--volume", type=float, default=1/3)
    p.add_argument("--filter-radius", type=float, default=.06)
    p.add_argument("--alpha-max", type=float, default=25000)
    p.add_argument("--objective", choices=["dissipation", "pressure-power"], default="dissipation")
    p.add_argument("--iterations", nargs="+", type=int, default=[25, 40, 60])
    p.add_argument("--move", type=float, default=.15)
    p.add_argument("--stationarity-tolerance", type=float, default=1e-3,
                   help="normalised projected KKT design-step norm; not a CFD accuracy threshold")
    p.add_argument("--max-seconds", type=float, default=600)
    p.add_argument("--snapshot-every", type=int, default=5)
    p.add_argument("--threshold", type=float, default=.5)
    p.add_argument("--no-plot", action="store_true",
                   help="analysis only; skip optional contour extraction and plotting")
    return p


if __name__ == "__main__":
    p = parser()
    args = p.parse_args()
    if (len(args.iterations) not in (1, 3) or not all(1 <= n <= 1000 for n in args.iterations) or
            not 0 < args.move <= 1 or not 0 < args.stationarity_tolerance < 1 or
            not 0 < args.max_seconds <= 7200 or not 1 <= args.snapshot_every <= 1000 or
            not 0 < args.threshold < 1):
        p.error("invalid iteration, update, time, output or threshold controls")
    try:
        result = run(args)
        if result["status"] in ("analysis-failed", "initial-analysis-failed", "interrupted"):
            raise SystemExit(2)
    except (ValueError, ArithmeticError, OSError, RuntimeError) as exc:
        p.exit(1, f"topology optimisation failed: {exc}\n")
