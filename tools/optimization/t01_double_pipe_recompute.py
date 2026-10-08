#!/usr/bin/env python3
"""Recompute the fixed T01 double-pipe inputs with the existing solvers.

The independent reference is constructed from four corners per rectangular
channel, not from the optimized density. Native solves run one Reynolds number
and one level at a time so the coarse result can be inspected before refinement.
No optimization, equation reconstruction or analytical replacement is performed.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import time
import zipfile

import numpy as np

from brinkman import Problem
from compare_sharp_designs import area, flow_metrics
from engineering_baselines import distance_to_path
from fidelity_metrics import leakage
from navier_stokes_brinkman import NavierStokesBrinkman
import native_flow
from topology_artifacts import contours, port_connectivity, signed_area


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = ROOT / "outputs/topology-fidelity/b-retries"
DEFAULT_OUTPUT = ROOT / "outputs/t01-agent"
DEFAULT_ARCHIVE = ROOT / "artifacts/current/laminar-foundation/t01-originals.zip"


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def rectangles(problem):
    width, height, aperture = (problem[k] for k in ("width", "height", "port_width"))
    return [[np.array([[0., centre-aperture/2], [width, centre-aperture/2],
                       [width, centre+aperture/2], [0., centre+aperture/2]])]
            for centre in (height/4, 3*height/4)]


def boundary_text(groups):
    return ("# Independent four-corner rectangular fluid interiors; explicit interior mode.\n"
            + "\n\n".join("\n".join(f"{x:.17g} {y:.17g}" for x, y in group[0])
                             for group in groups) + "\n")


def prepare(args):
    started = time.monotonic()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    record = dict(schema="cartmesh2d-t01-double-pipe-inputs-v1", physicalAccuracyQualified=False,
                  scope="Frozen-version evidence; rerun after the general CFD repair.",
                  source=str(args.source.resolve()), cases=[], sourceHashes={})
    if args.archive:
        with zipfile.ZipFile(args.archive) as archive:
            entries = []
            for member in archive.namelist():
                if member.endswith("/"):
                    continue
                local = ROOT/member
                same = local.is_file() and local.read_bytes() == archive.read(member)
                entries.append(dict(member=member, unpackedFileMatches=same))
            if not all(entry["unpackedFileMatches"] for entry in entries):
                raise ValueError("unpacked original archive differs from a recorded member")
        record["archive"] = dict(path=str(args.archive.resolve()), sha256=native_flow.sha(args.archive),
                                 members=entries)
    executable_dir = output/"bin"
    executable_dir.mkdir(exist_ok=True)
    frozen_mesh = executable_dir/"mesh-cli"
    if frozen_mesh.exists() and native_flow.sha(frozen_mesh) != native_flow.sha(args.mesh_cli):
        raise ValueError("existing frozen mesher differs from the requested mesher")
    if not frozen_mesh.exists():
        shutil.copy2(args.mesh_cli, frozen_mesh)
    record["executables"] = {str(p.resolve()): native_flow.sha(p)
                             for p in (frozen_mesh, args.flow_cli)}
    for reynolds in (0, 50):
        source = args.source/f"double-pipe-re-{reynolds}"
        original = read(source/"summary.json")
        problem = original["problem"]
        if problem["case"] != "double-pipe" or problem["reynolds"] != reynolds:
            raise ValueError("wrong T01 problem")
        model = NavierStokesBrinkman(Problem(**problem))
        groups = rectangles(problem)
        design_dir = output/f"independent-re{reynolds}"
        design_dir.mkdir(exist_ok=True)
        (design_dir/"fluid.xy").write_text(boundary_text(groups))
        regions = []
        for index, group in enumerate(groups):
            path = design_dir/f"fluid-component-{index}.xy"
            path.write_text(boundary_text([group]))
            regions.append(dict(file=path.name, holes=0, fluidArea=abs(signed_area(group[0])),
                                boundarySha256=native_flow.sha(path)))
        extraction = dict(schema="cartmesh2d-topology-extraction-v1", construction="independent rectangles",
                          fluidRegion="interior", components=2, holes=0, regions=regions, vertices=8,
                          fluidArea=area(groups), designArea=problem["width"]*problem["height"],
                          boundarySha256=native_flow.sha(design_dir/"fluid.xy"),
                          topologyQualified=False, solverQualityQualified=False,
                          portConnectivity=port_connectivity(groups, problem))
        write(design_dir/"extraction.json", extraction)
        write(design_dir/"summary.json", dict(problem=problem, construction=extraction["construction"]))
        entry = dict(reynolds=reynolds, problem=problem, originalStatus=original["status"],
                     originalOptimizationSeconds=original["elapsedSeconds"],
                     independentGeometry=extraction, fields=[])
        fields = {}
        for label, filename in (("candidate", "final.npz"), ("geometricSeed", "reference-geometricSeed.npz"),
                                ("engineering", "reference-engineering.npz")):
            path = source/filename
            with np.load(path) as stored:
                rho, design = stored["rho"].copy(), stored["design"].ravel().copy()
                q, beta = float(stored["q"]), float(stored["beta"])
                objective = float(stored["objective"])
            fields[label] = rho
            began = time.monotonic()
            evaluated = model.evaluate(design, q, beta, original["objective"])
            elapsed = time.monotonic()-began
            raw_groups = contours(rho, problem["width"], problem["height"])
            geometry = dict(components=len(raw_groups), holes=sum(len(g)-1 for g in raw_groups),
                            vertices=sum(len(loop) for g in raw_groups for loop in g),
                            rawSharpArea=area(raw_groups),
                            portConnectivity=port_connectivity(raw_groups, problem))
            if label != "engineering":
                ordered = sorted(raw_groups, key=lambda g: np.mean(g[0][:, 1]))
                geometry["maximumDistanceToIndependentRectangle"] = max(
                    float(np.max(distance_to_path(g[0], np.vstack([r[0], r[0][0]]))))
                    for g, r in zip(ordered, groups)) if len(ordered) == 2 else None
                geometry["areaDifferenceFromIndependent"] = area(raw_groups)-area(groups)
            controls = native_flow.matched_controls(problem, reynolds, .02)
            normalization = controls["nu"]*.02**2/problem["viscosity"]
            entry["fields"].append(dict(label=label, file=str(path.resolve()), sha256=native_flow.sha(path),
                q=q, beta=beta, storedObjective=objective, porousObjective=evaluated.objective,
                objectiveDifference=evaluated.objective-objective,
                recomputedRhoMaximumDifference=float(np.max(np.abs(evaluated.rho.reshape(rho.shape)-rho))),
                densityVolumeFraction=float(rho.mean()), geometry=geometry, porousSolveSeconds=elapsed,
                porousLinearResidual=evaluated.linear_residual, porousContinuity=evaluated.continuity,
                porousAdjointResidual=evaluated.adjoint_residual, newtonHistory=model.newton_history,
                porousLeakage=leakage(model, evaluated, q), physicalNormalization=normalization,
                scaledPorousTotalPressurePower=evaluated.objective*normalization,
                scaledPorousBoundaryPressureDrop=evaluated.objective*normalization/(model.inflow*.02),
                normalizationNote="J_B converted to the native speed/viscosity scale; not a sharp-wall CFD result."))
        entry["candidateVsSeedRhoMaximumDifference"] = float(np.max(np.abs(fields["candidate"]-fields["geometricSeed"])))
        entry["candidateVsSeedRhoUnequalValues"] = int(np.count_nonzero(fields["candidate"] != fields["geometricSeed"]))
        record["cases"].append(entry)
    for source in (Path(__file__), Path(native_flow.__file__), ROOT/"tools/optimization/brinkman.py",
                   ROOT/"tools/optimization/navier_stokes_brinkman.py", ROOT/"tools/optimization/topology_artifacts.py"):
        record["sourceHashes"][str(source.relative_to(ROOT))] = native_flow.sha(source)
    record["preparationSeconds"] = time.monotonic()-started
    write(output/"inputs.json", record)
    print(json.dumps(dict(output=str(output), preparationSeconds=record["preparationSeconds"],
                         rhoDifferences=[e["candidateVsSeedRhoMaximumDifference"] for e in record["cases"]])))


def run_native(args):
    prepared = read(args.output/"inputs.json")
    mesh_cli = args.output/"bin/mesh-cli"
    for executable in (mesh_cli, args.flow_cli):
        if prepared["executables"].get(str(executable.resolve())) != native_flow.sha(executable):
            raise ValueError("native executable differs from input manifest")
    directory = args.output/f"independent-re{args.reynolds}"
    target = args.output/f"independent-re{args.reynolds}-level-{args.level}"
    controls = argparse.Namespace(directory=directory, output=target, mesh_cli=mesh_cli,
        flow_cli=args.flow_cli, levels=[args.level], timeout=args.timeout, iterations=args.iterations,
        speed=.02, nu=1., reynolds=args.reynolds, tolerance=1e-8, small_alpha=.1,
        padding_fraction=1/30, mesh_only=False, component=None, convergence="strict",
        velocity_relaxation=.2, steady_acceleration="anderson")
    began = time.monotonic()
    native = native_flow.run(controls)
    elapsed = time.monotonic()-began
    metric = flow_metrics(native)
    write(target/"run-cost.json", dict(totalElapsedSeconds=elapsed, metrics=metric,
        scope="Geometry reading, both native meshing/solves, export and Python postprocessing; cold start, no coarse-map."))
    print(json.dumps(dict(reynolds=args.reynolds, level=args.level, metrics=metric, totalElapsedSeconds=elapsed)))
    if metric is None:
        raise SystemExit(1)


def collect_native(directory, label, reynolds, level, prepared):
    summary = read(directory/"summary.json")
    metric = flow_metrics(summary)
    if metric is None:
        raise ValueError(f"native case is not converged: {directory}")
    problem = next(e["problem"] for e in prepared["cases"] if e["reynolds"] == reynolds)
    expected = set(prepared["executables"].values())
    components, mesh_seconds, flow_seconds = [], 0., 0.
    leaves = summary.get("components", [summary])
    recomputed_power = 0.
    for leaf in leaves:
        if set(leaf["executables"].values()) != expected:
            raise ValueError("cannot compare different native executables")
        case = next(c for c in leaf["cases"] if c["status"] == "flow-converged")
        mesh = Path(case["mesh"]["path"])
        if native_flow.sha(mesh) != case["mesh"]["sha256"]:
            raise ValueError("native mesh differs from its solve provenance")
        flow = case["flow"]
        required = dict(convergenceMode="strict", tolerance=1e-8, velocityRelaxation=.2,
            pressureCorrectionPasses=1, pressureRelaxation=.25, coupling="simple",
            steadyAcceleration="anderson", andersonHistory=4,
            andersonStart=10, pressurePreconditioner="aggregation", adaptiveLinear=True,
            strictLinearFinal=True, convection="limited-linear", viscousStress="laplacian",
            momentumInertia=0 if reynolds == 0 else 1, speed=.02,
            nu=native_flow.matched_controls(problem, reynolds, .02)["nu"])
        if any(flow.get(key) != value for key, value in required.items()):
            raise ValueError(f"different native numerical controls: {directory}")
        measured = native_flow.pressure_metrics(mesh.parent/"flow", mesh.parent/"flow.boundaries",
                        problem, flow["speed"], flow["nu"], flow["momentumInertia"])
        recomputed_power += measured["totalPressurePower"]
        sizing = read(mesh.parent/"mesh.sizing.json")
        resolution = read(mesh.parent/"mesh.resolution.json")
        quality = read(mesh.parent/"mesh.construction-quality.json")
        for run in case["runs"]:
            if run["returncode"] != 0 or run["timedOut"]:
                raise ValueError("reported converged case contains a failed native process")
            if "--mesh" in run["command"]:
                flow_seconds += run["elapsedSeconds"]
            else:
                mesh_seconds += run["elapsedSeconds"]
        residuals = {k:flow[k] for k in ("continuity", "globalRelativeImbalance", "momentumResidual",
                                        "velocityChange", "pressureChange")}
        with (mesh.parent/"flow.cells.csv").open() as stream:
            field = list(csv.DictReader(stream))
        extrema = dict(minU=min(float(row["u"]) for row in field), maxU=max(float(row["u"]) for row in field),
                       maxAbsV=max(abs(float(row["v"])) for row in field),
                       minPressure=min(float(row["p"]) for row in field),
                       maxPressure=max(float(row["p"]) for row in field))
        components.append(dict(cells=case["mesh"]["cells"], faces=case["mesh"]["faces"],
            area=case["mesh"]["area"], iterations=flow["iterations"], residuals=residuals,
            fieldExtrema=extrema,
            meshPath=str(mesh), meshSha256=native_flow.sha(mesh),
            flowSummarySha256=native_flow.sha(mesh.parent/"flow.json"),
            boundarySha256=native_flow.sha(mesh.parent/"flow.boundaries"),
            meshDomain=sizing["domain"], levelHistogram=sizing["level_histogram"],
            referenceLength=resolution["reference_length"], actualResolution=resolution["actual"],
            constructionCounts=quality["counts"], constructionQuality=quality["quality"],
            meshTopologyValid=quality["valid"], solverQualityPassed=leaf["solverQualityPassed"]))
    end_to_end = read(directory/"run-cost.json")["totalElapsedSeconds"] if (directory/"run-cost.json").exists() else None
    return dict(reynolds=reynolds, label=label, level=level, **metric,
        source=str(directory/"summary.json"), summarySha256=native_flow.sha(directory/"summary.json"),
        components=components, effectiveH=math.sqrt(metric["area"]/metric["cells"]),
        meshSeconds=mesh_seconds, flowSeconds=flow_seconds,
        nativeProcessSeconds=mesh_seconds+flow_seconds, measuredEndToEndSeconds=end_to_end,
        postprocessedPowerDifference=recomputed_power-metric["totalPressurePower"],
        nativeFlowConverged=True, physicalAccuracyQualified=False)


def make_plots(report, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    colors = dict(candidate="#b04492", engineering="#157d9c", independent="#222222")
    labels = dict(candidate="Original candidate", engineering="Equal-area engineering", independent="Independent rectangles")
    figure, axes = plt.subplots(3, 2, figsize=(12, 12), constrained_layout=True)
    for column, reynolds in enumerate((0, 50)):
        ranking = next(item for item in report["comparisons"] if item["reynolds"] == reynolds)
        native = [row for row in report["native"] if row["reynolds"] == reynolds]
        porous = {row["label"]:row for row in report["inputs"]["cases"][column]["fields"]}
        x = np.arange(3)
        for offset, label in ((-.18, "candidate"), (.18, "engineering")):
            values = [porous[label]["porousObjective"]/porous["geometricSeed"]["porousObjective"]]
            values.extend(next(r["normalizedToIndependent"] for r in native if r["label"] == label and r["level"] == level)
                          for level in (6, 7))
            bars = axes[0, column].bar(x+offset, values, .34, color=colors[label], label=labels[label])
            axes[0, column].bar_label(bars, fmt="%.3f", padding=3, fontsize=9)
        axes[0, column].axhline(1, color=colors["independent"], linestyle="--", linewidth=1)
        axes[0, column].set(xticks=x, xticklabels=["Porous MAC 36 x 24", "Native level 6", "Native level 7"],
                            ylabel="J / straight reference in the same model", title=f"Re = {reynolds}: lower is better")
        axes[0, column].set_ylim(0, max(2.6 if reynolds else 1.6, axes[0, column].get_ylim()[1]))
        axes[0, column].text(.03, .97, "Observed reversal" if ranking["reversalObservedBothLevels"] else "Same observed order",
                             transform=axes[0, column].transAxes, va="top", fontsize=10)
        for label in ("engineering", "candidate", "independent"):
            selected = sorted((r for r in native if r["label"] == label), key=lambda r:r["level"])
            axes[1, column].plot([r["cells"] for r in selected], [r["fluxWeightedPressureDrop"] for r in selected],
                "o-" if label != "independent" else "x--", color=colors[label], label=labels[label],
                linewidth=1.5, markersize=7)
            axes[2, column].plot([r["cells"] for r in selected], [r["nativeProcessSeconds"] for r in selected],
                "o-" if label != "independent" else "x--", color=colors[label], label=labels[label],
                linewidth=1.5, markersize=7)
        axes[1, column].set(xscale="log", xlabel="Actual fluid cells (both components summed)",
                           ylabel="Kinematic pressure drop [m²/s²]")
        axes[2, column].set(xscale="log", yscale="log", xlabel="Actual fluid cells (both components summed)",
                           ylabel="Native meshing + solve + export [wall s]")
        for ax in axes[:, column]:
            ax.grid(axis="y", alpha=.2)
        axes[1, column].ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    axes[0, 0].legend(loc="lower right", fontsize=8)
    axes[1, 0].legend(fontsize=8)
    axes[2, 0].legend(fontsize=8)
    figure.suptitle("T01 fixed-input recomputation — frozen current solver\n"
                    "Two-grid observations; physical accuracy and repaired-solver qualification pending", fontsize=14)
    figure.savefig(str(args.artifact)+".png", dpi=170)
    plt.close(figure)
    figure, axes = plt.subplots(2, 3, figsize=(12, 6), constrained_layout=True)
    for row, reynolds in enumerate((0, 50)):
        for column, label in enumerate(("candidate", "engineering", "independent")):
            directory = (args.output/f"independent-re{reynolds}" if label == "independent" else
                         args.comparison_root/f"t01-re{reynolds}-{args.comparison_suffix}"/
                         ("candidate" if label == "candidate" else "baseline"))
            extraction = read(directory/"extraction.json")
            ax = axes[row, column]
            ax.set_facecolor("#ebedf0")
            for region in extraction["regions"]:
                vertices = np.loadtxt(directory/region["file"], comments="#")
                ax.add_patch(Polygon(vertices, closed=True, facecolor="#bde7f0", edgecolor=colors[label], linewidth=1.4))
            ax.set(xlim=(-.04, 1.54), ylim=(-.03, 1.03), aspect="equal", xlabel="x", ylabel="y",
                   title=f"Re {reynolds} — {labels[label]}\n{extraction['vertices']} vertices; area {extraction['fluidArea']:.9f}")
    figure.suptitle("Actual sharp boundaries: independent four-corner rectangles reproduce the candidate geometry", fontsize=12)
    figure.savefig(str(args.artifact)+"-geometry.png", dpi=170)
    plt.close(figure)


def summarize(args):
    inputs = read(args.output/"inputs.json")
    preserved_sha = native_flow.sha(args.archive)
    if preserved_sha != inputs["archive"]["sha256"]:
        raise ValueError("preserved input archive differs from the input-preparation archive")
    preserved_path = args.archive.resolve()
    report = dict(schema="cartmesh2d-t01-double-pipe-recompute-v1", inputs=inputs, native=[], comparisons=[],
        physicalAccuracyQualified=False, repairedSolverQualified=False,
        preservedOriginalArchive=dict(path=str(preserved_path.relative_to(ROOT)) if preserved_path.is_relative_to(ROOT) else str(preserved_path),
                                      sha256=preserved_sha, bytes=args.archive.stat().st_size),
        scope="Fixed original T01 material fields and independent rectangle CFD; current frozen solver only.",
        metricDefinitions=dict(pressureDrop="Flux-weighted inlet-minus-outlet static kinematic pressure [m²/s²].",
            power="Pressure power per fluid density and unit out-of-plane depth [m4/s3]; not W.",
            dimensionlessJ="Total kinematic pressure power / (nu*S²/mu), S=.02, mu=1. Equal prescribed ports cancel kinetic power.",
            normalization="Native J divided by independently solved rectangles at that Re and level; porous J divided by the existing geometric seed at that Re.",
            cost="NativeProcessSeconds includes all recorded native meshing, solves and export for both components. All solves cold-start. Legacy Python orchestration overhead was not measured; the independent run additionally records complete orchestration elapsed time. Porous recomputation is a subset of preparation time and is not added again. Interpreter startup and report/plot rendering are outside the recorded orchestration intervals.",
            mesh="Same meshing controls and level, different bounding boxes and actual cell sizes. Report actual cells, domains and wall resolution; equal level is not equal cell density."),
        limitations=["Two grids provide observed sensitivity, not an uncertainty bound or physical-accuracy qualification.",
            "General incompressible CFD repair is pending; all native conclusions here refer to the explicitly frozen binary.",
            "Porous scores retain the original 36 x 24 MAC discretization; no MAC spatial convergence is established.",
            "Engineering porous fields have density volume 1/3; sharp contours are shifted to exact candidate area .5. This projection is part of the reported model-to-sharp gap.",
            "Recorded wall times were acquired while other research jobs used the shared machine; they are total measured costs, not an isolated speedup benchmark."],
        reportGeneratorSha256=native_flow.sha(__file__))
    for reynolds in (0, 50):
        for level, suffix in ((6, ""), (7, "-fine")):
            comparison = args.comparison_root/f"t01-re{reynolds}-{args.comparison_suffix}{suffix}"
            for label, directory in (("candidate", comparison/f"candidate-level-{level}"),
                    ("engineering", comparison/f"baseline-level-{level}"),
                    ("independent", args.output/f"independent-re{reynolds}-level-{level}")):
                report["native"].append(collect_native(directory, label, reynolds, level, inputs))
        rows = [r for r in report["native"] if r["reynolds"] == reynolds]
        source_case = next(e for e in inputs["cases"] if e["reynolds"] == reynolds)
        porous = {f["label"]:f for f in source_case["fields"]}
        paired = []
        for level in (6, 7):
            by_label = {r["label"]:r for r in rows if r["level"] == level}
            straight = by_label["independent"]["dimensionlessTotalPower"]
            for label, row in by_label.items():
                row["normalizedToIndependent"] = row["dimensionlessTotalPower"]/straight
                row["candidatePhysicalImprovementVsIndependent"] = 1-row["normalizedToIndependent"] if label == "candidate" else None
                row["porousGapRelativeToNative"] = ((porous[label]["porousObjective"]-row["dimensionlessTotalPower"])/row["dimensionlessTotalPower"]
                                                      if label in porous else None)
            candidate, engineering = by_label["candidate"], by_label["engineering"]
            paired.append(dict(level=level, nativeCandidateMinusEngineering=candidate["dimensionlessTotalPower"]-engineering["dimensionlessTotalPower"],
                nativeCandidateReductionVsEngineering=1-candidate["dimensionlessTotalPower"]/engineering["dimensionlessTotalPower"],
                candidateMinusIndependentRelative=candidate["normalizedToIndependent"]-1))
        changes = {}
        for label in ("candidate", "engineering", "independent"):
            selected = sorted((r for r in rows if r["label"] == label), key=lambda r:r["level"])
            change = selected[1]["dimensionlessTotalPower"]-selected[0]["dimensionlessTotalPower"]
            changes[label] = dict(signedJChange=change, relativeChange=change/selected[1]["dimensionlessTotalPower"],
                cells=[r["cells"] for r in selected], effectiveH=[r["effectiveH"] for r in selected],
                totalNativeProcessSeconds=sum(r["nativeProcessSeconds"] for r in selected))
        porous_difference = porous["candidate"]["porousObjective"]-porous["engineering"]["porousObjective"]
        fine_difference = paired[-1]["nativeCandidateMinusEngineering"]
        report["comparisons"].append(dict(reynolds=reynolds,
            porousCandidateMinusEngineering=porous_difference,
            porousCandidateReductionVsEngineering=1-porous["candidate"]["porousObjective"]/porous["engineering"]["porousObjective"],
            pairedLevels=paired, observedGridChanges=changes,
            reversalObservedBothLevels=all(porous_difference*p["nativeCandidateMinusEngineering"] < 0 for p in paired),
            nativeOrderSameBothLevels=paired[0]["nativeCandidateMinusEngineering"]*fine_difference > 0,
            absoluteFineSeparationMinusObservedGridChanges=abs(fine_difference)-abs(changes["candidate"]["signedJChange"])-abs(changes["engineering"]["signedJChange"]),
            physicalRankingResolved=False, candidateTopologyBenefitEstablished=False,
            topologyConclusion="Candidate and geometric-seed density fields differ only at roundoff, and both sharp boundaries match two independent rectangles. No topology gain over straight pipes is demonstrated."))
    report["totalRecordedNativeProcessSeconds"] = sum(row["nativeProcessSeconds"] for row in report["native"])
    report["newIndependentNativeProcessSeconds"] = sum(row["nativeProcessSeconds"] for row in report["native"] if row["label"] == "independent")
    report["newIndependentEndToEndSeconds"] = sum(row["measuredEndToEndSeconds"] for row in report["native"] if row["label"] == "independent")
    report["newPreparationSeconds"] = inputs["preparationSeconds"]
    report["newPreparationAndIndependentSeconds"] = report["newPreparationSeconds"]+report["newIndependentEndToEndSeconds"]
    report["newPorousRecomputeSeconds"] = sum(f["porousSolveSeconds"] for c in inputs["cases"] for f in c["fields"])
    args.artifact.parent.mkdir(parents=True, exist_ok=True)
    write(args.output/"summary.json", report)
    write(Path(str(args.artifact)+".json"), report)
    columns = ["reynolds", "label", "level", "cells", "area", "effectiveH", "inletFlux", "fluxWeightedPressureDrop",
               "pressurePower", "dimensionlessTotalPower", "normalizedToIndependent", "porousGapRelativeToNative",
               "meshSeconds", "flowSeconds", "nativeProcessSeconds", "measuredEndToEndSeconds", "nativeFlowConverged"]
    with Path(str(args.artifact)+".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(report["native"])
    porous_columns = ["reynolds", "label", "cells", "porousObjective", "normalizedToSeed",
                      "densityVolumeFraction", "rawSharpArea", "scaledPorousBoundaryPressureDrop",
                      "scaledPorousTotalPressurePower", "porousSolveSeconds", "porousLinearResidual",
                      "porousContinuity", "porousAdjointResidual"]
    with Path(str(args.artifact)+"-porous.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, porous_columns, extrasaction="ignore")
        writer.writeheader()
        for case in inputs["cases"]:
            seed = next(f["porousObjective"] for f in case["fields"] if f["label"] == "geometricSeed")
            for field in case["fields"]:
                writer.writerow(dict(field, reynolds=case["reynolds"],
                    cells=case["problem"]["nx"]*case["problem"]["ny"],
                    normalizedToSeed=field["porousObjective"]/seed,
                    rawSharpArea=field["geometry"]["rawSharpArea"]))
    make_plots(report, args)
    print(json.dumps(report["comparisons"], indent=2))


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "native", "report"))
    p.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    p.add_argument("--mesh-cli", type=Path, default=ROOT/"build/cartmesh2d_cli")
    p.add_argument("--flow-cli", type=Path, required=True)
    p.add_argument("--reynolds", type=int, choices=(0, 50), default=0)
    p.add_argument("--level", type=int, default=6)
    p.add_argument("--iterations", type=int, default=16000)
    p.add_argument("--timeout", type=float, default=1800)
    p.add_argument("--artifact", type=Path, default=ROOT/"artifacts/current/t01-frozen-double-pipe")
    p.add_argument("--comparison-root", type=Path, default=ROOT/"outputs/laminar-repair")
    p.add_argument("--comparison-suffix", default="engineering-before",
                   help="read t01-re{0,50}-SUFFIX and SUFFIX-fine from comparison-root")
    return p


if __name__ == "__main__":
    arguments = parser().parse_args()
    dict(prepare=prepare, native=run_native, report=summarize)[arguments.action](arguments)
