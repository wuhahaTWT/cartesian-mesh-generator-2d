#!/usr/bin/env python3
"""Traceable connectivity/initialisation study from actual saved density fields.

The figure shows contours, not a CFD field. Native comparisons and all supplied
failed mesh probes are retained separately in its JSON evidence sidecar.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from topology_artifacts import contours, port_connectivity, signed_area
from compare_sharp_designs import assess, flow_metrics


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path):
    return dict(path=str(path.resolve()), sha256=sha(path))


def native_record(directory):
    path = directory/"summary.json"
    report = json.loads(path.read_text())
    if "components" in report:
        return dict(source=record(path), components=[native_record(directory/f"component-{index}")
                    for index in range(len(report["components"]))])
    result = {key:report[key] for key in ("meshAccepted", "independentTopologyPassed", "solverQualityPassed",
               "nativeFlowConverged", "independentFlowPassed", "externalCheckMesh", "executables")}
    result.update(source=record(path), cases=[])
    if "steadyContinuation" in report:
        result["steadyContinuation"] = report["steadyContinuation"]
    for case in report["cases"]:
        entry = {key:case[key] for key in ("level", "status", "mesh", "metrics", "independentFlowIssues") if key in case}
        entry["runs"] = case["runs"]
        if "flow" in case:
            entry["flow"] = {key:case["flow"].get(key) for key in ("iterations", "converged", "tolerance")}
            entry["flowSources"] = [record(directory/f"level-{case['level']}"/name) for name in
                                    ("flow.json", "flow.cells.csv", "flow.faces.csv", "flow.boundaries", "independent-flow.json")
                                    if (directory/f"level-{case['level']}"/name).exists()]
        failed = directory/f"level-{case['level']}"/"mesh.failed.solver-quality.json"
        if failed.exists():
            quality = json.loads(failed.read_text())
            entry["failedQuality"] = dict(source=record(failed), valid=quality["valid"],
                    policy=quality["policy"], metrics=quality["metrics"], issueCount=quality["issue_count"])
        accepted_quality = directory/f"level-{case['level']}"/"openfoam/solver_quality.json"
        if accepted_quality.exists():
            quality = json.loads(accepted_quality.read_text())
            entry["acceptedQuality"] = dict(source=record(accepted_quality), valid=quality["valid"],
                    policy=quality["policy"], metrics=quality["metrics"], issueCount=quality["issue_count"])
        diagnostic = directory/f"level-{case['level']}"/"independent-flow-diagnostic.json"
        if diagnostic.exists():
            audit = json.loads(diagnostic.read_text())
            entry["diagnosticOnlyAudit"] = dict(source=record(diagnostic), valid=audit["valid"], issues=audit["issues"],
                    momentumResidual=audit.get("momentumAudit", {}).get("cellResidual", {}).get("maxNormalized"),
                    scope="Extra diagnostic of a rejected run; cannot replace native convergence or qualify a design comparison.")
        result["cases"].append(entry)
    return result


def grouped_sensitivity(comparisons):
    """Never combine different geometry or native control settings into a trend."""
    groups = {}
    for report in comparisons:
        def executables(native):
            if "components" in native:
                return [entry for child in native["components"] for entry in executables(child)]
            return [(Path(path).name, value) for path, value in native.get("executables", {}).items()]
        hashes = sorted(set(entry for native in report.get("nativeEvidence", {}).values()
                            for entry in executables(native)))
        identity = dict(controls={key:value for key,value in report["controls"].items() if key != "levels"},
                        executables=hashes,
                        boundaries={key:report["designs"][key]["extraction"]["boundarySha256"]
                                    for key in ("baseline", "candidate")},
                        ports={key:report["problem"][key] for key in ("width", "height", "port_width", "case")})
        key = json.dumps(identity, sort_keys=True)
        group = groups.setdefault(key, dict(identity=identity, rows={}, sources=[], analyticReferences=[]))
        group["sources"].append(report["source"])
        group["analyticReferences"].append(report.get("analyticBaseline"))
        for row in report["rows"]:
            if row["level"] in group["rows"]:
                raise ValueError("duplicate grid in one sensitivity group; do not silently choose a preferred run")
            group["rows"][row["level"]] = row
    result = []
    for group in groups.values():
        group["rows"] = [group["rows"][level] for level in sorted(group["rows"])]
        group["allAttemptedGrids"] = assess(group["rows"])
        valid = [row for row in group["rows"] if row.get("baseline") and row.get("candidate")]
        group["acceptedPairsOnly"] = dict(levels=[row["level"] for row in valid], assessment=assess(valid),
                scope="Sensitivity among fully audited pairs only; rejected grids remain in allAttemptedGrids and invalidate any all-grid success claim.")
        references = group.pop("analyticReferences")
        if references and all(reference is not None for reference in references):
            drops = {reference["kinematicPressureDrop"] for reference in references}
            if len(drops) != 1:
                raise ValueError("analytic baselines differ within one sensitivity group")
            identity = group["identity"]
            throughput = 4*identity["controls"]["speed"]*identity["ports"]["port_width"]/3
            group["analyticReferenceComparison"] = analytic_reference_comparison(group["rows"], drops.pop(), throughput)
        result.append(group)
    return result


def analytic_reference_comparison(rows, reference_drop, throughput):
    """Separate check against an exact rectangular-channel reference.

    Does not change native-pair acceptance, rescue rejected flow fields, or
    establish a discretisation error bound for the optimized candidate.
    """
    values = []
    for row in rows:
        candidate = row.get("candidate")
        if candidate is None:
            continue
        if abs(candidate["inletFlux"]-throughput) > 1e-10*throughput:
            raise ValueError("candidate throughput differs from analytic straight-pipe reference")
        drop = candidate["fluxWeightedPressureDrop"]
        if not np.isfinite(drop) or drop <= 0:
            raise ValueError("candidate requires a positive finite pressure drop")
        values.append(dict(level=row["level"], candidateDrop=drop, relativeReduction=1-drop/reference_drop))
    result = dict(referencePressureDrop=reference_drop, throughput=throughput, candidates=values,
                  rankingStableOnAvailableCandidateGrids=False, physicalAccuracyQualified=False,
                  scope="Analytic planar Poiseuille reference versus converged/audited native candidates; not a native paired-grid comparison. All native baseline failures remain rejected.")
    if len(values) >= 2:
        previous, final = values[-2:]
        change = abs(final["candidateDrop"]-previous["candidateDrop"])
        margin = reference_drop-final["candidateDrop"]-change
        result.update(observedCandidateGridChange=change, gainBeyondObservedChange=margin,
                      rankingStableOnAvailableCandidateGrids=bool(margin > 0 and previous["candidateDrop"] < reference_drop),
                      interpretation="Observed two-grid sensitivity only, not a rigorous error estimate or mesh-independence proof.")
    return result


def straight_pipe_baseline(directory, report):
    """Analytic pressure reference only when the actual baseline is two rectangles."""
    p = report["problem"]
    if p["case"] != "double-pipe":
        return None
    with np.load(directory/"baseline/final.npz") as saved:
        groups = contours(saved["rho"], p["width"], p["height"])
    if len(groups) != 2 or any(len(group) != 1 for group in groups):
        return None
    tolerance = 1e-11+1e-9*max(p["width"], p["height"])
    for group, centre in zip(sorted(groups, key=lambda g:np.mean(g[0][:, 1])), [p["height"]/4, 3*p["height"]/4]):
        loop = group[0]
        lower = np.array([0, centre-p["port_width"]/2])
        upper = np.array([p["width"], centre+p["port_width"]/2])
        if (np.max(np.abs(loop.min(axis=0)-lower)) > tolerance or
                np.max(np.abs(loop.max(axis=0)-upper)) > tolerance):
            return None
        for a, b in zip(loop, np.roll(loop, -1, axis=0)):
            if not any(abs(a[axis]-value) <= tolerance and abs(b[axis]-value) <= tolerance
                       for axis in (0, 1) for value in (lower[axis], upper[axis])):
                return None
    controls = report["controls"]
    drop = 8*controls["nu"]*controls["speed"]*p["width"]/p["port_width"]**2
    return dict(kinematicPressureDrop=drop, formula="8 * nu * peak_inlet_speed * length / channel_width**2",
                scope="Fully developed planar Poiseuille solution; actual baseline rectangles verified. No candidate accuracy claim or new acceptance threshold.",
                observedErrors=[dict(level=row["level"], relativeError=abs(row["baseline"]["fluxWeightedPressureDrop"]/drop-1))
                                for row in report["rows"] if row.get("baseline")])


def render(runs, comparisons, probes, output, continuations=()):
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2]/"outputs/topology-plot-cache"))
    os.environ.setdefault("XDG_CACHE_HOME", os.environ["MPLCONFIGDIR"])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as PlotPath
    from matplotlib.patches import PathPatch

    evidence = dict(schema="cartmesh2d-connectivity-study-v1", runs=[], nativeComparisons=[], meshProbes=[],
                    physicalAccuracyQualified=False,
                    meaning="Connected-component port reachability, not flow transfer or mixing.",
                    limits=["Porous objectives across different alpha_max values are not the same objective.",
                            "Raw isocontour area differs from projected sharp CFD area; native comparisons match the latter.",
                            "Iteration limits are not optimization convergence; two grids are not a rigorous error estimate."])
    panels = []
    baseline_hash = None
    for directory in comparisons:
        path = directory/"summary.json"
        data = json.loads(path.read_text())
        if data.get("status") == "running":
            raise ValueError("study evidence requires completed comparisons; a running solve is not a result")
        baseline = data["designs"]["baseline"]
        if baseline_hash is None:
            baseline_hash = baseline["extraction"]["boundarySha256"]
            panels.append((directory/"baseline/final.npz", data["problem"], "Reference: two straight pipes", "Equal sharp-area baseline"))
        elif baseline["extraction"]["boundarySha256"] != baseline_hash:
            raise ValueError("study comparisons must share the exact same sharp baseline")
        item = {key:data[key] for key in ("problem", "controls", "areaTarget", "targetArea", "status", "assessment", "rows", "issues", "designs", "sourceHashes")}
        item["source"] = record(path)
        item["analyticBaseline"] = straight_pipe_baseline(directory, data)
        item["nativeEvidence"] = {str(Path(row[label+"Output"])):native_record(Path(row[label+"Output"]))
            for row in data["rows"] for label in ("baseline", "candidate") if label+"Output" in row}
        evidence["nativeComparisons"].append(item)
    shared_problem, shared_controls = None, None
    for directory in runs:
        path, field = directory/"summary.json", directory/"final.npz"
        data = json.loads(path.read_text())
        problem = {k:v for k,v in data["problem"].items() if k != "alpha_max"}
        controls = {k:v for k,v in data["controls"].items() if k != "initialization"}
        if shared_problem is not None and (problem != shared_problem or controls != shared_controls):
            raise ValueError("sensitivity runs differ beyond initialisation and solid resistance")
        shared_problem, shared_controls = problem, controls
        item = {key:data[key] for key in ("problem", "controls", "initialization", "status", "optimizationConverged",
                    "final", "finalParameters", "sameModelReferences", "sourceHashes", "elapsedSeconds")}
        item.update(source=record(path), field=record(field))
        with np.load(field) as saved:
            groups = contours(saved["rho"], float(saved["width"]), float(saved["height"]))
        item["rawContour"] = dict(components=len(groups), holes=sum(len(g)-1 for g in groups),
                                 fluidArea=sum(abs(signed_area(g[0]))-sum(abs(signed_area(h)) for h in g[1:]) for g in groups),
                                 portConnectivity=port_connectivity(groups, data["problem"]))
        evidence["runs"].append(item)
        name = data["initialization"]["kind"]
        title = f"{name.capitalize()} initial field | solid resistance {data['problem']['alpha_max']:g}"
        origin = "Merged connection supplied in initial seed" if name == "merged" else "No interior connection supplied"
        panels.append((field, data["problem"], title, origin))
    for directory in probes:
        evidence["meshProbes"].append(native_record(directory))
    evidence["nativeSensitivity"] = grouped_sensitivity(evidence["nativeComparisons"])
    evidence["continuedCandidates"] = []
    for directory in continuations:
        report = json.loads((directory/"summary.json").read_text())
        if "steadyContinuation" not in report or len(report["cases"]) != 1:
            raise ValueError("continued candidate must record a single same-mesh continuation")
        item = dict(nativeEvidence=native_record(directory), accepted=False)
        metric = flow_metrics(report)
        if metric is not None:
            continuation = report["steadyContinuation"]
            original = continuation.get("originalNativeDirectory", continuation["source"])
            source_comparison = record(Path(original).parent/"summary.json")
            groups = [group for group in evidence["nativeSensitivity"] if source_comparison in group["sources"]]
            if len(groups) != 1:
                raise ValueError("include the original failed comparison for this continued candidate")
            group = groups[0]
            if report["sourceSha256"] != group["identity"]["boundaries"]["candidate"]:
                raise ValueError("continuation is not the original candidate boundary")
            analytic = group.get("analyticReferenceComparison")
            if analytic is None:
                raise ValueError("continued-candidate comparison requires a verified analytic straight-pipe reference")
            # A separate evidence item; the original native-pair rows stay failed.
            level = report["cases"][0]["level"]
            rows = [dict(row) for row in group["rows"]]
            old = next(row for row in rows if row["level"] == level)
            if old.get("candidate") is not None:
                raise ValueError("continuation cannot silently replace an already accepted candidate")
            old["candidate"] = metric
            item.update(accepted=True, metrics=metric, sourceComparison=source_comparison,
                        analyticReferenceComparison=analytic_reference_comparison(rows,
                            analytic["referencePressureDrop"], analytic["throughput"]),
                        scope="Additional explicitly recorded iteration budget and same-mesh initial guess. Same equations/tolerance; not a cold-run cost comparison. Original native-pair failures remain unchanged.")
        evidence["continuedCandidates"].append(item)
    if not panels:
        raise ValueError("no saved fields to render")
    fig, axes = plt.subplots((len(panels)+1)//2, 2, figsize=(12, 3.45*((len(panels)+1)//2)),
                             constrained_layout=True, squeeze=False)
    for axis, (field, problem, title, origin) in zip(axes.flat, panels):
        with np.load(field) as saved:
            width, height = float(saved["width"]), float(saved["height"])
            groups = contours(saved["rho"], width, height)
        total = 0.
        for index, group in enumerate(groups):
            # Orient holes opposite outer boundaries for matplotlib's nonzero rule.
            vertices, codes = [], []
            for loop_index, loop in enumerate(group):
                want_positive = loop_index == 0
                if (signed_area(loop) > 0) != want_positive:
                    loop = loop[::-1]
                total += abs(signed_area(loop))*(1 if want_positive else -1)
                vertices.extend(np.vstack([loop, loop[0]]))
                codes.extend([PlotPath.MOVETO]+[PlotPath.LINETO]*(len(loop)-1)+[PlotPath.CLOSEPOLY])
            axis.add_patch(PathPatch(PlotPath(vertices, codes), facecolor=["#91d0e3", "#f3bd7a"][index%2],
                                     edgecolor="#183949", linewidth=1.05))
        graph = port_connectivity(groups, problem)
        axis.set(xlim=(-.025, width+.025), ylim=(-.025, height+.025), aspect="equal",
                 xlabel="x", ylabel="y", title=title)
        axis.set_facecolor("#e6e9ed")
        axis.text(.5, -.26, f"{origin}\n{len(groups)} fluid region(s); contour area {total:.6f}; all ports covered: {graph['allPortsCovered']}",
                  transform=axis.transAxes, ha="center", va="top", fontsize=9)
    for axis in list(axes.flat)[len(panels):]:
        axis.set_visible(False)
    fig.suptitle("Flow topology: separate pipes or a shared trunk?\nActual un-smoothed density contours; colours denote connected fluid regions", fontsize=14)
    fig.savefig(output, dpi=170)
    plt.close(fig)
    evidence["figure"] = record(output)
    output.with_suffix(".json").write_text(json.dumps(evidence, indent=2, allow_nan=False)+"\n")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", type=Path, required=True)
    parser.add_argument("--comparisons", nargs="+", type=Path, required=True)
    parser.add_argument("--mesh-probes", nargs="*", type=Path, default=[])
    parser.add_argument("--continued-candidates", nargs="*", type=Path, default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(render(args.runs, args.comparisons, args.mesh_probes, args.output, args.continued_candidates))
