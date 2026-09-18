#!/usr/bin/env python3
"""Render independently verified native 2-D flow fields and benchmark evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools/verification"))
from verify_native_flow import read_cm2d  # noqa: E402


@dataclass
class CaseData:
    record: dict[str, Any]
    label: str
    case: str
    mesh_path: Path
    prefix: Path
    mesh: Any
    cells: list[dict[str, float]]
    residuals: list[dict[str, float]]
    native: dict[str, Any]
    status: str
    polygons: list[list[tuple[float, float]]]


def finite(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{label} is not finite")
    return number


def command_option(command: Any, option: str) -> str | None:
    if not isinstance(command, list):
        return None
    try:
        index = command.index(option)
    except ValueError:
        return None
    return str(command[index + 1]) if index + 1 < len(command) else None


def path_from_summary(value: Any, summary_path: Path) -> Path | None:
    if not isinstance(value, str) or not value:
        return None
    path = Path(value)
    return (path if path.is_absolute() else summary_path.parent / path).resolve()


def case_paths(record: dict[str, Any], summary_path: Path) -> tuple[Path | None, Path | None]:
    stage = record.get("stage") if isinstance(record.get("stage"), dict) else {}
    mesh = record.get("mesh") or command_option(stage.get("command"), "--mesh")
    prefix = record.get("prefix") or command_option(stage.get("command"), "--output")
    return path_from_summary(mesh, summary_path), path_from_summary(prefix, summary_path)


def read_cells(path: Path, expected: int) -> list[dict[str, float]]:
    required = ("cell", "x", "y", "area", "u", "v", "p", "speed")
    indexed: dict[int, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        missing = [name for name in required if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        for row in reader:
            cell_id = int(row["cell"])
            if cell_id in indexed or not 0 <= cell_id < expected:
                raise ValueError(f"{path}: duplicate or out-of-range cell {cell_id}")
            indexed[cell_id] = {name: finite(row[name], f"cell {cell_id} {name}")
                                for name in required[1:]}
    if set(indexed) != set(range(expected)):
        raise ValueError(f"{path}: cell ids do not match CM2D")
    return [indexed[index] for index in range(expected)]


def read_residuals(path: Path) -> list[dict[str, float]]:
    columns = ("iteration", "momentumResidual", "continuity", "velocityChange", "pressureChange")
    result: list[dict[str, float]] = []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        missing = [name for name in columns if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        previous = 0
        for row in reader:
            iteration = int(row["iteration"])
            if iteration <= previous:
                raise ValueError(f"{path}: iterations are not strictly increasing")
            previous = iteration
            parsed = {name: finite(row[name], f"iteration {iteration} {name}")
                      for name in columns[1:]}
            if any(value < 0.0 for value in parsed.values()):
                raise ValueError(f"{path}: residual/change is negative")
            parsed["iteration"] = float(iteration)
            result.append(parsed)
    if not result:
        raise ValueError(f"{path}: no residual rows")
    return result


def native_status(native: dict[str, Any], record: dict[str, Any]) -> str:
    status, converged = native.get("status"), native.get("converged")
    if status == "converged" and converged is True:
        return "converged"
    if status == "iteration_limit" and converged is False:
        return "iteration limit — not converged"
    stage = record.get("stage", {})
    if isinstance(stage, dict) and stage.get("timedOut"):
        return "failed — timed out"
    return "failed — no valid convergence state"


def load_cases(summary: dict[str, Any], summary_path: Path) -> tuple[list[CaseData], list[str]]:
    loaded: list[CaseData] = []
    skipped: list[str] = []
    for record in summary.get("cases", []):
        label = str(record.get("label", "unnamed"))
        case = str(record.get("case", "unknown"))
        mesh_path, prefix = case_paths(record, summary_path)
        if mesh_path is None or prefix is None:
            skipped.append(f"{label}: summary has no mesh/output prefix")
            continue
        paths = {
            "mesh": mesh_path,
            "cells": Path(str(prefix) + ".cells.csv"),
            "residuals": Path(str(prefix) + ".residuals.csv"),
            "native": Path(str(prefix) + ".json"),
        }
        missing = [str(path) for path in paths.values() if not path.is_file()]
        if missing:
            skipped.append(f"{label}: missing {', '.join(missing)}")
            continue
        try:
            mesh = read_cm2d(paths["mesh"])
            cells = read_cells(paths["cells"], len(mesh.cells))
            residuals = read_residuals(paths["residuals"])
            native = json.loads(paths["native"].read_text(encoding="utf-8-sig"))
            if native.get("format") != "cartmesh2d-flow-summary-v1":
                raise ValueError("unsupported native flow JSON")
            polygons = [[mesh.vertices[vertex] for vertex in cell.vertices] for cell in mesh.cells]
            loaded.append(CaseData(record, label, case, mesh_path, prefix, mesh, cells,
                                   residuals, native, native_status(native, record), polygons))
        except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
            skipped.append(f"{label}: {exc}")
    return loaded, skipped


def select_finest(cases: list[CaseData], case_name: str) -> CaseData | None:
    selected = [case for case in cases if case.case == case_name]
    if not selected:
        return None
    return min(selected, key=lambda item: float(item.record.get("meshMeasurement", {}).get(
        "characteristicH", math.inf)))


def configure_axis(axis: Any, title: str) -> None:
    axis.set_facecolor("white")
    axis.set_title(title, fontsize=10.5, loc="left")
    axis.set_xlabel("x [m]")
    axis.set_ylabel("y [m]")
    axis.set_aspect("equal", adjustable="box")
    axis.tick_params(labelsize=8)


def draw_boundaries(axis: Any, case: CaseData) -> None:
    for patch, colour, width in ((2, "#52616b", 0.55), (1, "#111820", 1.15)):
        segments = []
        for edge in case.mesh.edges:
            if edge.neighbour < 0 and edge.patch == patch:
                segments.append((case.mesh.vertices[edge.v0], case.mesh.vertices[edge.v1]))
        if segments:
            from matplotlib.collections import LineCollection
            axis.add_collection(LineCollection(segments, colors=colour, linewidths=width,
                                               zorder=4, rasterized=True))


def render_fields(cases: list[CaseData], output: Path, dpi: int, source_note: str) -> Path:
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize, TwoSlopeNorm

    selected = [case for name in ("channel", "cavity", "external")
                if (case := select_finest(cases, name)) is not None]
    if not selected:
        raise ValueError("summary has no loadable flow field")
    figure, axes = plt.subplots(len(selected), 2, figsize=(12, 3.8 * len(selected)),
                                squeeze=False)
    for row, case in enumerate(selected):
        fields = (("speed", "|U| [m/s]", "viridis"),
                  ("p", "kinematic pressure p [m²/s²]", "coolwarm"))
        for column, (field, label, colourmap) in enumerate(fields):
            axis = axes[row, column]
            values = np.asarray([cell[field] for cell in case.cells], dtype=float)
            if field == "speed":
                norm = Normalize(vmin=0.0, vmax=max(float(values.max()), 1e-15))
            elif float(values.min()) < 0.0 < float(values.max()):
                norm = TwoSlopeNorm(vmin=float(values.min()), vcenter=0.0, vmax=float(values.max()))
            else:
                norm = Normalize(vmin=float(values.min()), vmax=max(float(values.max()),
                                                                     float(values.min()) + 1e-15))
            collection = PolyCollection(case.polygons, array=values, cmap=colourmap, norm=norm,
                                        edgecolors="none", rasterized=True)
            axis.add_collection(collection)
            axis.autoscale_view()
            draw_boundaries(axis, case)
            verification = "verification PASS" if case.record.get("valid") is True else "verification FAIL"
            configure_axis(axis, f"{case.label} · {len(case.mesh.cells):,} cells\n"
                                 f"{case.status}\n{verification} · {label}")
            figure.colorbar(collection, ax=axis, shrink=0.82, label=label)
    figure.suptitle("Native 2-D steady laminar flow on final CM2D polygons", fontsize=14)
    figure.text(0.5, 0.002, source_note + "\n"
                "External solid interiors remain white and are bounded in black. "
                "Iteration-limit and failed runs are not labelled converged.",
                ha="center", fontsize=8.5, color="#455a64")
    figure.tight_layout(rect=(0.0, 0.055, 1.0, 0.96))
    target = output / "native-flow-fields.png"
    figure.savefig(target, dpi=dpi, facecolor="white")
    plt.close(figure)
    return target


def render_residuals(cases: list[CaseData], output: Path, dpi: int, source_note: str) -> Path:
    import matplotlib.pyplot as plt

    if not cases:
        raise ValueError("summary has no residual history")
    figure, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    metrics = (("momentumResidual", "momentum residual"),
               ("velocityChange", "velocity change"),
               ("pressureChange", "pressure change"),
               ("continuity", "scaled local continuity"))
    for axis, (key, label) in zip(axes.flat, metrics):
        for case in cases:
            x = [row["iteration"] for row in case.residuals]
            y = [max(row[key], 1e-300) for row in case.residuals]
            axis.semilogy(x, y, linewidth=1.7, label=f"{case.label} · {case.status}")
        axis.set_title(label, loc="left", fontsize=11)
        axis.set_xlabel("SIMPLE iteration")
        axis.set_ylabel(label)
        axis.grid(True, which="both", color="#dce3e7", linewidth=0.6)
        axis.legend(fontsize=7.5)
    figure.suptitle("Native SIMPLE iteration histories from residual CSV files", fontsize=14)
    figure.text(0.5, 0.005, source_note + "\n"
                "These histories show the recorded stopping state; they do not by themselves certify accuracy.",
                ha="center", fontsize=8.5, color="#455a64")
    figure.tight_layout(rect=(0.0, 0.06, 1.0, 0.96))
    target = output / "native-flow-residuals.png"
    figure.savefig(target, dpi=dpi, facecolor="white")
    plt.close(figure)
    return target


def empty_benchmark(axis: Any, title: str, message: str) -> None:
    axis.set_title(title, loc="left", fontsize=11)
    axis.text(0.5, 0.5, message, transform=axis.transAxes, ha="center", va="center",
              color="#607d8b")
    axis.set_axis_off()


def channel_mid_profile(channel: CaseData) -> tuple[list[dict[str, float]], float, int]:
    bounds = channel.record.get("meshMeasurement", {}).get("bounds")
    xmin, ymin, xmax, ymax = map(float, bounds)
    mid = 0.5 * (xmin + xmax)
    nearest_distance = min(abs(cell["x"] - mid) for cell in channel.cells)
    x_tolerance = max(1e-12, 1e-10 * max(1.0, xmax - xmin))
    near = [cell for cell in channel.cells
            if abs(cell["x"] - mid) <= nearest_distance + x_tolerance]
    assert near, "nearest channel mid-plane band must contain cells"

    y_tolerance = max(1e-12, 1e-10 * max(1.0, ymax - ymin))
    groups: list[list[dict[str, float]]] = []
    for cell in sorted(near, key=lambda item: item["y"]):
        if not groups or abs(cell["y"] - groups[-1][0]["y"]) > y_tolerance:
            groups.append([cell])
        else:
            groups[-1].append(cell)
    profile = [{"y": sum(cell["y"] for cell in group) / len(group),
                "u": sum(cell["u"] for cell in group) / len(group)}
               for group in groups]
    assert profile, "channel mid-plane profile must contain averaged samples"
    return profile, nearest_distance, len(near)


def render_benchmarks(cases: list[CaseData], output: Path, dpi: int, source_note: str) -> Path:
    import matplotlib.pyplot as plt
    import numpy as np

    channel = select_finest([case for case in cases if case.record.get("valid") is True], "channel")
    cavity = select_finest([case for case in cases if case.record.get("valid") is True], "cavity")
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.4))

    if channel is None:
        empty_benchmark(axes[0], "Plane Poiseuille channel", "No verified channel case")
    else:
        bounds = channel.record.get("meshMeasurement", {}).get("bounds")
        xmin, ymin, xmax, ymax = map(float, bounds)
        profile, nearest_distance, band_cells = channel_mid_profile(channel)
        eta = np.linspace(0.0, 1.0, 300)
        axes[0].plot(4.0 * eta * (1.0 - eta), eta, color="#263238", linewidth=2,
                     label=r"analytic $4\eta(1-\eta)$")
        axes[0].plot([sample["u"] / channel.native["speed"] for sample in profile],
                     [(sample["y"] - ymin) / (ymax - ymin) for sample in profile], "o-",
                     color="#00796b", markersize=3, linewidth=1.2,
                     label=f"native nearest-midplane band; y-averaged ({len(profile)} points)")
        benchmark = channel.record.get("benchmark", {})
        axes[0].set_title(f"Channel · {channel.label}\nvelocity L2 rel. "
                          f"{float(benchmark.get('velocityL2Relative', math.nan)):.3g} · "
                          f"|x-xmid|={nearest_distance:.3g} m · {band_cells} cells",
                          loc="left", fontsize=11)
        axes[0].set_xlabel("u / inlet maximum speed")
        axes[0].set_ylabel("y / H")
        axes[0].grid(True, color="#dce3e7", linewidth=0.6)
        axes[0].legend(fontsize=8)

    samples = cavity.record.get("benchmark", {}).get("samples", []) if cavity else []
    for axis, field, title, x_label, y_label in (
        (axes[1], "u", "Cavity vertical centreline", "u / lid speed", "y / H"),
        (axes[2], "v", "Cavity horizontal centreline", "x / W", "v / lid speed"),
    ):
        chosen = [sample for sample in samples if sample.get("field") == field]
        if not chosen:
            empty_benchmark(axis, title, "No verified Re=100 Ghia comparison")
            continue
        coordinate = [float(sample["coordinate"]) for sample in chosen]
        reference = [float(sample["reference"]) for sample in chosen]
        actual = [float(sample["actual"]) for sample in chosen]
        if field == "u":
            axis.plot(reference, coordinate, "s", color="#37474f", markersize=4, label="Ghia et al. (1982)")
            axis.plot(actual, coordinate, "o-", color="#c62828", markersize=3, linewidth=1.2,
                      label="native anchored interpolation")
        else:
            axis.plot(coordinate, reference, "s", color="#37474f", markersize=4, label="Ghia et al. (1982)")
            axis.plot(coordinate, actual, "o-", color="#1565c0", markersize=3, linewidth=1.2,
                      label="native anchored interpolation")
        benchmark = cavity.record.get("benchmark", {})
        axis.set_title(f"{title} · {cavity.label}\nRMSE "
                       f"{float(benchmark.get('centrelineRmse', math.nan)):.3g}", loc="left", fontsize=11)
        axis.set_xlabel(x_label)
        axis.set_ylabel(y_label)
        axis.grid(True, color="#dce3e7", linewidth=0.6)
        axis.legend(fontsize=8)
    figure.suptitle("Case-specific comparisons on verified native flow outputs", fontsize=14)
    figure.text(0.5, 0.005, source_note, ha="center", fontsize=8.5, color="#455a64")
    figure.tight_layout(rect=(0.0, 0.04, 1.0, 0.96))
    target = output / "native-flow-benchmarks.png"
    figure.savefig(target, dpi=dpi, facecolor="white")
    plt.close(figure)
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True,
                        help="cartmesh2d-native-flow-verification-v1 summary")
    parser.add_argument("--output", type=Path, required=True, help="directory for three PNG figures")
    parser.add_argument("--dpi", type=int, default=180)
    arguments = parser.parse_args()
    if arguments.dpi < 72 or arguments.dpi > 600:
        raise ValueError("--dpi must be between 72 and 600")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "font.size": 9.5, "axes.edgecolor": "#455a64", "axes.labelcolor": "#263238",
        "xtick.color": "#455a64", "ytick.color": "#455a64",
    })

    summary_path = arguments.summary.resolve()
    summary = json.loads(summary_path.read_text(encoding="utf-8-sig"))
    if summary.get("format") != "cartmesh2d-native-flow-verification-v1":
        raise ValueError("unsupported native flow verification summary")
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cases, skipped = load_cases(summary, summary_path)
    if not cases:
        raise ValueError("summary has no loadable real CM2D/CSV flow case")

    summary_state = "PASS" if summary.get("valid") is True else "FAIL"
    source_note = f"Source verification summary: {summary_state}"
    if skipped:
        source_note += f" · {len(skipped)} unavailable case{'s' if len(skipped) != 1 else ''} skipped"
    figures = {
        "fields": str(render_fields(cases, output, arguments.dpi, source_note)),
        "residuals": str(render_residuals(cases, output, arguments.dpi, source_note)),
        "benchmarks": str(render_benchmarks(cases, output, arguments.dpi, source_note)),
    }
    report = {
        "format": "cartmesh2d-native-flow-render-v1",
        "sourceSummary": str(summary_path),
        "sourceSummaryValid": summary.get("valid") is True,
        "loadedCases": [{"label": case.label, "case": case.case, "status": case.status,
                         "verificationValid": case.record.get("valid") is True,
                         "mesh": str(case.mesh_path), "prefix": str(case.prefix)} for case in cases],
        "skipped": skipped,
        "figures": figures,
    }
    report_path = output / "render-summary.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({**figures, "report": str(report_path)}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
