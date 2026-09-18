#!/usr/bin/env python3
"""Render native finite-volume fields and fixed-sequence convergence evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools/verification"))
from verify_native_fv import read_cm2d  # noqa: E402


def load_cells(path: Path) -> dict[int, dict[str, float]]:
    result: dict[int, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            cell_id = int(row["cell"])
            result[cell_id] = {name: float(row[name]) for name in ("value", "exact", "error", "source")}
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--problem", choices=("constant", "linear", "sine", "diffusion"), default="sine")
    parser.add_argument("--dpi", type=int, default=180)
    arguments = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import TwoSlopeNorm

    summary = json.loads(arguments.summary.read_text(encoding="utf-8"))
    if summary.get("format") != "cartmesh2d-native-fv-verification-v1":
        raise ValueError("unsupported native FV summary format")
    output = (arguments.output_dir or arguments.summary.parent / "figures").resolve()
    output.mkdir(parents=True, exist_ok=True)
    meshes = {item["name"]: item for item in summary["meshes"]}
    selected_cases = [case for case in summary["cases"] if case["problem"] == arguments.problem]
    selected_cases.sort(key=lambda case: float(meshes[case["grid"]]["convergence_h"]), reverse=True)
    if not selected_cases:
        raise ValueError(f"summary has no {arguments.problem} cases")

    loaded = []
    for case in selected_cases:
        mesh = read_cm2d(Path(meshes[case["grid"]]["path"]))
        values = load_cells(Path(case["artifacts"]["cells"]))
        if set(values) != set(range(len(mesh.cells))):
            raise ValueError(f"{case['grid']}: cell CSV ids do not match CM2D")
        polygons = [[mesh.vertices[vertex] for vertex in cell.vertices] for cell in mesh.cells]
        loaded.append((case, mesh, values, polygons))

    all_values = [row["value"] for _, _, values, _ in loaded for row in values.values()]
    field_min, field_max = min(all_values), max(all_values)
    manufactured = arguments.problem != "diffusion"
    all_errors = ([row["error"] for _, _, values, _ in loaded for row in values.values()]
                  if manufactured else [])
    error_limit = max(abs(min(all_errors)), abs(max(all_errors)), 1e-16) if all_errors else None
    columns = 2 if manufactured else 1
    fig, axes = plt.subplots(len(loaded), columns, figsize=((7 if columns == 1 else 12), 4 * len(loaded)),
                             squeeze=False, constrained_layout=True)
    for row_index, (case, mesh, values, polygons) in enumerate(loaded):
        field = [values[cell.id]["value"] for cell in mesh.cells]
        field_collection = PolyCollection(polygons, array=field, cmap="viridis", edgecolors="#263238",
                                          linewidths=0.08, clim=(field_min, field_max), rasterized=True)
        panels = [(axes[row_index, 0], field_collection, "numerical value")]
        if manufactured:
            error = [values[cell.id]["error"] for cell in mesh.cells]
            error_collection = PolyCollection(
                polygons, array=error, cmap="coolwarm", edgecolors="#263238", linewidths=0.08,
                norm=TwoSlopeNorm(vmin=-error_limit, vcenter=0.0, vmax=error_limit), rasterized=True)
            panels.append((axes[row_index, 1], error_collection, "value - exact"))
        for axis, collection, label in panels:
            axis.add_collection(collection)
            axis.autoscale_view()
            axis.set_aspect("equal")
            axis.set_xlabel("x")
            axis.set_ylabel("y")
            axis.set_title(f"{case['grid']} · {len(mesh.cells)} cells · {label}")
            fig.colorbar(collection, ax=axis, shrink=0.86)
    title = (f"Native 2-D FVM {arguments.problem} manufactured solution on final CM2D polygons"
             if manufactured else "Native 2-D FVM diffusion field\nfinal CM2D polygons; no analytic exact field")
    fig.suptitle(title, fontsize=14)
    fields_path = output / f"native-fv-{arguments.problem}-fields.png"
    fig.savefig(fields_path, dpi=arguments.dpi)
    plt.close(fig)

    convergence_path = None
    if manufactured and len(selected_cases) >= 3:
        hs = [float(meshes[case["grid"]]["convergence_h"]) for case in selected_cases]
        l2 = [float(case["independent"]["l2Error"]) for case in selected_cases]
        linf = [float(case["independent"]["linfError"]) for case in selected_cases]
        fig, axis = plt.subplots(figsize=(8.5, 6), constrained_layout=True)
        axis.loglog(hs, l2, "o-", linewidth=2, label="area-weighted L2")
        axis.loglog(hs, linf, "s--", linewidth=2, label="Linf")
        axis.invert_xaxis()
        axis.grid(True, which="both", alpha=0.3)
        axis.set_xlabel("requested wall size h")
        axis.set_ylabel("error")
        axis.set_title(f"Fixed geometry/domain {arguments.problem} convergence")
        if arguments.problem == "sine":
            for item in summary.get("convergence", {}).get("orders", []):
                coarse_index = next(i for i, case in enumerate(selected_cases) if case["grid"] == item["coarse"])
                fine_index = next(i for i, case in enumerate(selected_cases) if case["grid"] == item["fine"])
                x = math.sqrt(hs[coarse_index] * hs[fine_index])
                y = math.sqrt(l2[coarse_index] * l2[fine_index])
                axis.annotate(f"observed p={item['observed_order']:.3f}", (x, y), xytext=(8, 8),
                              textcoords="offset points", fontsize=10)
        axis.legend()
        axis.text(0.02, 0.02, "Observed orders apply only to this fixed sequence; no general order claim.",
                  transform=axis.transAxes, fontsize=9, color="#37474f")
        convergence_path = output / f"native-fv-{arguments.problem}-convergence.png"
        fig.savefig(convergence_path, dpi=arguments.dpi)
        plt.close(fig)

    result = {"problem": arguments.problem, "fields": str(fields_path),
              "convergence": str(convergence_path) if convergence_path else None,
              "source_summary": str(arguments.summary.resolve())}
    (output / f"native-fv-{arguments.problem}-render.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
