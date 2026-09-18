#!/usr/bin/env python3
"""Show a verified channel's real fields, refinement interface and profile."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

from render_native_flow import load_cases, channel_mid_profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--zoom", type=float, nargs=4, required=True,
                        metavar=("XMIN", "XMAX", "YMIN", "YMAX"))
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text())
    summary["cases"] = [c for c in summary["cases"] if c["label"] == args.label]
    if len(summary["cases"]) != 1 or not summary["cases"][0].get("valid"):
        raise ValueError("requires one independently verified case")
    cases, skipped = load_cases(summary, args.summary.resolve())
    if skipped or len(cases) != 1 or cases[0].case != "channel":
        raise ValueError("requires loadable real channel mesh and fields")
    case = cases[0]
    left, right, bottom, top = args.zoom
    if not np.isfinite(args.zoom).all() or not left < right or not bottom < top:
        raise ValueError("invalid zoom bounds")

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#596774"})
    fig, axes = plt.subplots(2, 2, figsize=(12, 7),
                             gridspec_kw={"height_ratios": [1, 1.6]})
    for axis, field, title, cmap in (
        (axes[0, 0], "speed", "Velocity magnitude [m/s]", "viridis"),
        (axes[0, 1], "p", "Kinematic pressure [m²/s²]", "coolwarm"),
    ):
        collection = PolyCollection(case.polygons, array=np.array([c[field] for c in case.cells]),
                                    cmap=cmap, edgecolors="none", rasterized=True)
        axis.add_collection(collection)
        axis.autoscale_view()
        axis.set_aspect("equal")
        axis.set(xlabel="x [m]", ylabel="y [m]", title=title)
        fig.colorbar(collection, ax=axis, shrink=.8, pad=.02)

    zoom = [p for p in case.polygons
            if max(v[0] for v in p) >= left and min(v[0] for v in p) <= right
            and max(v[1] for v in p) >= bottom and min(v[1] for v in p) <= top]
    if not zoom:
        raise ValueError("zoom does not intersect the actual mesh")
    axis = axes[1, 0]
    axis.add_collection(PolyCollection(zoom, facecolors="white", edgecolors="#71818c", linewidths=.45))
    axis.set(xlim=(left, right), ylim=(bottom, top), xlabel="x [m]", ylabel="y [m]",
             title="Actual final mesh: local refinement interface")
    axis.set_aspect("equal")

    profile, _, _ = channel_mid_profile(case)
    bounds = case.record["meshMeasurement"]["bounds"]
    ymin, ymax = bounds[1], bounds[3]
    eta = np.linspace(0, 1, 400)
    axis = axes[1, 1]
    axis.plot(4 * eta * (1 - eta), eta, color="#293741", lw=2, label="Analytic Poiseuille")
    axis.plot([p["u"] / case.native["speed"] for p in profile],
              [(p["y"] - ymin) / (ymax - ymin) for p in profile], ".", ms=2,
              color="#00796b", label="Computed midplane band")
    error = case.record["benchmark"]["velocityL2Relative"]
    axis.set(xlabel="u / inlet maximum speed", ylabel="y / H",
             title=f"Velocity profile · relative L2 error {error:.3g}")
    axis.grid(color="#dce3e7", lw=.5)
    axis.legend(fontsize=9)
    fig.suptitle(f"Native 2D laminar channel · {len(case.cells):,} cells · {case.status}", fontsize=15)
    fig.text(.5, .018, "Real CM2D polygons and CSV fields · independent geometry, continuity and analytic checks passed",
             ha="center", fontsize=9, color="#455a64")
    fig.tight_layout(rect=(0, .04, 1, .95))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160, facecolor="white")
    plt.close(fig)
    print(json.dumps({"image": str(args.output), "cells": len(case.cells),
                      "visibleZoomPolygons": len(zoom), "label": case.label}))


if __name__ == "__main__":
    main()
