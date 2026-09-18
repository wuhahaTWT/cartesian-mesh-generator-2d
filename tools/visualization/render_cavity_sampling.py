#!/usr/bin/env python3
"""Visualize real cavity fields and a separately labelled sampling calibration."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "verification"))
from verify_native_flow import affine_sample, idw, read_cm2d


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit = json.loads(args.validation.read_text())
    comparison = json.loads(args.comparison.read_text())
    cases = sorted((c for c in audit["cases"] if c["scheme"] == "limited-linear"),
                   key=lambda c: c["verification"]["counts"]["cells"])
    if len(cases) != 3 or any(not c["verification"]["valid"] for c in cases):
        raise ValueError("Expected three individually verified cavity results")
    case = cases[-1]["verification"]
    mesh_path = Path(case["mesh"])
    csv_path = Path(case["prefix"] + ".cells.csv")
    if hashlib.sha256(mesh_path.read_bytes()).hexdigest() != case["meshSha256"]:
        raise ValueError("Mesh hash changed")
    if hashlib.sha256(csv_path.read_bytes()).hexdigest() != cases[-1]["hashes"]["cells.csv"]:
        raise ValueError("Field hash changed")
    mesh = read_cm2d(mesh_path)
    with csv_path.open() as stream:
        cells = sorted(csv.DictReader(stream), key=lambda row: int(row["cell"]))
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    polygons = [[mesh.vertices[i] for i in cell.vertices] for cell in mesh.cells]
    collection = PolyCollection(polygons, array=np.array([float(c["speed"]) for c in cells]),
                                cmap="viridis", edgecolor="#cbd5dc70", linewidth=.1)
    axes[0, 0].add_collection(collection)
    axes[0, 0].set(xlim=(0, 1), ylim=(0, 1), xlabel="x", ylabel="y",
                   title="Actual Re100 cavity | 3,844 cells, unchanged solver")
    axes[0, 0].set_aspect("equal")
    fig.colorbar(collection, ax=axes[0, 0], label="Speed / lid speed")

    n = 30
    f = lambda x, y: .7 + 1.3*x - .8*y
    synthetic = [dict(x=(i+.5)/n, y=(j+.5)/n, u=f((i+.5)/n, (j+.5)/n))
                 for j in range(n) for i in range(n)]
    walls = [dict(x=.5, y=y, u=f(.5, y)) for y in (0., 1.)]
    axis = np.linspace(.9, 1., 151)
    axes[0, 1].plot(axis, [idw(synthetic, .5, y, "u", boundary=walls)-f(.5, y)
                         for y in axis], color="#b7623d", label="Legacy IDW8")
    axes[0, 1].plot(axis, [affine_sample(synthetic, .5, y, "u", boundary=walls)-f(.5, y)
                         for y in axis], color="#207ca8", label="Affine fit")
    axes[0, 1].set(xlabel="y", ylabel="Sample minus exact value",
                   title="Calibration only: a known linear field, 30 x 30 samples")

    counts = [c["verification"]["counts"]["cells"] for c in cases]
    axes[1, 0].loglog(counts, [c["legacyRmse"] for c in cases], "o--",
                     color="#b7623d", label="Ghia / legacy IDW8")
    axes[1, 0].loglog(counts, [c["newRmse"] for c in cases], "o-",
                     color="#207ca8", label="Ghia / affine fit")
    axes[1, 0].set(xlabel="Final cells", ylabel="Centreline RMSE / lid speed",
                   title="Same saved fields | refinement gate still fails")

    samples = case["benchmark"]["samples"]
    old = case["benchmark"]["legacyIdw"]["samples"]
    u = sorted((s for s in samples if s["field"] == "u"), key=lambda s: s["coordinate"])
    legacy_u = sorted((s for s in old if s["field"] == "u"), key=lambda s: s["coordinate"])
    axes[1, 1].plot([s["actual"] for s in legacy_u], [s["coordinate"] for s in legacy_u],
                     ".--", color="#b7623d", label="Same FVM / IDW8")
    axes[1, 1].plot([s["actual"] for s in u], [s["coordinate"] for s in u],
                     ".-", color="#207ca8", label="Same FVM / affine")
    axes[1, 1].scatter([s["reference"] for s in u], [s["coordinate"] for s in u],
                       facecolors="none", edgecolors="#242424", s=35, label="Ghia table")
    axes[1, 1].set(xlabel="u / lid speed", ylabel="y / H", ylim=(0, 1),
                   title="Finest vertical centreline | original reference retained")
    for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.grid(alpha=.18, which="both"); ax.legend(frameon=False, fontsize=9)
    diffs = sorted((c for c in comparison["cases"] if c["scheme"] == "limited-linear"),
                   key=lambda c: c["cells"])
    values = " -> ".join(f"{c['errorsVsFD']['FD129']['combinedRms']:.3g}" for c in diffs)
    fig.suptitle("Cavity sampling audit | no changes to the flow solution", fontsize=16)
    fig.supxlabel("Separate 101-point FD129 comparison: " + values +
                  ". FD is a finite-grid diagnostic, not exact truth.", fontsize=10)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=165)


if __name__ == "__main__":
    main()
