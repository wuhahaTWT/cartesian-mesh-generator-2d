#!/usr/bin/env python3
"""Render the actual OpenFOAM cell polygons and manufactured-solution error."""

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "verification"))
from check_openfoam2d import read_boundary, read_faces, read_labels, read_points  # noqa: E402
from openfoam_harmonic_mms import exact_value, scalar_internal  # noqa: E402


def cell_polygons(case: Path):
    mesh = case / "constant" / "polyMesh"
    points = read_points(mesh / "points")
    faces = read_faces(mesh / "faces")
    owner = read_labels(mesh / "owner")
    patches = read_boundary(mesh / "boundary")
    front = next(patch for patch in patches if patch["name"] == "frontAndBack")
    per_cell = {}
    for face_id in range(front["startFace"], front["startFace"] + front["nFaces"]):
        face = faces[face_id]
        mean_z = sum(points[index][2] for index in face) / len(face)
        cell = owner[face_id]
        candidate = (mean_z, [(points[index][0], points[index][1]) for index in face])
        if cell not in per_cell or candidate[0] < per_cell[cell][0]:
            per_cell[cell] = candidate
    return [per_cell[cell][1] for cell in range(len(per_cell))]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--field", choices=("constant", "linear", "quadratic"), required=True)
    parser.add_argument("--time", default="30")
    parser.add_argument("--reports", nargs="*", type=Path, default=[])
    args = parser.parse_args()
    polygons = cell_polygons(args.case)
    cx = scalar_internal(args.case / "0" / "Cx")
    cy = scalar_internal(args.case / "0" / "Cy")
    numerical = scalar_internal(args.case / args.time / "T")
    exact = [exact_value(x, y, args.field) for x, y in zip(cx, cy)]
    error = [abs(value - reference) for value, reference in zip(numerical, exact)]
    if not (len(polygons) == len(numerical) == len(exact)):
        raise ValueError("mesh and field cell counts differ")

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), constrained_layout=True)
    solution = PolyCollection(polygons, array=numerical, cmap="coolwarm",
                              edgecolors=(0, 0, 0, 0.22), linewidths=0.12)
    axes[0].add_collection(solution)
    axes[0].autoscale_view()
    axes[0].set_aspect("equal")
    axes[0].set_title(f"OpenFOAM solution T ({len(polygons)} cells)")
    fig.colorbar(solution, ax=axes[0], shrink=0.82)

    positive = [value for value in error if value > 0.0]
    floor = min(positive, default=1e-16)
    log_error = [math.log10(max(value, floor)) for value in error]
    errors = PolyCollection(polygons, array=log_error, cmap="magma",
                            edgecolors=(0, 0, 0, 0.24), linewidths=0.12)
    axes[1].add_collection(errors)
    axes[1].autoscale_view()
    axes[1].set_aspect("equal")
    axes[1].set_title("Cell error log10(|T - exact|)")
    fig.colorbar(errors, ax=axes[1], shrink=0.82)

    reports = [json.loads(path.read_text(encoding="utf-8")) for path in args.reports]
    if reports:
        cells = [report["cell_count"] for report in reports]
        axes[2].loglog(cells, [report["l1"] for report in reports], "o-", label="L1")
        axes[2].loglog(cells, [report["l2"] for report in reports], "s-", label="L2")
        axes[2].loglog(cells, [report["linf"] for report in reports], "^-", label="Linf")
        axes[2].set_xlabel("solver cells")
        axes[2].set_ylabel("error norm")
        axes[2].grid(True, which="both", alpha=0.3)
        axes[2].legend()
        axes[2].set_title("Global-refinement convergence")
    else:
        axes[2].axis("off")
    for axis in axes[:2]:
        axis.set_xlabel("x")
        axis.set_ylabel("y")
    fig.suptitle("cartmesh2d V1c: harmonic manufactured solution on final polyMesh")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=190)


if __name__ == "__main__":
    main()
