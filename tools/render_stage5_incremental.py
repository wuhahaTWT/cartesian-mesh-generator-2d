#!/usr/bin/env python3
"""渲染阶段 5 修改前、修改后和实际重构区域三联图。"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meshio
import numpy as np


def flattened_cells(mesh: meshio.Mesh) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    centers = []
    for block in mesh.cells:
        centers.append(np.mean(mesh.points[np.asarray(block.data)], axis=1))
    data = {
        name: np.concatenate([np.asarray(array, dtype=float) for array in arrays])
        for name, arrays in mesh.cell_data.items()
    }
    return np.concatenate(centers), data


def slice_mask(centers: np.ndarray, target: float) -> np.ndarray:
    distance = np.abs(centers[:, 2] - target)
    threshold = np.quantile(distance, 0.08)
    return distance <= threshold + 1.0e-15


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old-mesh", type=Path, required=True)
    parser.add_argument("--new-mesh", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--slice-z", type=float, default=0.5)
    args = parser.parse_args()

    old_centers, old_data = flattened_cells(meshio.read(args.old_mesh))
    new_centers, new_data = flattened_cells(meshio.read(args.new_mesh))
    old_mask = slice_mask(old_centers, args.slice_z)
    new_mask = slice_mask(new_centers, args.slice_z)

    figure, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    panels = (
        (axes[0], old_centers[old_mask], old_data["fluid_volume_fraction"][old_mask],
         "Before: fluid volume fraction", "viridis", 0.0, 1.0),
        (axes[1], new_centers[new_mask], new_data["fluid_volume_fraction"][new_mask],
         "After: fluid volume fraction", "viridis", 0.0, 1.0),
        (axes[2], new_centers[new_mask], new_data["stage5_rebuilt"][new_mask],
         "Incremental work: rebuilt=1", "coolwarm", 0.0, 1.0),
    )
    for axis, centers, values, title, color_map, minimum, maximum in panels:
        scatter = axis.scatter(centers[:, 0], centers[:, 1], c=values, s=16,
                               cmap=color_map, vmin=minimum, vmax=maximum,
                               marker="s", linewidths=0)
        axis.set_title(title)
        axis.set_aspect("equal")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        figure.colorbar(scatter, ax=axis, shrink=0.78)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=170)
    plt.close(figure)


if __name__ == "__main__":
    main()
