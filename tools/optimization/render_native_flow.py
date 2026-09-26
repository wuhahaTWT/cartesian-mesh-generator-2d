#!/usr/bin/env python3
"""Plot actual accepted Cut-cell polygons and native CSV fields, without interpolation."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"tools"/"verification"))
from verify_native_flow import read_cm2d


def leaf_reports(report):
    if "components" in report:
        return [leaf for child in report["components"] for leaf in leaf_reports(child)]
    return [report]


def load_case(directory):
    report = json.loads((directory/"summary.json").read_text())
    polygons, speed, pressure, boundaries, sources, iterations = [], [], [], [], [], []
    for leaf in leaf_reports(report):
        if not leaf["independentTopologyPassed"] or not leaf["independentFlowPassed"]:
            raise ValueError("preview requires an accepted mesh and independently audited flow")
        path = Path(leaf["acceptedMesh"]["path"])
        if hashlib.sha256(path.read_bytes()).hexdigest() != leaf["acceptedMesh"]["sha256"]:
            raise ValueError("mesh changed after verification")
        mesh = read_cm2d(path)
        flow = path.parent/"flow.cells.csv"
        with flow.open() as stream:
            rows = {int(row["cell"]):row for row in csv.DictReader(stream)}
        if set(rows) != {cell.id for cell in mesh.cells}:
            raise ValueError("field IDs do not match the verified mesh")
        for cell in mesh.cells:
            row = rows[cell.id]
            velocity = np.hypot(float(row["u"]), float(row["v"]))
            value = float(row["p"])
            if not np.isfinite([velocity, value]).all():
                raise ValueError("nonfinite native field")
            polygons.append([mesh.vertices[index] for index in cell.vertices])
            speed.append(velocity)
            pressure.append(value)
        boundaries.extend([[mesh.vertices[edge.v0], mesh.vertices[edge.v1]]
                           for edge in mesh.edges if edge.neighbour < 0])
        sources.append(dict(mesh=str(path), meshSha256=leaf["acceptedMesh"]["sha256"],
                            field=str(flow), fieldSha256=hashlib.sha256(flow.read_bytes()).hexdigest()))
        accepted = next(case for case in leaf["cases"] if case["status"] == "flow-audited")
        iterations.append(accepted["flow"]["iterations"])
    return polygons, speed, pressure, boundaries, sources, iterations


def render(directories, output, labels=None):
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT/"outputs"/"topology-plot-cache"))
    os.environ.setdefault("XDG_CACHE_HOME", str(ROOT/"outputs"/"topology-plot-cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection, LineCollection
    labels = labels or [path.name for path in directories]
    if len(labels) != len(directories):
        raise ValueError("provide one label per result directory")
    fig, axes = plt.subplots(len(directories), 3, figsize=(14, 3.9*len(directories)),
                             squeeze=False, constrained_layout=True)
    records = []
    for row, (directory, label) in enumerate(zip(directories, labels)):
        polygons, speed, pressure, boundary, sources, iterations = load_case(directory)
        points = np.concatenate([np.array(polygon) for polygon in polygons])
        xmin, ymin = points.min(axis=0); xmax, ymax = points.max(axis=0)
        for col, (values, title, cmap, unit) in enumerate([
                (None, f"{len(polygons):,} real Cut-cells", None, None),
                (speed, "Native velocity magnitude", "viridis", "m/s"),
                (pressure, "Native kinematic pressure", "coolwarm", "m²/s²")]):
            axis = axes[row, col]
            if values is None:
                collection = PolyCollection(polygons, facecolors="#edf4f8", edgecolors="#536b7b", linewidths=.23)
            else:
                collection = PolyCollection(polygons, array=np.array(values), cmap=cmap, edgecolors="none")
                fig.colorbar(collection, ax=axis, shrink=.72, label=unit)
            axis.add_collection(collection)
            axis.add_collection(LineCollection(boundary, colors="#192c3f", linewidths=.6))
            axis.set(xlim=(xmin-.02, xmax+.02), ylim=(min(0,ymin)-.02, max(1,ymax)+.02),
                     aspect="equal", xlabel="x (m)", ylabel="y (m)", title=f"{label} | {title}")
            axis.set_facecolor("#e5e8eb")
        records.append(dict(label=label, cells=len(polygons), sources=sources, iterations=iterations))
    fig.suptitle("CartMesh2D | Density topology → extracted wall → native Cut-cell CFD\n"
                 "All extracted regions retained; low-Re flow converged and independently audited; physical accuracy not qualified",
                 fontsize=12)
    fig.savefig(output, dpi=160)
    plt.close(fig)
    output.with_suffix(output.suffix+".json").write_text(json.dumps(
        dict(schema="cartmesh2d-topology-native-preview-v1", cases=records,
             fieldRendering="One constant scalar per actual native polygon; no interpolated or fabricated field."),
        indent=2, allow_nan=False)+"\n")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--labels", nargs="+")
    args = parser.parse_args()
    print(render(args.directories, args.output, args.labels))
