"""Artifacts derived from actual topology analysis fields, with exact provenance."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np


def contours(rho, width, height, threshold=.5):
    """Piecewise-linear contour of centre samples with edge extension.

    No cell-centre deletion, polyline smoothing, dropped components or filled
    holes. contourpy returns closed outer loops and their holes separately.
    Native geometry/quality checks remain authoritative on the extracted XY.
    """
    import contourpy
    if (rho.ndim != 2 or not np.isfinite(rho).all() or not 0 < threshold < 1 or
            width <= 0 or height <= 0):
        raise ValueError("invalid field or extraction threshold")
    ny, nx = rho.shape
    x = np.r_[0, (np.arange(nx)+.5)*width/nx, width]
    y = np.r_[0, (np.arange(ny)+.5)*height/ny, height]
    field = np.pad(rho, 1, mode="edge")
    generator = contourpy.contour_generator(x=x, y=y, z=field, fill_type="OuterOffset",
                                            corner_mask=False)
    points, offsets = generator.filled(threshold, max(1.0, float(np.max(rho)))+1)
    groups = []
    for vertices, boundaries in zip(points, offsets):
        loops = []
        for a, b in zip(boundaries[:-1], boundaries[1:]):
            loop = vertices[a:b]
            if np.array_equal(loop[0], loop[-1]):
                loop = loop[:-1]
            # Only identical consecutive points are omitted; all coordinates
            # and all nonzero line segments from the contour are preserved.
            keep = np.r_[True, np.any(loop[1:] != loop[:-1], axis=1)]
            loop = loop[keep]
            if len(loop) < 3 or abs(signed_area(loop)) <= 0:
                raise ValueError("zero-area contour component; extraction rejected")
            loops.append(loop)
        groups.append(loops)
    if not groups:
        raise ValueError("threshold produced no fluid region")
    return groups


def signed_area(loop):
    return float(np.sum(loop[:, 0]*np.roll(loop[:, 1], -1)-
                        loop[:, 1]*np.roll(loop[:, 0], -1))/2)


def port_connectivity(groups, problem):
    """Geometric port reachability, not a mixing/flow-transfer matrix."""
    width, height = problem["width"], problem["height"]
    centres = [height/4, 3*height/4] if problem["case"] == "double-pipe" else [height/4]
    right = centres if problem["case"] == "double-pipe" else [3*height/4]
    ports = [(f"inlet_{i}", 0., y) for i, y in enumerate(centres)]
    ports += [(f"outlet_{i}", width, y) for i, y in enumerate(right)]
    tolerance = 1e-11+1e-9*max(width, height)
    coverage = np.zeros((len(groups), len(ports)))
    for component, group in enumerate(groups):
        loop = group[0]
        for index, (_, x, y) in enumerate(ports):
            intervals = []
            lo, hi = y-problem["port_width"]/2, y+problem["port_width"]/2
            for a, b in zip(loop, np.roll(loop, -1, axis=0)):
                if abs(a[0]-x) <= tolerance and abs(b[0]-x) <= tolerance:
                    start, end = max(lo, min(a[1], b[1])), min(hi, max(a[1], b[1]))
                    if end > start:
                        intervals.append((start, end))
            # Union avoids double counting shared endpoint representations.
            end = -float("inf")
            for start, stop in sorted(intervals):
                coverage[component, index] += max(0, stop-max(start, end))
                end = max(end, stop)
    present = coverage > tolerance
    count = len(centres)
    matrix = [[bool(np.any(present[:, i] & present[:, count+j])) for j in range(len(right))]
              for i in range(count)]
    total = np.sum(coverage, axis=0)
    return dict(ports=[name for name, _, _ in ports], coverageLengths=coverage.tolist(),
                allPortsCovered=bool(np.all(np.abs(total-problem["port_width"]) <= tolerance)),
                componentPorts=[[ports[i][0] for i in range(len(ports)) if present[c, i]] for c in range(len(groups))],
                inletToOutletReachability=matrix,
                meaning="Paths within the same fluid component; not measured flow splitting or molecular mixing.")


def extract(root, threshold=.5, problem=None):
    root = Path(root)
    with np.load(root/"final.npz") as f:
        rho, width, height = f["rho"], float(f["width"]), float(f["height"])
    groups = contours(rho, width, height, threshold)
    if problem is None and (root/"summary.json").exists():
        problem = json.loads((root/"summary.json").read_text()).get("problem")
    def write_boundary(path, selected):
        with path.open("w") as f:
            f.write("# Fluid interior of rho isocontour; use explicit interior mode.\n")
            f.write(f"# threshold={threshold:.17g}; no smoothing or component removal\n")
            for group in selected:
                for loop in group:
                    f.writelines(f"{x:.17g} {y:.17g}\n" for x, y in loop)
                    f.write("\n")
    path = root/"fluid.xy"
    write_boundary(path, groups)
    regions = []
    for index, group in enumerate(groups):
        component = root/f"fluid-component-{index}.xy"
        write_boundary(component, [group])
        regions.append(dict(file=component.name, holes=len(group)-1,
                            fluidArea=abs(signed_area(group[0]))-sum(abs(signed_area(h)) for h in group[1:]),
                            boundarySha256=hashlib.sha256(component.read_bytes()).hexdigest()))
    fluid_area = sum(abs(signed_area(group[0]))-sum(abs(signed_area(h)) for h in group[1:])
                     for group in groups)
    records = dict(schema="cartmesh2d-topology-extraction-v1", threshold=threshold,
                   field="filtered projected rho at cell centres with edge extension",
                   fluidRegion="interior", components=len(groups),
                   holes=sum(len(group)-1 for group in groups), regions=regions,
                   vertices=sum(len(loop) for group in groups for loop in group),
                   fluidArea=fluid_area, designArea=width*height,
                   extractedVolumeFraction=fluid_area/(width*height),
                   sourceFieldSha256=hashlib.sha256((root/"final.npz").read_bytes()).hexdigest(),
                   boundarySha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                   topologyQualified=False, solverQualityQualified=False,
                   note="Contour construction only. Sharp-wall extraction changes the porous model; native checks required.")
    if problem is not None:
        records["portConnectivity"] = port_connectivity(groups, problem)
    (root/"extraction.json").write_text(json.dumps(records, indent=2)+"\n")
    # Standalone vector preview of exactly the exported polygon coordinates.
    paths = []
    for group in groups:
        commands = []
        for loop in group:
            commands.append("M "+" L ".join(f"{x:.12g},{height-y:.12g}" for x, y in loop)+" Z")
        paths.append(f'<path d="{" ".join(commands)}" fill="#e1f5ff" fill-rule="evenodd" stroke="#156e96" stroke-width="{width/700}"/>')
    (root/"fluid.svg").write_text(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">'
                                  f'<rect width="{width}" height="{height}" fill="#24364b"/>'+
                                  "".join(paths)+"</svg>\n")
    return records


def render(root):
    root = Path(root).resolve()
    os.environ.setdefault("MPLCONFIGDIR", str(root/"plot-cache"))
    os.environ.setdefault("XDG_CACHE_HOME", str(root/"plot-cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    report = json.loads((root/"summary.json").read_text())
    with np.load(root/"final.npz") as f:
        rho, u, v, pressure = (f[k] for k in ("rho", "u", "v", "p"))
        width, height = float(f["width"]), float(f["height"])
    with np.load(root/"reference-uniformPorous.npz") as f:
        initial = f["rho"].copy()
    ny, nx = rho.shape
    xs, ys = (np.arange(nx)+.5)*width/nx, (np.arange(ny)+.5)*height/ny
    extent = (0, width, 0, height)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.3), constrained_layout=True)
    for axis, field, title in ((axes[0, 0], initial, "Initial porous design (same final model)"),
                               (axes[0, 1], rho, "Optimised fluid fraction / extracted wall")):
        im = axis.imshow(field, origin="lower", extent=extent, vmin=0, vmax=1, cmap="gray")
        axis.set(title=title, xlabel="x", ylabel="y", aspect="equal")
    threshold = report.get("extraction", {}).get("threshold", .5)
    groups = contours(rho, width, height, threshold)
    for group in groups:
        for loop in group:
            closed = np.vstack((loop, loop[0]))
            axes[0, 1].plot(closed[:, 0], closed[:, 1], color="#ea714b", lw=.85)
    fig.colorbar(im, ax=axes[0, :2], shrink=.7, label="0 = solid resistance, 1 = fluid")
    speed = np.hypot(u, v)
    image = axes[0, 2].imshow(speed, origin="lower", extent=extent, cmap="viridis")
    axes[0, 2].streamplot(xs, ys, u, v, density=.8, color="white", linewidth=.4, arrowsize=.6)
    axes[0, 2].set(title="Computed velocity (Brinkman analysis)", xlabel="x", ylabel="y")
    fig.colorbar(image, ax=axes[0, 2], shrink=.7, label="Nondimensional speed")
    shown = np.ma.masked_where(rho < threshold, pressure)
    image = axes[1, 0].imshow(shown, origin="lower", extent=extent, cmap="coolwarm")
    axes[1, 0].set_facecolor("#e6e8ea")
    axes[1, 0].set(title="Brinkman pressure (rho above threshold)", xlabel="x", ylabel="y")
    fig.colorbar(image, ax=axes[1, 0], shrink=.7, label="Nondimensional pressure")
    for stage in report["stages"]:
        rows = [r for r in report["history"] if r["stage"] == stage["stage"]]
        x = [r["iteration"] for r in rows]
        values = [r["objective"] for r in rows]
        if rows:
            x.append(x[-1]+1)
            values.append(stage["final"]["objective"])
        axes[1, 1].plot(x, values, "-", lw=1.7, label=f"q={stage['q']}, beta={stage['beta']}")
    axes[1, 1].set(title="Objective within each continuation stage", xlabel="Accepted step", ylabel="Dissipation / selected objective")
    axes[1, 1].legend(fontsize=8); axes[1, 1].grid(alpha=.2)
    rows = report["history"]
    axes[1, 2].plot([r["iteration"] for r in rows], [r["grayness"] for r in rows], label="Grayness")
    axes[1, 2].plot([r["iteration"] for r in rows], [r["volumeFraction"] for r in rows], label="Fluid fraction")
    axes[1, 2].axhline(report["problem"]["volume_fraction"], color="black", ls="--", lw=.8, label="Volume limit")
    axes[1, 2].set(title="Constraint and intermediate material", xlabel="Accepted step", ylim=(0, 1))
    axes[1, 2].legend(fontsize=8); axes[1, 2].grid(alpha=.2)
    fig.suptitle(f"CartMesh2D fluid topology prototype | {report['problem']['case']} | {nx} x {ny}\n"
                 f"Stokes-Brinkman model; {report['status']}; physical accuracy not qualified", fontsize=13)
    fig.savefig(root/"overview.png", dpi=160)
    plt.close(fig)
    return root/"overview.png"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--threshold", type=float, default=.5)
    args = parser.parse_args()
    extraction = extract(args.directory, args.threshold)
    path = args.directory/"summary.json"
    report = json.loads(path.read_text())
    report["extraction"] = extraction
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    render(args.directory)
    print(json.dumps(extraction, indent=2))
