"""Plot the structured SIMPLE run outputs. Requires NumPy and Matplotlib."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np


def field(directory):
    values = np.genfromtxt(directory / "field.csv", delimiter=",", names=True)
    return values.reshape(int(values["j"].max()) + 1, int(values["i"].max()) + 1)


def centreline(data, name, *, horizontal, position):
    positions = data["y"][:, 0] if horizontal else data["x"][0, :]
    values = data[name] if horizontal else data[name].T
    index = np.searchsorted(positions, position)
    weight = (position - positions[index - 1]) / (positions[index] - positions[index - 1])
    return (1 - weight) * values[index - 1] + weight * values[index]


def save(fig, directory, name):
    fig.savefig(directory / name, dpi=205, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def square(ax):
    ax.add_patch(Rectangle((4.5, 4.5), 1, 1, color="#263444", zorder=5))
    ax.set_aspect("equal")
    ax.set(xlabel="x/D", ylabel="y/D")


def edges(data):
    x = data["x"][0, :]
    y = data["y"][:, 0]
    dx = data["dx"][0, :]
    dy = data["dy"][:, 0]
    return np.r_[x - dx / 2, x[-1] + dx[-1] / 2], np.r_[y - dy / 2, y[-1] + dy[-1] / 2]


def contours(data, directory, name, re, n):
    x, y = edges(data)
    fluid = data["fluid"] == 1
    velocity = np.hypot(data["u"], data["v"])
    omega = data["vorticity"]
    limit = np.quantile(np.abs(omega[fluid]), .985)
    fig, axes = plt.subplots(2, 2, figsize=(11, 5.9), layout="constrained")
    fields = [
        (velocity, "Speed  |u| / U", "viridis", {}),
        (data["p"], "Pressure  p / (ρU²)", "cividis", {}),
        (data["theta"], "Temperature  θ", "inferno", {"vmin": 0, "vmax": 1}),
        (omega, "Vorticity  ωD / U", "RdBu_r", {"vmin": -limit, "vmax": limit}),
    ]
    for ax, (values, title, color, controls) in zip(axes.flat, fields):
        mesh = ax.pcolormesh(x, y, np.ma.masked_where(~fluid, values), cmap=color, **controls)
        square(ax)
        ax.set_title(title, loc="left")
        fig.colorbar(mesh, ax=ax, pad=.025, fraction=.04, extend="both" if values is omega else "neither")
    fig.suptitle(f"Heated square cylinder  |  Re = {re}, N = {n}")
    save(fig, directory, name)


def grid_figure(coarse, fine, directory):
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.3), layout="constrained")
    for ax, data, label in zip(axes, (coarse, fine), ("N = 10", "N = 40")):
        x, y = edges(data)
        ax.vlines(x, 3.3, 6.7, color="#5e768f", linewidth=.35)
        ax.hlines(y, 3.3, 7.7, color="#5e768f", linewidth=.35)
        square(ax)
        ax.set(xlim=(3.3, 7.7), ylim=(3.3, 6.7), title=label)
    save(fig, directory, "mesh.png")


def wake(coarse, fine, faster, directory):
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.4), layout="constrained")
    for data, label in [(coarse, "Re 20, N 10"), (fine, "Re 20, N 40"), (faster, "Re 40, N 20")]:
        x = data["x"][0, :]
        for ax, name, wall in zip(axes, ("u", "theta"), (0, 1)):
            values = centreline(data, name, horizontal=True, position=5)
            downstream = x > 5.5
            ax.plot(np.r_[0, x[downstream] - 5.5], np.r_[wall, values[downstream]], label=label)
    axes[0].axhline(0, color=".45", linewidth=.8)
    axes[0].set(ylabel="u / U", xlim=(0, 5))
    axes[1].set(ylabel="θ", xlim=(0, 12))
    for ax in axes:
        ax.set_xlabel("Distance from rear wall / D")
        ax.legend()
    save(fig, directory, "wake.png")


def verification(cavity, outlet, reference, directory):
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.3), layout="constrained")
    for ax, name, horizontal, zero, one in zip(axes[:2], ("u", "v"), (False, True), (0, 0), (1, 0)):
        pos = cavity["x"][0, :] if horizontal else cavity["y"][:, 0]
        vals = centreline(cavity, name, horizontal=horizontal, position=.5)
        ax.plot(np.r_[0, pos, 1], np.r_[zero, vals, one], label="Structured SIMPLE")
        rows = reference[reference["component"] == name]
        ax.plot(rows["position"], rows["value"], "o", mfc="none", markersize=4, label="Ghia et al. (1982)")
        ax.set(xlabel="x" if horizontal else "y", ylabel=name, title="Lid-driven cavity, Re = 100")
        ax.legend()
    exact = np.linspace(0, 1, 200)
    axes[2].plot(exact, 6 * exact * (1 - exact), label="Poiseuille")
    axes[2].plot(outlet["y"], outlet["u"], "o", mfc="none", markersize=3, label="Outlet")
    axes[2].set(xlabel="y/H", ylabel="u / U", title="Plane channel, Re = 20")
    axes[2].legend()
    save(fig, directory, "verification.png")


def convergence(root, directory):
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.4), layout="constrained")
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for index, n in enumerate((10, 20, 40)):
        history = np.genfromtxt(root / f"re20-n{n}" / "residuals.csv", delimiter=",", names=True,
                                dtype=None, encoding="utf-8")
        flow = history[history["phase"] == "flow"]
        thermal = history[history["phase"] == "thermal"]
        for key, style in (("u", "-"), ("v", "--"), ("continuity", ":")):
            axes[0].semilogy(flow["iteration"], flow[key], style, color=colors[index],
                            label=f"N {n}" if key == "u" else None)
        if thermal.size:
            axes[1].semilogy(thermal["iteration"], thermal["theta"], color=colors[index], label=f"N {n}")
    for ax, xlabel in zip(axes, ("SIMPLE iteration", "Energy iteration")):
        ax.axhline(1e-6, color=".45", linewidth=.9, label="Tolerance")
        ax.set(xlabel=xlabel, ylabel="Normalised residual")
        ax.legend()
    axes[0].set_title("Flow: u (solid), v (dashed), mass (dotted)", loc="left")
    axes[1].set_title("Temperature", loc="left")
    save(fig, directory, "residuals.png")


def native_profiles(structured, native, directory):
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.15), layout="constrained")
    for data, label in [(structured, "Structured central"), (native, "Native FV")]:
        x = data["x"][0, :]
        fluid = (x < 4.5) | (x > 5.5)
        for ax, name, ylabel in zip(axes, ("u", "p", "theta"), ("u / U", "p / (ρU²)", "θ")):
            values = centreline(data, name, horizontal=True, position=5)
            ax.plot(x, np.where(fluid, values, np.nan), label=label)
            ax.axvspan(4.5, 5.5, color=".94", zorder=0)
            ax.set(xlabel="x/D", ylabel=ylabel, xlim=(0, 15))
            ax.legend()
    save(fig, directory, "native-profiles.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reference", type=Path,
                        default=Path(__file__).resolve().parents[1] / "examples/square-cylinder-heat/ghia-re100.csv")
    args = parser.parse_args()
    directory = args.output / "figures"
    directory.mkdir(exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update({"font.size": 9, "legend.frameon": False, "axes.titleweight": "bold"})
    coarse, fine, faster = [field(args.output / name) for name in ("re20-n10", "re20-n40", "re40-n20")]
    contours(fine, directory, "re20-fields.png", 20, 40)
    contours(faster, directory, "re40-fields.png", 40, 20)
    grid_figure(coarse, fine, directory)
    wake(coarse, fine, faster, directory)
    cavity = field(args.output / "cavity-re100")
    outlet = np.genfromtxt(args.output / "channel-re20" / "outlet.csv", delimiter=",", names=True)
    reference = np.genfromtxt(args.reference, delimiter=",", names=True, dtype=None, encoding="utf-8")
    verification(cavity, outlet, reference, directory)
    convergence(args.output, directory)
    native = args.output / "native-re20-n20"
    if (native / "field.csv").exists():
        native_profiles(field(args.output / "re20-n20"), field(native), directory)


if __name__ == "__main__":
    main()
