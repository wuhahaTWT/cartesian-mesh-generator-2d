"""Plot the structured SIMPLE run outputs. Requires NumPy and Matplotlib."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager
from matplotlib.patches import Rectangle
import numpy as np

PRESSURE_LIMITS = {"vmin": -.6, "vmax": 1.1}


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
    ax.set(xlabel="流向坐标 x/D", ylabel="横向坐标 y/D")


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
        (velocity, "速度模  |u| / U", "viridis", {}, "neither"),
        (data["p"], "压力  p / (ρU²)", "RdBu_r", PRESSURE_LIMITS, "both"),
        (data["theta"], "温度  θ", "inferno", {"vmin": 0, "vmax": 1}, "neither"),
        (omega, "涡量  ωD / U", "RdBu_r", {"vmin": -limit, "vmax": limit}, "both"),
    ]
    for ax, (values, title, color, controls, extend) in zip(axes.flat, fields):
        mesh = ax.pcolormesh(x, y, np.ma.masked_where(~fluid, values), cmap=color, **controls)
        square(ax)
        ax.set_title(title, loc="left")
        fig.colorbar(mesh, ax=ax, pad=.025, fraction=.04, extend=extend)
    fig.suptitle(f"恒温方柱绕流与换热  ·  Re = {re}，N = {n}")
    save(fig, directory, name)


def local_fields(data, directory):
    xe, ye = edges(data)
    x, y = data["x"][0], data["y"][:, 0]
    # streamline integrates on a uniform display grid, with the square masked.
    x_line, y_line = np.linspace(4, 7, 181), np.linspace(3.5, 6.5, 181)
    xx, yy = np.meshgrid(x_line, y_line)
    solid = (xx >= 4.5) & (xx <= 5.5) & (yy >= 4.5) & (yy <= 5.5)

    def velocity_on_display_grid(name):
        along_x = np.array([np.interp(x_line, x, row) for row in data[name]])
        values = np.array([np.interp(y_line, y, col) for col in along_x.T]).T
        return np.ma.masked_where(solid, values)

    u_line = velocity_on_display_grid("u")
    v_line = velocity_on_display_grid("v")
    fluid = data["fluid"] == 1
    mid = len(y_line) // 2
    for name, label, color, limits in [
        ("p", "压力等值线与流线", "RdBu_r", PRESSURE_LIMITS),
        ("u", "流向速度", "PuOr_r", {"vmin": -.25, "vmax": 1.35}),
        ("theta", "无量纲温度", "inferno", {"vmin": 0, "vmax": 1}),
    ]:
        file_name, unit = {"p": ("pressure", "p / (ρU²)"), "u": ("velocity", "u / U"),
                           "theta": ("temperature", "θ")}[name]
        fig, ax = plt.subplots(figsize=(5.3, 4.55), layout="constrained")
        values = np.ma.masked_where(~fluid, data[name])
        mesh = ax.pcolormesh(xe, ye, values, cmap=color, rasterized=True, **limits)
        fig.colorbar(mesh, ax=ax, label=unit, extend="neither" if name == "theta" else "both")
        if name == "p":
            ax.contour(x, y, values, levels=[-.6, -.4, -.2, 0, .2, .4, .6, .8, 1.],
                       colors="#334155", linewidths=.55, alpha=.55)
            ax.streamplot(x_line, y_line, u_line, v_line, color="#193952", linewidth=.6,
                          density=1.15, arrowsize=.68)
            # Explicit seeds keep the two recirculating streamlines visible.
            ends = np.flatnonzero(u_line[mid, :] < 0)
            if ends.size:
                seed_x = x_line[ends[(ends.size - 1) // 2]]
                ax.streamplot(x_line, y_line, u_line, v_line, color="#8c2238", linewidth=1.05,
                              start_points=[[seed_x, 4.88], [seed_x, 5.12]], density=2, arrowsize=.8)
        elif name == "u":
            ax.contour(x_line, y_line, u_line, levels=[0.], colors="#24152b", linewidths=.9)
        square(ax)
        ax.set(xlim=(4, 7), ylim=(3.5, 6.5), xlabel="流向坐标 x/D", ylabel="横向坐标 y/D",
               title=f"{label}  ·  Re = 20，N = 80")
        save(fig, directory, f"{file_name}-detail.png")


def grid_figure(coarse, fine, directory):
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.3), layout="constrained")
    for ax, data, label in zip(axes, (coarse, fine), ("N = 10", "N = 80")):
        x, y = edges(data)
        ax.vlines(x, 3.3, 6.7, color="#5e768f", linewidth=.35)
        ax.hlines(y, 3.3, 7.7, color="#5e768f", linewidth=.35)
        square(ax)
        ax.set(xlim=(3.3, 7.7), ylim=(3.3, 6.7), title=label)
    save(fig, directory, "mesh.png")


def wake(coarse, fine, faster, directory):
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.4), layout="constrained")
    for data, label in [(coarse, "Re = 20，N = 10"), (fine, "Re = 20，N = 80"), (faster, "Re = 40，N = 20")]:
        x = data["x"][0, :]
        for ax, name, wall in zip(axes, ("u", "theta"), (0, 1)):
            values = centreline(data, name, horizontal=True, position=5)
            downstream = x > 5.5
            ax.plot(np.r_[0, x[downstream] - 5.5], np.r_[wall, values[downstream]], label=label)
    axes[0].axhline(0, color=".45", linewidth=.8)
    axes[0].set(ylabel="u / U", xlim=(0, 5))
    axes[1].set(ylabel="θ", xlim=(0, 12))
    for ax in axes:
        ax.set_xlabel("距柱后壁的距离 / D")
        ax.legend()
    save(fig, directory, "wake.png")


def verification(cavity, outlet, reference, directory):
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.3), layout="constrained")
    for ax, name, horizontal, zero, one in zip(axes[:2], ("u", "v"), (False, True), (0, 0), (1, 0)):
        pos = cavity["x"][0, :] if horizontal else cavity["y"][:, 0]
        vals = centreline(cavity, name, horizontal=horizontal, position=.5)
        ax.plot(np.r_[0, pos, 1], np.r_[zero, vals, one], label="结构化计算")
        rows = reference[reference["component"] == name]
        ax.plot(rows["position"], rows["value"], "o", mfc="none", markersize=4, label="Ghia 等（1982）")
        ax.set(xlabel="流向坐标 x" if horizontal else "横向坐标 y", ylabel=name, title="顶盖驱动方腔 · Re = 100")
        ax.legend()
    exact = np.linspace(0, 1, 200)
    axes[2].plot(exact, 6 * exact * (1 - exact), label="抛物线解析解")
    axes[2].plot(outlet["y"], outlet["u"], "o", mfc="none", markersize=3, label="出口速度")
    axes[2].set(xlabel="横向坐标 y/H", ylabel="u / U", title="平板通道 · Re = 20")
    axes[2].legend()
    save(fig, directory, "verification.png")


def convergence(root, directory):
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.65), layout="constrained")
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for index, n in enumerate((10, 20, 40, 80)):
        history = np.genfromtxt(root / f"re20-n{n}" / "residuals.csv", delimiter=",", names=True,
                                dtype=None, encoding="utf-8")
        flow = history[history["phase"] == "flow"]
        thermal = history[history["phase"] == "thermal"]
        for key, style in (("u", "-"), ("v", "--")):
            axes[0].semilogy(flow["iteration"], flow[key], style, color=colors[index],
                            label=f"N = {n}" if key == "u" else None)
        axes[1].semilogy(flow["iteration"], flow["continuity"], color=colors[index], label=f"N = {n}")
        if n == 80:
            axes[1].semilogy(flow["iteration"], flow["corrected_continuity"], ":", color=colors[index],
                            label="N = 80 校正后", alpha=.8)
        if thermal.size:
            axes[2].semilogy(thermal["iteration"], thermal["theta"], color=colors[index], label=f"N = {n}")
    for ax, xlabel in zip(axes, ("SIMPLE 迭代步", "SIMPLE 迭代步", "能量方程迭代步")):
        ax.axhline(1e-6, color=".45", linewidth=.9, label="外迭代容差")
        ax.set(xlabel=xlabel, ylabel="归一化残差")
        ax.legend()
    axes[0].set_title("动量（u 实线，v 虚线）", loc="left")
    axes[1].set_title("连续性（预测值为实线）", loc="left")
    axes[2].set_title("能量", loc="left")
    save(fig, directory, "residuals.png")


def native_profiles(structured, native, directory):
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.15), layout="constrained")
    for data, label in [(structured, "结构化计算"), (native, "项目有限体积")]:
        x = data["x"][0, :]
        fluid = (x < 4.5) | (x > 5.5)
        for ax, name, ylabel in zip(axes, ("u", "p", "theta"), ("u / U", "p / (ρU²)", "θ")):
            values = centreline(data, name, horizontal=True, position=5)
            ax.plot(x, np.where(fluid, values, np.nan), label=label)
            ax.axvspan(4.5, 5.5, color=".94", zorder=0)
            ax.set(xlabel="流向坐标 x/D", ylabel=ylabel, xlim=(0, 15))
            ax.legend()
    save(fig, directory, "native-profiles.png")


def alternating_pressure(data):
    p = data["p"][6:-6, 6:-6]
    return (p[:-1, :-1] - p[:-1, 1:] - p[1:, :-1] + p[1:, 1:]) / 4


def checkerboard_demo(reference, direct, directory):
    fig, axes = plt.subplots(2, 2, figsize=(7, 6), layout="constrained")
    cases = (reference, direct)
    centred = [data["p"] - data["p"].mean() for data in cases]
    alternating = [alternating_pressure(data) for data in cases]
    limits = [max(np.quantile(np.abs(values), .99) for values in row) for row in (centred, alternating)]
    for i, (data, title) in enumerate(zip(cases, ("Rhie–Chow 动量插值", "面速度直接线性插值"))):
        x, y = edges(data)
        pressure = axes[0, i].pcolormesh(x, y, centred[i], cmap="RdBu_r", vmin=-limits[0], vmax=limits[0])
        local = axes[1, i].pcolormesh(data["x"][0, 6:-6], data["y"][6:-6, 0], alternating[i],
                                      cmap="RdBu_r", vmin=-limits[1], vmax=limits[1])
        axes[0, i].set_title(title)
        for ax in axes[:, i]:
            ax.set(xlabel="流向坐标 x", ylabel="横向坐标 y", aspect="equal")
    fig.colorbar(pressure, ax=axes[0], label="减去域平均的压力", extend="both")
    fig.colorbar(local, ax=axes[1], label="内部 2×2 压力交替差", extend="both")
    save(fig, directory, "checkerboard-pressure.png")


def surface_profiles(root, directory):
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.7), layout="constrained")
    faces = ("front", "top", "rear")
    for n, style, color in [(40, "--", "#9b5b19"), (80, "-", "#17658b")]:
        data = np.genfromtxt(root / f"re20-n{n}" / "surface.csv", delimiter=",", names=True,
                             dtype=None, encoding="utf-8")
        for name, ax in zip(("nu", "cp"), axes):
            for face in faces:
                on_face = (data["face"] == face) & (data["s_over_d"] <= 2)
                ax.plot(data["s_over_d"][on_face], data[name][on_face], style, color=color,
                        label=f"N = {n}" if face == "front" else None)
    for ax, ylabel in zip(axes, ("局部努塞尔数 Nu", "表面压力系数 Cp")):
        for s in (.5, 1.5):
            ax.axvline(s, color=".5", linewidth=.8)
        for s, label in zip((.25, 1., 1.75), ("迎风面", "上表面", "背风面")):
            ax.text(s, 1.02, label, transform=ax.get_xaxis_transform(), ha="center")
        ax.set(xlim=(0, 2), xlabel="从前驻点起的弧长 s/D", ylabel=ylabel)
        ax.legend(loc="best")
    save(fig, directory, "surface-profiles.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reference", type=Path,
                        default=Path(__file__).resolve().parents[1] / "examples/square-cylinder-heat/ghia-re100.csv")
    args = parser.parse_args()
    directory = args.output / "figures"
    directory.mkdir(exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    available = {font.name for font in fontManager.ttflist}
    cjk = next((name for name in ("Microsoft YaHei", "Noto Sans CJK SC", "PingFang SC", "SimHei")
                if name in available), None)
    if cjk is None:
        raise RuntimeError("绘制中文图表需要微软雅黑、Noto Sans CJK SC、苹方或黑体")
    plt.rcParams.update({"font.family": cjk, "font.size": 9, "legend.frameon": False, "axes.titleweight": "bold"})
    coarse, fine, faster = [field(args.output / name) for name in ("re20-n10", "re20-n80", "re40-n20")]
    contours(fine, directory, "re20-fields.png", 20, 80)
    contours(faster, directory, "re40-fields.png", 40, 20)
    local_fields(fine, directory)
    grid_figure(coarse, fine, directory)
    wake(coarse, fine, faster, directory)
    cavity = field(args.output / "cavity-re100")
    outlet = np.genfromtxt(args.output / "channel-re20" / "outlet.csv", delimiter=",", names=True)
    reference = np.genfromtxt(args.reference, delimiter=",", names=True, dtype=None, encoding="utf-8")
    verification(cavity, outlet, reference, directory)
    checkerboard_demo(cavity, field(args.output / "cavity-direct-interpolation"), directory)
    surface_profiles(args.output, directory)
    convergence(args.output, directory)
    native = args.output / "native-re20-n20"
    if (native / "field.csv").exists():
        native_profiles(field(args.output / "re20-n20"), field(native), directory)


if __name__ == "__main__":
    main()
