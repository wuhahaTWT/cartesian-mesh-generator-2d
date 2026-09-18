#!/usr/bin/env python3
"""Independent Re=100 lid-driven-cavity streamfunction/vorticity probe.

This is a diagnostic reference only.  It is deliberately separate from the
native finite-volume solver: a uniform finite-difference grid, centered
vorticity transport, Thom wall vorticity, and a SciPy DST-I Poisson solve.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import shlex
import sys
import time
from pathlib import Path

import numpy as np
from scipy import __version__ as scipy_version
from scipy.fft import dstn, idstn

SCRIPT_PATH = Path(__file__).resolve()
NU = 0.01
LID_U = 1.0
TOL_RESIDUAL = 1.0e-8
TOL_CENTERLINE = 1.0e-8
MAX_ITER = 200000

GHIA_U = (
    (1.0000, 1.00000), (0.9766, 0.84123), (0.9688, 0.78871),
    (0.9609, 0.73722), (0.9531, 0.68717), (0.8516, 0.23151),
    (0.7344, 0.00332), (0.6172, -0.13641), (0.5000, -0.20581),
    (0.4531, -0.21090), (0.2813, -0.15662), (0.1719, -0.10150),
    (0.1016, -0.06434), (0.0703, -0.04775), (0.0625, -0.04192),
    (0.0547, -0.03717), (0.0000, 0.00000),
)
GHIA_V = (
    (1.0000, 0.00000), (0.9688, -0.05906), (0.9609, -0.07391),
    (0.9531, -0.08864), (0.9453, -0.10313), (0.9063, -0.16914),
    (0.8594, -0.22445), (0.8047, -0.24533), (0.5000, 0.05454),
    (0.2344, 0.17527), (0.2266, 0.17507), (0.1563, 0.16077),
    (0.0938, 0.12317), (0.0781, 0.10890), (0.0703, 0.10091),
    (0.0625, 0.09233), (0.0000, 0.00000),
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def eigenvalues(n: int, h: float) -> np.ndarray:
    k = np.arange(1, n + 1, dtype=float)
    s2 = np.sin(np.pi * k / (2.0 * (n + 1.0))) ** 2
    return -4.0 / h**2 * (s2[:, None] + s2[None, :])


def poisson(omega: np.ndarray, lambdas: np.ndarray) -> np.ndarray:
    # Δψ=-ω, with orthonormal DST-I in both axes.
    return idstn(-dstn(omega, type=1, norm="ortho") / lambdas,
                 type=1, norm="ortho")


def poisson_validation() -> dict:
    n = 31
    h = 1.0 / (n + 1.0)
    ij = np.arange(1, n + 1, dtype=float)
    x = ij * h
    X, Y = np.meshgrid(x, x, indexing="ij")
    exact = np.sin(np.pi * X) * np.sin(np.pi * Y)
    lam11 = -4.0 / h**2 * (2.0 * np.sin(np.pi / (2.0 * (n + 1.0))) ** 2)
    omega = -lam11 * exact
    got = poisson(omega, eigenvalues(n, h))
    return {
        "nInterior": n,
        "mode": "sin(pi*x)sin(pi*y)",
        "discreteEigenvalue11": float(lam11),
        "maxAbsPsiError": float(np.max(np.abs(got - exact))),
        "centerPsi": float(got[n // 2, n // 2]),
        "expectedCenterPsi": 1.0,
        "passedRoundoff": bool(np.max(np.abs(got - exact)) < 5e-13),
    }


def apply_boundary_vorticity(omega: np.ndarray, psi: np.ndarray, h: float) -> np.ndarray:
    """Return full-grid omega, including Thom values on all four walls."""
    full = np.zeros_like(psi)
    full[1:-1, 1:-1] = omega
    inv_h2 = 1.0 / h**2
    # psi=0 on every wall; adjacent interior value is the only extrapolation.
    full[0, 1:-1] = -2.0 * psi[1, 1:-1] * inv_h2
    full[-1, 1:-1] = -2.0 * psi[-2, 1:-1] * inv_h2 - 2.0 * LID_U / h
    full[1:-1, 0] = -2.0 * psi[1:-1, 1] * inv_h2
    full[1:-1, -1] = -2.0 * psi[1:-1, -2] * inv_h2
    # Corner values are diagnostic only; use the adjacent wall formula.
    full[0, 0] = -2.0 * psi[1, 1] * inv_h2
    full[0, -1] = -2.0 * psi[1, -2] * inv_h2
    full[-1, 0] = -2.0 * psi[-2, 1] * inv_h2 - 2.0 * LID_U / h
    full[-1, -1] = -2.0 * psi[-2, -2] * inv_h2 - 2.0 * LID_U / h
    return full


def velocities(psi: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    # Array axis 0 is y and axis 1 is x: u=psi_y, v=-psi_x.
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)
    u[1:-1, 1:-1] = (psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2.0 * h)
    v[1:-1, 1:-1] = -(psi[1:-1, 2:] - psi[1:-1, :-2]) / (2.0 * h)
    u[-1, :] = LID_U
    return u, v


def rhs(omega: np.ndarray, h: float, lambdas: np.ndarray):
    psi_full = np.zeros((omega.shape[0] + 2, omega.shape[1] + 2), dtype=float)
    psi_full[1:-1, 1:-1] = poisson(omega, lambdas)
    omega_full = apply_boundary_vorticity(omega, psi_full, h)
    u, v = velocities(psi_full, h)
    wx = (omega_full[1:-1, 2:] - omega_full[1:-1, :-2]) / (2.0 * h)
    wy = (omega_full[2:, 1:-1] - omega_full[:-2, 1:-1]) / (2.0 * h)
    lap = (omega_full[1:-1, 2:] + omega_full[1:-1, :-2]
           + omega_full[2:, 1:-1] + omega_full[:-2, 1:-1]
           - 4.0 * omega) / h**2
    transport = u[1:-1, 1:-1] * wx + v[1:-1, 1:-1] * wy
    return -transport + NU * lap, psi_full, omega_full, u, v


def centerline(psi: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
    mid = psi.shape[0] // 2
    return np.concatenate((u[:, mid], v[mid, :]))


def divergence_metrics(u: np.ndarray, v: np.ndarray, h: float) -> dict:
    # Away from walls, the mixed derivatives cancel algebraically.  The
    # first interior layer is reported separately because wall velocities use
    # one-sided physical boundary values rather than streamfunction derivatives.
    div = ((u[1:-1, 2:] - u[1:-1, :-2])
           + (v[2:, 1:-1] - v[:-2, 1:-1])) / (2.0 * h)
    deep = div[1:-1, 1:-1]
    adjacent = np.concatenate((div[0, :], div[-1, :], div[:, 0], div[:, -1]))
    def stats(a: np.ndarray) -> dict:
        return {"maxAbs": float(np.max(np.abs(a))),
                "rms": float(np.sqrt(np.mean(a * a))),
                "count": int(a.size)}
    return {"interiorDistanceGe2": stats(deep),
            "boundaryAdjacent": stats(adjacent),
            "definition": "central ux+vy; deep excludes one-cell wall layer"}


def comparison(coords: np.ndarray, values: np.ndarray, table: tuple[tuple[float, float], ...]) -> dict:
    ref_x = np.asarray([p[0] for p in table], dtype=float)
    ref_y = np.asarray([p[1] for p in table], dtype=float)
    order = np.argsort(ref_x)
    sampled = np.interp(ref_x[order], coords, values)
    err = sampled - ref_y[order]
    return {"count": int(err.size), "bias": float(np.mean(err)),
            "maxAbs": float(np.max(np.abs(err))),
            "rmse": float(np.sqrt(np.mean(err * err)))}


def run_case(N: int, output_root: Path, time_limit: float) -> dict:
    start = time.monotonic()
    h = 1.0 / (N - 1.0)
    dt = min(0.2 * h, 0.18 * h * h / NU)
    n = N - 2
    lambdas = eigenvalues(n, h)
    omega = np.zeros((n, n), dtype=float)
    previous_centerline = None
    last_res = math.inf
    last_change = math.inf
    iterations = 0
    timed_out = False
    while iterations < MAX_ITER:
        if time.monotonic() - start >= time_limit:
            timed_out = True
            break
        k1, _, _, _, _ = rhs(omega, h, lambdas)
        stage = omega + dt * k1
        k2, _, _, _, _ = rhs(stage, h, lambdas)
        new_omega = omega + 0.5 * dt * (k1 + k2)
        new_residual, psi_new, omega_full, u, v = rhs(new_omega, h, lambdas)
        current = centerline(psi_new, u, v)
        if previous_centerline is not None:
            last_change = float(np.max(np.abs(current - previous_centerline)))
        # Use the residual of the accepted new state.  The RK2 stage residual
        # k2 is not the residual of new_omega and can cross the tolerance
        # first (the old diagnostic stopped early at N129 for this reason).
        last_res = float(np.max(np.abs(new_residual)))
        omega = new_omega
        previous_centerline = current
        iterations += 1
        if iterations >= 10 and last_res <= TOL_RESIDUAL and last_change <= TOL_CENTERLINE:
            break
    _, psi, omega_full, u, v = rhs(omega, h, lambdas)
    elapsed = time.monotonic() - start
    residual, _, _, _, _ = rhs(omega, h, lambdas)
    last_res = float(np.max(np.abs(residual)))
    if iterations == 0:
        last_change = 0.0
    divergence = divergence_metrics(u, v, h)

    npy = {}
    for name, array in (("omega", omega_full), ("psi", psi), ("u", u), ("v", v)):
        path = output_root / f"N{N}_{name}.npy"
        np.save(path, array)
        npy[name] = {"path": str(path.relative_to(output_root)),
                     "shape": list(array.shape), "sha256": sha256(path)}
    xs = np.linspace(0.0, 1.0, N)
    ys = xs.copy()
    mid = N // 2
    csv_path = output_root / f"N{N}_centerlines.csv"
    with csv_path.open("w") as f:
        f.write("coordinate,u_centerline,v_centerline\n")
        for i in range(N):
            f.write(f"{ys[i]:.17g},{u[i, mid]:.17g},{v[mid, i]:.17g}\n")
    ghia_u = comparison(ys, u[:, mid], GHIA_U)
    ghia_v = comparison(xs, v[mid, :], GHIA_V)
    both = np.array([ghia_u["rmse"], ghia_v["rmse"]])
    reported_change = None if not math.isfinite(last_change) else last_change
    return {
        "N": N, "h": h, "dt": dt, "nu": NU,
        "iterations": iterations, "elapsedSeconds": elapsed,
        "timeLimitSeconds": time_limit, "timedOut": timed_out,
        "converged": (not timed_out and last_res <= TOL_RESIDUAL and last_change <= TOL_CENTERLINE),
        "residualInf": last_res, "centerlineChangeInf": reported_change,
        "residualTolerance": TOL_RESIDUAL, "centerlineTolerance": TOL_CENTERLINE,
        "psiCenter": float(psi[mid, mid]), "uCenter": float(u[mid, mid]),
        "vCenter": float(v[mid, mid]), "psiSignAndCenterU": bool(psi[mid, mid] < 0 and u[mid, mid] < 0),
        "divergence": divergence,
        "ghia": {"u": ghia_u, "v": ghia_v,
                 "combinedRMSE": float(np.sqrt(np.mean(both * both)))},
        "arrays": npy,
        "centerlineCsv": str(csv_path.relative_to(output_root)),
        "centerlineCsvSha256": sha256(csv_path),
        "boundaryVorticity": "Thom: stationary -2*psi_adj/h^2; moving top adds -2/h",
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Independent centered-FD/DST-I Re=100 cavity diagnostic; not the native FVM solver."
    )
    parser.add_argument(
        "--output-root", default="outputs/native-flow/cavity-vorticity",
        help="new or empty directory for arrays, CSV files, and summary.json",
    )
    parser.add_argument(
        "--grid-points", nargs="+", type=int, default=[33, 65, 129],
        help="odd square node counts (default: 33 65 129)",
    )
    parser.add_argument(
        "--timeout", type=float, default=180.0,
        help="per-grid wall-clock limit in seconds (default: 180)",
    )
    args = parser.parse_args(argv)
    if not args.grid_points or any(n < 5 or n % 2 == 0 for n in args.grid_points):
        parser.error("--grid-points must contain odd values >= 5")
    if len(set(args.grid_points)) != len(args.grid_points):
        parser.error("--grid-points values must be unique")
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be a finite positive number")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    output_root = Path(args.output_root).expanduser()
    if output_root.exists() and any(output_root.iterdir()):
        raise RuntimeError(f"refusing non-empty output directory: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)

    poisson_check = poisson_validation()
    if not poisson_check["passedRoundoff"]:
        raise RuntimeError(f"DST-I Poisson validation failed: {poisson_check}")

    summary = {
        "format": "cartmesh2d-independent-ghia-reference-v3",
        "diagnosticOnly": True,
        "ghiaComparisonOnly": True,
        "productionSolver": "not used; this tool is independent centeredFD/DST-I and does not call native FVM output",
        "method": {
            "grid": "uniform square N by N node finite difference",
            "poisson": "SciPy orthonormal DST-I solving discrete Delta psi=-omega",
            "velocity": "u=psi_y, v=-psi_x; central differences interior",
            "vorticityTransport": "centered second-order advection and Laplacian; explicit Heun/RK2",
            "dt": "min(0.2*h, 0.18*h^2/nu)", "nu": NU,
            "boundaryVorticity": "Thom stationary walls -2*psi_adj/h^2, moving top additionally -2/h",
            "convergence": "accepted-new-state omega equation residual infinity norm and centerline infinity change",
            "timeLimitPerCaseSeconds": args.timeout,
            "maxIterations": MAX_ITER,
        },
        "parameters": {
            "outputRoot": str(output_root),
            "gridPoints": args.grid_points,
            "timeoutSecondsPerGrid": args.timeout,
            "nu": NU,
            "lidSpeed": LID_U,
            "residualTolerance": TOL_RESIDUAL,
            "centerlineTolerance": TOL_CENTERLINE,
            "maxIterations": MAX_ITER,
            "timeStepRule": "min(0.2*h, 0.18*h^2/nu)",
        },
        "command": [sys.executable, str(SCRIPT_PATH), *(sys.argv[1:] if argv is None else argv)],
        "commandString": shlex.join([sys.executable, str(SCRIPT_PATH), *(sys.argv[1:] if argv is None else argv)]),
        "scriptSHA256": sha256(SCRIPT_PATH),
        "software": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy_version},
        "poissonValidation": poisson_check,
        "ghiaSource": "GHIA_U/GHIA_V constants currently used by verify_native_flow.py; 17 points, linear interpolation",
        "ghiaQualification": "comparison only; no pass/fail gate",
        "cases": {},
    }
    for N in args.grid_points:
        summary["cases"][f"N{N}"] = run_case(N, output_root, args.timeout)
    all_converged = all(case["converged"] for case in summary["cases"].values())
    summary["success"] = bool(summary["poissonValidation"]["passedRoundoff"] and all_converged)
    path = output_root / "summary.json"
    path.write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"success": summary["success"], "poisson": summary["poissonValidation"], "cases": {
        k: {x: v[x] for x in ("iterations", "elapsedSeconds", "converged", "residualInf", "centerlineChangeInf", "uCenter", "vCenter", "timedOut")}
        for k, v in summary["cases"].items()}}, indent=2))
    return 0 if summary["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
