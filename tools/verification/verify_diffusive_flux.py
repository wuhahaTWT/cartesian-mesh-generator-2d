#!/usr/bin/env python3
"""Check native multicomponent fluxes using a separate Stefan-Maxwell solve.

Binary collision/Soret/enthalpy data still come from Python Cantera. The linear
flux solve is independent of its multicomponent D matrix; this is constitutive
and interface evidence, not experimental transport or flame validation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import cantera as ct
import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stefan_maxwell(gas, driving):
    """Mass-flux unknowns, no division by mole fractions of absent species.

    d_k = sum_{j != k} X_k X_j (V_j - V_k)/D_kj; sum(j_k) = 0.
    Substituting V_k = j_k/(rho Y_k) cancels X_k/Y_k = Wbar/W_k.
    """
    x, mw, mean, rho = gas.X, gas.molecular_weights, gas.mean_molecular_weight, gas.density
    binary = gas.binary_diff_coeffs
    n = gas.n_species
    matrix = np.zeros((n, n))
    for k in range(n):
        for j in range(n):
            if j != k:
                matrix[k, j] = x[k] * mean / (rho * mw[j] * binary[k, j])
                matrix[k, k] -= x[j] * mean / (rho * mw[k] * binary[k, j])
    constraint = int(np.argmax(x))
    matrix[constraint, :] = 1
    rhs = np.array(driving, copy=True)
    rhs[constraint] = 0
    result = np.linalg.solve(matrix, rhs)
    residual = matrix @ result - rhs
    scale = np.abs(matrix) @ np.abs(result) + np.abs(rhs)
    error = float(np.max(np.abs(residual) / np.maximum(scale, np.finfo(float).tiny)))
    return result, error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--mechanism-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    probe_digest = digest(args.probe)
    rng = np.random.default_rng(93271)
    cases = []
    # No spatial/experimental threshold: bound native arithmetic, the independent
    # dense linear solve, and roundoff from the conservative energy inversion.
    tolerance = 2e-9
    for mechanism, composition in [
        ("h2o2.yaml", "H2:2,O2:1,N2:3.76"),
        ("gri30.yaml", "CH4:1,O2:2,N2:7.52"),
        ("h2o2.yaml", "H2:1,N2:1"),
        ("h2o2.yaml", "H2O:1"),
    ]:
        path = (args.mechanism_root / mechanism).resolve()
        gas = ct.Solution(str(path), transport_model="multicomponent")
        for case in ("concentration", "pressure", "soret", "combined", "balanced_pressure", "all_species"):
            temperature, pressure = (950, 101325) if mechanism == "h2o2.yaml" else (1300, 7 * 101325)
            gas.TPX = temperature, pressure, composition
            if case == "all_species":
                gas.TPX = temperature, pressure, np.exp(rng.uniform(-12, 0, gas.n_species))
            x = gas.X.copy()
            gradient = x * rng.uniform(-200, 200, gas.n_species)
            gradient -= x * gradient.sum()
            pivot = int(np.argmax(x))
            gradient[pivot] = -sum(gradient[k] for k in range(gas.n_species) if k != pivot)
            grad_t, grad_p = 0.0, 0.0
            if case == "pressure":
                gradient[:] = 0
                grad_p = 120.0
            elif case == "soret":
                gradient[:] = 0
                grad_t = 2e5
            elif case in ("combined", "all_species"):
                grad_t, grad_p = -1.5e5, 80.0
            elif case == "balanced_pressure":
                grad_p = 120.0
                gradient = -(gas.X - gas.Y) * grad_p
                gradient[pivot] = -sum(gradient[k] for k in range(gas.n_species) if k != pivot)
            drive = gradient + (gas.X - gas.Y) * grad_p
            mass_flux, linear_error = stefan_maxwell(gas, drive)
            thermal = -gas.thermal_diff_coeffs * grad_t / gas.T
            # Exact kinetic-theory limits, not the backend's artificial traces.
            thermal[x == 0] = 0
            if np.count_nonzero(x) == 1:
                thermal[:] = 0
            expected = mass_flux + thermal
            h = gas.partial_molar_enthalpies / gas.molecular_weights
            conduction = -gas.thermal_conductivity * grad_t
            enthalpy = h @ expected
            # Scales use the separately evaluated driving contributions, so a
            # near-zero equilibrium/cancelled net flux is not its own divisor.
            composition_flux, _ = stefan_maxwell(gas, gradient)
            pressure_flux, _ = stefan_maxwell(gas, (gas.X - gas.Y) * grad_p)
            flux_scale = float(np.sum(abs(composition_flux) + abs(pressure_flux) + abs(thermal)))
            heat_scale = float(abs(conduction) + np.sum(abs(h) * (abs(composition_flux) + abs(pressure_flux) + abs(thermal))))
            payload = " ".join(format(float(v), ".17g") for v in [temperature, pressure, grad_t, grad_p, *x, *gradient]) + "\n"
            index = len(cases)
            (args.output / f"case-{index:02d}.in").write_text(payload)
            if digest(args.probe) != probe_digest:
                raise RuntimeError("probe executable changed during validation; keep output and rerun after build completion")
            process = subprocess.run([str(args.probe.resolve()), str(path)], input=payload, text=True, capture_output=True, check=False)
            (args.output / f"case-{index:02d}.stdout").write_text(process.stdout)
            (args.output / f"case-{index:02d}.stderr").write_text(process.stderr)
            record = {"mechanism": mechanism, "mechanismSha256": digest(path), "composition": composition,
                      "case": case, "temperature": temperature, "pressure": pressure, "X": x.tolist(),
                      "gradX": gradient.tolist(), "gradT": grad_t, "gradLogP": grad_p,
                      "independentSpeciesFlux": expected.tolist(), "independentConduction": float(conduction),
                      "independentEnthalpy": float(enthalpy), "linearSolveRelative": linear_error,
                      "returncode": process.returncode, "passed": False}
            if process.returncode == 0:
                native = json.loads(process.stdout)
                j = np.array(native["species"])
                species_error = float(np.max(abs(j - expected)) / max(flux_scale, 1e-300))
                energy_error = float(abs(native["energy"] - (conduction + enthalpy)) / max(heat_scale, 1e-300))
                mass_error = abs(float(j.sum())) / max(flux_scale, 1e-300)
                record.update(native=native, speciesRelative=species_error, energyRelative=float(energy_error), massRelative=mass_error)
                record["passed"] = bool(max(species_error, energy_error, mass_error, linear_error) < tolerance)
                if composition == "H2:1,N2:1" and case in ("concentration", "pressure"):
                    k, other = gas.species_index("H2"), gas.species_index("N2")
                    binary_j = -gas.density * gas.molecular_weights[k] * gas.molecular_weights[other] / gas.mean_molecular_weight**2 * gas.binary_diff_coeffs[k, other] * drive[k]
                    binary_error = float(abs(native["species"][k] - binary_j) / max(flux_scale, 1e-300))
                    record["binaryAnalyticRelative"] = float(binary_error)
                    record["passed"] = record["passed"] and binary_error < tolerance
            cases.append(record)
    if digest(args.probe) != probe_digest:
        raise RuntimeError("probe executable changed during validation; results are not a single-build qualification")
    report = {"passed": all(c["passed"] for c in cases), "cases": cases, "tolerance": tolerance,
              "cantera": ct.__version__, "numpy": np.__version__, "probeSha256": probe_digest,
              "oracle": "Independent Stefan-Maxwell mass-flux linear system and binary analytic limit; collision, thermal diffusion and enthalpy data from Python Cantera",
              "physicalQualification": False}
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"passed": report["passed"], "cases": len(cases),
                      "failed": [i for i, c in enumerate(cases) if not c["passed"]]}, indent=2))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
