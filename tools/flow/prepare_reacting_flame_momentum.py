#!/usr/bin/env python3
"""Prepare an optional momentum-compatible planar-flame initial condition.

Use the initial T/Y of a native fixture. Adjust pressure,
density, velocity and total energy consistently with constant axial momentum
flux, including variable viscosity. This is a diagnostic initial condition,
not a coupled steady flame, a restart correction, or a change of solver defaults.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import cantera as ct
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(reference_directory, source_run, fixture_index, output):
    reference_path = reference_directory / "reference.json"
    reference = json.loads(reference_path.read_text())
    run_path = source_run / "report.json"
    run = json.loads(run_path.read_text())
    selected = [case for case in run["cases"] if case["fixture_index"] == fixture_index]
    if len(selected) != 1:
        raise ValueError("source run must contain the selected fixture exactly once")
    case = selected[0]
    fixture = reference["fixtures"][fixture_index]
    field_path = Path(case["directory"]) / "field.json"
    field = json.loads(field_path.read_text())
    if field.get("restart") is not None or field["initial"]["time"] != 0:
        raise ValueError("this prepares new initial data, not a restart state")
    nx, ny = field["nx"], field["ny"]
    u = np.asarray(field["initial"]["U"])[:nx]
    if not np.array_equal(np.tile(u, (ny, 1)), field["initial"]["U"]):
        raise ValueError("source initial data must be identical across strip rows")
    lines = Path(fixture["path"]).read_text().splitlines()
    if not np.array_equal(np.array([[float(v) for v in row.split()] for row in lines[4:4 + nx]]), u):
        raise ValueError("source initial state differs from its fixture")
    inlet = np.asarray([float(v) for v in lines[4 + nx].split()])
    if np.any(u[:, 2] != 0) or inlet[2] != 0:
        raise ValueError("planar axial flow with zero transverse velocity required")
    edges = np.asarray(field["xEdges"])
    if len(edges) != nx + 1 or nx < 3 or not np.all(np.diff(edges) > 0):
        raise ValueError("invalid source strip geometry")
    x = (edges[:-1] + edges[1:]) / 2
    temperature = np.asarray(field["initial"]["temperature"])[:nx]
    fractions = u[:, 4:] / u[:, 0, None]
    gas = ct.Solution(reference["mechanism"], transport_model="multicomponent")
    if gas.thermo_model != "ideal-gas":
        raise ValueError("this preparation requires an ideal-gas mixture")
    mass_flux = float(inlet[1])
    if not np.isfinite(mass_flux) or mass_flux <= 0:
        raise ValueError("positive axial inlet mass flux required")
    momentum_flux = reference["pressure_Pa"] + inlet[1] ** 2 / inlet[0]
    thermal_eos, viscosity = [], []
    for t, q, y in zip(temperature, u, fractions):
        gas.set_unnormalized_mass_fractions(y)
        gas.TD = t, q[0]
        thermal_eos.append(gas.P / gas.density)
        viscosity.append(gas.viscosity)
    thermal_eos, viscosity = np.asarray(thermal_eos), np.asarray(viscosity)
    pressure = np.asarray(field["initial"]["pressure"])[:nx].copy()
    iterations = []
    for iteration in range(50):
        velocity = mass_flux * thermal_eos / pressure
        stress = (4 / 3) * viscosity * np.gradient(velocity, x, edge_order=2)
        discriminant = (momentum_flux + stress) ** 2 - 4 * mass_flux ** 2 * thermal_eos
        if not np.all(np.isfinite(discriminant)) or not np.all(discriminant > 0):
            raise ValueError("no separated low-Mach momentum root at the source T/Y")
        updated = ((momentum_flux + stress) + np.sqrt(discriminant)) / 2
        delta = float(np.max(abs(updated - pressure)))
        iterations.append({"iteration": iteration + 1, "maximum_pressure_update_Pa": delta})
        pressure = updated
        # Input fixed-point termination at a pressure-scaled floating-point
        # budget, not a physical residual or combustion-accuracy threshold.
        if delta <= 64 * np.finfo(float).eps * max(abs(pressure)):
            break
    else:
        raise ValueError("initial-pressure fixed point did not converge within 50 iterations")
    prepared = []
    for t, p, y in zip(temperature, pressure, fractions):
        gas.set_unnormalized_mass_fractions(y)
        gas.TP = t, p
        rho = gas.density
        prepared.append([rho, mass_flux, 0., rho * gas.int_energy_mass + mass_flux ** 2 / (2 * rho), *(rho * y)])
    prepared = np.asarray(prepared)
    if not np.all(np.isfinite(prepared)) or not np.all(prepared[:, 0] > 0) or not np.all(prepared[:, 4:] >= 0):
        raise ValueError("invalid prepared conservative state")
    velocity = prepared[:, 1] / prepared[:, 0]
    stress = (4 / 3) * viscosity * np.gradient(velocity, x, edge_order=2)
    momentum_defect = pressure + mass_flux * velocity - stress - momentum_flux
    widths = np.diff(edges)
    before, after = (np.sum(widths[:, None] * values, axis=0) for values in (u, prepared))
    output.mkdir(parents=True, exist_ok=False)
    path = output / "grid.fixture"
    with path.open("w") as stream:
        stream.write("\n".join(lines[:4]) + "\n")
        for row in prepared:
            stream.write(" ".join(format(v, ".17g") for v in row) + "\n")
        stream.write(lines[4 + nx] + "\nEND\n")
    profile_file = reference.get("profile_file", "profile-2.json")
    shutil.copyfile(reference_directory / profile_file, output / profile_file)
    reference.update(
        scope="Optional momentum-compatible initial data at original T/Y; not a new coupled steady BVP or flame qualification.",
        source_reference={"path": str(reference_path), "sha256": sha(reference_path)},
        fixtures=[{**fixture, "path": str(path), "sha256": sha(path), "original_fixture_index": fixture_index}])
    reference["files_sha256"] = {p.name: sha(p) for p in (path, output / profile_file)}
    (output / "reference.json").write_text(json.dumps(reference, indent=2) + "\n")
    result = {
        "source_run": {"path": str(run_path), "sha256": sha(run_path)},
        "source_kind": "initial_evaluation" if run.get("initial_evaluation_only", False) else "evolution_initial_state",
        "source_field": {"path": str(field_path), "sha256": sha(field_path)},
        "source_fixture": {"path": fixture["path"], "sha256": fixture["sha256"]},
        "script_sha256": sha(__file__), "output_fixture_sha256": sha(path), "cantera": ct.__version__,
        "method": "Hold cell T/Y and the original inlet state; solve p + m^2/rho - tau_xx = inlet momentum flux with rho=p/(Rmix*T), tau_xx=(4/3)*mu*du/dx, then reconstruct conservative initial data.",
        "gradient": "Nonuniform second-order cell-centre finite difference; not the native discrete momentum flux.",
        "fixed_point_stop": {"quantity": "maximum pressure update", "unit": "Pa", "relative_roundoff_budget": float(64 * np.finfo(float).eps), "maximum_iterations": 50},
        "iterations": iterations,
        "maximum_preparation_momentum_flux_defect_Pa": float(np.max(abs(momentum_defect))),
        "pressure_range_Pa": [float(min(pressure)), float(max(pressure))],
        "maximum_pressure_change_Pa": float(np.max(abs(pressure - np.asarray(field["initial"]["pressure"])[:nx]))),
        "maximum_mass_fraction_change": float(np.max(abs(prepared[:, 4:] / prepared[:, 0, None] - fractions))),
        "integrated_conserved_initial_change_per_height": (after - before).tolist(),
        "relative_initial_mass_change": float((after[0] - before[0]) / before[0]),
        "native_solver_boundary_or_defaults_changed": False,
        "qualified_flame": False,
        "limits": ["This creates new initial rho/velocity/energy and changes initial integrals; it must not be applied to an accepted restart.",
                   "It does not solve steady chemistry, species, or energy balance. Full native evolution is still required to assess the prepared state."]}
    (output / "preparation.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--source-run", type=Path, required=True)
    parser.add_argument("--grid", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.reference.resolve(), args.source_run.resolve(), args.grid, args.output.resolve())
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
