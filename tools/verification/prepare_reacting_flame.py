#!/usr/bin/env python3
"""Solve an independent planar flame and prepare conservative native fixtures.

Cantera's steady 1D BVP is independent of the native transient FV algorithm but
shares kinetics and transport data. The reference's floating-point species
defects are retained and explicitly bounded during initial-data preparation.
No native time step calls this preparation or changes a species by clipping.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import cantera as ct
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare_species(raw, absolute_tolerance):
    raw = np.asarray(raw)
    allowance = 10 * absolute_tolerance + 64 * (raw.shape[1] + 1) * np.finfo(float).eps
    if not np.all(np.isfinite(raw)):
        raise ValueError("nonfinite reference species")
    prepared = np.maximum(raw, 0)
    for row in prepared:
        dependent = int(np.argmax(row))
        row[dependent] = 1 - float(np.sum(np.delete(row, dependent), dtype=np.longdouble))
    maximum = float(np.max(abs(prepared - raw)))
    if np.min(prepared) < 0 or maximum > allowance:
        raise ValueError(f"reference preparation exceeds numerical allowance: {maximum} > {allowance}")
    return prepared, allowance, maximum


def integrate_profile(edges, grid, conserved):
    result = []
    for a, b in zip(edges[:-1], edges[1:]):
        points = np.r_[a, grid[(grid > a) & (grid < b)], b]
        values = np.array([np.interp(points, grid, column) for column in conserved.T]).T
        result.append(np.trapezoid(values, points, axis=0) / (b - a))
    return np.array(result)


def fixture_grid(grid, temperature, spacing):
    centre = grid[np.argmax(np.gradient(temperature, grid))]
    # Resolve the H2 reference reaction layer and its downstream relaxation;
    # retain the complete upstream/downstream domain with graded spacing.
    left, right = max(grid[0], centre - .0006), min(grid[-1], centre + .002)
    core = np.linspace(left, right, int(np.ceil((right - left) / spacing)) + 1)
    before, after = [left], [right]
    dx = spacing
    while before[-1] > grid[0]:
        before.append(max(grid[0], before[-1] - dx)); dx = min(dx * 1.15, 2.5e-4)
    dx = spacing
    while after[-1] < grid[-1]:
        after.append(min(grid[-1], after[-1] + dx)); dx = min(dx * 1.15, 2.5e-4)
    return np.r_[before[:0:-1], core, after[1:]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--composition", default="H2:1.6,O2:1,N2:3.76")
    parser.add_argument("--temperature", type=float, default=400)
    parser.add_argument("--pressure", type=float, default=101325)
    parser.add_argument("--width", type=float, default=.02)
    parser.add_argument("--spacings", nargs="+", type=float, default=[40e-6, 20e-6, 10e-6])
    args = parser.parse_args()
    if any(not np.isfinite(x) or x <= 0 for x in [args.temperature, args.pressure, args.width, *args.spacings]):
        raise ValueError("all dimensional inputs must be positive and finite")
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    gas = ct.Solution(str(args.mechanism.resolve()))
    gas.TPX = args.temperature, args.pressure, args.composition
    flame = ct.FreeFlame(gas, width=args.width)
    flame.transport_model = "multicomponent"
    flame.soret_enabled = True
    absolute_tolerance = 1e-15
    flame.flame.set_steady_tolerances(default=(1e-9, absolute_tolerance))
    flame.flame.set_transient_tolerances(default=(1e-7, absolute_tolerance))
    flame.set_max_grid_points(flame.flame, 1500)
    stages = []
    for stage, criterion in enumerate([.06, .03, .015]):
        flame.set_refine_criteria(ratio=3, slope=criterion, curve=2 * criterion, prune=0)
        flame.solve(loglevel=1, auto=(stage == 0))
        flame.save(args.output / f"reference-{stage}.yaml", name="flame")
        raw = {"grid": flame.grid.tolist(), "T": flame.T.tolist(), "Y": flame.Y.T.tolist(),
               "u": flame.velocity.tolist(), "rho": flame.density.tolist(), "p": flame.P, "species": gas.species_names}
        (args.output / f"profile-{stage}.json").write_text(json.dumps(raw) + "\n")
        stages.append({"points": len(flame.grid), "speed_m_per_s": float(flame.velocity[0]),
                       "peak_temperature_K": float(max(flame.T)), "minimum_raw_mass_fraction": float(np.min(flame.Y)),
                       "maximum_raw_species_sum_error": float(max(abs(np.sum(flame.Y, axis=0) - 1))),
                       "refine_slope": criterion, "refine_curve": 2 * criterion})
    grid, temperature, velocity, raw_y = flame.grid, flame.T, flame.velocity, flame.Y.T
    prepared_y, allowance, maximum = prepare_species(raw_y, absolute_tolerance)
    # A deliberately invalid input must not be silently repaired by this
    # verification-only import path.
    invalid = raw_y.copy(); invalid[0, 0] = -1e-4
    try:
        prepare_species(invalid, absolute_tolerance)
    except ValueError:
        pass
    else:
        raise AssertionError("physical-size negative reference was repaired")
    elements = np.array([[gas.n_atoms(k, m) * gas.atomic_weight(m) / gas.molecular_weights[k]
                         for m in range(gas.n_elements)] for k in range(gas.n_species)])
    conserved = []
    for i in range(len(grid)):
        gas.set_unnormalized_mass_fractions(raw_y[i])
        gas.TP = temperature[i], flame.P
        rho = gas.density
        # Preserve rho, momentum and formation-inclusive total energy of each
        # reference point; only the documented species preparation is applied.
        conserved.append([rho, rho * velocity[i], 0, rho * (gas.int_energy_mass + .5 * velocity[i] ** 2),
                          *(rho * prepared_y[i])])
    conserved = np.array(conserved)
    np.savez(args.output / "prepared-reference.npz", grid=grid, conserved=conserved,
             raw_y=raw_y, prepared_y=prepared_y, temperature=temperature, velocity=velocity)
    fixtures = []
    for number, spacing in enumerate(args.spacings):
        edges = fixture_grid(grid, temperature, spacing)
        columns = integrate_profile(edges, grid, conserved)
        path = args.output / f"grid-{number}.fixture"
        with path.open("w") as out:
            out.write(f"CM2D_FLAME_FIXTURE 1\n{len(columns)} 2 {gas.n_species} 0.001\n" + " ".join(gas.species_names) + "\n")
            out.write(" ".join(format(v, ".17g") for v in edges) + "\n")
            for row in columns:
                out.write(" ".join(format(v, ".17g") for v in row) + "\n")
            out.write(" ".join(format(v, ".17g") for v in conserved[0]) + "\nEND\n")
        fixtures.append({"path": str(path), "sha256": sha(path), "spacing_m": spacing,
                         "columns": len(columns), "rows": 2, "height_m": .001,
                         "minimum_dx_m": float(min(np.diff(edges)))})
    report = {"scope": "independent 1D BVP reference and native conservative initial fixtures; not native flame qualification",
              "cantera": ct.__version__, "mechanism": str(args.mechanism.resolve()), "mechanism_sha256": sha(args.mechanism),
              "composition": args.composition, "temperature_K": args.temperature, "pressure_Pa": args.pressure,
              "transport": "multicomponent with Soret", "stages": stages, "fixtures": fixtures,
              "thermal_thickness_m": float((temperature[-1] - temperature[0]) / max(np.gradient(temperature, grid))),
              "projection": {"maximum_absolute_mass_fraction_change": maximum, "allowance": allowance,
                             "allowance_basis": "10 times BVP absolute species tolerance plus 64*(N+1)*machine epsilon for closure arithmetic",
                             "maximum_element_mass_fraction_change": float(np.max(abs((prepared_y - raw_y) @ elements))),
                             "rho_momentum_total_energy_preserved_per_reference_point": True,
                             "native_time_step_projection": False, "large_negative_input_rejected": True},
              "elapsed_seconds_including_reference_and_mapping": time.perf_counter() - started,
              "primary_equations": "https://www.cantera.org/3.2/reference/onedim/governing-equations.html",
              "reference_model_difference": "Cantera free flame has constant thermodynamic pressure and no longitudinal kinetic-energy/viscous-work balance; native planar flow retains both"}
    report["files_sha256"] = {p.name: sha(p) for p in args.output.iterdir() if p.is_file()}
    (args.output / "reference.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"reference_stages": stages, "projection": report["projection"], "fixtures": fixtures}, indent=2))


if __name__ == "__main__":
    main()
