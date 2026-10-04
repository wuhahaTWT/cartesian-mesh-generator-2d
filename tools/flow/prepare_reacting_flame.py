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
import subprocess
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


def fixed_grids(directory, indices):
    """Read only geometry from a hash-verified reference, never its state."""
    path = directory / "reference.json"
    reference = json.loads(path.read_text())
    if any(sha(directory / name) != digest for name, digest in reference["files_sha256"].items()):
        raise ValueError("fixed-grid reference file hash mismatch")
    indices = list(range(len(reference["fixtures"]))) if indices is None else indices
    if not indices or len(set(indices)) != len(indices) or any(i < 0 or i >= len(reference["fixtures"]) for i in indices):
        raise ValueError("invalid or repeated fixed-grid fixture index")
    grids = []
    for index in indices:
        fixture = reference["fixtures"][index]
        source = Path(fixture["path"])
        if sha(source) != fixture["sha256"]:
            raise ValueError("fixed-grid fixture hash mismatch")
        lines = source.read_text().splitlines()
        if lines[0] != "CM2D_FLAME_FIXTURE 1":
            raise ValueError("unsupported fixed-grid fixture")
        nx, ny, ns, height = lines[1].split()
        nx, ny, ns, height = int(nx), int(ny), int(ns), float(height)
        edges = np.asarray([float(v) for v in lines[3].split()])
        if (nx < 1 or ny < 1 or ns < 1 or not np.isfinite(height) or height <= 0
                or edges.shape != (nx + 1,) or not np.all(np.isfinite(edges))
                or not np.all(np.diff(edges) > 0) or len(lines) != nx + 6 or lines[-1] != "END"):
            raise ValueError("invalid fixed-grid geometry or fixture layout")
        grids.append({"edges": edges, "rows": ny, "height_m": height,
                      "spacing_m": fixture["spacing_m"], "source_fixture_index": index,
                      "source_fixture": str(source), "source_fixture_sha256": sha(source)})
    return grids, {"path": str(path), "sha256": sha(path), "indices": indices,
                   "reuse": "Only x edges, row count and height; all conservative initial states are newly integrated from the final BVP profile."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--composition", default="H2:1.6,O2:1,N2:3.76")
    parser.add_argument("--temperature", type=float, default=400)
    parser.add_argument("--pressure", type=float, default=101325)
    parser.add_argument("--width", type=float, default=.02)
    mesh = parser.add_mutually_exclusive_group()
    mesh.add_argument("--spacings", nargs="+", type=float)
    mesh.add_argument("--fixed-grid-reference", type=Path,
                      help="reuse exact geometry from an existing reference to isolate BVP input sensitivity")
    parser.add_argument("--fixed-grid-index", type=int, action="append",
                        help="select an existing fixture; repeat for several (requires --fixed-grid-reference)")
    parser.add_argument("--refine-slopes", nargs="+", type=float, default=[.06, .03, .015],
                        help="decreasing BVP refinement slopes; curve is twice each slope, original defaults retained")
    parser.add_argument("--maximum-reference-points", type=int, default=1500)
    parser.add_argument("--closure-probe", type=Path,
                        help="optional independent BVP executable enforcing sum(Y)=1 during a fixed-grid solve; all final import gates remain unchanged")
    parser.add_argument("--conservative-reference", action="store_true",
                        help="explicit total-species/total-enthalpy flux BVP discretization (requires --closure-probe); no change to native equations or import gates")
    args = parser.parse_args()
    spacings = args.spacings if args.spacings is not None else [40e-6, 20e-6, 10e-6]
    if any(not np.isfinite(x) or x <= 0 for x in [args.temperature, args.pressure, args.width, *spacings]):
        raise ValueError("all dimensional inputs must be positive and finite")
    if (any(not np.isfinite(x) or not 0 < x <= .5 for x in args.refine_slopes)
            or any(b >= a for a, b in zip(args.refine_slopes, args.refine_slopes[1:]))
            or args.maximum_reference_points < 6):
        raise ValueError("refinement slopes must decrease within (0, 0.5], with at least six reference points allowed")
    if args.fixed_grid_index is not None and args.fixed_grid_reference is None:
        raise ValueError("fixed-grid indices require a reference")
    if args.conservative_reference and args.closure_probe is None:
        raise ValueError("conservative reference requires an explicit closure probe")
    templates, fixed_source = (None, None) if args.fixed_grid_reference is None else fixed_grids(
        args.fixed_grid_reference.resolve(), args.fixed_grid_index)
    preparer_hash = sha(__file__)
    mechanism_hash = sha(args.mechanism)
    probe_hash = None if args.closure_probe is None else sha(args.closure_probe)
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
    flame.set_max_grid_points(flame.flame, args.maximum_reference_points)
    stages = []
    for stage, criterion in enumerate(args.refine_slopes):
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
    profile_file = f"profile-{len(stages) - 1}.json"
    closure_reference = None
    if args.closure_probe is not None:
        probe = args.closure_probe.resolve()
        dependent = int(np.argmax(np.mean(raw_y, axis=0)))
        guess = args.output / "closure-guess.txt"
        with guess.open("w") as stream:
            stream.write(f"CM2D_BVP_GUESS 1\n{len(grid)} {gas.n_species} {dependent} {flame.P:.17g} {flame.inlet.T:.17g} {flame.inlet.mdot:.17g} {flame.fixed_temperature:.17g}\n")
            stream.write(" ".join(gas.species_names) + "\n")
            stream.write(" ".join(format(v, ".17g") for v in flame.inlet.Y) + "\n")
            for x, t, u, y in zip(grid, temperature, velocity, raw_y):
                stream.write(" ".join(format(v, ".17g") for v in [x, t, u, *y]) + "\n")
            stream.write("END\n")
        profile_file = "profile-conservative.json" if args.conservative_reference else "profile-closed.json"
        command = [str(probe), str(args.mechanism.resolve()), str(guess.resolve()),
                   str((args.output / profile_file).resolve()), "1"]
        if args.conservative_reference:
            command.append("1")
        with (args.output / "closure.stdout").open("w") as stdout, (args.output / "closure.stderr").open("w") as stderr:
            process = subprocess.run(command, stdout=stdout, stderr=stderr)
        if process.returncode != 0 or sha(probe) != probe_hash:
            raise ValueError("independent closure BVP failed or executable changed; preserved stdout/stderr")
        closed = json.loads((args.output / profile_file).read_text())
        if (not closed["constraint_enabled"] or closed["species"] != gas.species_names
                or not np.array_equal(closed["grid"], grid) or closed["p"] != flame.P
                or closed["transport"] != "multicomponent" or not closed["soret"] or not closed["energy"]
                or closed["steady_relative_tolerance"] != 1e-9 or closed["steady_absolute_tolerance"] != absolute_tolerance
                or closed.get("conservative_flux_form", False) != args.conservative_reference):
            raise ValueError("independent closure BVP metadata differs from the requested model/grid/controls")
        closure_reference = {
            "command": command, "probe_sha256": probe_hash, "dependent_species": closed["dependent_species"],
            "maximum_original_species_residual_per_s": float(np.max(abs(np.asarray(closed["original_species_residual_per_s_interior"])))),
            "maximum_temperature_change_from_unconstrained_BVP_K": float(max(abs(np.asarray(closed["T"]) - temperature))),
            "speed_change_from_unconstrained_BVP_m_per_s": float(closed["u"][0] - velocity[0]),
            "maximum_mass_fraction_change_from_unconstrained_BVP": float(np.max(abs(np.asarray(closed["Y"]) - raw_y))),
            "method": "Re-solve all original thermochemistry and transport with one redundant interior species equation replaced by algebraic mass closure. Original species-equation residuals are separately output; no solved profile is normalized or clipped."}
        if args.conservative_reference:
            closure_reference.update(
                conservative_flux_form=True,
                maximum_selected_species_residual_per_s=float(np.max(abs(np.asarray(closed["selected_species_residual_per_s_interior"])))),
                method="Re-solve complete thermochemistry and multicomponent/Soret transport using interval total-species and total-enthalpy flux divergences on nodal dual volumes, with algebraic mass closure. Original and selected equation residuals are both retained. Pseudo-time is BVP stabilization only; no solved profile is normalized or clipped.")
        temperature, velocity, raw_y = np.asarray(closed["T"]), np.asarray(closed["u"]), np.asarray(closed["Y"])
        stages.append({**stages[-1], "fixed_grid_algebraic_mass_closure": True,
                       "speed_m_per_s": float(velocity[0]), "peak_temperature_K": float(max(temperature)),
                       "minimum_raw_mass_fraction": float(np.min(raw_y)),
                       "maximum_raw_species_sum_error": float(max(abs(raw_y.sum(axis=1) - 1)))})
        if args.conservative_reference:
            stages[-1]["conservative_flux_form"] = True
    prepared_y, allowance, maximum = prepare_species(raw_y, absolute_tolerance)
    # A deliberately invalid input must not be silently repaired by this
    # optional Python input-preparation path.
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
    if templates is None:
        templates = [{"edges": fixture_grid(grid, temperature, spacing), "spacing_m": spacing,
                      "rows": 2, "height_m": .001} for spacing in spacings]
    for number, template in enumerate(templates):
        edges = template["edges"]
        if edges[0] != grid[0] or edges[-1] != grid[-1]:
            raise ValueError("fixed native grid must cover exactly the final BVP domain; extrapolation is not permitted")
        columns = integrate_profile(edges, grid, conserved)
        path = args.output / f"grid-{number}.fixture"
        with path.open("w") as out:
            out.write(f"CM2D_FLAME_FIXTURE 1\n{len(columns)} {template['rows']} {gas.n_species} {template['height_m']:.17g}\n" + " ".join(gas.species_names) + "\n")
            out.write(" ".join(format(v, ".17g") for v in edges) + "\n")
            for row in columns:
                out.write(" ".join(format(v, ".17g") for v in row) + "\n")
            out.write(" ".join(format(v, ".17g") for v in conserved[0]) + "\nEND\n")
        fixtures.append({**{k: v for k, v in template.items() if k != "edges"},
                         "path": str(path), "sha256": sha(path), "columns": len(columns),
                         "minimum_dx_m": float(min(np.diff(edges)))})
    report = {"scope": "independent 1D BVP reference and native conservative initial fixtures; not native flame qualification",
              "cantera": ct.__version__, "mechanism": str(args.mechanism.resolve()), "mechanism_sha256": mechanism_hash,
              "composition": args.composition, "temperature_K": args.temperature, "pressure_Pa": args.pressure,
              "transport": "multicomponent with Soret", "stages": stages, "fixtures": fixtures,
              "profile_file": profile_file, "fixed_grid_reference": fixed_source,
              "closure_reference": closure_reference,
              "preparer_sha256": preparer_hash, "maximum_reference_points": args.maximum_reference_points,
              "thermal_thickness_m": float((temperature[-1] - temperature[0]) / max(np.gradient(temperature, grid))),
              "projection": {"maximum_absolute_mass_fraction_change": maximum, "allowance": allowance,
                             "allowance_basis": "10 times BVP absolute species tolerance plus 64*(N+1)*machine epsilon for closure arithmetic",
                             "maximum_element_mass_fraction_change": float(np.max(abs((prepared_y - raw_y) @ elements))),
                             "rho_momentum_total_energy_preserved_per_reference_point": True,
                             "native_time_step_projection": False, "large_negative_input_rejected": True},
              "elapsed_seconds_including_reference_and_mapping": time.perf_counter() - started,
              "primary_equations": "https://www.cantera.org/3.2/reference/onedim/governing-equations.html",
              "reference_model_difference": "Cantera free flame has constant thermodynamic pressure and no longitudinal kinetic-energy/viscous-work balance; native planar flow retains both"}
    if sha(__file__) != preparer_hash or sha(args.mechanism) != mechanism_hash:
        raise ValueError("preparer or mechanism changed during reference calculation")
    report["files_sha256"] = {p.name: sha(p) for p in args.output.iterdir() if p.is_file()}
    (args.output / "reference.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"reference_stages": stages, "projection": report["projection"], "fixtures": fixtures}, indent=2))


if __name__ == "__main__":
    main()
