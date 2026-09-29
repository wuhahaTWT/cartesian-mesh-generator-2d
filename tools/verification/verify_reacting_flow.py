#!/usr/bin/env python3
"""Independently read native reacting fields, audit budgets and render real cells.

This checks geometry, conservative-state thermodynamics and numerical budgets.
The Python Cantera EOS shares a backend with the native run; no experimental
accuracy or resolved-flame qualification is inferred from this audit.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import cantera as ct
import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compare", type=Path, help="another audited run on the same mesh, for a descriptive time-step comparison")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    path = args.case / "field.json"
    source_files = ["field.json", "steps.jsonl", "resolved-mechanism.yaml", "accepted.checkpoint"]
    source_hashes = {name: sha(args.case / name) for name in source_files}
    r = json.loads(path.read_text())
    assert r["complete"] and not r["qualifiedFlame"]
    gas = ct.Solution(str(args.case / "resolved-mechanism.yaml"))
    assert gas.species_names == r["species"]
    atoms = np.array([[gas.n_atoms(k, m) for m in range(gas.n_elements)] for k in range(gas.n_species)])
    element_weights = atoms * gas.atomic_weights[None, :] / gas.molecular_weights[:, None]
    assert np.array_equal(atoms.ravel(), r["atomCounts"])
    assert np.allclose(gas.molecular_weights, r["molecularWeights"], rtol=0, atol=1e-12)
    cells, faces = r["mesh"]["cells"], r["mesh"]["faces"]
    area = np.array([c["area"] for c in cells])
    assert np.all(area > 0)
    polygons = []
    incidence = [[] for _ in faces]
    geometry_error = 0.0
    for i, c in enumerate(cells):
        assert len(c["faces"]) >= 3 and len(set(c["faces"])) == len(c["faces"])
        vertices = []
        closure = np.zeros(2)
        for f in c["faces"]:
            assert isinstance(f, int) and 0 <= f < len(faces)
            incidence[f].append(i)
            face = faces[f]
            sign = 1 if face["owner"] == i else -1
            assert sign == 1 or face["neighbour"] == i
            s = np.array(face["S"]) * sign
            closure += s
            vertices.append(np.array(face["centre"]) + np.array([s[1], -s[0]]) / 2)
        points = np.array(vertices)
        measured = .5 * np.sum(points[:, 0] * np.roll(points[:, 1], -1) - points[:, 1] * np.roll(points[:, 0], -1))
        geometry_error = max(geometry_error, abs(measured / c["area"] - 1))
        assert np.linalg.norm(closure) < 1e-12 * np.sqrt(c["area"])
        polygons.append(points)
    for f, uses in zip(faces, incidence):
        expected = [f["owner"]] + ([] if f["neighbour"] is None else [f["neighbour"]])
        assert len(set(expected)) == len(expected) and sorted(uses) == sorted(expected)
    assert geometry_error < 1e-10
    time_previous, steps_previous = -1.0, -1
    thermo_error, closure_error, velocity_error = 0.0, 0.0, 0.0
    integrals, thermo_scale, sound_scale = [], 0.0, 0.0
    for frame in r["frames"]:
        assert frame["time"] > time_previous and frame["steps"] > steps_previous
        time_previous, steps_previous = frame["time"], frame["steps"]
        u = np.array(frame["U"])
        assert u.shape == (len(cells), gas.n_species + 4) and np.all(np.isfinite(u))
        assert np.all(u[:, 0] > 0) and np.all(u[:, 4:] >= 0)
        y = u[:, 4:] / u[:, 0, None]
        closure_error = max(closure_error, float(np.max(abs(y.sum(axis=1) - 1))))
        for i in range(len(cells)):
            gas.set_unnormalized_mass_fractions(y[i])
            gas.TD = frame["temperature"][i], u[i, 0]
            assert max(sp.thermo.min_temp for sp in gas.species()) <= gas.T <= min(sp.thermo.max_temp for sp in gas.species())
            kinetic = (u[i, 1] ** 2 + u[i, 2] ** 2) / (2 * u[i, 0])
            energy = u[i, 0] * gas.int_energy_mass + kinetic
            scale = max(abs(u[i, 3]), u[i, 0] * gas.cv_mass * gas.T, kinetic)
            thermo_error = max(thermo_error, abs(energy - u[i, 3]) / scale,
                               abs(gas.P - frame["pressure"][i]) / gas.P)
            velocity_error = max(velocity_error, float(np.max(abs(u[i, 1:3] / u[i, 0] - frame["velocity"][i]))))
            if frame["steps"] == 0:
                thermo_scale += area[i] * scale
                sound_scale = max(sound_scale, gas.sound_speed)
        integrals.append(np.sum(u.astype(np.longdouble) * area[:, None], axis=0))
    assert r["frames"][-1]["time"] == r["targetTime"]
    assert closure_error < 1e-10 and thermo_error < 2e-9 and velocity_error < 1e-10
    before, after = integrals[0], integrals[-1]
    boundary = np.array(r["boundaryImpulse"])
    # Closed static adiabatic walls: exactly no mass/species or total energy
    # crosses a boundary. Momentum may be exchanged through wall traction.
    assert boundary[0] == 0 and boundary[3] == 0 and np.all(boundary[4:] == 0)
    mass_drift = float(abs(after[0] - before[0]) / before[0])
    energy_drift = float(abs(after[3] - before[3]) / thermo_scale)
    element_drift = np.asarray((after[4:] - before[4:]) @ element_weights / before[0], dtype=float)
    momentum_error = float(np.max(abs(after[1:3] - before[1:3] + boundary[1:3])) / (before[0] * sound_scale))
    # Normalized numerical conservation gates, inherited from the chemical
    # substep allowance. They are NOT tolerances on flame temperature or speed.
    assert max(mass_drift, energy_drift, max(abs(element_drift)), momentum_error) < 1e-8
    steps = [json.loads(line) for line in (args.case / "steps.jsonl").read_text().splitlines()]
    assert len(steps) == r["frames"][-1]["steps"]
    clock = 0.0
    for step in steps:
        assert step["accepted"] and step["dt"] > 0 and step["previousTime"] == clock
        assert step["courant"] <= r.get("controls", {}).get("courant", .35) * (1 + 1e-12)
        if "rejectedReasons" in step:
            assert len(step["rejectedReasons"]) == step["rejections"] and step["failure"] == ""
        clock += step["dt"]
    assert clock == r["targetTime"]
    assert sum(s["rejections"] for s in steps) == r["rejections"]
    fuel = gas.species_index("H2")
    fuel_consumed = float(1 - after[4 + fuel] / before[4 + fuel])
    first, last = r["frames"][0], r["frames"][-1]
    peak_speed = float(np.max(np.linalg.norm(last["velocity"], axis=1)))
    assert fuel_consumed > 0 and max(last["temperature"]) > max(first["temperature"]) and peak_speed > 0
    report = {"passed": True, "qualified_flame": False,
              "scope": "short-time detailed-chemistry/flow coupling on a coarse native 2D verification mesh",
              "cells": len(cells), "faces": len(faces), "time": clock, "steps": len(steps),
              "initial_temperature_range": [min(first["temperature"]), max(first["temperature"])],
              "final_temperature_range": [min(last["temperature"]), max(last["temperature"])],
              "final_pressure_range": [min(last["pressure"]), max(last["pressure"])],
              "maximum_speed_m_per_s": peak_speed, "hydrogen_mass_consumed_fraction": fuel_consumed,
              "mass_drift_relative": mass_drift, "total_energy_drift_relative": energy_drift,
              "element_drift_per_initial_mass": element_drift.tolist(), "momentum_balance_scaled": momentum_error,
              "energy_scale_J_per_m": float(thermo_scale), "initial_maximum_sound_speed_m_per_s": float(sound_scale),
              "geometric_area_relative_error": geometry_error, "thermo_interface_relative_error": thermo_error,
              "species_sum_absolute_error": closure_error, "velocity_roundtrip_absolute_error": velocity_error,
              "cantera_reference": ct.__version__, "raw_field_sha256": sha(path),
              "source_files": source_hashes}
    if args.compare:
        other_path = args.compare / "field.json"
        other_hash = sha(other_path)
        other = json.loads(other_path.read_text())
        assert other["complete"] and r["targetTime"] == other["targetTime"]
        assert r["mesh"] == other["mesh"] and r["species"] == other["species"]
        assert r["frames"][0] == other["frames"][0]
        assert source_hashes["resolved-mechanism.yaml"] == sha(args.compare / "resolved-mechanism.yaml")
        assert {k:v for k,v in r["controls"].items() if k != "courant"} == {k:v for k,v in other["controls"].items() if k != "courant"}
        reference = other["frames"][-1]
        area_l2 = lambda delta: float(np.sqrt(np.sum(area * delta ** 2) / np.sum(area)))
        temperature_delta = np.array(last["temperature"]) - reference["temperature"]
        pressure_delta = np.array(last["pressure"]) - reference["pressure"]
        velocity_delta = np.linalg.norm(np.array(last["velocity"]) - reference["velocity"], axis=1)
        u, v = np.array(last["U"]), np.array(reference["U"])
        species_delta = abs(u[:, 4:] / u[:, 0, None] - v[:, 4:] / v[:, 0, None])
        report["time_step_comparison"] = {
            "scope": "two time-step settings on one coarse mesh; not a convergence-order or spatial-accuracy qualification",
            "case_courant": r["controls"]["courant"], "other_courant": other["controls"]["courant"],
            "other_raw_sha256": other_hash, "other_steps": reference["steps"],
            "temperature_max_difference_K": float(np.max(abs(temperature_delta))),
            "temperature_area_l2_difference_K": area_l2(temperature_delta),
            "pressure_max_difference_Pa": float(np.max(abs(pressure_delta))),
            "pressure_area_l2_difference_Pa": area_l2(pressure_delta),
            "velocity_max_difference_m_per_s": float(np.max(velocity_delta)),
            "mass_fraction_max_absolute_difference": float(np.max(species_delta)),
            "mass_fraction_area_mean_l1_difference": float(np.sum(area[:, None] * species_delta) / np.sum(area))}
        assert sha(other_path) == other_hash
    assert {name: sha(args.case / name) for name in source_files} == source_hashes
    (args.output / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / "build" / "matplotlib-cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    fig, axes = plt.subplots(2, 2, figsize=(11, 10), layout="constrained")
    common = Normalize(min(first["temperature"]), max(last["temperature"]))
    centres = np.array([c["centre"] for c in cells]) * 1000
    panels = [(first["temperature"], "Initial temperature", "K", "inferno", common),
              (last["temperature"], f"Temperature at {clock * 1e6:.0f} microseconds", "K", "inferno", common),
              (np.array(last["pressure"]) / 1000, "Pressure", "kPa", "viridis", None),
              (np.linalg.norm(last["velocity"], axis=1), "Velocity magnitude and direction", "m/s", "cividis", None)]
    for ax, (values, title, unit, cmap, norm) in zip(axes.flat, panels):
        collection = PolyCollection([p * 1000 for p in polygons], array=np.array(values), cmap=cmap, norm=norm,
                                    edgecolors=(1, 1, 1, .25), linewidths=.5)
        ax.add_collection(collection); ax.autoscale_view(); ax.set_aspect("equal")
        ax.set(xlabel="x (mm)", ylabel="y (mm)", title=title)
        fig.colorbar(collection, ax=ax, label=unit, shrink=.83)
    v = np.array(last["velocity"])
    axes[1, 1].quiver(centres[:, 0], centres[:, 1], v[:, 0], v[:, 1], angles="xy", scale_units="xy", scale=max(peak_speed / 4, 1), color="white", width=.004)
    fig.suptitle(f"Native detailed H2/air reacting flow\n{len(cells)} cells | {len(steps)} accepted steps | H2 consumed: {100 * fuel_consumed:.3f}%", fontsize=14)
    fig.supxlabel("Adiabatic no-slip enclosure. Early ignition / coupling test; no resolved-flame or experimental qualification.", fontsize=9)
    fig.savefig(args.output / "reacting-flow.png", dpi=160)
    plt.close(fig)
    print(json.dumps({k:report[k] for k in ["passed", "steps", "final_temperature_range", "maximum_speed_m_per_s", "hydrogen_mass_consumed_fraction", "mass_drift_relative", "total_energy_drift_relative"]}, indent=2))


if __name__ == "__main__":
    main()
