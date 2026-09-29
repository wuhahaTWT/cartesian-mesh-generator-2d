#!/usr/bin/env python3
"""Run native planar-flame fixtures and independently audit their actual fields.

Steady reference profiles are initial data, not accepted native checkpoints.
Report short-time drift and spatial residuals separately from flame qualification.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import traceback

import cantera as ct
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(path, reference, fixture, audit_output=None):
    r = json.loads((path / "field.json").read_text())
    gas = ct.Solution(str(path / "resolved-mechanism.yaml"))
    area = np.array(r["areas"])
    faces, nx, ny = r["faces"], r["nx"], r["ny"]
    assert len(area) == nx * ny and np.all(area > 0)
    assert len(faces) == nx * (ny + 1) + (nx + 1) * ny
    expected_area = np.tile(np.diff(r["xEdges"]) * r["height"] / ny, ny)
    assert max(abs(area / expected_area - 1)) < 1e-10
    atoms = np.array([[gas.n_atoms(k, m) for m in range(gas.n_elements)] for k in range(gas.n_species)])
    elements = atoms * gas.atomic_weights[None, :] / gas.molecular_weights[:, None]
    thermo_error, chemical_error, flux_error, closure_error = 0., 0., 0., 0.
    integrals, energy_scale, sound_scale = [], 0., 0.
    state_metrics = []
    for label in ["initial", "final"]:
        state = r[label]
        u = np.array(state["U"])
        assert u.shape == (nx * ny, gas.n_species + 4) and np.all(np.isfinite(u))
        assert np.all(u[:, 0] > 0) and np.all(u[:, 4:] >= 0)
        y = u[:, 4:] / u[:, 0, None]
        closure_error = max(closure_error, float(np.max(abs(y.sum(axis=1) - 1))))
        transport = np.array(state["transportDerivative"])
        chemistry = np.array(state["chemistryDerivative"])
        derivative = np.array(state["derivative"])
        assert np.array_equal(derivative, transport + chemistry)
        fluxes = np.array(state["faceFlux"])
        balance, activity, normals = np.zeros_like(u), np.zeros_like(u), np.zeros((len(area), 2))
        boundary = np.zeros(u.shape[1])
        for face, flux in zip(faces, fluxes):
            a, b = face["owner"], face["neighbour"]
            assert 0 <= a < len(area) and (b is None or 0 <= b < len(area) and b != a)
            balance[a] += flux; activity[a] += abs(flux); normals[a] += face["S"]
            if b is not None:
                balance[b] -= flux; activity[b] += abs(flux); normals[b] -= face["S"]
            else:
                boundary += flux
        assert np.max(np.linalg.norm(normals, axis=1) / np.sqrt(area)) < 1e-10
        flux_error = max(flux_error, float(np.max(abs(balance + area[:, None] * transport)
            / np.maximum(activity + abs(area[:, None] * transport), 1e-300))))
        assert np.array_equal(boundary, state["boundaryFlux"])
        assert np.max(abs(fluxes[:, 4:].sum(axis=1) - fluxes[:, 0])
                      / np.maximum(abs(fluxes[:, 4:]).sum(axis=1) + abs(fluxes[:, 0]), 1e-300)) < 1e-10
        tdot = []
        for i, q in enumerate(u):
            gas.set_unnormalized_mass_fractions(y[i]); gas.TD = state["temperature"][i], q[0]
            kinetic = (q[1] ** 2 + q[2] ** 2) / (2 * q[0])
            escale = max(abs(q[3]), q[0] * gas.cv_mass * gas.T, kinetic)
            thermo_error = max(thermo_error, abs(q[0] * gas.int_energy_mass + kinetic - q[3]) / escale,
                               abs(gas.P - state["pressure"][i]) / gas.P)
            expected = gas.net_production_rates * gas.molecular_weights
            gross = (gas.creation_rates + gas.destruction_rates) * gas.molecular_weights
            chemical_error = max(chemical_error, float(np.max(abs(chemistry[i, 4:] - expected) / np.maximum(gross, 1e-300))))
            assert np.all(chemistry[i, :4] == 0)
            du = derivative[i]
            de = du[3] - q[1] / q[0] * du[1] - q[2] / q[0] * du[2] + kinetic / q[0] * du[0]
            species_energy = gas.partial_molar_int_energies / gas.molecular_weights
            tdot.append((de - species_energy @ du[4:]) / (q[0] * gas.cv_mass))
            if label == "initial":
                energy_scale += area[i] * escale; sound_scale = max(sound_scale, gas.sound_speed)
        integrals.append(np.sum(u.astype(np.longdouble) * area[:, None], axis=0))
        fuel = gas.species_index("H2")
        # This is source evaluation on the reference-initialized profile,
        # not an independently predicted flame speed from a steady native run.
        consumption = -float(np.sum(area * chemistry[:, 4 + fuel]) / r["height"])
        scale = np.sum(area[:, None] * abs(chemistry[:, 4:]))
        state_metrics.append({"species_residual_L1_over_chemical_activity": float(np.sum(area[:, None] * abs(derivative[:, 4:])) / scale),
                              "maximum_temperature_derivative_K_per_s": float(np.max(np.abs(tdot))),
                              "hydrogen_consumption_kg_per_m2_s": consumption,
                              "pressure_max_relative_to_reference": float(np.max(abs(np.array(state["pressure"]) / reference["pressure_Pa"] - 1)))})
    assert closure_error < 1e-10 and thermo_error < 2e-9 and chemical_error < 2e-8 and flux_error < 2e-11
    # Verify the native initial state is the supplied conservative average;
    # the solver must not silently repair imported Y. Algebraic N-1 stage
    # closure is reported explicitly and remains in the physical budgets.
    lines = Path(fixture["path"]).read_text().splitlines()
    columns = np.array([[float(v) for v in line.split()] for line in lines[4:4 + nx]])
    assert np.array_equal(np.tile(columns, (ny, 1)), r["initial"]["U"])
    before, after = integrals
    impulse = np.array(r["boundaryImpulse"])
    physical_budget = after - before + impulse
    mass_error = float(abs(physical_budget[0]) / before[0])
    energy_error = float(abs(physical_budget[3]) / energy_scale)
    momentum_error = float(np.max(abs(physical_budget[1:3])) / (before[0] * sound_scale))
    element_error = float(np.max(abs(physical_budget[4:] @ elements / before[0])))
    assert max(mass_error, energy_error, momentum_error, element_error) < 1e-8
    first, last = r["initial"], r["final"]
    steps = [json.loads(line) for line in (path / "steps.jsonl").read_text().splitlines()]
    clock, accepted = 0., 0
    for step in steps:
        assert step["time"] == clock
        if step["accepted"]:
            assert step["dt"] > 0 and step["courant"] <= .35 * (1 + 1e-12)
            clock += step["dt"]; accepted += 1
        else:
            assert step is steps[-1]
    assert clock == last["time"] and accepted == last["steps"]
    closure_report = None
    if "massClosure" in r:
        closure = r["massClosure"]
        accepted_log = [s for s in steps if s["accepted"]]
        assert all(0 <= s["maximumMassClosureFraction"] <= 64 * (gas.n_species + 1) * np.finfo(float).eps
                   for s in accepted_log)
        assert max((s["maximumMassClosureFraction"] for s in accepted_log), default=0.) == closure["maximumFraction"]
        # Python 3.12+ sum() uses compensated accumulation; the native log
        # accumulator uses ordered IEEE additions. Reproduce those additions
        # for the exact serialization check, without relaxing its threshold.
        absolute = 0.
        for s in accepted_log:
            absolute += s["absoluteMassClosureIntegral"]
        assert absolute == closure["absoluteIntegral"]
        for kind in ["transport", "chemistry"]:
            total = np.zeros(gas.n_species + 4)
            for s in accepted_log:
                change = np.array(s[kind + "MassClosureChange"])
                assert change.shape == total.shape and np.all(np.isfinite(change)) and np.all(change[:4] == 0)
                total += change
            assert np.array_equal(total, closure[kind])
        # Unsubtracted physical budgets above still include all closure
        # changes. Separately expose the total absolute roundoff operation.
        assert closure["absoluteIntegral"] >= 0
        closure_report = {"maximum_stage_fraction": closure["maximumFraction"],
                          "absolute_integral_per_initial_mass": float(closure["absoluteIntegral"] / before[0]),
                          "net_element_change_per_initial_mass": ((np.array(closure["transport"])[4:]
                              + np.array(closure["chemistry"])[4:]) @ elements / float(before[0])).tolist()}
    report = {"numerical_checks_passed": True, "physical_endpoint_reached": r["complete"] and accepted > 0,
              "qualified_flame": False, "cells": nx * ny, "spacing_m": fixture["spacing_m"],
              "time_s": clock, "accepted_steps": accepted, "rejections": r["rejections"], "failure": r["failure"],
              "thermo_relative_error": thermo_error, "source_error_over_gross_activity": chemical_error,
              "face_assembly_relative_error": flux_error, "species_sum_absolute_error": closure_error,
              "mass_budget_per_initial_mass": mass_error, "energy_budget_scaled": energy_error,
              "momentum_budget_scaled": momentum_error, "element_budget_per_initial_mass": element_error,
              "initial_residual": state_metrics[0], "final_residual": state_metrics[1],
              "maximum_temperature_drift_K": float(np.max(abs(np.array(last["temperature"]) - first["temperature"]))),
              "field_sha256": sha(path / "field.json")}
    if closure_report is not None:
        report["mass_closure_diagnostics"] = closure_report
    (Path(audit_output or path) / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    return report, r


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--duration", type=float, default=1e-6)
    parser.add_argument("--jobs", type=int, default=3, help="independent native cases in parallel; 1 for sequential execution")
    parser.add_argument("--audit-native-output", type=Path,
                        help="audit an already completed verifier run into a NEW output directory, without running the solver")
    args = parser.parse_args()
    verifier_hash = sha(Path(__file__))
    assert np.isfinite(args.duration) and args.duration > 0
    assert 1 <= args.jobs <= 3
    args.output.mkdir(parents=True, exist_ok=False)
    reference = json.loads((args.reference / "reference.json").read_text())
    for name, digest in reference["files_sha256"].items():
        assert sha(args.reference / name) == digest
    assert sha(reference["mechanism"]) == reference["mechanism_sha256"]
    probe_hash = sha(args.probe)
    previous = None
    if args.audit_native_output:
        previous = json.loads((args.audit_native_output / "report.json").read_text())
        assert previous["probe_sha256"] == probe_hash and previous["duration_s"] == args.duration
        assert previous["reference_report_sha256"] == sha(args.reference / "reference.json")
        assert len(previous["cases"]) == len(reference["fixtures"])
    def run_case(item):
        number, fixture = item
        assert sha(fixture["path"]) == fixture["sha256"] and sha(args.probe) == probe_hash
        audit_path = args.output / f"native-{number}"
        path = (args.audit_native_output or args.output) / f"native-{number}"
        if previous is None:
            start = time.perf_counter()
            process = subprocess.run([str(args.probe.resolve()), reference["mechanism"], fixture["path"],
                                      str(path), format(args.duration, ".17g")], capture_output=True, text=True)
            elapsed, returncode, failure = time.perf_counter() - start, process.returncode, process.stderr
            (args.output / f"native-{number}.stdout").write_text(process.stdout)
            (args.output / f"native-{number}.stderr").write_text(process.stderr)
        else:
            saved = previous["cases"][number]
            assert Path(saved["directory"]).resolve() == path.resolve()
            elapsed, returncode, failure = saved["elapsed_seconds_including_all_IO"], saved["returncode"], saved.get("failure", "")
            audit_path.mkdir()
        raw_hashes = {str(p): sha(p) for p in (path.iterdir() if path.is_dir() else [])
                      if p.is_file() and p.name not in {"audit.json", "audit-error.txt"}}
        case = {"returncode": returncode, "elapsed_seconds_including_all_IO": elapsed, "directory": str(path), "raw_sha256": raw_hashes}
        field = None
        if (path / "field.json").exists():
            try:
                checked, field = audit(path, reference, fixture, audit_path); case.update(checked)
            except Exception as error:
                (audit_path / "audit-error.txt").write_text(traceback.format_exc())
                case.update(numerical_checks_passed=False, physical_endpoint_reached=False,
                            failure=f"Independent audit failed: {type(error).__name__}: {error}; see audit-error.txt")
        else:
            case.update(numerical_checks_passed=False, physical_endpoint_reached=False, failure=failure)
        assert all(sha(p) == digest for p, digest in raw_hashes.items())
        print(json.dumps(case), flush=True)
        return case, field
    wall_start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        results = list(executor.map(run_case, enumerate(reference["fixtures"])))
    cases = [case for case, _ in results]
    fields = [field for _, field in results if field is not None]
    assert sha(args.probe) == probe_hash
    assert sha(Path(__file__)) == verifier_hash
    report = {"qualified_flame": False, "duration_s": args.duration, "probe_sha256": probe_hash,
              "execution_mode": "audit_existing" if previous is not None else "run_and_audit",
              "verifier_sha256": verifier_hash,
              "reference_report_sha256": sha(args.reference / "reference.json"), "cases": cases,
              "parallel_jobs": args.jobs, "case_wall_seconds_including_audits": time.perf_counter() - wall_start,
              "all_numerical_checks_passed": all(c["numerical_checks_passed"] for c in cases),
              "all_endpoints_reached": all(c["physical_endpoint_reached"] for c in cases),
              "all_native_runs_succeeded": all(c["returncode"] == 0 for c in cases),
              "limits": ["reference and native share chemistry/transport backend", "constant-pressure BVP versus full compressible FV",
                         "short evolution from a reference initial condition does not predict or qualify flame speed",
                         "reference roundoff preparation recorded separately; no native species clipping", "no experiment or actual combustor qualification"]}
    if previous is not None:
        report["audit_source_report"] = {"path": str(args.audit_native_output / "report.json"),
                                         "sha256": sha(args.audit_native_output / "report.json")}
        report["audit_wall_seconds"] = report.pop("case_wall_seconds_including_audits")
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    if len(fields) == len(cases):
        os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / "build" / "matplotlib-cache"))
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        profile = json.loads((args.reference / "profile-2.json").read_text())
        grid, temperature = np.array(profile["grid"]), np.array(profile["T"])
        centre = grid[np.argmax(np.gradient(temperature, grid))]
        fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
        axes[0, 0].plot((grid-centre)*1000, temperature, "k--", label="Cantera steady reference")
        fuel = profile["species"].index("H2")
        axes[0, 1].plot((grid-centre)*1000, np.array(profile["Y"])[:, fuel], "k--", label="Reference H2")
        for case, field in zip(cases, fields):
            nx = field["nx"]; x = (np.array(field["xEdges"])[1:] + np.array(field["xEdges"])[:-1])/2
            u = np.array(field["final"]["U"])[:nx]; label = f"FV {case['spacing_m']*1e6:.0f} um"
            axes[0, 0].plot((x-centre)*1000, field["final"]["temperature"][:nx], label=label)
            axes[0, 1].plot((x-centre)*1000, u[:, 4+fuel]/u[:, 0], label=label)
            axes[1, 0].plot((x-centre)*1000, np.array(field["final"]["temperature"][:nx])-field["initial"]["temperature"][:nx], label=label)
        h = [c["spacing_m"]*1e6 for c in cases]
        axes[1, 1].loglog(h, [c["initial_residual"]["species_residual_L1_over_chemical_activity"] for c in cases], "o-")
        for ax, title, unit in [(axes[0, 0], "Temperature", "K"), (axes[0, 1], "Hydrogen mass fraction", "Y(H2)"),
                                (axes[1, 0], "Temperature change during native evolution", "K")]:
            ax.set(title=title, xlabel="x relative to reference flame (mm)", ylabel=unit, xlim=(-.6, 2)); ax.legend(fontsize=8); ax.grid(alpha=.2)
        axes[1, 1].set(title="Initial species equation imbalance", xlabel="Core spacing (um)", ylabel="L1 residual / L1 reaction activity")
        axes[1, 1].grid(alpha=.2)
        fig.suptitle(f"Detailed H2/air planar flame | native evolution requested: {args.duration*1e6:g} us")
        fig.supxlabel("Reference-initialized short-time verification; no independently predicted flame-speed or experimental qualification.", fontsize=8)
        fig.savefig(args.output / "reacting-flame.png", dpi=160); plt.close(fig)
    if not report["all_numerical_checks_passed"] or not report["all_endpoints_reached"] or not report["all_native_runs_succeeded"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
