#!/usr/bin/env python3
"""Run native planar-flame fixtures and independently audit their actual fields.

Steady reference profiles are initial data, not accepted native checkpoints.
Report short-time drift and spatial residuals separately from flame qualification.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
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


def read_checkpoint(path, cells, variables):
    """Read the serialized state independently; native loading checks binding."""
    def scalar(token):
        value = float(token)
        assert np.isfinite(value) and (value != 0 or Decimal(token) == 0)
        return value

    with Path(path).open("rb") as stream:
        assert stream.readline() == b"CM2D_REACTING_CHECKPOINT 1\n"
        size_text = stream.readline().strip()
        assert size_text.isdigit()
        size = int(size_text)
        assert 0 < size <= Path(path).stat().st_size
        binding = stream.read(size)
        assert len(binding) == size
        clock = stream.readline().decode("ascii").split()
        assert len(clock) == 3 and clock[0] == "STATE" and clock[2].isdigit()
        time_s, count = scalar(clock[1]), int(clock[2])
        assert np.isfinite(time_s) and time_s > 0 and count > 0
        rows = [[scalar(v) for v in stream.readline().decode("ascii").split()] for _ in range(cells)]
        u = np.array(rows)
        assert u.shape == (cells, variables) and np.all(np.isfinite(u))
        assert stream.read().split() == [b"END"]
    return {"time": time_s, "steps": count, "U": rows}


def audit_fields(r, mechanism, steps, reference, fixture, checkpoint=None):
    """Audit serialized native data; callers identify its raw file or sample."""
    gas = ct.Solution(str(mechanism))
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
    # Verify exact import, including the clock. A restart starts a new BDF
    # session and new budget interval, while retaining the accepted physical
    # time and total step count. Native checkpoint binding is also enforced.
    if r.get("restart") is None:
        assert checkpoint is None and r["initial"]["time"] == 0 and r["initial"]["steps"] == 0
        lines = Path(fixture["path"]).read_text().splitlines()
        columns = np.array([[float(v) for v in line.split()] for line in lines[4:4 + nx]])
        assert np.array_equal(np.tile(columns, (ny, 1)), r["initial"]["U"])
    else:
        assert checkpoint is not None and r["restart"]["restoresIntegratorHistory"] is False
        assert r["restart"]["budgetOrigin"] == "initial-state"
        assert all(r["initial"][key] == checkpoint[key] for key in ("time", "steps", "U"))
    before, after = integrals
    impulse = np.array(r["boundaryImpulse"])
    physical_budget = after - before + impulse
    mass_error = float(abs(physical_budget[0]) / before[0])
    energy_error = float(abs(physical_budget[3]) / energy_scale)
    momentum_error = float(np.max(abs(physical_budget[1:3])) / (before[0] * sound_scale))
    element_error = float(np.max(abs(physical_budget[4:] @ elements / before[0])))
    assert max(mass_error, energy_error, momentum_error, element_error) < 1e-8
    first, last = r["initial"], r["final"]
    end_time = r.get("endTime", r["duration"])
    assert np.isfinite(r["duration"]) and r["duration"] >= 0 and np.isfinite(end_time)
    assert end_time == first["time"] + r["duration"]
    assert first["time"] <= last["time"] <= end_time
    assert r["complete"] == (last["time"] == end_time)
    implicit = r.get("integration", {}).get("method") == "CVODES-BDF"
    clock, accepted = first["time"], 0
    for step in steps:
        assert step["time"] == clock
        if step["accepted"]:
            assert step["dt"] > 0
            if implicit:
                assert step["internalSteps"] == accepted + 1
                assert step["dt"] <= r["integration"]["maximumStep"] * (1 + 1e-12)
            else:
                assert step["courant"] <= .35 * (1 + 1e-12)
            clock += step["dt"]; accepted += 1
        else:
            assert step is steps[-1]
    assert clock == last["time"] and first["steps"] + accepted == last["steps"]
    if implicit and "acceptedByBdfOrder" in r["integration"]:
        orders = [0] * 6
        for step in steps:
            if step["accepted"]:
                assert 1 <= step["bdfOrder"] <= r["integration"]["maximumOrder"]
                orders[step["bdfOrder"]] += 1
        assert orders == r["integration"]["acceptedByBdfOrder"] and sum(orders) == accepted
    if implicit and r["integration"].get("errorControlCorrection") == "fullPredictorCorrection":
        errors = [s["localErrorNorm"] for s in steps if s["accepted"]]
        # Backend's weighted local-error acceptance invariant only. It does
        # not establish an independently measured physical time-step error.
        allowance = 1 + 64 * np.finfo(float).eps * nx * ny * (gas.n_species + 3)
        assert all(np.isfinite(e) and 0 <= e <= allowance for e in errors)
        assert max(errors, default=0.) == r["integration"]["maximumLocalErrorNorm"]
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
              "initial_time_s": first["time"], "initial_accepted_steps": first["steps"],
              "total_accepted_steps": last["steps"], "restarted": checkpoint is not None,
              "thermo_relative_error": thermo_error, "source_error_over_gross_activity": chemical_error,
              "face_assembly_relative_error": flux_error, "species_sum_absolute_error": closure_error,
              "mass_budget_per_initial_mass": mass_error, "energy_budget_scaled": energy_error,
              "momentum_budget_scaled": momentum_error, "element_budget_per_initial_mass": element_error,
              "initial_residual": state_metrics[0], "final_residual": state_metrics[1],
              "maximum_temperature_drift_K": float(np.max(abs(np.array(last["temperature"]) - first["temperature"])))}
    if closure_report is not None:
        report["mass_closure_diagnostics"] = closure_report
    if implicit:
        report["integration"] = r["integration"]
        constraint = np.array(r["constraintChange"])
        assert constraint.shape == impulse.shape and np.all(np.isfinite(constraint)) and np.all(constraint[:4] == 0)
        # A BDF quadrature audit is separate from the unsubtracted physical
        # element/energy budgets above. Same 1e-8 mass-normalized gate.
        equation_defect = float(np.max(abs(physical_budget[4:] - np.array(r["chemistryChange"])[4:] - constraint[4:])) / before[0])
        assert equation_defect < 1e-8
        report["bdf_species_quadrature_defect_per_mass"] = equation_defect
        report["constraint_change_per_initial_mass"] = (constraint / float(before[0])).tolist()
    return report


def audit_sample_counters(field, sample, accepted_rhs, last_sample):
    """Separate accepted-step statistics from a stopped trial's final snapshot.

    A stopped backend can evaluate more RHS/Jacobian candidates without
    accepting a state. Its final sample may contain those terminal counters,
    or the last regular sample may predate them. Neither is a new time step.
    """
    saved, terminal = sample["integration"], field["integration"]
    rhs = saved["rhsCalls"]
    assert type(rhs) is int and type(accepted_rhs) is int and rhs >= 0 and accepted_rhs >= 0
    stopped_final = (last_sample and not field["complete"] and bool(field["failure"])
                     and sample["state"] == field["final"])
    if rhs != accepted_rhs:
        assert stopped_final and rhs > accepted_rhs and saved == terminal
    if not last_sample or field["integration"].get("samplingFailure"):
        return
    if not stopped_final:
        assert saved == terminal
        return
    assert type(terminal["rhsCalls"]) is int and terminal["rhsCalls"] >= rhs
    monotone_counts = {"rhsCalls", "rejectedRhsCalls", "errorTestFailures", "linearSetups", "jacobianEvaluations",
                       "nonlinearIterations", "nonlinearConvergenceFailures", "dampedNewtonUpdates",
                       "continuedDampedNewtonUpdates", "reflectedNewtonUpdates", "reflectionAttempts"}
    assert saved.keys() == terminal.keys()
    for key, before in saved.items():
        after = terminal[key]
        if key in monotone_counts:
            assert type(before) is int and type(after) is int and 0 <= before <= after
        elif key == "maximumReflectedSpeciesScaledChange":
            assert np.isfinite(before) and np.isfinite(after) and 0 <= before <= after
        elif key == "minimumNewtonFraction":
            assert np.isfinite(before) and np.isfinite(after) and 0 <= after <= before <= 1
        elif key == "lastDampedTrialFailure":
            assert isinstance(before, str) and isinstance(after, str)
        elif key == "canceled":
            assert type(before) is bool and type(after) is bool and (not before or after)
        else:
            # Controls, sample counts, accepted-order counts and accepted LTE
            # statistics cannot change while retaining this accepted state.
            assert before == after


def audit_samples(path, field, steps, reference, fixture, audit_output, checkpoint=None):
    every = field.get("integration", {}).get("sampleEveryAcceptedSteps", 0)
    if not every:
        assert not (path / "samples.jsonl").exists()
        return None
    samples = [json.loads(line) for line in (path / "samples.jsonl").read_text().splitlines()]
    accepted = [step for step in steps if step["accepted"]]
    first_step = field["initial"]["steps"]
    expected = [first_step] + [first_step + i for i in range(1, len(accepted) + 1) if i % every == 0]
    if expected[-1] != first_step + len(accepted):
        expected.append(first_step + len(accepted))
    sampling_failure = field["integration"].get("samplingFailure", "")
    assert samples
    assert [sample["state"]["steps"] for sample in samples] == (expected[:len(samples)] if sampling_failure else expected)
    assert len(samples) == field["integration"]["sampleResidualEvaluations"]
    assert samples[0]["kind"] == "initial" and samples[0]["state"] == field["initial"]
    if not sampling_failure:
        assert samples[-1]["state"] == field["final"]
    else:
        assert sampling_failure in field["failure"]
    initial_temperature = np.array(field["initial"]["temperature"]).reshape(field["ny"], field["nx"]).mean(axis=0)
    marker_temperature = .5 * (float(initial_temperature.min()) + float(initial_temperature.max()))
    x = .5 * (np.array(field["xEdges"][:-1]) + field["xEdges"][1:])
    rows = []
    for index, sample in enumerate(samples):
        n = sample["state"]["steps"] - first_step
        assert sample["kind"] == ("accepted" if n else "initial")
        assert sample["integration"]["sampleResidualEvaluations"] == index + 1
        audit_sample_counters(field, sample, accepted[n - 1]["rhsCalls"] if n else 0, index == len(samples)-1)
        virtual = dict(field)
        virtual["final"] = sample["state"]
        for key in ("boundaryImpulse", "chemistryChange", "constraintChange", "integration"):
            virtual[key] = sample[key]
        virtual["complete"] = field["complete"] and sample["state"]["time"] == field.get("endTime", field["duration"])
        virtual["failure"] = field["failure"] if index == len(samples) - 1 else ""
        checked = audit_fields(virtual, path / "resolved-mechanism.yaml", accepted[:n], reference, fixture, checkpoint)
        temperature = np.array(sample["state"]["temperature"]).reshape(field["ny"], field["nx"]).mean(axis=0)
        markers = []
        for j, (a, b) in enumerate(zip(temperature[:-1], temperature[1:])):
            if a < marker_temperature <= b or b < marker_temperature <= a:
                markers.append(float(x[j] + (marker_temperature - a) / (b - a) * (x[j + 1] - x[j])))
        rows.append({"sample_index": index, "kind": sample["kind"], "temperature_marker_positions_m": markers,
                     "audit": checked})
    result = {"all_samples_passed": True, "sampling_complete": not sampling_failure, "sampling_failure": sampling_failure,
              "sample_count": len(samples), "every_accepted_steps": every,
              "samples_sha256": sha(path / "samples.jsonl"), "scaffold_field_sha256": sha(path / "field.json"),
              "marker_temperature_K": marker_temperature,
              "marker_method": "all crossings of the fixed initial midpoint temperature in the row-mean cell-centre profile, with linear interpolation; a position diagnostic, not an independently qualified flame speed",
              "samples": rows}
    if not field["complete"] and field["failure"] and samples[-1]["state"] == field["final"]:
        last_accepted_rhs = accepted[-1]["rhsCalls"] if accepted else 0
        result["terminal_attempt_statistics"] = {
            "last_accepted_step_rhs_calls": last_accepted_rhs,
            "terminal_rhs_calls": field["integration"]["rhsCalls"],
            "rhs_calls_after_last_accepted_step": field["integration"]["rhsCalls"]-last_accepted_rhs,
            "last_sample_rhs_calls": samples[-1]["integration"]["rhsCalls"],
            "last_sample_counters_include_stopped_trials": samples[-1]["integration"]["rhsCalls"] != last_accepted_rhs}
    target = audit_output / "sample-audits.json"
    target.write_text(json.dumps(result, indent=2) + "\n")
    return {"all_samples_passed": True, "sampling_complete": not sampling_failure, "sampling_failure": sampling_failure,
            "sample_count": len(samples), "sample_audits": str(target),
            "sample_audits_sha256": sha(target), "samples_sha256": result["samples_sha256"]}


def audit(path, reference, fixture, audit_output=None):
    r = json.loads((path / "field.json").read_text())
    steps = [json.loads(line) for line in (path / "steps.jsonl").read_text().splitlines()]
    checkpoint = None
    if r.get("restart") is not None:
        assert r["restart"]["checkpointFile"] == "restart.checkpoint"
        checkpoint = read_checkpoint(path / "restart.checkpoint", r["cells"], len(r["initial"]["U"][0]))
    report = audit_fields(r, path / "resolved-mechanism.yaml", steps, reference, fixture, checkpoint)
    if checkpoint is not None:
        report["restart_checkpoint_sha256"] = sha(path / "restart.checkpoint")
    report["field_sha256"] = sha(path / "field.json")
    sampling = audit_samples(path, r, steps, reference, fixture, Path(audit_output or path), checkpoint)
    if sampling is not None:
        report["accepted_state_samples"] = sampling
    (Path(audit_output or path) / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    return report, r


def check_initial_evaluation(field, case, directory):
    """A zero-time diagnostic is never a physical endpoint or checkpoint."""
    assert field["duration"] == field["endTime"] == field["initial"]["time"] == field["final"]["time"] == 0
    assert field["complete"] and not field["failure"] and field.get("restart") is None
    assert field["initial"] == field["final"] and field["initial"]["steps"] == 0
    assert case["accepted_steps"] == case["total_accepted_steps"] == 0
    assert case["numerical_checks_passed"] and not case["physical_endpoint_reached"]
    assert not (directory / "accepted.checkpoint").exists()


def plot_sample_history(case, field, output, plt):
    sampling = case.get("accepted_state_samples")
    if not sampling:
        return
    assert sha(sampling["sample_audits"]) == sampling["sample_audits_sha256"]
    history = json.loads(Path(sampling["sample_audits"]).read_text())
    samples_path = Path(case["directory"]) / "samples.jsonl"
    assert sha(samples_path) == history["samples_sha256"]
    rows = history["samples"]
    times = np.array([row["audit"]["time_s"] for row in rows])
    selected = sorted({int(np.argmin(abs(times - target))) for target in np.linspace(times[0], times[-1], 4)})
    snapshots = {}
    with samples_path.open() as stream:
        for index, line in enumerate(stream):
            if index in selected:
                snapshots[index] = json.loads(line)["state"]
    x = .5 * (np.array(field["xEdges"][:-1]) + field["xEdges"][1:])
    initial_t = np.array(field["initial"]["temperature"]).reshape(field["ny"], field["nx"]).mean(axis=0)
    origin = x[np.argmax(np.gradient(initial_t, x))]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), layout="constrained")
    for index in selected:
        state = snapshots[index]
        temperature = np.array(state["temperature"]).reshape(field["ny"], field["nx"]).mean(axis=0)
        axes[0, 0].plot((x - origin) * 1000, temperature, label=f"{state['time'] * 1e6:.3g} us")
    axes[0, 0].set(title="Actual accepted temperature profiles", xlabel="Distance from initial flame (mm)",
                   ylabel="Temperature (K)", xlim=(-.6, 2))
    axes[0, 0].legend(fontsize=8)
    positions = np.array([row["temperature_marker_positions_m"][0]
                          if len(row["temperature_marker_positions_m"]) == 1 else np.nan for row in rows])
    axes[0, 1].plot(times * 1e6, (positions - positions[0]) * 1e6, "o-", markersize=3)
    axes[0, 1].set(title=f"Fixed {history['marker_temperature_K']:.1f} K position marker",
                   ylabel="Displacement from initial position (um)")
    axes[1, 0].plot(times * 1e6, [row["audit"]["final_residual"]["hydrogen_consumption_kg_per_m2_s"] for row in rows], "o-", markersize=3)
    axes[1, 0].set(title="Actual integrated H2 consumption", ylabel="kg / m² / s")
    axes[1, 0].ticklabel_format(axis="y", style="plain", useOffset=False)
    axes[1, 1].plot(times * 1e6, [row["audit"]["final_residual"]["species_residual_L1_over_chemical_activity"] * 100 for row in rows], "o-", markersize=3)
    axes[1, 1].set(title="Species equation imbalance", ylabel="L1 residual / L1 reaction activity (%)")
    for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.set_xlabel("Actual accepted time (us)")
    for ax in axes.flat:
        ax.grid(alpha=.2)
    suffix = "" if history["sampling_complete"] else " | sampling stopped early"
    fig.suptitle(f"Detailed H2/air | {field['cells']} cells | actual end {field['final']['time'] * 1e6:g} / requested {field.get('endTime', field['duration']) * 1e6:g} us{suffix}")
    fig.supxlabel("Recorded accepted states and budgets; position marker is not a flame-speed qualification; multiple crossings remain in the JSON", fontsize=8)
    fig.savefig(output / f"native-{case['fixture_index']}-history.png", dpi=160); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--duration", type=float, default=1e-6, help="additional physical time from fixture or restart checkpoint")
    mode.add_argument("--evaluate-initial", action="store_true",
                      help="evaluate and audit the original initial state with zero advancement; never qualifies a physical endpoint")
    parser.add_argument("--restart-checkpoint", type=Path,
                        help="resume one selected fixture from an accepted checkpoint; BDF history is rebuilt")
    parser.add_argument("--rtol", type=float,
                        help="implicit relative local integration control; defaults remain unchanged")
    parser.add_argument("--conserved-atol", type=float,
                        help="implicit absolute local tolerance in the conserved-variable scales")
    parser.add_argument("--species-atol", type=float,
                        help="explicit implicit-probe local integration control; physical audit gates are unchanged")
    parser.add_argument("--bdf-order", type=int, choices=range(1, 6))
    parser.add_argument("--jacobian-order", type=int, choices=(1, 2))
    parser.add_argument("--newton-iterations", type=int)
    parser.add_argument("--continue-damped", type=int, choices=(0, 1))
    parser.add_argument("--reflect-species", type=int, choices=(0, 1))
    parser.add_argument("--sample-every", type=int, help="save and audit every Nth accepted state plus initial/final states")
    parser.add_argument("--max-samples", type=int, help="explicit diagnostic record budget; exhaustion preserves the accepted state and fails")
    parser.add_argument("--jobs", type=int, default=3, help="independent native cases in parallel; 1 for sequential execution")
    parser.add_argument("--grid", type=int, action="append",
                        help="reference fixture index; repeat to select cases (default: every fixture)")
    parser.add_argument("--audit-native-output", type=Path,
                        help="audit an already completed verifier run into a NEW output directory, without running the solver")
    args = parser.parse_args()
    verifier_hash = sha(Path(__file__))
    if args.evaluate_initial:
        assert args.restart_checkpoint is None and args.sample_every is None and args.max_samples is None
        args.duration = 0.
    else:
        assert np.isfinite(args.duration) and args.duration > 0
    assert args.rtol is None or np.isfinite(args.rtol) and 0 < args.rtol < 1
    assert args.conserved_atol is None or np.isfinite(args.conserved_atol) and args.conserved_atol > 0
    assert args.species_atol is None or np.isfinite(args.species_atol) and args.species_atol > 0
    assert args.newton_iterations is None or 0 < args.newton_iterations <= np.iinfo(np.int32).max
    assert args.sample_every is None or 0 < args.sample_every <= np.iinfo(np.int32).max
    assert args.max_samples is None or 0 < args.max_samples <= np.iinfo(np.int32).max
    requested_controls = {key: value for key, value in {
        "relativeTolerance": args.rtol, "absoluteConservedTolerance": args.conserved_atol,
        "absoluteSpeciesFraction": args.species_atol, "maximumOrder": args.bdf_order,
        "jacobianAdvectionOrder": args.jacobian_order, "maximumNonlinearIterations": args.newton_iterations,
        "continueDampedNewton": None if args.continue_damped is None else bool(args.continue_damped),
        "reflectSpeciesNewton": None if args.reflect_species is None else bool(args.reflect_species),
        "sampleEveryAcceptedSteps": args.sample_every, "maximumSamples": args.max_samples
    }.items() if value is not None}
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
        assert bool(previous.get("initial_evaluation_only", False)) == args.evaluate_initial
        assert previous["probe_sha256"] == probe_hash and previous["duration_s"] == args.duration
        assert previous["reference_report_sha256"] == sha(args.reference / "reference.json")
    indices = args.grid if args.grid is not None else (previous or {}).get("fixture_indices", list(range(len(reference["fixtures"]))))
    assert indices and len(set(indices)) == len(indices) and all(0 <= i < len(reference["fixtures"]) for i in indices)
    restart = None if previous is None else previous.get("restart")
    if args.restart_checkpoint is not None:
        assert len(indices) == 1
        requested_restart = {"path": str(args.restart_checkpoint.resolve()), "sha256": sha(args.restart_checkpoint)}
        assert previous is None or restart == requested_restart
        restart = requested_restart
    saved_cases = {} if previous is None else {c.get("fixture_index", i): c for i, c in enumerate(previous["cases"])}
    assert previous is None or all(i in saved_cases for i in indices)
    def run_case(item):
        number, fixture = item
        assert sha(fixture["path"]) == fixture["sha256"] and sha(args.probe) == probe_hash
        audit_path = args.output / f"native-{number}"
        path = (args.audit_native_output or args.output) / f"native-{number}"
        if previous is None:
            start = time.perf_counter()
            command = [str(args.probe.resolve()), reference["mechanism"], fixture["path"], str(path), format(args.duration, ".17g")]
            if restart is not None:
                assert sha(restart["path"]) == restart["sha256"]
                command += ["--restart", restart["path"]]
            for option in ("rtol", "conserved_atol", "species_atol"):
                if getattr(args, option) is not None:
                    command += ["--" + option.replace("_", "-"), format(getattr(args, option), ".17g")]
            for option in ("bdf_order", "jacobian_order", "newton_iterations", "continue_damped", "reflect_species", "sample_every", "max_samples"):
                if getattr(args, option) is not None:
                    command += ["--" + option.replace("_", "-"), str(getattr(args, option))]
            process = subprocess.run(command, capture_output=True, text=True)
            if restart is not None:
                assert sha(restart["path"]) == restart["sha256"]
            elapsed, returncode, failure = time.perf_counter() - start, process.returncode, process.stderr
            (args.output / f"native-{number}.stdout").write_text(process.stdout)
            (args.output / f"native-{number}.stderr").write_text(process.stderr)
        else:
            saved = saved_cases[number]
            assert Path(saved["directory"]).resolve() == path.resolve()
            assert all(sha(p) == digest for p, digest in saved.get("raw_sha256", {}).items())
            elapsed, returncode, failure = saved["elapsed_seconds_including_all_IO"], saved["returncode"], saved.get("failure", "")
            audit_path.mkdir()
        raw_hashes = {str(p): sha(p) for p in (path.iterdir() if path.is_dir() else [])
                      if p.is_file() and p.name not in {"audit.json", "audit-error.txt", "sample-audits.json"}}
        case = {"fixture_index": number, "returncode": returncode, "elapsed_seconds_including_all_IO": elapsed, "directory": str(path), "raw_sha256": raw_hashes}
        field = None
        if (path / "field.json").exists():
            try:
                checked, field = audit(path, reference, fixture, audit_path); case.update(checked)
                if args.evaluate_initial:
                    check_initial_evaluation(field, case, path)
                    case["initial_evaluation_passed"] = True
                assert field["duration"] == args.duration
                assert (field.get("restart") is not None) == (restart is not None)
                if restart is not None:
                    assert field["restart"]["sourcePath"] == restart["path"]
                    assert sha(path / "restart.checkpoint") == restart["sha256"]
                for key, value in requested_controls.items():
                    assert field["integration"][key] == value
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
        results = list(executor.map(run_case, ((i, reference["fixtures"][i]) for i in indices)))
    cases = [case for case, _ in results]
    fields = [field for _, field in results if field is not None]
    assert sha(args.probe) == probe_hash
    assert sha(Path(__file__)) == verifier_hash
    report = {"qualified_flame": False, "duration_s": args.duration, "restart": restart, "probe_sha256": probe_hash,
              "requested_species_absolute_tolerance": args.species_atol,
              "requested_implicit_controls": requested_controls,
              "execution_mode": "audit_existing" if previous is not None else ("evaluate_initial" if args.evaluate_initial else "run_and_audit"),
              "verifier_sha256": verifier_hash,
              "reference_report_sha256": sha(args.reference / "reference.json"), "fixture_indices": indices, "cases": cases,
              "parallel_jobs": args.jobs, "case_wall_seconds_including_audits": time.perf_counter() - wall_start,
              "all_numerical_checks_passed": all(c["numerical_checks_passed"] for c in cases),
              "all_endpoints_reached": all(c["physical_endpoint_reached"] for c in cases),
              "all_native_runs_succeeded": all(c["returncode"] == 0 for c in cases),
              "limits": ["reference and native share chemistry/transport backend", "constant-pressure BVP versus full compressible FV",
                         "short evolution from a reference initial condition does not predict or qualify flame speed",
                         "reference roundoff preparation recorded separately; no native species clipping", "no experiment or actual combustor qualification"]}
    if args.evaluate_initial:
        report.update(initial_evaluation_only=True,
                      all_initial_evaluations_passed=all(c.get("initial_evaluation_passed", False) for c in cases))
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
        profile = json.loads((args.reference / reference.get("profile_file", "profile-2.json")).read_text())
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
        change_title = "Initial readback temperature difference" if args.evaluate_initial else "Temperature change during native evolution"
        for ax, title, unit in [(axes[0, 0], "Temperature", "K"), (axes[0, 1], "Hydrogen mass fraction", "Y(H2)"),
                                (axes[1, 0], change_title, "K")]:
            ax.set(title=title, xlabel="x relative to reference flame (mm)", ylabel=unit, xlim=(-.6, 2)); ax.legend(fontsize=8); ax.grid(alpha=.2)
        axes[1, 1].set(title="Initial species equation imbalance", xlabel="Core spacing (um)", ylabel="L1 residual / L1 reaction activity")
        axes[1, 1].grid(alpha=.2)
        fig.suptitle("Detailed H2/air | initial-state evaluation | zero accepted steps" if args.evaluate_initial else
                     f"Detailed H2/air planar flame | native evolution requested: {args.duration*1e6:g} us")
        fig.supxlabel("Initial-data and residual audit only; no physical advancement or accepted checkpoint." if args.evaluate_initial else
                      "Reference-initialized short-time verification; no independently predicted flame-speed or experimental qualification.", fontsize=8)
        fig.savefig(args.output / "reacting-flame.png", dpi=160); plt.close(fig)
        for case, field in results:
            plot_sample_history(case, field, args.output, plt)
    requested_operation_passed = report.get("all_initial_evaluations_passed", report["all_endpoints_reached"])
    if not report["all_numerical_checks_passed"] or not requested_operation_passed or not report["all_native_runs_succeeded"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
