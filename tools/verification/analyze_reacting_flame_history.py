#!/usr/bin/env python3
"""Diagnose a completed, independently audited planar-flame history.

Isotherm motion, shape change and velocity relative to explicit upstream
stations are reported separately. No steady-flame-speed qualification is
inferred. Multiple crossings and plateaus remain visible rather than being
silently converted into a selected flame position.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sample_history(path, initial, final, expected_count):
    """Stream complete records and require both audited endpoints and count."""
    count, previous_time, last_state = 0, None, None
    with Path(path).open() as stream:
        for line in stream:
            sample = json.loads(line)
            state = sample["state"]
            assert np.isfinite(state["time"])
            assert state == initial if count == 0 else state["time"] > previous_time
            count += 1
            assert count <= expected_count
            previous_time, last_state = state["time"], state
            yield sample
    assert count == expected_count and last_state == final


def marker_positions(x, temperature, level):
    x, temperature = np.asarray(x), np.asarray(temperature)
    assert x.ndim == 1 and x.shape == temperature.shape and len(x) >= 2
    assert np.all(np.isfinite(x)) and np.all(np.diff(x) > 0)
    assert np.all(np.isfinite(temperature)) and np.isfinite(level)
    crossings, plateaus = set(), []
    for j, (a, b) in enumerate(zip(temperature[:-1] - level, temperature[1:] - level)):
        if a == 0 and b == 0:
            if plateaus and plateaus[-1][1] == x[j]:
                plateaus[-1][1] = float(x[j + 1])
            else:
                plateaus.append([float(x[j]), float(x[j + 1])])
        elif a == 0:
            crossings.add(float(x[j]))
        elif b == 0:
            crossings.add(float(x[j + 1]))
        elif (a < 0 < b) or (b < 0 < a):
            crossings.add(float(x[j] - a * (x[j + 1] - x[j]) / (b - a)))
    positions = sorted(crossings)
    unique = positions[0] if len(positions) == 1 and not plateaus else None
    return {"crossings_m": positions, "plateaus_m": plateaus, "unique_position_m": unique}


def fit_motion(times, positions):
    times = np.asarray(times, dtype=float)
    assert times.ndim == 1 and len(times) == len(positions) and np.all(np.isfinite(times))
    assert np.all(np.diff(times) > 0)
    if len(times) < 3 or any(p is None for p in positions):
        return {"available": False, "reason": "Need at least three saved times with one isolated crossing at every time."}
    positions = np.asarray(positions, dtype=float)
    assert np.all(np.isfinite(positions))
    t, displacement = times - times[0], positions - positions[0]
    velocity, intercept = np.polyfit(t, displacement, 1)
    deviation = displacement - (intercept + velocity * t)
    interval = np.diff(positions) / np.diff(times)
    return {"available": True, "samples": len(times), "first_time_s": float(times[0]), "last_time_s": float(times[-1]),
            "lab_velocity_m_per_s": float(velocity),
            "endpoint_secant_velocity_m_per_s": float(displacement[-1] / t[-1]),
            "maximum_linear_fit_position_deviation_m": float(np.max(abs(deviation))),
            "interval_lab_velocity_min_m_per_s": float(interval.min()),
            "interval_lab_velocity_max_m_per_s": float(interval.max()),
            "method": "Unweighted least-squares line through actual saved positions and times; endpoint secant and interval range also retained."}


def translated_shape(x, initial, current, displacement, widths):
    if displacement is None:
        return None
    x, initial, current, widths = map(np.asarray, (x, initial, current, widths))
    assert x.shape == initial.shape == current.shape == widths.shape
    assert np.all(np.diff(x) > 0) and np.all(widths > 0) and np.isfinite(displacement)
    assert all(np.all(np.isfinite(a)) for a in (x, initial, current, widths))
    source_x = x - displacement
    overlap = (source_x >= x[0]) & (source_x <= x[-1])
    if not np.any(overlap):
        return None
    difference = current[overlap] - np.interp(source_x[overlap], x, initial)
    return {"maximum_K": float(np.max(abs(difference))),
            "width_weighted_RMS_K": float(np.sqrt(np.average(difference * difference, weights=widths[overlap]))),
            "evaluated_columns": int(overlap.sum()), "excluded_columns": int((~overlap).sum()),
            "evaluated_width_fraction": float(widths[overlap].sum() / widths.sum()),
            "first_evaluated_centre_m": float(x[overlap][0]), "last_evaluated_centre_m": float(x[overlap][-1])}


def read_audited_history(report_path, fixture_index, reference_directory):
    report = json.loads(report_path.read_text())
    assert all(report[k] for k in ("all_numerical_checks_passed", "all_endpoints_reached", "all_native_runs_succeeded"))
    assert not report.get("initial_evaluation_only", False)
    cases = [c for c in report["cases"] if c["fixture_index"] == fixture_index]
    assert len(cases) == 1
    case = cases[0]
    assert case["returncode"] == 0 and case["numerical_checks_passed"] and case["physical_endpoint_reached"]
    directory = Path(case["directory"])
    tracked = {str(report_path): sha(report_path), **case["raw_sha256"]}
    assert case["raw_sha256"] and all(sha(p) == h for p, h in tracked.items())
    field_path = directory / "field.json"
    field = json.loads(field_path.read_text())
    assert sha(field_path) == case["field_sha256"]
    assert field["complete"] and not field["failure"]
    assert field["duration"] == report["duration_s"] > 0
    if "endTime" not in field:
        # Histories predating restart support serialized only a cold duration.
        assert field["initial"]["time"] == field["initial"]["steps"] == 0 and field.get("restart") is None
    end_time = field.get("endTime", field["duration"])
    assert field["final"]["time"] == end_time == case["time_s"] > field["initial"]["time"]
    sampling = case["accepted_state_samples"]
    assert sampling["sampling_complete"] and sampling["all_samples_passed"]
    history_path = Path(sampling["sample_audits"])
    tracked[str(history_path)] = sampling["sample_audits_sha256"]
    assert sha(history_path) == sampling["sample_audits_sha256"]
    history = json.loads(history_path.read_text())
    assert history["all_samples_passed"] and history["sampling_complete"]
    assert history["scaffold_field_sha256"] == sha(field_path)
    samples_path = directory / "samples.jsonl"
    assert sha(samples_path) == history["samples_sha256"] == sampling["samples_sha256"]
    assert len(history["samples"]) == sampling["sample_count"]
    samples = sample_history(samples_path, field["initial"], field["final"], sampling["sample_count"])
    reference_path = reference_directory / "reference.json"
    assert sha(reference_path) == report["reference_report_sha256"]
    reference = json.loads(reference_path.read_text())
    tracked[str(reference_path)] = sha(reference_path)
    tracked.update({str(reference_directory / p): h for p, h in reference["files_sha256"].items()})
    fixture = reference["fixtures"][fixture_index]
    tracked[fixture["path"]] = fixture["sha256"]
    tracked[reference["mechanism"]] = reference["mechanism_sha256"]
    assert all(sha(p) == h for p, h in tracked.items())
    lines = Path(fixture["path"]).read_text().splitlines()
    inlet = np.array([float(v) for v in lines[4 + field["nx"]].split()])
    assert inlet.shape == (len(field["initial"]["U"][0]),) and inlet[0] > 0
    return report, case, field, history, samples, reference, inlet, tracked


def analyze(field, history, samples, reference, inlet, levels, alignment_level, stations, window):
    nx, ny = field["nx"], field["ny"]
    edges, area = np.asarray(field["xEdges"]), np.asarray(field["areas"])
    x, widths = .5 * (edges[:-1] + edges[1:]), np.diff(edges)
    assert nx >= 3 and ny >= 1 and field["height"] > 0 and area.shape == (nx * ny,)
    assert all(x[0] <= station <= x[-1] for station in stations)
    cells = [int(np.argmin(abs(x - station))) for station in stations]
    assert len(set(cells)) == len(cells), "Requested stations must resolve to different columns."
    inflow = [i for i, f in enumerate(field["faces"]) if f["neighbour"] is None and f["S"][0] < 0]
    outflow = [i for i, f in enumerate(field["faces"]) if f["neighbour"] is None and f["S"][0] > 0]
    assert len(inflow) == len(outflow) == ny, "Requires the planar-strip geometry of the flame probe."
    def row_mean(values):
        a = np.asarray(values)
        return a.reshape((ny, nx) + a.shape[1:]).mean(axis=0)
    initial_temperature = row_mean(field["initial"]["temperature"])
    initial_position = marker_positions(x, initial_temperature, alignment_level)["unique_position_m"]
    frames = []
    for index, sample in enumerate(samples):
        audited = history["samples"][index]
        state, checked = sample["state"], audited["audit"]
        assert audited["sample_index"] == index and checked["numerical_checks_passed"]
        if "total_accepted_steps" not in checked:
            assert field["initial"]["time"] == field["initial"]["steps"] == 0 and field.get("restart") is None
        assert checked["time_s"] == state["time"]
        assert checked.get("total_accepted_steps", checked["accepted_steps"]) == state["steps"]
        t = row_mean(state["temperature"])
        u = row_mean(state["U"])
        pressure = row_mean(state["pressure"])
        flux = np.asarray(state["faceFlux"])
        markers = [{"temperature_K": level, **marker_positions(x, t, level)} for level in levels]
        position = marker_positions(x, t, alignment_level)["unique_position_m"]
        displacement = None if position is None or initial_position is None else position - initial_position
        points = []
        for station, cell in zip(stations, cells):
            points.append({"requested_x_m": station, "actual_x_m": float(x[cell]),
                           "velocity_m_per_s": float(u[cell, 1] / u[cell, 0]),
                           "density_kg_per_m3": float(u[cell, 0]), "temperature_K": float(t[cell]),
                           "temperature_difference_from_inlet_K": float(t[cell] - reference["temperature_K"]),
                           "pressure_Pa": float(pressure[cell]),
                           "maximum_mass_fraction_difference_from_inlet": float(np.max(abs(u[cell, 4:] / u[cell, 0] - inlet[4:] / inlet[0]))),
                           "left_of_all_unique_markers": all(m["unique_position_m"] is not None and x[cell] < m["unique_position_m"] for m in markers)})
        frames.append({"time_s": state["time"], "total_accepted_steps": state["steps"], "markers": markers,
                       "alignment_displacement_m": displacement,
                       "maximum_temperature_drift_K": checked["maximum_temperature_drift_K"],
                       "row_mean_shape_after_translation": translated_shape(x, initial_temperature, t, displacement, widths),
                       "maximum_temperature_spread_across_rows_K": float(np.ptp(np.array(state["temperature"]).reshape(ny, nx), axis=0).max()),
                       "stations": points,
                       "inlet_mass_flux_kg_per_m2_s": float(-flux[inflow, 0].sum() / field["height"]),
                       "outlet_mass_flux_kg_per_m2_s": float(flux[outflow, 0].sum() / field["height"]),
                       "native_residual": checked["final_residual"]})
    assert len(frames) == len(history["samples"])
    chosen = [f for f in frames if f["time_s"] >= frames[-1]["time_s"] - window]
    times = np.array([f["time_s"] for f in chosen])
    late_motion = []
    for i, level in enumerate(levels):
        motion = {"temperature_K": level, **fit_motion(times, [f["markers"][i]["unique_position_m"] for f in chosen])}
        if motion["available"]:
            relative = []
            for j, station in enumerate(stations):
                velocity = [f["stations"][j]["velocity_m_per_s"] for f in chosen]
                mean = float(np.trapezoid(velocity, times) / (times[-1] - times[0]))
                relative.append({"requested_station_m": station, "actual_station_m": float(x[cells[j]]),
                                 "time_mean_station_velocity_m_per_s": mean,
                                 "station_relative_marker_speed_m_per_s": mean - motion["lab_velocity_m_per_s"],
                                 "station_left_of_all_unique_markers_in_window": all(f["stations"][j]["left_of_all_unique_markers"] for f in chosen),
                                 "maximum_station_temperature_difference_from_inlet_K": max(abs(f["stations"][j]["temperature_difference_from_inlet_K"]) for f in chosen),
                                 "maximum_station_mass_fraction_difference_from_inlet": max(f["stations"][j]["maximum_mass_fraction_difference_from_inlet"] for f in chosen)})
            motion["relative_to_stations"] = relative
        late_motion.append(motion)
    return {"analysis_complete": True, "qualified_flame_speed": False, "mature_combustion_goal_complete": False,
            "cells": field["cells"], "initial_time_s": frames[0]["time_s"], "final_time_s": frames[-1]["time_s"],
            "frame_count": len(frames), "marker_temperatures_K": levels, "alignment_temperature_K": alignment_level,
            "frames": frames, "late_window": {"requested_duration_s": window, "actual_first_s": float(times[0]),
                "actual_last_s": float(times[-1]), "samples": len(times), "marker_motion": late_motion},
            "reference_BVP_speed_m_per_s": reference["stages"][-1]["speed_m_per_s"],
            "inlet_mass_flux_initial_to_final_fraction": (frames[-1]["inlet_mass_flux_kg_per_m2_s"] / frames[0]["inlet_mass_flux_kg_per_m2_s"] - 1)
                if frames[0]["inlet_mass_flux_kg_per_m2_s"] else None,
            "limits": ["Planar-strip row-mean diagnostic; transverse temperature spread is reported, not discarded.",
                       "All isotherm crossings and exact plateaus retained; ambiguous or missing markers have no velocity estimate.",
                       "Station velocity minus lab marker velocity is a diagnostic, not a qualified burning velocity. Station temperature/composition departures from inlet are explicit.",
                       "Shape alignment interpolates initial row-mean temperature only on its overlap; it never changes a state, residual or accepted checkpoint.",
                       "Requested late window selects actual stored times; unequal sampling is explicit and no steady-state threshold is imposed.",
                       "Reference speed shares chemistry/transport backend and is not experimental validation; grid/time/domain errors remain separate."]}


def plot(result, output):
    os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / "build/matplotlib-cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    frames = result["frames"]
    times = np.array([f["time_s"] for f in frames])
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.5), layout="constrained")
    axes[0, 0].plot(times * 1e6, [f["maximum_temperature_drift_K"] for f in frames], label="Fixed-coordinate maximum")
    axes[0, 0].plot(times * 1e6, [f["row_mean_shape_after_translation"]["maximum_K"] if f["row_mean_shape_after_translation"] else np.nan for f in frames], label="Shape after translation (overlap)")
    axes[0, 0].set(title="Temperature evolution", ylabel="Difference from initial (K)")
    for i, level in enumerate(result["marker_temperatures_K"]):
        positions = [f["markers"][i]["unique_position_m"] for f in frames]
        values = np.array([np.nan if p is None else p for p in positions])
        axes[0, 1].plot(times * 1e6, (values - values[0]) * 1e6, label=f"{level:g} K")
        axes[1, 0].plot(.5 * (times[:-1] + times[1:]) * 1e6, np.diff(values) / np.diff(times), label=f"{level:g} K")
    axes[0, 1].set(title="Isotherm positions in fixed coordinates", ylabel="Displacement (um)")
    axes[1, 0].set(title="Interval marker velocities (not burning speeds)", ylabel="Lab velocity (m/s)")
    axes[1, 1].plot(times * 1e6, [100 * f["native_residual"]["species_residual_L1_over_chemical_activity"] for f in frames], label="Original native species residual")
    axes[1, 1].set(title="Species equation imbalance", ylabel="L1 residual / net chemical activity (%)")
    for ax in axes.flat:
        ax.set_xlabel("Physical time (us)"); ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle(f"Audited native flame history | {result['cells']} cells | {result['frame_count']} frames")
    fig.supxlabel("Detailed-physics history diagnostic; no steady flame-speed, grid-independence or experimental qualification.", fontsize=8)
    fig.savefig(output / "history.png", dpi=160); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--fixture-index", type=int, default=0)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--marker", type=float, action="append", required=True, help="isotherm in K; repeat for multiple markers")
    parser.add_argument("--align-temperature", type=float, required=True, help="single isotherm used only for shape-translation diagnosis (K)")
    parser.add_argument("--station", type=float, action="append", required=True, help="upstream x in m; nearest actual column is reported")
    parser.add_argument("--late-window", type=float, required=True, help="diagnostic end-window length in seconds, not a convergence criterion")
    args = parser.parse_args()
    assert all(np.isfinite(v) for v in [*args.marker, *args.station, args.align_temperature, args.late_window])
    assert all(v > 0 for v in [*args.marker, args.align_temperature, args.late_window])
    assert len(set(args.marker)) == len(args.marker) and len(set(args.station)) == len(args.station)
    assert not args.output.exists()
    source_hash = sha(__file__)
    report, case, field, history, samples, reference, inlet, tracked = read_audited_history(args.report, args.fixture_index, args.reference)
    result = analyze(field, history, samples, reference, inlet, args.marker, args.align_temperature, args.station, args.late_window)
    result.update(analyzer_sha256=source_hash, source_report=str(args.report), fixture_index=args.fixture_index,
                  source_sha256=tracked, native_probe_sha256=report["probe_sha256"], original_verifier_sha256=report["verifier_sha256"])
    assert sha(__file__) == source_hash and all(sha(p) == h for p, h in tracked.items())
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    plot(result, args.output)
    print(json.dumps({"frames": result["frame_count"], "final_time_s": result["final_time_s"],
                      "late_window": result["late_window"], "qualified_flame_speed": False}, indent=2))


if __name__ == "__main__":
    main()
