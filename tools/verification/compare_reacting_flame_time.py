#!/usr/bin/env python3
"""Measure temporal-control sensitivity of completed native flame fields.

Reports are supplied in order of decreasing tolerance, on the same grid and
physical endpoint. This measures field differences; it does not infer BDF order
from tolerance ratios or impose an experimental flame-accuracy qualification.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import cantera as ct
import numpy as np


TOLERANCES = ("relativeTolerance", "absoluteConservedTolerance", "absoluteSpeciesFraction")
CONTROLS = ("method", "errorControlCorrection", "maximumOrder", "maximumNonlinearIterations",
            "continueDampedNewton", "reflectSpeciesNewton", "spatialOrder",
            "jacobianAdvectionOrder", "maximumStep", "maximumResidualEvaluations")
GEOMETRY = ("cells", "nx", "ny", "height", "xEdges", "areas", "faces")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_case(path, fixture_index):
    report = json.loads(path.read_text())
    case = next(c for c in report["cases"] if c["fixture_index"] == fixture_index)
    assert case["returncode"] == 0 and case["numerical_checks_passed"] and case["physical_endpoint_reached"]
    directory = Path(case["directory"])
    field_path, mechanism = directory / "field.json", directory / "resolved-mechanism.yaml"
    assert case["raw_sha256"] and all(sha(p) == h for p, h in case["raw_sha256"].items())
    assert sha(field_path) == case["field_sha256"]
    field = json.loads(field_path.read_text())
    assert field["complete"] and not field["failure"]
    assert field["duration"] == field["final"]["time"] == case["time_s"] == report["duration_s"] > 0
    assert field["integration"] == case["integration"]
    return {"path": path, "report": report, "case": case, "field": field,
            "mechanism": mechanism, "mechanism_sha256": sha(mechanism)}


def compare(a, b, species, fuel):
    area = np.array(a["areas"])
    first, second = a["final"], b["final"]
    u, v = np.array(first["U"]), np.array(second["U"])
    dy = v[:, 4:] / v[:, 0, None] - u[:, 4:] / u[:, 0, None]
    dt = np.array(second["temperature"]) - first["temperature"]
    dp = (np.array(second["pressure"]) - first["pressure"]) / first["pressure"]
    dv = v[:, 1:3] / v[:, 0, None] - u[:, 1:3] / u[:, 0, None]
    source = np.array(first["chemistryDerivative"])
    refined_source = np.array(second["chemistryDerivative"])
    consumption = -float(area @ source[:, 4 + fuel] / a["height"])
    refined_consumption = -float(area @ refined_source[:, 4 + fuel] / a["height"])
    return {"maximum_temperature_difference_K": float(np.max(abs(dt))),
            "volume_RMS_temperature_difference_K": float(np.sqrt(area @ (dt * dt) / area.sum())),
            "maximum_mass_fraction_difference": float(np.max(abs(dy))),
            "maximum_mass_fraction_difference_by_species": dict(zip(species, np.max(abs(dy), axis=0).tolist())),
            "volume_RMS_mass_fraction_difference_by_species": dict(zip(species, np.sqrt(area @ (dy * dy) / area.sum()).tolist())),
            "maximum_pressure_relative_difference": float(np.max(abs(dp))),
            "maximum_velocity_difference_m_per_s": float(np.max(np.linalg.norm(dv, axis=1))),
            "fuel_consumption_kg_per_m2_s": [consumption, refined_consumption],
            "fuel_consumption_relative_difference": (abs(refined_consumption - consumption)
                / max(abs(consumption), abs(refined_consumption))) if consumption or refined_consumption else 0.}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, action="append", required=True,
                        help="completed verifier report, repeated from largest to smallest tolerances")
    parser.add_argument("--fixture-index", type=int, default=0)
    parser.add_argument("--fuel", default="H2")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert len(args.report) >= 3 and len(set(p.resolve() for p in args.report)) == len(args.report)
    data = [load_case(p, args.fixture_index) for p in args.report]
    base = data[0]
    gas = ct.Solution(str(base["mechanism"]))
    species, fuel = gas.species_names, gas.species_index(args.fuel)
    controls = base["field"]["integration"]
    previous_tolerances = None
    records = []
    for item in data:
        field, case, report = item["field"], item["case"], item["report"]
        assert item["mechanism_sha256"] == base["mechanism_sha256"]
        assert report["reference_report_sha256"] == base["report"]["reference_report_sha256"]
        assert all(field[key] == base["field"][key] for key in GEOMETRY)
        assert field["initial"] == base["field"]["initial"]
        assert field["final"]["time"] == base["field"]["final"]["time"]
        assert all(field["integration"][key] == controls[key] for key in CONTROLS)
        tolerances = np.array([field["integration"][key] for key in TOLERANCES])
        assert np.all(np.isfinite(tolerances)) and np.all(tolerances > 0)
        if previous_tolerances is not None:
            assert np.all(tolerances < previous_tolerances), "all three local tolerances must decrease"
        previous_tolerances = tolerances
        records.append({"report": str(item["path"]), "report_sha256": sha(item["path"]),
                        "field_sha256": sha(Path(case["directory"]) / "field.json"),
                        "probe_sha256": report["probe_sha256"],
                        "controls": {key: field["integration"][key] for key in TOLERANCES},
                        "accepted_steps": case["accepted_steps"], "rhs_calls": field["integration"]["rhsCalls"],
                        "elapsed_seconds_including_all_IO": case["elapsed_seconds_including_all_IO"],
                        "final_residual": case["final_residual"],
                        "fuel_consumption_kg_per_m2_s": -float(np.array(field["areas"])
                            @ np.array(field["final"]["chemistryDerivative"])[:, 4 + fuel] / field["height"]),
                        "maximum_temperature_drift_K": case["maximum_temperature_drift_K"]})
    pairs = [{"coarser": i, "finer": i + 1,
              **compare(data[i]["field"], data[i + 1]["field"], species, fuel)} for i in range(len(data) - 1)]
    reduction = {}
    for key in ("maximum_temperature_difference_K", "volume_RMS_temperature_difference_K",
                "maximum_mass_fraction_difference", "fuel_consumption_relative_difference"):
        values = [p[key] for p in pairs]
        reduction[key] = {"successive_differences": values,
                          "strictly_decreasing": all(y < x for x, y in zip(values, values[1:])),
                          "last_ratio": values[-1] / values[-2] if values[-2] else None}
    result = {"comparison_complete": True, "physical_time_s": base["field"]["final"]["time"],
              "cells": base["field"]["cells"], "fixture_index": args.fixture_index, "fuel": args.fuel,
              "mechanism_sha256": base["mechanism_sha256"], "runs": records, "successive_pairs": pairs,
              "coarsest_to_finest": compare(base["field"], data[-1]["field"], species, fuel),
              "difference_trends": reduction, "comparator_sha256": sha(__file__),
              "same_probe_binary": len(set(r["probe_sha256"] for r in records)) == 1,
              "limits": ["same serialized grid, initial residual and state, mechanism, endpoint and non-tolerance controls checked",
                         "different probe binaries require separately reviewed source provenance",
                         "adaptive tolerance refinement is not fixed-step convergence or BDF-order measurement",
                         "no experimental threshold imposed on these measured field differences",
                         "spatial error, long-time propagation and actual combustor qualification remain separate"],
              "qualified_flame": False, "general_performance_qualified": False}
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / "build/matplotlib-cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    field = base["field"]; nx = field["nx"]
    x = (np.array(field["xEdges"][:-1]) + field["xEdges"][1:]) / 2
    initial = np.array(field["initial"]["temperature"])[:nx]
    origin = x[np.argmax(np.gradient(initial, x))]; x = (x - origin) * 1000
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), layout="constrained")
    labels = [f"rtol = {record['controls']['relativeTolerance']:.2g}" for record in records]
    axes[0, 0].plot(x, initial, "k--", label="Initial reference")
    for label, item in zip(labels, data):
        axes[0, 0].plot(x, item["field"]["final"]["temperature"][:nx], label=label)
    for i in range(len(data) - 1):
        delta = np.array(data[i]["field"]["final"]["temperature"]) - data[i + 1]["field"]["final"]["temperature"]
        axes[0, 1].plot(x, delta[:nx], label=f"Level {i} minus {i + 1}")
    for ax, title, ylabel in [(axes[0, 0], "Native profiles near the flame", "Temperature (K)"),
                              (axes[0, 1], "Whole-domain temperature difference", "Temperature difference (K)")]:
        ax.set(title=title, xlabel="Distance from initial flame (mm)", ylabel=ylabel)
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    axes[0, 0].set_xlim(-.6, 2)
    axes[0, 1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    pair_labels = [f"{i} to {i + 1}" for i in range(len(pairs))]
    for key, label in [("maximum_temperature_difference_K", "Maximum"),
                       ("volume_RMS_temperature_difference_K", "Volume RMS")]:
        axes[1, 0].plot(pair_labels, [p[key] for p in pairs], "o-", label=label)
    axes[1, 0].set(title="Measured temperature sensitivity", xlabel="Successive tolerance levels", ylabel="Difference (K)")
    axes[1, 0].grid(alpha=.2); axes[1, 0].legend(fontsize=8)
    consumption = np.array([record["fuel_consumption_kg_per_m2_s"] for record in records])
    if consumption[-1] != 0:
        shown = (consumption - consumption[-1]) / abs(consumption[-1]) * 1e9
        ylabel = "Difference from finest control (parts per billion)"
    else:
        shown = consumption; ylabel = "kg / m² / s"
    axes[1, 1].plot(range(len(data)), shown, "o-")
    axes[1, 1].set(title=f"Integrated {args.fuel} consumption sensitivity", xlabel="Successively smaller tolerances", ylabel=ylabel, xticks=range(len(data)))
    axes[1, 1].grid(alpha=.2)
    axes[1, 1].ticklabel_format(axis="y", style="plain", useOffset=False)
    fig.suptitle(f"Detailed {args.fuel} flame | {field['cells']} cells | {result['physical_time_s'] * 1e6:g} microseconds")
    fig.supxlabel("Full chemistry, multicomponent transport and Soret; temporal sensitivity only, no spatial or experimental qualification", fontsize=8)
    fig.savefig(args.output / "time-refinement.png", dpi=160); plt.close(fig)
    print(json.dumps({"comparison_complete": True, "difference_trends": reduction,
                      "coarsest_to_finest": result["coarsest_to_finest"]}))


if __name__ == "__main__":
    main()
