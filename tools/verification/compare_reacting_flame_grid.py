#!/usr/bin/env python3
"""Compare completed planar-flame grids at the same physical time and controls.

Fine cell averages are restricted by exact rectangle overlaps before converting
conserved variables to temperature. Initial projection differences are reported
separately; neither a reference initial profile nor a small drift qualifies a
steady native flame speed.
"""
import argparse
import json
import os
from pathlib import Path

import cantera as ct
import numpy as np

from compare_reacting_flame_time import CONTROLS, TOLERANCES, load_case, sha


def restriction(source_edges, target_edges):
    source, target = np.array(source_edges), np.array(target_edges)
    assert source[0] == target[0] and source[-1] == target[-1]
    assert np.all(np.diff(source) > 0) and np.all(np.diff(target) > 0)
    overlap = np.maximum(0., np.minimum(target[1:, None], source[None, 1:])
                         - np.maximum(target[:-1, None], source[None, :-1]))
    weights = overlap / np.diff(target)[:, None]
    # Exact overlap identities, allowing floating subtraction/accumulation
    # roundoff scaled by the number of participating cells. Not a CFD gate.
    allowance = 64 * np.finfo(float).eps * max(len(source), len(target))
    assert np.max(abs(weights.sum(axis=1) - 1)) <= allowance
    assert np.max(abs(np.diff(target) @ weights - np.diff(source))) / (source[-1] - source[0]) <= allowance
    return weights


def restrict_state(state, field, weights):
    u = np.array(state["U"]).reshape(field["ny"], field["nx"], -1)
    return np.einsum("ij,rjk->rik", weights, u).reshape(-1, u.shape[-1])


def temperature(gas, u, guesses):
    values = []
    for q, guess in zip(u, guesses):
        assert np.all(np.isfinite(q)) and q[0] > 0 and np.all(q[4:] >= 0)
        gas.set_unnormalized_mass_fractions(q[4:] / q[0])
        gas.TD = guess, q[0]
        gas.UV = q[3] / q[0] - (q[1] ** 2 + q[2] ** 2) / (2 * q[0] ** 2), 1 / q[0]
        values.append(gas.T)
    return np.array(values)


def norms(values, area):
    return {"maximum_K": float(np.max(abs(values))),
            "volume_RMS_K": float(np.sqrt(area @ (values * values) / area.sum()))}


def pair(coarse, fine, gas):
    weights = restriction(fine["xEdges"], coarse["xEdges"])
    area = np.array(coarse["areas"])
    differences = {}
    for label in ("initial", "final"):
        restricted = restrict_state(fine[label], fine, weights)
        projected_t = temperature(gas, restricted, coarse[label]["temperature"])
        differences[label] = projected_t - coarse[label]["temperature"]
    coarse_y = np.array(coarse["final"]["U"])
    coarse_y = coarse_y[:, 4:] / coarse_y[:, 0, None]
    restricted_y = restricted[:, 4:] / restricted[:, 0, None]
    return {"initial_projection_temperature_difference": norms(differences["initial"], area),
            "final_temperature_difference_after_restriction": norms(differences["final"], area),
            "temperature_evolution_difference_after_restriction": norms(differences["final"] - differences["initial"], area),
            "maximum_final_mass_fraction_difference_after_restriction": float(np.max(abs(coarse_y - restricted_y)))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", nargs=2, action="append", required=True, metavar=("REPORT", "FIXTURE_INDEX"),
                        help="repeat in coarse-to-fine order; multiple cases may come from the same report")
    parser.add_argument("--fuel", default="H2")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert len(args.case) >= 3
    data = [load_case(Path(p), int(index)) for p, index in args.case]
    base = data[0]; gas = ct.Solution(str(base["mechanism"])); fuel = gas.species_index(args.fuel)
    records = []
    previous_spacing = float("inf")
    for item in data:
        field, case, report = item["field"], item["case"], item["report"]
        assert item["mechanism_sha256"] == base["mechanism_sha256"]
        assert report["reference_report_sha256"] == base["report"]["reference_report_sha256"]
        assert field["final"]["time"] == base["field"]["final"]["time"]
        assert field["initial"]["time"] == base["field"]["initial"]["time"] == 0
        assert field["ny"] == base["field"]["ny"] and field["height"] == base["field"]["height"]
        assert all(field["integration"][key] == base["field"]["integration"][key] for key in CONTROLS + TOLERANCES)
        assert 0 < case["spacing_m"] < previous_spacing
        previous_spacing = case["spacing_m"]
        consumption = {}
        for label in ("initial", "final"):
            consumption[label] = -float(np.array(field["areas"])
                @ np.array(field[label]["chemistryDerivative"])[:, 4 + fuel] / field["height"])
        records.append({"report": str(item["path"]), "report_sha256": sha(item["path"]),
                        "fixture_index": case["fixture_index"], "field_sha256": case["field_sha256"],
                        "probe_sha256": report["probe_sha256"], "spacing_m": case["spacing_m"], "cells": case["cells"],
                        "accepted_steps": case["accepted_steps"], "rhs_calls": field["integration"]["rhsCalls"],
                        "initial_residual": case["initial_residual"], "final_residual": case["final_residual"],
                        "maximum_temperature_drift_K": case["maximum_temperature_drift_K"],
                        "fuel_consumption_kg_per_m2_s": consumption})
    pairs = []
    for i in range(len(data) - 1):
        consumed = [r["fuel_consumption_kg_per_m2_s"]["final"] for r in records[i:i + 2]]
        pairs.append({"coarser": i, "finer": i + 1, **pair(data[i]["field"], data[i + 1]["field"], gas),
                      "fuel_consumption_relative_difference": (abs(consumed[0] - consumed[1])
                          / max(abs(v) for v in consumed)) if any(consumed) else 0.})
    result = {"comparison_complete": True, "physical_time_s": base["field"]["final"]["time"],
              "fuel": args.fuel, "mechanism_sha256": base["mechanism_sha256"],
              "runs": records, "successive_pairs": pairs,
              "comparator_sha256": sha(__file__),
              "case_loader_sha256": sha(Path(__file__).with_name("compare_reacting_flame_time.py")),
              "projection": "piecewise-constant fine conserved averages integrated over coarse rectangles; EOS after restriction; initial projection discrepancy and evolution discrepancy shown separately",
              "limits": ["same mechanism, reference provenance, endpoint, y grid and integration controls checked",
                         "probe source provenance must be reviewed separately when binaries differ",
                         "non-nested grid restriction contributes projection error; not a formal spatial-order measurement",
                         "only one endpoint; no long-time or steady native flame-speed qualification",
                         "temporal refinement so far covers the coarsest grid, not every grid",
                         "reference and native share chemistry/transport data; no experimental qualification"],
              "qualified_flame": False, "mesh_independence_qualified": False}
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / "build/matplotlib-cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x0 = (np.array(base["field"]["xEdges"][:-1]) + base["field"]["xEdges"][1:]) / 2
    t0 = np.array(base["field"]["initial"]["temperature"])[:base["field"]["nx"]]
    origin = x0[np.argmax(np.gradient(t0, x0))]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), layout="constrained")
    for record, item in zip(records, data):
        field = item["field"]; nx = field["nx"]
        x = ((np.array(field["xEdges"][:-1]) + field["xEdges"][1:]) / 2 - origin) * 1000
        label = f"{record['spacing_m'] * 1e6:g} um ({field['cells']} cells)"
        axes[0, 0].plot(x, field["final"]["temperature"][:nx], label=label)
        change = np.array(field["final"]["temperature"])[:nx] - field["initial"]["temperature"][:nx]
        axes[0, 1].plot(x, change, label=label)
    for ax, title, ylabel in [(axes[0, 0], "Native temperature profiles", "Temperature (K)"),
                              (axes[0, 1], "Evolution from each conservative initial profile", "Temperature change (K)")]:
        ax.set(title=title, xlabel="Distance from initial flame (mm)", ylabel=ylabel, xlim=(-.6, 2))
        ax.legend(fontsize=8); ax.grid(alpha=.2)
    spacing = [r["spacing_m"] * 1e6 for r in records]
    axes[1, 0].plot(spacing, [r["fuel_consumption_kg_per_m2_s"]["final"] for r in records], "o-")
    axes[1, 0].set(title=f"Actual integrated {args.fuel} consumption", ylabel="kg / m² / s")
    for label in ("initial", "final"):
        axes[1, 1].plot(spacing, [r[label + "_residual"]["species_residual_L1_over_chemical_activity"] * 100 for r in records], "o-", label=label)
    axes[1, 1].set(title="Species equation imbalance", ylabel="L1 residual / L1 reaction activity (%)")
    axes[1, 1].legend(fontsize=8)
    for ax in axes[1]:
        ax.set(xlabel="Core grid spacing (um)", xticks=spacing); ax.invert_xaxis(); ax.grid(alpha=.2)
    fig.suptitle(f"Detailed {args.fuel} flame | native grid refinement at {result['physical_time_s'] * 1e6:g} microseconds")
    fig.supxlabel("Full chemistry, multicomponent transport and Soret; short transient from reference profiles, no steady flame-speed qualification", fontsize=8)
    fig.savefig(args.output / "grid-refinement.png", dpi=160); plt.close(fig)
    print(json.dumps({"comparison_complete": True, "successive_pairs": pairs,
                      "final_residuals": [r["final_residual"] for r in records]}))


if __name__ == "__main__":
    main()
