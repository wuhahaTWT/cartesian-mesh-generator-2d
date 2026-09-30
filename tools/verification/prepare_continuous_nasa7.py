#!/usr/bin/env python3
"""Prepare NASA7 with continuous h and s; never replace source data.

Keep every reaction, transport parameter, low-temperature coefficient and cp(T)
coefficient. Set only the two high-temperature integration constants so h(T) and
s(T) join their low-temperature values. This changes thermochemistry (including
reverse rates) and therefore requires its own physical sensitivity assessment.
The variant is opt-in input preparation, not a hidden solver normalization.
"""
import argparse
import copy
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

import cantera as ct
import numpy as np
from ruamel.yaml import YAML


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def primitive(coefficients, temperature, entropy=False, constant=True):
    """High-precision NASA7 h/R or s/R from the exact input doubles."""
    a = [Decimal.from_float(float(v)) for v in coefficients]
    t = Decimal.from_float(float(temperature))
    if entropy:
        value = a[0] * t.ln() + sum(a[j] * t ** j / Decimal(j) for j in range(1, 5))
        return value + (a[6] if constant else 0)
    value = sum(a[j] * t ** (j + 1) / Decimal(j + 1) for j in range(5))
    return value + (a[5] if constant else 0)


def join(coefficients):
    c = np.asarray(coefficients, dtype=float).copy()
    assert c.shape == (15,) and np.isfinite(c).all()
    temperature, high, low = c[0], c[1:8], c[8:15]
    with localcontext() as context:
        context.prec = 60
        high[5] = float(primitive(low, temperature) - primitive(high, temperature, constant=False))
        high[6] = float(primitive(low, temperature, entropy=True)
                        - primitive(high, temperature, entropy=True, constant=False))
    return c


def resolved_yaml(gas):
    # Cantera emits these generated scalar header fields even with header=False.
    # Remove only their complete root-level lines for repeatable source bytes.
    skip = ("generator:", "date:", "git-commit:")
    return "\n".join(line for line in gas.write_yaml(precision=17, header=False).splitlines()
                     if not line.startswith(skip)) + "\n"


def boundary_metrics(species, weight):
    t = species.thermo.coeffs[0]
    left = float(t)
    right = float(np.nextafter(t, np.inf))
    thermo = species.thermo
    h_left, h_right = thermo.h(left) / weight, thermo.h(right) / weight
    cv = thermo.cp(left) / weight - ct.gas_constant / weight
    e_left = h_left - ct.gas_constant * left / weight
    energy_allowance = 64 * np.finfo(float).eps * max(abs(e_left), abs(cv * left), 1.)
    return {"temperature_K": left, "right_temperature_K": right,
            "enthalpy_jump_J_per_kg": h_right - h_left,
            "entropy_jump_J_per_kg_K": (thermo.s(right) - thermo.s(left)) / weight,
            "cp_jump_J_per_kg_K": (thermo.cp(right) - thermo.cp(left)) / weight,
            "energy_join_defect_over_native_EOS_allowance": abs(h_right - h_left) / energy_allowance}


def prepare(mechanism, phase, output):
    mechanism, output = Path(mechanism).resolve(), Path(output).resolve()
    original_hash = sha(mechanism)
    gas = ct.Solution(str(mechanism), phase) if phase else ct.Solution(str(mechanism))
    assert gas.thermo_model == "ideal-gas", "preparation covers ideal-gas NASA7 mechanisms only"
    output.mkdir(parents=True, exist_ok=False)
    source = output / "original-resolved.yaml"
    source.write_text(resolved_yaml(gas))
    configured_transport_model = gas.transport_model
    input_reactions = [dict(reaction.input_data) for reaction in gas.reactions()]
    # Change the expanded YAML tree directly. Re-exporting modified Cantera
    # objects can round Ea/R through the YAML unit conversion a second time.
    # Both variants must load exactly the same forward reaction definitions.
    gas = ct.Solution(str(source), gas.name, transport_model=None)
    yaml = YAML(typ="safe")
    original_document = yaml.load(source.read_text())
    document = copy.deepcopy(original_document)
    entries = {item["name"]: item for item in document["species"]}
    reactions = [dict(reaction.input_data) for reaction in gas.reactions()]
    source_species = gas.species()
    changes = []
    for k, original in enumerate(source_species):
        assert isinstance(original.thermo, ct.NasaPoly2), f"unsupported thermo model for {original.name}; no fallback"
        coefficients = original.thermo.coeffs.copy()
        updated = join(coefficients)
        assert np.array_equal(updated[:6], coefficients[:6]) and np.array_equal(updated[8:], coefficients[8:])
        species = ct.Species.from_dict(original.input_data)
        species.thermo = ct.NasaPoly2(original.thermo.min_temp, original.thermo.max_temp,
                                    original.thermo.reference_pressure, updated)
        data = entries[original.name]["thermo"]
        assert data["model"] == "NASA7" and len(data["data"]) == 2
        assert np.array_equal(data["data"][0], coefficients[8:]) and np.array_equal(data["data"][1], coefficients[1:8])
        data["data"][1] = updated[1:8].tolist()
        weight = gas.molecular_weights[k]
        changes.append({"species": original.name, "original_coefficients": coefficients.tolist(),
                        "continuous_coefficients": updated.tolist(),
                        "high_enthalpy_offset_J_per_kg": float((updated[6] - coefficients[6]) * ct.gas_constant / weight),
                        "high_entropy_offset_J_per_kg_K": float((updated[7] - coefficients[7]) * ct.gas_constant / weight),
                        "original_boundary": boundary_metrics(original, weight),
                        "continuous_boundary": boundary_metrics(species, weight)})
    derived = output / "continuous.yaml"
    with derived.open("w") as stream:
        stream.write("# Explicit continuous NASA7 variant; see report.json for source and all coefficient changes.\n")
        yaml.dump(document, stream)
    reloaded = ct.Solution(str(derived), gas.name, transport_model=None)
    assert [dict(reaction.input_data) for reaction in reloaded.reactions()] == reactions
    assert reloaded.species_names == gas.species_names and reloaded.element_names == gas.element_names
    assert np.array_equal(reloaded.molecular_weights, gas.molecular_weights)
    for k, original in enumerate(source_species):
        prepared = reloaded.species(k)
        left, right = dict(original.input_data), dict(prepared.input_data)
        left.pop("thermo"); right.pop("thermo")
        assert left == right, "non-thermodynamic species data changed"
        assert np.array_equal(prepared.thermo.coeffs, changes[k]["continuous_coefficients"])
        assert changes[k]["continuous_boundary"]["energy_join_defect_over_native_EOS_allowance"] <= 1
    assert sha(mechanism) == original_hash
    report = {"qualified_mechanism": False, "preparation_complete": True,
              "method": "Preserve low-temperature h/s anchors and all NASA7 cp coefficients; integrate cp and cp/T across the join by changing only high-region a5 and a6.",
              "formula_reference": "https://cantera.org/stable/reference/thermo/species-thermo.html#the-nasa-7-coefficient-polynomial-parameterization",
              "cantera_version": ct.__version__, "script_sha256": sha(__file__),
              "source": {"path": str(mechanism), "sha256": original_hash, "phase": gas.name},
              "configured_transport_model": configured_transport_model,
              "source_resolution": {"reaction_definitions_compared_to": "original-resolved.yaml",
                                    "input_to_resolved_numeric_round_trip_reaction_indices": [i for i, (a, b) in enumerate(zip(input_reactions, reactions)) if a != b],
                                    "note": "Cantera source expansion can round activation-energy unit conversions; both exported variants share identical forward definitions. Original-source sensitivity must include this round trip."},
              "species": gas.n_species, "reactions": gas.n_reactions,
              "unchanged": ["species/elements/order", "all forward reaction definitions", "transport parameters", "temperature ranges", "low-region NASA7 coefficients", "all cp(T) coefficients"],
              "changed": ["high-region enthalpy and entropy integration constants", "corresponding equilibrium constants and reverse reaction rates"],
              "checks": {"coefficient_and_reaction_identity": True,
                         "join_allowance": "64*double_epsilon*max(abs(e),abs(cv*T),1), J/kg; existing native energy inversion roundoff allowance, not a physical accuracy gate",
                         "maximum_join_defect_over_native_EOS_allowance": max(c["continuous_boundary"]["energy_join_defect_over_native_EOS_allowance"] for c in changes)},
              "changes": changes,
              "files_sha256": {p.name: sha(p) for p in (source, derived)},
              "limitations": ["NASA7 two-region ideal-gas preparation only", "original cp(T) joins, including jumps in cp or its derivatives, are preserved", "original mechanism and solver defaults are untouched", "no kinetics, flame, experiment or combustor qualification from continuity alone"]}
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--phase")
    parser.add_argument("--output", type=Path, required=True, help="new directory; source file is never modified")
    args = parser.parse_args()
    report = prepare(args.mechanism, args.phase, args.output)
    print(json.dumps({"preparation_complete": True, "qualified_mechanism": False,
                      "species": report["species"], "reactions": report["reactions"], "output": str(args.output)}))


if __name__ == "__main__":
    main()
