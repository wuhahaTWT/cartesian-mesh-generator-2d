#!/usr/bin/env python3
"""Check a prepared NASA7 variant and measure its H2/CH4 reactor sensitivity.

Independent Gauss quadrature checks h/s against cp and cp/T. Native round trips
and closed-reactor trajectories check the implementation. Both Cantera interfaces
share underlying physical data; this does not qualify a combustion mechanism.
"""
import argparse
import json
from pathlib import Path
import subprocess
import time
import traceback

import cantera as ct
import numpy as np

from prepare_continuous_nasa7 import sha


def relative(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
    return float(np.max(abs(a - b)) / max(1., float(np.max(abs(b)))))


def integral_check(gas):
    nodes, weights = np.polynomial.legendre.leggauss(32)
    rows = []
    for species in gas.species():
        model = species.thermo
        lo, mid, hi = model.min_temp, model.coeffs[0], model.max_temp
        temperatures = sorted(set([lo, .5 * (lo + mid), np.nextafter(mid, -np.inf), mid,
                                   np.nextafter(mid, np.inf), .5 * (mid + hi), hi]))
        h_error, s_error = 0., 0.
        for t in temperatures:
            h, s = model.h(lo), model.s(lo)
            cuts = [lo] + ([mid] if lo < mid < t else []) + [t]
            for left, right in zip(cuts[:-1], cuts[1:]):
                if left == right:
                    continue
                x = .5 * (right - left) * nodes + .5 * (right + left)
                cp = np.array([model.cp(float(a)) for a in x])
                assert np.all(cp > ct.gas_constant), "nonpositive ideal-gas cv in quadrature samples"
                h += .5 * (right - left) * float(weights @ cp)
                s += .5 * (right - left) * float(weights @ (cp / x))
            h_error = max(h_error, abs(model.h(t) - h) / max(abs(h), abs(model.cp(t) * t), 1.))
            s_error = max(s_error, abs(model.s(t) - s) / max(abs(s), 1.))
        rows.append({"species": species.name, "evaluated_temperatures_K": temperatures,
                     "enthalpy_integral_scaled_error": h_error, "entropy_integral_scaled_error": s_error})
    return rows


def verify(prepared, probe, output, fuel):
    prepared, probe, output = Path(prepared).resolve(), Path(probe).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    preparation = json.loads((prepared / "report.json").read_text())
    assert preparation["preparation_complete"] and not preparation["qualified_mechanism"]
    assert all(sha(prepared / name) == digest for name, digest in preparation["files_sha256"].items())
    original = Path(preparation["source"]["path"])
    assert sha(original) == preparation["source"]["sha256"]
    variant = prepared / "continuous.yaml"
    phase = preparation["source"]["phase"]
    x = {"H2": 1.6, "O2": 1., "N2": 3.76} if fuel == "H2" else {"CH4": 1., "O2": 2., "N2": 7.52}
    reactor_temperature, duration = (1100., .001) if fuel == "H2" else (1400., .002)
    report = {"passed": False, "qualified_mechanism": False, "fuel": fuel,
              "scope": "Thermodynamic consistency, native interface and short closed-reactor sensitivity only",
              "preparation_report": str(prepared / "report.json"), "preparation_report_sha256": sha(prepared / "report.json"),
              "probe_sha256": sha(probe), "verifier_sha256": sha(__file__),
              "existing_verification_gates": {"interface_normalized": 5e-10, "thermo_formula_normalized": 5e-11,
                                              "reactor_temperature_relative": 2e-6, "reactor_Y_absolute": 2e-6},
              "native_cases": [], "limits": ["backend and mechanism data are shared", "no flame, experiment or combustor qualification", "new thermochemistry requires its own source-bound flame reference and physical sensitivity validation"]}
    started = time.perf_counter()
    try:
        source_gas = ct.Solution(str(original), transport_model=None)
        assert source_gas.name == phase, "native probe currently selects the first source phase"
        resolved_gas = ct.Solution(str(prepared / "original-resolved.yaml"), phase, transport_model=None)
        gas = ct.Solution(str(variant), phase, transport_model=None)
        original_integrals, corrected_integrals = integral_check(source_gas), integral_check(gas)
        report["cp_integral_checks"] = {"original": original_integrals, "continuous": corrected_integrals}
        assert max(max(r["enthalpy_integral_scaled_error"], r["entropy_integral_scaled_error"]) for r in corrected_integrals) < 5e-11
        rates = []
        for t in [400., 999.9999, 1000., 1000.0001, 1800., 3000.]:
            for p in [101325., 1013250.]:
                for g in (source_gas, resolved_gas, gas):
                    g.TPX = t, p, x
                assert np.array_equal(resolved_gas.forward_rate_constants, gas.forward_rate_constants)
                assert source_gas.cp_mass == gas.cp_mass
                rates.append({"temperature_K": t, "pressure_Pa": p,
                              "source_expansion_forward_rate_normalized_change": relative(resolved_gas.forward_rate_constants, source_gas.forward_rate_constants),
                              "reverse_rate_normalized_change": relative(gas.reverse_rate_constants, source_gas.reverse_rate_constants),
                              "internal_energy_change_J_per_kg": gas.int_energy_mass - source_gas.int_energy_mass})
        report["fixed_state_sensitivity"] = rates
        states = {}
        for label, mechanism in [("original", original), ("continuous", variant)]:
            gas = ct.Solution(str(mechanism), phase, transport_model=None)
            for index, (temperature, dt) in enumerate([(999.9999, 0.), (1000., 0.), (1000.0001, 0.), (reactor_temperature, duration)]):
                name = f"{label}-{index}"
                command = [str(probe), str(mechanism), format(temperature, ".17g"), "101325", format(dt, ".17g")]
                command.extend(f"{k}={v}" for k, v in x.items())
                start = time.perf_counter()
                proc = subprocess.run(command, capture_output=True, text=True, timeout=90)
                stdout, stderr = output / (name + ".stdout"), output / (name + ".stderr")
                stdout.write_text(proc.stdout); stderr.write_text(proc.stderr)
                row = {"name": name, "command": command, "returncode": proc.returncode,
                       "elapsed_seconds_including_IO": time.perf_counter() - start,
                       "stdout_sha256": sha(stdout), "stderr_sha256": sha(stderr)}
                report["native_cases"].append(row)
                assert proc.returncode == 0, proc.stderr
                native = json.loads(proc.stdout)
                gas.TPX = temperature, 101325., x
                q = native["initial"]
                errors = {"temperature": relative(q["temperature"], temperature), "pressure": relative(q["pressure"], gas.P),
                          "cp": relative(q["cp"], gas.cp_mass), "internalEnergyDensity": relative(q["internalEnergyDensity"], gas.density * gas.int_energy_mass),
                          "speciesEnthalpies": relative(q["speciesEnthalpies"], gas.partial_molar_enthalpies / gas.molecular_weights)}
                row.update(temperature_round_trip_difference_K=q["temperature"] - temperature,
                           interface_errors=errors, requested_temperature_K=temperature)
                if label == "continuous":
                    assert max(errors.values()) < 5e-10
                if dt:
                    assert native["time"] == dt and native["internalSteps"] > 0
                    assert native["elementDrift"] <= 1e-8 and native["energyDrift"] <= 1e-8
                    states[label] = native["final"]
                    reactor = ct.Reactor(gas, clone=True)
                    network = ct.ReactorNet([reactor])
                    network.rtol, network.atol, network.max_steps = 1e-8, 1e-18, 20000
                    network.advance(dt)
                    t_error = abs(native["final"]["temperature"] - reactor.T) / reactor.T
                    y_error = float(np.max(abs(np.array(native["final"]["Y"]) - reactor.phase.Y)))
                    assert t_error < 2e-6 and y_error < 2e-6
                    row.update(reactor_temperature_relative_error=t_error, reactor_mass_fraction_absolute_error=y_error,
                               final_temperature_K=native["final"]["temperature"], element_drift=native["elementDrift"], energy_drift=native["energyDrift"])
        report["reactor_sensitivity"] = {"duration_s": duration, "initial_temperature_K": reactor_temperature,
                                         "final_temperature_change_K": states["continuous"]["temperature"] - states["original"]["temperature"],
                                         "maximum_final_mass_fraction_change": float(np.max(abs(np.array(states["continuous"]["Y"]) - states["original"]["Y"])))}
        assert sha(probe) == report["probe_sha256"] and sha(__file__) == report["verifier_sha256"]
        report["passed"] = True
    except Exception:
        report["failure"] = traceback.format_exc()
        raise
    finally:
        report["elapsed_seconds"] = time.perf_counter() - started
        (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fuel", choices=("H2", "CH4"), required=True, help="module-check mixture; not a target combustor specification")
    args = parser.parse_args()
    report = verify(args.prepared, args.probe, args.output, args.fuel)
    print(json.dumps({"passed": report["passed"], "qualified_mechanism": False, "native_cases": len(report["native_cases"]),
                      "reactor_sensitivity": report["reactor_sensitivity"]}))


if __name__ == "__main__":
    main()
