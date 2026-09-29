#!/usr/bin/env python3
"""Audit the native detailed-gas interface, not physical flame accuracy.

NASA7 cp/internal energy are evaluated directly from input coefficients. Reaction
trajectories and full transport matrices are compared against the separately
installed Cantera Python interface with the SAME mechanism and local tolerances.
That comparison shares a chemistry backend and is explicitly not an independent
mechanism or experimental validation. All tolerances below are dimensionless
except the reported temperatures/thermophysical data (SI).
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import time

import cantera as ct
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def difference(actual, expected):
    """Max difference / max(1, max|reference|); scale and units are in report."""
    a, b = np.asarray(actual), np.asarray(expected)
    assert a.shape == b.shape, (a.shape, b.shape)
    assert np.all(np.isfinite(a)) and np.all(np.isfinite(b)), "nonfinite reference data"
    return float(np.max(np.abs(a - b)) / max(1., float(np.max(np.abs(b)))))


def nasa7(gas, temperature, fractions):
    cp_species, u_species = [], []
    for species, weight in zip(gas.species(), gas.molecular_weights):
        data = species.input_data['thermo']
        assert data['model'] == 'NASA7', 'independent formula only covers NASA7'
        ranges = data['temperature-ranges']
        assert ranges[0] <= temperature <= ranges[-1]
        interval = next(i for i in range(len(ranges) - 1) if temperature <= ranges[i + 1])
        a = data['data'][interval]
        t = temperature
        cp_r = a[0] + a[1]*t + a[2]*t*t + a[3]*t**3 + a[4]*t**4
        h_rt = a[0] + a[1]*t/2 + a[2]*t*t/3 + a[3]*t**3/4 + a[4]*t**4/5 + a[5]/t
        cp_species.append(ct.gas_constant * cp_r / weight)
        u_species.append(ct.gas_constant * t * (h_rt - 1) / weight)
    return (math.fsum(y*c for y, c in zip(fractions, cp_species)),
            math.fsum(y*u for y, u in zip(fractions, u_species)))


def compare_state(native, gas):
    expected = {
        'temperature': gas.T, 'pressure': gas.P, 'density': gas.density,
        'internalEnergyDensity': gas.density*gas.int_energy_mass,
        'cp': gas.cp_mass, 'cv': gas.cv_mass, 'molecularWeight': gas.mean_molecular_weight,
        'viscosity': gas.viscosity, 'conductivity': gas.thermal_conductivity,
        'Y': gas.Y, 'speciesEnthalpies': gas.partial_molar_enthalpies/gas.molecular_weights,
        'massProductionRates': gas.net_production_rates*gas.molecular_weights,
        'massRateActivities': (gas.creation_rates+gas.destruction_rates)*gas.molecular_weights,
        'multicomponentDiffusion': gas.multi_diff_coeffs.ravel(order='F'),
        'thermalDiffusion': gas.thermal_diff_coeffs,
        'enthalpyReleaseRate': -float(np.dot(gas.partial_molar_enthalpies, gas.net_production_rates)),
    }
    return {name: difference(native[name], value) for name, value in expected.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--mechanism-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    assert ct.__version__ == '3.2.0', 'reference version must match pinned backend'
    cases = [
        ('hydrogen-ignition', 'h2o2.yaml', 1100., 101325., .001, {'H2': 2, 'O2': 1, 'N2': 3.76}),
        ('methane-reaction', 'gri30.yaml', 1400., 101325., .002, {'CH4': 1, 'O2': 2, 'N2': 7.52}),
        ('methane-high-pressure-thermo', 'gri30.yaml', 800., 1013250., 0., {'CH4': 1, 'O2': 2, 'N2': 7.52}),
        ('water-negative-formation-energy', 'h2o2.yaml', 800., 101325., 0., {'H2O': 1}),
    ]
    report = {
        'schema': 'cartmesh2d-detailed-gas-verification-v1', 'qualified_reacting_flow': False,
        'scope': 'native thermochemistry source module only; no spatial flow or experimental validation',
        'reference': {'cantera': ct.__version__, 'numpy': np.__version__, 'platform': platform.platform()},
        'probe_sha256': digest(args.probe), 'cases': [],
        'tolerances': {
            'state_and_transport_interface_normalized': 5e-10,
            'independent_nasa7_normalized': 5e-11,
            'integrated_temperature_normalized': 2e-6,
            'integrated_mass_fraction_absolute': 2e-6,
            'basis': 'floating point interface and local integration error allowances, not physical accuracy',
        },
    }
    start = time.monotonic()
    try:
        for name, filename, temperature, pressure, duration, composition in cases:
            mechanism = (args.mechanism_root / filename).resolve()
            command = [str(args.probe.resolve()), str(mechanism), str(temperature), str(pressure), str(duration)]
            command.extend(f'{species}={value}' for species, value in composition.items())
            proc = subprocess.run(command, capture_output=True, text=True, timeout=60)
            (args.output / (name + '.stdout')).write_text(proc.stdout)
            (args.output / (name + '.stderr')).write_text(proc.stderr)
            assert proc.returncode == 0, proc.stderr
            native = json.loads(proc.stdout)
            gas = ct.Solution(str(mechanism), transport_model='multicomponent')
            gas.TPX = temperature, pressure, composition
            comparisons = compare_state(native['initial'], gas)
            assert max(comparisons.values()) < 5e-10, (name, comparisons)
            nasa_cp, nasa_u = nasa7(gas, native['initial']['temperature'], native['initial']['Y'])
            independent = {
                'cp': difference(native['initial']['cp'], nasa_cp),
                'internalEnergy': difference(native['initial']['internalEnergyDensity']/native['initial']['density'], nasa_u),
            }
            assert max(independent.values()) < 5e-11, (name, independent)
            item = {'name': name, 'command': command, 'mechanism_sha256': digest(mechanism),
                    'species': gas.n_species, 'reactions': gas.n_reactions,
                    'initial_interface_errors': comparisons, 'independent_nasa7_errors': independent}
            if duration:
                reactor = ct.Reactor(gas, clone=True)
                network = ct.ReactorNet([reactor])
                network.rtol, network.atol = 1e-8, 1e-18
                network.max_steps = 20000
                network.advance(duration)
                final = native['final']
                temp_error = abs(final['temperature'] - reactor.T) / reactor.T
                y_error = float(np.max(np.abs(np.array(final['Y']) - reactor.phase.Y)))
                assert temp_error < 2e-6 and y_error < 2e-6, (name, temp_error, y_error)
                assert min(final['Y']) >= 0, 'accepted negative species'
                final_cp, final_u = nasa7(gas, final['temperature'], final['Y'])
                assert difference(final['cp'], final_cp) < 5e-11
                assert difference(final['internalEnergyDensity']/final['density'], final_u) < 5e-11
                item.update({'initial_temperature_K': temperature, 'final_temperature_K': final['temperature'],
                             'reference_temperature_K': reactor.T, 'temperature_relative_difference': temp_error,
                             'maximum_species_absolute_difference': y_error,
                             'element_mass_fraction_drift': native['elementDrift'],
                             'normalized_energy_drift': native['energyDrift'],
                             'time_s': native['time'], 'native_internal_steps': native['internalSteps']})
            report['cases'].append(item)
        report['passed'] = True
    except Exception as error:
        report['passed'] = False
        report['failure'] = str(error)
        raise
    finally:
        report['elapsed_seconds'] = time.monotonic() - start
        (args.output / 'summary.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'passed': True, 'cases': len(cases), 'summary': str(args.output/'summary.json')}))


if __name__ == '__main__':
    main()
