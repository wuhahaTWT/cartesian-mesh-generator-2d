#!/usr/bin/env python3
"""Read back Re100 cavity runs against a second published numerical reference.

The original Ghia checks remain visible and participate in allChecksPassed.
A published finite-grid/extrapolated result is not an analytic exact solution.
This tool never changes fields, solver controls, or existing acceptance limits.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
import sys
import time
from typing import Any

import verify_native_flow as native

REFERENCE = Path(__file__).with_name('references') / 'cavity-marchi-2009-re100.json'


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def load_reference(path: Path) -> dict[str, Any]:
    reference = json.loads(path.read_text())
    require(reference.get('format') == 'cartmesh2d-published-cavity-reference-v1', 'Unknown reference format')
    require(reference.get('reynolds') == 100 and reference.get('lidSpeed') == 1
            and reference.get('domain') == [0, 0, 1, 1], 'Reference is not the unit Re100 cavity')
    require(reference.get('doi') == '10.1590/S1678-58782009000300004', 'Unexpected reference DOI')
    for component in ('u', 'v'):
        rows = reference.get(component)
        require(isinstance(rows, list) and len(rows) == 15, 'Reference must contain all 15 points per component')
        for i, row in enumerate(rows, 1):
            require(row.get('coordinate') == i / 16, 'Reference coordinate missing, reordered or duplicated')
            for key in ('value', 'estimatedError'):
                require(isinstance(row.get(key), (int, float)) and math.isfinite(row[key]),
                        'Nonfinite reference data')
            require(row['estimatedError'] >= 0, 'Negative reference error estimate')
    return reference


def compare_reference(cells: list[dict[str, float]], reference: dict[str, Any]) -> dict[str, Any]:
    """Use the same calibrated affine point sampler as the existing Ghia gate."""
    samples, errors, uncertainties = [], [], []
    for component in ('u', 'v'):
        for row in reference[component]:
            coordinate = row['coordinate']
            x, y = (.5, coordinate) if component == 'u' else (coordinate, .5)
            walls = ([dict(x=x, y=0., u=0.), dict(x=x, y=1., u=1.)] if component == 'u'
                     else [dict(x=0., y=y, v=0.), dict(x=1., y=y, v=0.)])
            actual = native.affine_sample(cells, x, y, component, boundary=walls)
            error = actual - row['value']
            samples.append(dict(field=component, coordinate=coordinate, actual=actual,
                                reference=row['value'], referenceEstimatedError=row['estimatedError'], difference=error))
            errors.append(error)
            uncertainties.append(row['estimatedError'])
    rms = lambda values: math.sqrt(math.fsum(value * value for value in values) / len(values))
    return dict(reference='Marchi, Suero & Araki 2009, Tables 6 and 7, Re100',
                samplingMethod='distance-weighted affine fit with Dirichlet wall anchors; v2',
                samples=samples, centrelineRmse=rms(errors), centrelineMaxError=max(map(abs, errors)),
                referenceEstimatedErrorRms=rms(uncertainties),
                referenceEstimatedErrorMax=max(uncertainties),
                estimatedRmseRange=[rms([max(0., abs(e) - u) for e, u in zip(errors, uncertainties)]),
                                    rms([abs(e) + u for e, u in zip(errors, uncertainties)])],
                uncertaintyScope='Uses the authors estimated discretization errors; not rigorous confidence bounds or exact truth')


def option(command: list[str], name: str) -> str:
    require(command.count(name) == 1, f'Run command must contain exactly one {name}')
    index = command.index(name)
    require(index + 1 < len(command), f'Missing value for {name}')
    return command[index + 1]


def audit_runs(path: Path, reference_path: Path = REFERENCE) -> dict[str, Any]:
    reference = load_reference(reference_path)
    runs = json.loads(path.read_text())
    require(isinstance(runs, list) and len(runs) >= 3, 'At least three actual runs are required')
    report: dict[str, Any] = dict(format='cartmesh2d-cavity-reference-audit-v1',
        scope='Independent readback; published reference comparison and existing Ghia checks reported separately',
        runsManifest=dict(path=str(path.resolve()), sha256=native.sha256_file(path)),
        reference=dict(path=str(reference_path.resolve()), sha256=native.sha256_file(reference_path), data=reference),
        verifierSha256=native.sha256_file(Path(native.__file__)), scriptSha256=native.sha256_file(Path(__file__)),
        cases=[], issues=[], independentAuditsPassed=False, allChecksPassed=False)
    signatures, prefixes = set(), set()
    for run in runs:
        require(run.get('returnCode') == 0 and run.get('timedOut') is False, 'Run failed or timed out')
        mesh, prefix = Path(run['mesh']).resolve(), Path(run['prefix']).resolve()
        require(prefix not in prefixes, 'Duplicate result prefix')
        prefixes.add(prefix)
        require(native.sha256_file(mesh) == run['meshSha256'], 'Mesh hash changed after run')
        command = run['command']
        require(isinstance(command, list) and command and all(isinstance(s, str) for s in command), 'Invalid command')
        require(native.sha256_file(Path(command[0])) == run['binarySha256'], 'Producer binary changed or unavailable')
        require(Path(option(command, '--mesh')).resolve() == mesh and
                Path(option(command, '--output')).resolve() == prefix, 'Command artifact paths mismatch')
        require(option(command, '--case') == 'cavity', 'Run is not cavity')
        require(not any(x in command for x in ('--time-step', '--steps', '--restart')), 'Steady refinement cannot include transient/restart')
        require(float(option(command, '--nu')) == .01 and float(option(command, '--speed')) == 1., 'Run is not Re100 unit-speed')
        convection = option(command, '--convection')
        method = option(command, '--pressure-preconditioner')
        tolerance = float(option(command, '--tolerance'))
        require(math.isfinite(tolerance) and tolerance > 0, 'Invalid run tolerance')
        controls = native.argument_parser().parse_args([])
        controls.max_iterations = int(option(command, '--max-iterations'))
        require(controls.max_iterations > 0, 'Invalid run iteration budget')
        result = native.verify_case(mesh, prefix, 'cavity', .01, 1., controls)
        summary = result['native']
        require(summary.get('case') == 'cavity' and summary.get('nu') == .01 and summary.get('speed') == 1., 'Result physics mismatch')
        require(not summary.get('temporalDiscretization'), 'Result is transient')
        require(summary.get('convection') == convection and summary.get('pressurePreconditioner') == method
                and summary.get('tolerance') == tolerance, 'Result controls do not match command')
        stress = option(command, '--viscous-stress') if '--viscous-stress' in command else 'symmetric'
        require(stress == 'symmetric' and summary.get('viscousStress') == stress,
                'Reference sequence requires matching symmetric stress')
        bounds = result['meshMeasurement']['bounds']
        require(len(bounds) == 4 and all(native.close(value, target, 1e-12, 1e-12)
                for value, target in zip(bounds, (0., 0., 1., 1.))),
                'Reference sequence requires unit square')
        signatures.add((run['binarySha256'], convection, method, tolerance, summary.get('viscousStress'),
                        summary.get('pressureDiscretization'), summary.get('pressureBoundaryReconstruction')))
        cells = native.load_csv(Path(str(prefix) + '.cells.csv'), ('cell', 'x', 'y', 'u', 'v', 'speed'))
        numeric = [{key: native.finite(row[key], key) for key in ('x', 'y', 'u', 'v', 'speed')} for row in cells]
        comparison = compare_reference(numeric, reference)
        report['cases'].append(dict(label=run['label'], run=run, independentAudit=result, published=comparison))
        if not result['valid']:
            report['issues'].append(f"{run['label']}: independent geometry/physics audit failed")
    require(len(signatures) == 1, 'Refinement sequence mixes binary or numerical controls')
    audits = [dict(c['independentAudit'], label=c['label']) for c in report['cases']]
    original = native.sequence_checks(audits)
    second = copy.deepcopy(audits)
    for item, case in zip(second, report['cases']):
        item['benchmark'] = case['published']
    published = native.sequence_checks(second)
    report.update(ghiaSequence=original, publishedSequence=published,
                  independentAuditsPassed=not report['issues'],
                  trendRule='Existing sequence_checks rule: fine RMSE <= 1.10 * coarse RMSE; unchanged for both references',
                  allChecksPassed=not report['issues'] and original['valid'] and published['valid'])
    return report


def generate_runs(root: Path, levels: list[int], timeout: float,
                  mesh_cli: Path, flow_cli: Path) -> Path:
    """Generate and solve an ordered, bounded series; preserve failed commands."""
    require(len(levels) >= 3 and sorted(set(levels)) == levels and all(4 <= n <= 7 for n in levels),
            'Need at least three distinct ascending levels from 4 through 7')
    require(math.isfinite(timeout) and timeout > 0, 'Timeout must be finite and positive')
    root = root.resolve()
    mesh_cli, flow_cli = mesh_cli.resolve(strict=True), flow_cli.resolve(strict=True)
    root.mkdir(parents=True, exist_ok=False)
    runs, generation = [], []
    for level in levels:
        label = f'cavity-l{level}'
        directory = root / label
        directory.mkdir()
        request = native.Request(label, 'cavity', level, 1. / (2 ** level - 2), .01, 1.)
        mesh_prefix = directory / 'mesh'
        mesh_command = native.mesh_command(mesh_cli, mesh_prefix, request, root)
        stage = native.run(mesh_command, root / 'logs' / f'mesh-{label}', timeout)
        generation.append(stage)
        native.write_json(root / 'mesh-generation.json', generation)
        require(stage['returncode'] == 0, f'{label}: mesher failed; see retained logs')
        mesh = Path(str(mesh_prefix) + '.solver.cm2d')
        prefix = directory / 'flow'
        command = [str(flow_cli), '--mesh', str(mesh), '--output', str(prefix),
                   '--case', 'cavity', '--nu', '.01', '--speed', '1',
                   '--convection', 'limited-linear', '--viscous-stress', 'symmetric',
                   '--pressure-preconditioner', 'aggregation', '--tolerance', '1e-8',
                   '--max-iterations', '20000', '--profile']
        started = time.monotonic()
        actual = native.run(command, root / 'logs' / f'flow-{label}', timeout)
        run = dict(label=label, command=command, mesh=str(mesh), meshSha256=native.sha256_file(mesh),
                   prefix=str(prefix), binarySha256=native.sha256_file(flow_cli),
                   returnCode=actual['returncode'], timedOut=actual['timedOut'],
                   timeoutSeconds=timeout, elapsedSeconds=time.monotonic() - started,
                   stdout=actual['stdout'], stderr=actual['stderr'])
        runs.append(run)
        native.write_json(root / 'runs.json', runs)
        print(json.dumps({key: run[key] for key in ('label', 'returnCode', 'timedOut', 'elapsedSeconds')}), flush=True)
        require(actual['returncode'] == 0, f'{label}: solve failed; see retained logs and manifest')
    return root / 'runs.json'


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--runs', type=Path, help='Read an existing manifest without launching solvers')
    mode.add_argument('--generate', type=Path, help='Generate and solve in a NEW output directory')
    parser.add_argument('--levels', nargs='+', type=int, default=[4, 5, 6])
    parser.add_argument('--timeout', type=float, default=360., help='Wall seconds per process')
    parser.add_argument('--mesh-cli', type=Path, default=Path('build/cartmesh2d_cli'))
    parser.add_argument('--flow-cli', type=Path, default=Path('build/cartmesh2d_flow_cli'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.generate and args.generate.exists():
        parser.error('Generation directory already exists; use a new directory to preserve evidence')
    try:
        runs = generate_runs(args.generate, args.levels, args.timeout, args.mesh_cli, args.flow_cli) if args.generate else args.runs
        report = audit_runs(runs)
    except (OSError, ValueError, KeyError, TypeError, ArithmeticError) as error:
        report = dict(format='cartmesh2d-cavity-reference-audit-v1', independentAuditsPassed=False,
                      allChecksPassed=False, issues=[str(error)])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({**{key: report[key] for key in ('independentAuditsPassed', 'allChecksPassed', 'issues')},
                      'ghiaSequencePassed': report.get('ghiaSequence', {}).get('valid'),
                      'publishedSequencePassed': report.get('publishedSequence', {}).get('valid')}))
    return 0 if report['allChecksPassed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
