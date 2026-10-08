#!/usr/bin/env python3
"""Observe an actual production Stokes SIMPLE step and attribute its response.

Linux/GNU ld research tool, linked against existing production archives. It
captures actual matrix and Rhie--Chow calls; all PDE operations are native C++.
No product/default control is changed. Results describe an observed mode, not
physical accuracy or a complete spectral bound.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mesh', type=Path)
    parser.add_argument('boundaries', type=Path)
    parser.add_argument('--output', type=Path, default=Path('outputs/cell-step'))
    parser.add_argument('--build-dir', type=Path, default=Path('build'))
    parser.add_argument('--target', type=int, default=3710,
                        help='zero-based cell to observe; default is retained annulus cell 3711 in one-based notation')
    parser.add_argument('--power-steps', type=int, default=180)
    args = parser.parse_args()
    if sys.platform != 'linux':
        parser.error('This read-only wrapper probe requires Linux with GNU ld --wrap.')
    if args.target < 0 or not 1 <= args.power_steps <= 100000:
        parser.error('target must be nonnegative and power-steps must be 1..100000')
    root = Path(__file__).resolve().parents[2]
    build = args.build_dir.resolve()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    result_path = out / 'result.json'
    result_path.write_text(json.dumps({'status': 'building_diagnostic',
                                      'target0': args.target, 'powerSteps': args.power_steps}, indent=2) + '\n')
    libraries = [build / 'libcartmesh2d_fv.a', build / 'libcartmesh2d.a']
    cache = (build / 'CMakeCache.txt').read_text()
    compiler = next(line.split('=', 1)[1] for line in cache.splitlines()
                    if line.startswith('CMAKE_CXX_COMPILER:'))
    symbols_text = subprocess.check_output(['nm', '-g', '--defined-only', str(libraries[0])], text=True)
    symbols = {}
    for name in ('momentum', 'rhieChowFlux'):
        matches = [line.split()[-1] for line in symbols_text.splitlines()
                   if ' T ' in line and name in line]
        if len(matches) != 1:
            raise RuntimeError(f'Expected one production {name} symbol, got {matches}')
        symbols[name] = matches[0]
    (out / 'symbols.hpp').write_text(''.join(
        f'#define {name.upper()}_SYMBOL "{symbol}"\n' for name, symbol in symbols.items()))
    build_commands = []
    for binary, source, wrapped in [('cell_step', 'cell_step.cpp', True),
                                    ('cell_step_plain', 'cell_step.cpp', False),
                                    ('cell_step_response', 'cell_step_response.cpp', True)]:
        command = [compiler, '-O3', '-DNDEBUG', '-std=gnu++20',
                   '-I' + str(root / 'include'), '-I' + str(root / 'src/fv'),
                   '-I' + str(out), str(root / 'tools/benchmarks' / source)]
        if wrapped:
            command += ['-DDIAG_WRAP'] + ['-Wl,--wrap=' + symbol for symbol in symbols.values()]
        command += [str(path) for path in libraries] + ['-o', str(out / binary)]
        compiled = subprocess.run(command, text=True, capture_output=True)
        if compiled.returncode:
            result_path.write_text(json.dumps({'status': 'failed', 'stage': 'compile',
                                               'command': command, 'returnCode': compiled.returncode,
                                               'stdout': compiled.stdout, 'stderr': compiled.stderr}, indent=2) + '\n')
        compiled.check_returncode()
        build_commands.append(command)
    before = {str(path): sha(path) for path in libraries}
    runs = []

    def run(name, command):
        started = time.monotonic()
        result = subprocess.run(command, cwd=root, text=True, capture_output=True)
        elapsed = time.monotonic() - started
        (out / (name + '.log')).write_text(result.stdout + result.stderr)
        runs.append({'name': name, 'command': command, 'returnCode': result.returncode,
                     'elapsedWallSecondsIncludingIO': elapsed, 'stdout': result.stdout, 'stderr': result.stderr})
        result_path.write_text(json.dumps({'status': 'failed' if result.returncode else 'running_diagnostic',
                                           'runs': runs}, indent=2) + '\n')
        result.check_returncode()

    mesh, boundary = args.mesh.resolve(), args.boundaries.resolve()
    for name, binary in [('trace', 'cell_step'), ('plain', 'cell_step_plain')]:
        run(name, [str(out / binary), str(mesh), str(boundary), str(out / name),
                   str(args.target), str(args.power_steps)])
    names = ['input-modal-cells.csv', 'input-modal-faces.csv',
             'output-modal-cells.csv', 'output-modal-faces.csv',
             'output-base-cells.csv', 'output-base-faces.csv',
             'output-delta-cells.csv', 'output-delta-faces.csv', 'power.csv']
    invariant = {name: (out / 'trace' / name).read_bytes() ==
                       (out / 'plain' / name).read_bytes() for name in names}
    if not all(invariant.values()):
        result_path.write_text(json.dumps({'status': 'failed', 'stage': 'wrapper_invariance',
                                           'comparison': invariant, 'runs': runs}, indent=2) + '\n')
        raise RuntimeError(f'Wrapper changed native outputs: {invariant}')
    run('response', [str(out / 'cell_step_response'), str(mesh), str(boundary),
                     str(out / 'trace'), str(out / 'response')])
    replay = {}
    for captured, original in [('replayed-modal', 'output-modal'), ('replayed-base', 'output-base')]:
        for kind in ('cells', 'faces'):
            name = captured + '-' + kind + '.csv'
            replay[name] = ((out / 'response' / name).read_bytes() ==
                            (out / 'trace' / (original + '-' + kind + '.csv')).read_bytes())
    if not all(replay.values()):
        result_path.write_text(json.dumps({'status': 'failed', 'stage': 'source_replay',
                                           'comparison': replay, 'runs': runs}, indent=2) + '\n')
        raise RuntimeError(f'Source replay differs from captured production step: {replay}')
    after = {str(path): sha(path) for path in libraries}
    if before != after:
        result_path.write_text(json.dumps({'status': 'failed', 'stage': 'production_archive_change',
                                           'before': before, 'after': after, 'runs': runs}, indent=2) + '\n')
        raise RuntimeError('Production archives changed during diagnostic')
    result = {
        'status': 'diagnostic_complete',
        'scope': 'Read-only native Stokes SIMPLE trace; closed fixed-velocity walls, constant viscosity. No product discretization/default change.',
        'sourceHead': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
        'inputs': [{'path': str(path), 'sha256': sha(path)} for path in (mesh, boundary)],
        'diagnosticSourceSha256': {str(path.relative_to(root)): sha(path) for path in
                                  [Path(__file__), root / 'tools/benchmarks/cell_step.cpp',
                                   root / 'tools/benchmarks/cell_step_response.cpp']},
        'buildCommands': build_commands, 'wrappedSymbols': symbols, 'runs': runs,
        'productionLibrariesSha256Before': before, 'productionLibrariesSha256After': after,
        'wrappedUnwrappedByteEqual': invariant, 'sourceReplayByteEqual': replay,
        'traceSummary': read_rows(out / 'trace/summary.csv'),
        'stageProjections': read_rows(out / 'trace/stage-projection.csv'),
        'traceAlignment': read_rows(out / 'trace/alignment.csv'),
        'sourceProjections': read_rows(out / 'response/source-projection.csv'),
        'sourceClosure': read_rows(out / 'response/response-closure.csv'),
        'projectionWeights': 'Cell area separately for velocity and pressure; face length for flux; each divided by its corresponding input squared norm.',
        'limitations': 'Each source is propagated through the same actual compact matrix and native RC/pressure correction. Their sum is numerical with finite Krylov stopping error. Large contribution is not proof of a unique instability cause or an ablation result. Mode fields have arbitrary normalization and are not physical solution fields.'
    }
    result_path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'result': str(out / 'result.json'), 'stageProjections': result['stageProjections'],
                      'sourceClosure': result['sourceClosure']}, indent=2))


if __name__ == '__main__':
    main()
