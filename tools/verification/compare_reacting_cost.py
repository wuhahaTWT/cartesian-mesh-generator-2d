#!/usr/bin/env python3
"""Compare complete native hot-spot runs at unchanged physics and controls.

Intended for query/cache optimizations: all saved fields, accepted step logs,
mechanisms and checkpoints must match exactly, excluding wall-clock timing.
Includes process startup, mechanism preparation, mesh, integration and output.
The prescribed 48-cell reacting fixture is not a general performance benchmark.
"""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-probe', type=Path, required=True)
    parser.add_argument('--candidate-probe', type=Path, required=True)
    parser.add_argument('--mechanism', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repetitions', type=int, default=3)
    args = parser.parse_args()
    if not 1 <= args.repetitions <= 10:
        raise ValueError('repetitions must be between 1 and 10')
    args.output.mkdir(parents=True, exist_ok=False)
    probes = {'baseline': args.baseline_probe.resolve(), 'candidate': args.candidate_probe.resolve()}
    hashes = {name: sha(path) for name, path in probes.items()}
    mechanism_hash = sha(args.mechanism)
    results, canonical = [], None
    started = time.perf_counter()
    for repetition in range(args.repetitions):
        order = ['baseline', 'candidate'] if repetition % 2 == 0 else ['candidate', 'baseline']
        for label in order:
            path = args.output / f'{label}-{repetition}'
            assert sha(probes[label]) == hashes[label] and sha(args.mechanism) == mechanism_hash
            begin = time.perf_counter()
            process = subprocess.run([str(probes[label]), str(args.mechanism.resolve()), str(path)],
                                     text=True, capture_output=True)
            elapsed = time.perf_counter() - begin
            (args.output / f'{label}-{repetition}.stdout').write_text(process.stdout)
            (args.output / f'{label}-{repetition}.stderr').write_text(process.stderr)
            record = {'label': label, 'repetition': repetition, 'returncode': process.returncode,
                      'complete_cost_seconds': elapsed, 'directory': str(path), 'exact_match': False}
            if process.returncode == 0:
                field = json.loads((path / 'field.json').read_text())
                record.update(complete=field['complete'], physical_time=field['frames'][-1]['time'],
                              steps=field['frames'][-1]['steps'], rejections=field['rejections'])
                field.pop('elapsedSeconds')
                content = (field, (path / 'steps.jsonl').read_bytes(),
                           (path / 'accepted.checkpoint').read_bytes(), (path / 'resolved-mechanism.yaml').read_bytes())
                if canonical is None:
                    canonical = content
                record['exact_match'] = content == canonical and field['complete']
                record['files_sha256'] = {name: sha(path / name) for name in
                    ['field.json', 'steps.jsonl', 'accepted.checkpoint', 'resolved-mechanism.yaml']}
            results.append(record)
            print(json.dumps(record), flush=True)
    assert all(sha(path) == hashes[name] for name, path in probes.items())
    assert sha(args.mechanism) == mechanism_hash
    medians = {name: statistics.median(r['complete_cost_seconds'] for r in results if r['label'] == name)
               for name in probes}
    report = {'passed': all(r['returncode'] == 0 and r['exact_match'] for r in results),
              'repetitions': args.repetitions, 'runs': results, 'probe_sha256': hashes,
              'mechanism_sha256': mechanism_hash, 'median_seconds': medians,
              'median_speedup': medians['baseline'] / medians['candidate'],
              'experiment_wall_seconds': time.perf_counter() - started,
              'scope': '48-cell warped sealed detailed H2/air hot spot, 40 us, default unchanged controls; timing includes startup through output',
              'field_comparison': 'exact parsed fields plus byte-identical step log, mechanism and checkpoint; wall time excluded',
              'general_performance_qualification': False, 'flame_physical_qualification': False}
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': report['passed'], 'median_seconds': medians, 'median_speedup': report['median_speedup']}))
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
