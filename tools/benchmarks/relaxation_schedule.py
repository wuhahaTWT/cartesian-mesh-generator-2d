#!/usr/bin/env python3
"""Experimental residual-feedback schedules through the unchanged native CLI.

Each segment reloads u/v/p AND owner-oriented face flux on the same mesh.
Restart overhead is included; fixed-factor segmented controls are required.
This is not an in-core adaptive solver or a physical-time checkpoint.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import time

import relaxation_study as study


def merit(summary, tolerance):
    values = [summary.get('momentumResidual', math.inf) / tolerance,
              summary.get('continuity', math.inf) / 1e-8,
              summary.get('globalRelativeImbalance', math.inf) / 1e-8]
    if 'fluxConsistency' in summary:
        values.append(summary['fluxConsistency'] / min(tolerance, 1e-8))
    return max(values)


def checkpoint(root, case, record, target):
    prefix = root / record['prefix']
    mesh = study.native_mesh.read_cm2d(root / case['mesh'])
    cells, faces = study.iterate.extract(mesh, study.iterate.rows(Path(str(prefix)+'.cells.csv')),
                                        study.iterate.rows(Path(str(prefix)+'.faces.csv')))
    for ext, header, data in [('.cells.csv', 'cell,x,y,u,v,p', cells),
                              ('.flux.csv', 'face,owner,neighbour,x,y,flux', faces)]:
        with Path(str(target)+ext).open('w') as stream:
            writer = csv.writer(stream)
            writer.writerow(header.split(','))
            writer.writerows(data)


def run(root, manifest, batch, case_name, method, policy, pair, chunk, budget):
    directory = root / 'schedules' / batch / f'{case_name}__{method}__{policy}'
    directory.mkdir(parents=True, exist_ok=True)
    settings = dict(case=case_name, method=method, policy=policy, initialPair=pair,
                    chunk=chunk, budget=budget, increase=1.25, decrease=.5, minimum=.05,
                    rejectGrowth=2., growthThreshold=.8,
                    rule='accept finite completed segments with merit <= 2*previous; reject otherwise; '
                         'increase both factors by 1.25 up to 1 if accepted merit ratio < .8; '
                         'halve both down to .05 on rejection; no reference field used by controller')
    study.freeze(directory / 'spec.json', settings)
    output = directory / 'record.json'
    if output.exists():
        return json.loads(output.read_text())
    # A restart during a schedule would make its elapsed time incomparable.
    if (directory / 'started.json').exists():
        raise RuntimeError(f'Interrupted schedule preserved: {directory}; use another batch')
    study.write(directory / 'started.json', dict(time=time.time()))
    started = time.perf_counter()
    steps = []; initial = None; previous = None; accepted = None
    used_iterations = 0; status = 'iteration-limit'
    while used_iterations < budget['iterations']:
        remaining = budget['seconds'] - (time.perf_counter()-started)
        if remaining <= 0:
            status = 'time-limit'; break
        index = len(steps)
        record = study.run_one(root, manifest, batch, case_name, method, pair, 0,
                               dict(iterations=min(chunk, budget['iterations']-used_iterations), seconds=remaining),
                               manifest['cpus'][0], initial=initial,
                               label=f'{case_name}__{method}__{policy}__s{index:03d}')
        iterations = record['summary'].get('iterations')
        # Failed native exceptions have no final iteration count. Charge a full
        # segment to the iteration budget and retain that uncertainty explicitly.
        used_iterations += iterations if iterations is not None else chunk
        value = merit(record['summary'], manifest['controls']['tolerance'])
        finite = math.isfinite(value)
        can_continue = record['status'] == 'iteration-limit' and not record['summary'].get('failureReason')
        accept = record['status'] == 'eligible' or (can_continue and finite and
                 (policy == 'fixed' or previous is None or value <= 2*previous))
        ratio = value/previous if finite and previous else None
        step = dict(run=record['prefix'], pair=list(pair), nativeStatus=record['status'],
                    merit=value if finite else None, ratio=ratio, accepted=accept,
                    chargedIterations=iterations if iterations is not None else chunk,
                    iterationCountKnown=iterations is not None, wallSeconds=record['wallSeconds'])
        steps.append(step)
        if record['status'] == 'eligible':
            accepted = record; status = 'eligible'; break
        if record['status'] in ['time-limit', 'different-solution', 'uncertified']:
            status = record['status']; break
        if accept:
            accepted = record
            target = directory / f'initial-{index:03d}'
            checkpoint(root, manifest['cases'][case_name], record, target)
            initial = target; previous = value
            if policy == 'adaptive' and ratio is not None and ratio < .8:
                pair = [min(1., a*1.25) for a in pair]
        elif policy == 'adaptive':
            next_pair = [max(.05, a*.5) for a in pair]
            if next_pair == pair:
                status = 'controller-stalled'; break
            pair = next_pair
        else:
            status = record['status'] if not can_continue else 'residual-growth'; break
        study.write(directory / 'progress.json', dict(steps=steps, chargedIterations=used_iterations))
    result = dict(**settings, status=status, steps=steps, chargedIterations=used_iterations,
                  wallSeconds=time.perf_counter()-started,
                  finalRun=None if accepted is None else accepted['prefix'],
                  agreement=None if accepted is None else accepted['agreement'],
                  costScope='all segment processes, field comparison, CSV extraction/reload and controller overhead',
                  controllerSha256=study.sha(__file__))
    study.write(output, result)
    print(json.dumps({k: result[k] for k in ['case', 'method', 'policy', 'status', 'wallSeconds', 'chargedIterations']}), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=study.ROOT / 'outputs/relaxation-study')
    parser.add_argument('--batch', required=True)
    parser.add_argument('--cases', nargs='+', required=True)
    parser.add_argument('--method', choices=study.METHODS, required=True)
    parser.add_argument('--policy', choices=['fixed', 'adaptive'], required=True)
    parser.add_argument('--pair', default='.6,.25')
    parser.add_argument('--chunk', type=int, default=50)
    parser.add_argument('--iterations', type=int, default=3000)
    parser.add_argument('--seconds', type=float, default=600)
    args = parser.parse_args()
    pair = list(map(float, args.pair.split(',')))
    if len(pair) != 2 or any(not 0 < a <= 1 for a in pair) or args.chunk < 10:
        parser.error('Need two factors in (0,1] and chunk >= 10')
    root = args.root.resolve(); manifest = study.load(root)
    for case in args.cases:
        run(root, manifest, args.batch, case, args.method, args.policy, pair.copy(), args.chunk,
            dict(iterations=args.iterations, seconds=args.seconds))


if __name__ == '__main__':
    main()
