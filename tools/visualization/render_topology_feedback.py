#!/usr/bin/env python3
"""Plot actual accepted Cut-cell fields and the online optimisation trajectory."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools/flow'))
from native_mesh import read_cm2d


def load(root):
    summary = json.loads((root/'summary.json').read_text())
    rows = json.loads((root/'evaluations.json').read_text())['evaluations']
    return summary, {r['id']:r for r in rows}


def native_field(row):
    report = json.loads(Path(row['nativeSummary']).read_text())
    polygons, speed, sources = [], [], []
    for leaf in report.get('components', [report]):
        accepted = next(c for c in leaf['cases'] if c['status'] == 'flow-audited')
        mesh_path = Path(accepted['mesh']['path'])
        csv_path = mesh_path.parent/'flow.cells.csv'
        mesh = read_cm2d(mesh_path)
        with csv_path.open() as stream:
            cells = {int(r['cell']):r for r in csv.DictReader(stream)}
        for c in mesh.cells:
            polygons.append([mesh.vertices[v] for v in c.vertices])
            speed.append(float(cells[c.id]['speed']))
        sources.extend([mesh_path, csv_path])
    return polygons, np.asarray(speed), sources


def evidence(path):
    return dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def run_index(root):
    summary, rows = load(root)
    report = {k:v for k, v in summary.items() if k != 'history'}
    report['directory'] = str(root.resolve())
    report['files'] = [evidence(root/p) for p in ('run.json', 'summary.json', 'evaluations.json')]
    report['rejectedEvaluations'] = [dict(id=r['id'], name=r['name'], status=r['status'],
        directory=r['directory'], issue=r.get('issue')) for r in rows.values() if r['high'] is None]
    report['acceptedDesigns'] = []
    for i in summary['acceptedIds']:
        r = rows[i]
        report['acceptedDesigns'].append(dict(id=i, name=r['name'], J_B=r['porous']['objective'],
            J_T=r['high'], extraction=r['sharp']['extraction'], feedback=r['feedback'],
            material=evidence(Path(r['directory'])/'material.npz'),
            boundary=evidence(Path(r['directory'])/'sharp'/'fluid.xy'),
            nativeSummary=evidence(Path(r['nativeSummary']))))
    report['feedbackSteps'] = [h for h in summary['history'] if h.get('feedback')]
    return report


def run(args):
    os.environ.setdefault('MPLCONFIGDIR', str(args.source/'plot-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    summary, rows = load(args.source)
    chosen = summary['acceptedIds']
    ids = list(dict.fromkeys([chosen[0], chosen[len(chosen)//2], chosen[-1]]))
    fig = plt.figure(figsize=(13.5, 8.5), constrained_layout=True)
    grid = fig.add_gridspec(2, max(3, len(ids)), height_ratios=[1.1, 1])
    fields = [native_field(rows[i]) for i in ids]
    vmax = max(float(values.max()) for _, values, _ in fields)
    sources = [args.source/'summary.json', args.source/'evaluations.json']
    for column, (i, (polygons, values, files)) in enumerate(zip(ids, fields)):
        axis = fig.add_subplot(grid[0, column])
        collection = PolyCollection(polygons, array=values, cmap='viridis', clim=(0, vmax),
                                    edgecolors='#26374a', linewidths=.08, rasterized=True)
        axis.add_collection(collection); axis.autoscale_view(); axis.set_aspect('equal')
        row = rows[i]
        extraction = row['sharp']['extraction']
        label = {'feedback-gradient':'Feedback update', 'initial':'Initial design',
                 'merge-shared-corridor':'Merge/split proposal'}.get(row['name'], 'Accepted material update')
        axis.set_title(f"#{i} {label}\nJ_T={row['high']:.5f}\n"
                       f"components={extraction['components']}; holes={extraction['holes']}", fontsize=10)
        axis.set_xlabel('x'); axis.set_ylabel('y')
        fig.colorbar(collection, ax=axis, label='Native speed', shrink=.7)
        sources.extend(files)
    axis = fig.add_subplot(grid[1, :2])
    evaluated = [r for r in rows.values() if r['high'] is not None]
    axis.scatter([r['id'] for r in evaluated], [r['high'] for r in evaluated],
                 color='#9ba8b8', label='Native candidate', s=26)
    corrected = [r for r in evaluated if r['feedback'] is not None]
    axis.scatter([r['id'] for r in corrected], [r['high'] for r in corrected],
                 color='#df7834', label='New feedback-gradient candidate', marker='s', s=38)
    axis.plot(chosen, [rows[i]['high'] for i in chosen], 'o-', color='#177f83',
              label='Accepted trajectory', linewidth=2)
    if args.baseline:
        base, base_rows = load(args.baseline)
        bid = base['acceptedIds']
        axis.plot(bid, [base_rows[i]['high'] for i in bid], 'x--', color='#9b5da1',
                  label='Same probes, gradient correction disabled')
        sources.extend([args.baseline/'summary.json', args.baseline/'evaluations.json'])
    axis.set_xlabel('Native evaluation index (separate runs)')
    axis.set_ylabel('Dimensionless total-pressure power J_T')
    axis.set_title('Actual native objectives; failed trials remain in the ledger')
    axis.grid(alpha=.2); axis.legend(fontsize=8)
    note = fig.add_subplot(grid[1, 2:]); note.axis('off')
    failed = sum(r['high'] is None for r in rows.values())
    message = [f"Sharp objective: {summary['initialObjective']:.5f} -> {summary['finalObjective']:.5f}",
        f"Observed reduction: {100*summary['relativeReduction']:.3f}%",
        f"Feedback-gradient accepted steps: {summary['feedbackDrivenAcceptedSteps']}",
        f"Accepted connectivity/holes changes: {summary['acceptedTopologyChanges']}",
        f"Native evaluations: {len(rows)}; rejected: {failed}",
        f"Elapsed: {summary['elapsedSeconds']:.1f} s", '',
        'Equal area, fixed ports, Re, mesh settings and solver gates.',
        'One native grid; no grid-independent or global-optimum claim.',
        'Real polygons and CSV values; no field interpolation.']
    note.text(0, .95, '\n'.join(message), va='top', fontsize=10, linespacing=1.8)
    fig.suptitle('Online Cut-cell feedback changes subsequent density designs', fontsize=16)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160); plt.close(fig)
    metadata = dict(source=str(args.source), baseline=str(args.baseline) if args.baseline else None,
        acceptedFieldIds=ids, files=[dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                                   for p in dict.fromkeys(sources)])
    args.output.with_suffix(args.output.suffix+'.json').write_text(json.dumps(metadata, indent=2)+'\n')
    index = dict(schema='cartmesh2d-topology-feedback-results-v1', primary=run_index(args.source),
                 physicalAccuracyQualified=False, preview=evidence(args.output),
                 relatedRuns=[run_index(p) for p in args.related])
    if args.baseline:
        index['baseline'] = run_index(args.baseline)
        current_run = json.loads((args.source/'run.json').read_text())
        baseline_run = json.loads((args.baseline/'run.json').read_text())
        same = all(current_run[k] == baseline_run[k] for k in ('source', 'problem', 'q', 'beta', 'nativeControls'))
        if not same:
            raise ValueError('baseline physical/source/native contracts differ')
        index['comparison'] = dict(samePhysicalAndNativeControls=True,
            additionalRelativeReduction=1-summary['finalObjective']/index['baseline']['finalObjective'],
            additionalNativeEvaluations=summary['nativeEvaluations']-index['baseline']['nativeEvaluations'],
            elapsedDifferenceSeconds=summary['elapsedSeconds']-index['baseline']['elapsedSeconds'],
            scope='Single recorded online run; no general speedup or physical-accuracy claim.')
    args.output.with_suffix('.json').write_text(json.dumps(index, indent=2, allow_nan=False)+'\n')
    print(args.output)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('--baseline', type=Path)
    p.add_argument('--related', type=Path, nargs='*', default=[])
    p.add_argument('--output', type=Path, required=True)
    run(p.parse_args())
