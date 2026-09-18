#!/usr/bin/env python3
"""Re-verify two physical-time ladders and plot real fields and errors.

The input root contains taylor-green/runs.json and circle-re20/runs.json.
No solver is run and no engineering-accuracy threshold is invented here.
"""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
import compare_transient_steps as comparison


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--circle-runs', type=Path, help='explicit combined runset retaining successful runs after a documented retry')
    p.add_argument('--output', type=Path, required=True, help='PNG and JSON prefix')
    a = p.parse_args()
    record = {'scope': 'fixed mesh, time-step comparison; no spatial/domain or shedding qualification', 'groups': {}}
    for name in ('taylor-green', 'circle-re20'):
        folder = a.root/name
        runset = a.circle_runs if name == 'circle-re20' and a.circle_runs else folder/'runs.json'
        compare = comparison.compare(runset, runset.parent/'comparison.json')
        source = json.loads(runset.read_text())
        runs = []
        for run in source['runs']:
            cmd = run['command']; prefix = Path(cmd[cmd.index('--output')+1]); mesh = Path(source['mesh'])
            audit_path = Path(str(prefix)+'.long-audit.json')
            audit = comparison.audit.verify(mesh, prefix, audit_path)
            runs.append({'execution': run, 'audit': audit, 'prefix': str(prefix),
                         'performance': json.loads(Path(str(prefix)+'.performance.json').read_text())})
        record['groups'][name] = {'comparison': compare, 'runs': runs,
            'source': str(runset), 'sourceSha256': comparison.audit.native.sha256_file(runset),
            'priorAttempts': source.get('priorAttempts', [])}
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout='constrained')
    tg = record['groups']['taylor-green']; circle = record['groups']['circle-re20']
    for r in tg['runs']:
        h = r['audit']['analyticDecay']['history']
        axes[0,1].plot([v['time'] for v in h], [100*v['relativeEnergyError'] for v in h], label=f"dt={r['execution']['dt']:g}")
    axes[0,1].set(xlabel='Physical time', ylabel='Energy error relative to analytic (%)', title='Taylor–Green decay: relative error accumulates')
    axes[0,1].axhline(0,color='0.6',linewidth=.6); axes[0,1].legend()
    for r in circle['runs']:
        h = r['audit']['history']
        axes[1,0].plot([v['time'] for v in h], [v['forceX'] for v in h], label=f"dt={r['execution']['dt']:g}")
    axes[1,0].set(xlabel='Physical time', ylabel='Fx / density / depth (log scale)', yscale='log', title='Abrupt cylinder startup: initial peaks are not qualified'); axes[1,0].legend()
    for label,g in [('Taylor–Green',tg),('Cylinder',circle)]:
        d = g['comparison']['adjacentDistances']
        axes[1,1].loglog([x['coarseDt'] for x in d],[x['vectorL2'] for x in d],'-o',label=label)
    axes[1,1].set(xlabel='Coarser dt of each adjacent pair', ylabel='Area-weighted velocity difference', title='Same mesh, same end time: temporal self-convergence'); axes[1,1].legend()
    fine = min(circle['runs'], key=lambda r:r['execution']['dt']); mesh = comparison.audit.native.read_cm2d(Path(circle['comparison']['mesh']))
    with Path(fine['prefix']+'.cells.csv').open() as stream: fields = {int(row['cell']):row for row in csv.DictReader(stream)}
    polygons = [[mesh.vertices[v] for v in cell.vertices] for cell in mesh.cells]
    values = np.array([float(fields[i]['speed']) for i in range(len(polygons))])
    collection = PolyCollection(polygons, array=values, cmap='viridis', edgecolors='#50565b', linewidths=.18)
    axes[0,0].add_collection(collection); axes[0,0].set(xlim=(-3,5),ylim=(-3,3),aspect='equal',xlabel='x',ylabel='y',title=f"Actual Cut-cell field at t={fine['audit']['time']:g} | {len(polygons):,} cells")
    fig.colorbar(collection,ax=axes[0,0],label='Speed')
    for ax in (axes[0,1],axes[1,0],axes[1,1]): ax.grid(True,alpha=.22)
    fig.suptitle('Native transient laminar verification | fixed Euler time steps\nDiscrete balances verified; spatial/domain independence remains unqualified',fontsize=13)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(str(a.output)+'.png',dpi=165);plt.close(fig)
    record['sourceFiles'] = {str(Path(__file__)):comparison.audit.native.sha256_file(Path(__file__)),
                            str(Path(comparison.__file__)):comparison.audit.native.sha256_file(Path(comparison.__file__)),
                            str(Path(comparison.audit.__file__)):comparison.audit.native.sha256_file(Path(comparison.audit.__file__))}
    Path(str(a.output)+'.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
