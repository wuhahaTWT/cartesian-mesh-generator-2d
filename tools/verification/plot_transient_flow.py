#!/usr/bin/env python3
"""Create compact evidence and a plot from the completed transient matrix.

Uses existing final meshes/CSV, verifies recorded hashes, and reports errors
without imposing a universal accuracy threshold. Does not run a solver.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
import verify_transient_flow as verifier


def rows(path):
    with path.open() as stream:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(stream)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='output prefix')
    args = parser.parse_args()
    groups = {}
    provenance = {}
    for folder in ['time', *[f'{kind}-l{level}' for kind in ('cartesian', 'warped') for level in (4,5,6)],
                   'startup-channel', 'startup-cavity', 'startup-external']:
        source = args.root / folder / 'runs.json'
        runset = json.loads(source.read_text())
        entries = []
        for run in runset['runs']:
            if not run['valid'] or run['returnCode'] != 0:
                raise ValueError(f'incomplete matrix: {source}')
            command = run['command']
            prefix = Path(command[command.index('--output')+1])
            mesh = Path(command[command.index('--mesh')+1])
            if verifier.native.sha256_file(mesh) != runset['meshSha256']:
                raise ValueError(f'changed mesh: {mesh}')
            report_path = Path(str(prefix)+'.audit.json')
            report = json.loads(report_path.read_text())
            if not report['valid']:
                raise ValueError(f'invalid report: {report_path}')
            for file, digest in report['sha256'].items():
                if verifier.native.sha256_file(Path(file)) != digest:
                    raise ValueError(f'changed artifact: {file}')
            for name in ('nu','speed','tolerance'):
                if report['controls'][name] != float(command[command.index('--'+name)+1]):
                    raise ValueError(f'invocation mismatch: {name}')
            provenance[str(source)] = verifier.native.sha256_file(source)
            provenance[str(report_path)] = verifier.native.sha256_file(report_path)
            entries.append(dict(run, prefix=str(prefix), mesh=str(mesh), binarySha256=runset['binarySha256'], meshSha256=runset['meshSha256'],
                maxFaceDeviation=max(report['momentumAudit']['maxFaceDeviation'].values()),
                momentumResidualDeviation=report['momentumAudit']['summaryDeviation']['momentumResidual']['absolute'],
                conservationDifference=report['momentumAudit']['cellBoundaryConservationDifference']['magnitude'],
                nativeDefinitionContinuity=report['independentContinuity']['nativeDefinitionContinuity']))
        groups[folder] = entries
    series = groups['time']
    distances = []
    for first, second in zip(series, series[1:]):
        a = rows(Path(first['prefix']+'.cells.csv')); b = rows(Path(second['prefix']+'.cells.csv'))
        distance = math.sqrt(math.fsum(x['area']*((x['u']-y['u'])**2+(x['v']-y['v'])**2) for x,y in zip(a,b)) /
                             math.fsum(x['area'] for x in a))
        distances.append(distance)
    orders = [math.log(a/b,2) for a,b in zip(distances,distances[1:])]
    fig, axes = plt.subplots(2,2,figsize=(12,9),layout='constrained')
    fine = series[-1]
    mesh = verifier.native.read_cm2d(Path(fine['mesh']))
    field = rows(Path(fine['prefix']+'.cells.csv'))
    polygons = [[mesh.vertices[v] for v in c.vertices] for c in mesh.cells]
    artist = PolyCollection(polygons,array=np.array([r['speed'] for r in field]),cmap='viridis',edgecolors='#ffffff50',linewidths=.18)
    axes[0,0].add_collection(artist); axes[0,0].autoscale_view(); axes[0,0].set_aspect('equal')
    selected=field[::3]
    axes[0,0].quiver([r['x'] for r in selected],[r['y'] for r in selected],[r['u'] for r in selected],[r['v'] for r in selected],color='white',scale=16,width=.002)
    axes[0,0].set(title='Actual 900-cell mesh | Taylor-Green at t = 0.2',xlabel='x',ylabel='y')
    fig.colorbar(artist,ax=axes[0,0],label='Velocity magnitude (m/s)',shrink=.8)
    for run in series:
        history=rows(Path(run['prefix']+'.time-history.csv'))
        axes[0,1].plot([0]+[r['time'] for r in history],[.25]+[r['kineticEnergy'] for r in history],label=f"dt = {run['dt']:g}",lw=1.2)
    ts=np.linspace(0,.2,201)
    axes[0,1].plot(ts,.25*np.exp(-4*.1*np.pi**2*ts),'k--',label='Analytic',lw=1)
    axes[0,1].set(title='Physical kinetic-energy decay',xlabel='Physical time (s)',ylabel='Energy per density and depth')
    axes[0,1].legend()
    dt=[r['dt'] for r in series]
    axes[1,0].loglog(dt,[math.hypot(r['analyticErrors']['uL2'],r['analyticErrors']['vL2']) for r in series],'o-',label='Velocity L2 / Uref')
    axes[1,0].loglog(dt,[r['analyticErrors']['pL2'] for r in series],'s-',label='Pressure L2 / Uref^2')
    axes[1,0].set(title='Fixed mesh, fixed final time',xlabel='Time step (s)',ylabel='Analytic error')
    axes[1,0].legend()
    for kind,marker in [('cartesian','o'),('warped','s')]:
        data=[groups[f'{kind}-l{level}'][0] for level in (4,5,6)]
        axes[1,1].loglog([d['counts']['cells'] for d in data],
            [math.hypot(d['analyticErrors']['uL2'],d['analyticErrors']['vL2']) for d in data],marker+'-',label=kind)
    axes[1,1].set(title='Spatial refinement | dt = 0.001, t = 0.05',xlabel='Final cell count',ylabel='Velocity L2 / Uref (includes time error)')
    axes[1,1].legend()
    for ax in axes.flat:
        if ax is not axes[0,0]: ax.grid(alpha=.22)
    fig.suptitle('Native 2D transient FVM | backward Euler, limited-linear convection',fontsize=14)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(str(args.output)+'.png',dpi=170)
    plt.close(fig)
    report={'scope':'CLI/core transient milestone; not engineering certification or desktop delivery',
            'groups':groups,'sourceReportsSha256':provenance,'timeSelfDistances':distances,'timeSelfObservedOrders':orders,
            'limits':['first-order fixed physical time step; no automatic CFL controller',
                      'spatial refinement includes finite dt error; no spatial-order certification',
                      'startup checks are discrete balance, not shedding frequency or force accuracy',
                      'no new external OpenFOAM or Linux/Windows execution',
                      'previous Ghia limited-linear trend failure remains unresolved']}
    validation_path = args.root/'validation.json'
    if validation_path.is_file():
        report['validation'] = json.loads(validation_path.read_text())
        for path, digest in report['validation']['sourceSha256'].items():
            if verifier.native.sha256_file(Path(path)) != digest:
                raise ValueError(f'changed validated source: {path}')
        expected = report['validation']['binarySha256']
        if any(r['binarySha256'] != expected for group in groups.values() for r in group):
            raise ValueError('matrix has mixed binaries')
        report['sourceReportsSha256'][str(validation_path)] = verifier.native.sha256_file(validation_path)
    verifier.native.write_json(Path(str(args.output)+'.json'),report)
    print(json.dumps({'cases':sum(len(v) for v in groups.values()),'timeSelfObservedOrders':orders}))


if __name__=='__main__':
    main()
