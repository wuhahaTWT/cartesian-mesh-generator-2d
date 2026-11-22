"""Plot actual native meshes and collect existing native/analytic diagnostics.

No PDE assembly or independent acceptance gate. Run the commands documented in
DEVELOPMENT_CN.md first; failed/coarse attempts remain in the same output folder.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.font_manager import FontProperties

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'outputs/shape-agglomeration'
OUT = ROOT / 'artifacts/current'
sys.path.insert(0, str(ROOT / 'tools/flow'))
import native_mesh as nm
from annulus_cross import errors


def read(path):
    return json.loads(path.read_text())


def reference(path):
    return dict(path=str(path.relative_to(ROOT)),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def flow(prefix):
    summary = read(BASE / (prefix + '.json'))
    performance = read(BASE / (prefix + '.performance.json'))
    with (BASE / (prefix + '.cells.csv')).open() as source:
        rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(source)]
    return dict(summary=reference(BASE / (prefix + '.json')),
                fields=reference(BASE / (prefix + '.fields.json')),
                controls={k: summary[k] for k in ['nu', 'speed', 'convection', 'coupling',
                    'velocityRelaxation', 'pressureRelaxation', 'pressureCorrectionPasses',
                    'convergenceMode', 'adaptiveLinear', 'tolerance', 'steadyAcceleration']},
                acceptance={k: summary[k] for k in ['converged', 'iterations', 'momentumResidual',
                    'continuity', 'velocityChange', 'pressureChange', 'globalRelativeImbalance',
                    'strictLinearFinal']},
                performance=performance, circularReferenceDiagnostics=errors(rows))


def main():
    font = Path('/Library/Fonts/Arial Unicode.ttf')
    if font.exists():
        plt.rcParams['font.family'] = FontProperties(fname=str(font)).get_name()
    plt.rcParams.update({'axes.unicode_minus': False, 'font.size': 11})
    fig, axes = plt.subplots(2, 2, figsize=(11, 10), layout='constrained')
    records = []
    for column, (prefix, label) in enumerate([
            ('legacy160', '旧处理：小单元合并＋凸划分'),
            ('final160', '新处理：按形状聚合，保留物面')]):
        path = BASE / (prefix + '.solver.cm2d')
        mesh = nm.read_cm2d(path)
        quality = read(BASE / (prefix + '.quality.json'))
        polygons = [[mesh.vertices[v] for v in c.vertices] for c in mesh.cells]
        walls = [(mesh.vertices[e.v0], mesh.vertices[e.v1]) for e in mesh.edges if e.neighbour < 0]
        for row in range(2):
            ax = axes[row, column]
            ax.add_collection(PolyCollection(polygons, facecolors='#e4eff5',
                edgecolors='#516d7b', linewidths=.25 if row == 0 else .8))
            ax.add_collection(LineCollection(walls, colors='#a83236', linewidths=1.1))
            ax.set_aspect('equal')
            ax.set_xlabel('x / m'); ax.set_ylabel('y / m')
            if row == 0:
                ax.set_xlim(-1.04, 1.04); ax.set_ylim(-1.04, 1.04)
                ax.set_title(f'{label}\n{len(mesh.cells):,} 个单元')
                box = plt.Rectangle((-.55, -.45), .30, .30, fill=False,
                                    edgecolor='#d78630', linewidth=1.5)
                ax.add_patch(box)
            else:
                ax.set_xlim(-.55, -.25); ax.set_ylim(-.45, -.15)
                m = quality['metrics']
                ax.set_title(f"最小角 {m['min_interior_angle_deg']:.2f}°；形状比 {m['max_cell_aspect']:.2f}\n"
                             f"最大内部非正交 {m['max_non_orthogonality_deg']:.2f}°")
        records.append(dict(mesh=reference(path), nativeQualityReport=quality))
    fig.suptitle('同一物面、同一 160×160 背景格：真实网格对照\n红线为原物面；原始 XY 未平滑', fontsize=16)
    fig.savefig(OUT / 'shape-agglomeration.png', dpi=160)
    plt.close(fig)

    generations = []
    for n in [128, 160, 192, 256]:
        prefix = BASE / f'final{n}'
        q = read(Path(str(prefix) + '.shape-quality.json'))
        mesh = Path(str(prefix) + ('.solver.cm2d' if q['valid'] else '.failed.solver.cm2d'))
        generations.append(dict(n=n, mesh=reference(mesh),
            nativeTargetReport=q, metrics=read(Path(str(prefix) + '.metrics.json')),
            nativeRun=read(Path(str(prefix) + '.shape-run.json'))))
    old_raw = BASE / 'legacy160.raw.cm2d'
    new_raw = BASE / 'final160.raw.cm2d'
    solved = BASE / 'endpoints160.solver.cm2d'
    final = BASE / 'final160.solver.cm2d'
    # Provenance checks, not new numerical acceptance tolerances.
    assert old_raw.read_bytes() == new_raw.read_bytes()
    assert solved.read_bytes() == final.read_bytes()
    evidence = dict(format='native-shape-agglomeration-v1', date='2026-10-31',
        baseCommit='0c7e82f', branch='codex/cfd-development',
        implementation='opt-in --shape-agglomeration; native exact unions before convex partition',
        defaultsChanged=False,
        targets=dict(minInteriorAngleDeg=20, maxNativeHydraulicAspect=10,
                     maxInternalNonOrthogonalityDeg=45, maximumAreaRatio=4,
                     maximumPatchCells=6, maximumPasses=64),
        geometry=reference(ROOT / 'outputs/annulus-cross/original-physical-walls.xy'),
        sourceFiles=[reference(ROOT / p) for p in [
            'src/quality/SolverTopology2D.cpp', 'include/cartmesh2d/quality/SolverTopology2D.hpp',
            'apps/cartmesh2d_cli.cpp', 'tools/benchmarks/annulus_mesh_cross.cpp']],
        consumerBinary=reference(ROOT / 'build/cartmesh2d_flow_cli'),
        identicalRaw160Meshes=True, solvedMeshMatchesFinalImplementation=True,
        generationCommands=read(BASE / 'generation-runs.json'), generations=generations,
        paired160=dict(meshes=records, old=flow('legacy160-simple'),
                       new=flow('endpoints160-simple')),
        nativeTests=dict(log='outputs/shape-agglomeration/ctest-final.log', passed=4,
            scope='shape unions, existing solver export, topology, small-cell agglomeration'),
        productCLI=dict(runs=read(BASE / 'cli-runs.json'), fine='outputs/shape-agglomeration/review-product256.log',
            fineCells=23136, coarseFailure='outputs/shape-agglomeration/review-product128-failed.log',
            otherInputs=['rectangle: passed', 'circle: passed',
                         'concave L: rejected; five physical 270-degree corners exceed retained concavity limit']),
        notQualified=['Only this annular CFD comparison; no mesh independence or general stability qualification.',
            'Circular analytic reference differs from polygon geometry/facet wall motion; errors include model approximation.',
            'Native shape targets are engineering hypotheses; all other legacy checks remain active.',
            'Docker daemon unavailable: no new external checkMesh or mature SIMPLE run.',
            'No desktop UI/runtime update, full regression, platform CI or installer build.',
            'The 128/192 grids fail the target and are retained as failed inputs; refinement is not monotonic.',
            'The 23136-cell SIMPLE trial was deliberately terminated for cost after a prefix, not counted as converged.'])
    archive = BASE / 'representative-fields.tar.gz'
    if archive.exists():
        evidence['representativeArchive'] = reference(archive)
    (OUT / 'shape-agglomeration.json').write_text(json.dumps(evidence, indent=2, ensure_ascii=False) + '\n')
    print(OUT / 'shape-agglomeration.png')


if __name__ == '__main__':
    main()
