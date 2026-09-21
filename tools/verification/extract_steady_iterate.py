#!/usr/bin/env python3
"""Extract a finite same-mesh steady iterate, including owner-oriented flux.

This does not certify convergence or create a physical-time checkpoint. The
solver must enforce its prescribed boundary fluxes and every acceptance gate.
Compressed source CSVs are supported, but the original solver mesh is required.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
import time
from pathlib import Path
from verify_native_flow import read_cm2d, measure


def source_file(prefix, suffix):
    path = Path(str(prefix) + suffix)
    if not path.exists():
        path = Path(str(path) + '.gz')
    return path


def rows(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as stream:
        return list(csv.DictReader(stream))


def extract(mesh, cell_rows, face_rows):
    geometry = measure(mesh, 1e-12, 1e-9)
    if geometry.issues:
        raise ValueError('Source mesh geometry is invalid: ' + '; '.join(geometry.issues))
    if len(cell_rows) != len(mesh.cells) or len(face_rows) != len(mesh.edges):
        raise ValueError('Source fields do not cover the original mesh')
    cells, faces = [], []
    for i, (point, row, area) in enumerate(zip(geometry.centroids, cell_rows, geometry.areas)):
        values = [float(row[k]) for k in ['cell', 'x', 'y', 'u', 'v', 'p']]
        if not all(map(math.isfinite, values)) or values[0] != i:
            raise ValueError('Nonfinite or unordered source cells')
        scale = max(abs(point[0]), abs(point[1]), math.sqrt(area))
        if any(abs(values[j + 1] - point[j]) > 32 * sys.float_info.epsilon * scale for j in [0, 1]):
            raise ValueError('Source cell coordinates differ from mesh')
        # Preserve the native CSV coordinates rather than introduce independent
        # polygon-moment rounding in the target-bound initial guess.
        cells.append([row[k] for k in ['cell', 'x', 'y', 'u', 'v', 'p']])
    for edge, row in zip(mesh.edges, face_rows):
        values = [float(row[k]) for k in ['face', 'owner', 'neighbour', 'flux']]
        if not all(map(math.isfinite, values)) or values[:3] != [edge.id, edge.owner, edge.neighbour]:
            raise ValueError('Source face incidence/order or flux is invalid')
        a, b = mesh.vertices[edge.v0], mesh.vertices[edge.v1]
        faces.append([edge.id, edge.owner, edge.neighbour, (a[0] + b[0]) / 2, (a[1] + b[1]) / 2, row['flux']])
    return cells, faces


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ['mesh', 'source-prefix', 'output']:
        parser.add_argument('--' + flag, type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    outputs = [Path(str(args.output) + suffix) for suffix in ['.cells.csv', '.flux.csv', '.json']]
    if any(path.exists() for path in outputs):
        parser.error('Output exists; preserve earlier evidence')
    summary_path = Path(str(args.source_prefix) + '.json')
    summary = json.loads(summary_path.read_text())
    if summary.get('format') != 'cartmesh2d-flow-summary-v1' or 'steadyFaceInterpolation' not in summary:
        parser.error('Source must be a native steady solve, not a physical-time checkpoint')
    cell_path = source_file(args.source_prefix, '.cells.csv')
    face_path = source_file(args.source_prefix, '.faces.csv')
    mesh = read_cm2d(args.mesh)
    if summary.get('cells') != len(mesh.cells):
        parser.error('Source summary does not match mesh size')
    cells, faces = extract(mesh, rows(cell_path), rows(face_path))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for path, header, records in zip(outputs[:2], ['cell,x,y,u,v,p', 'face,owner,neighbour,x,y,flux'], [cells, faces]):
        with path.open('x') as stream:
            writer = csv.writer(stream)
            writer.writerow(header.split(','))
            writer.writerows(records)
    report = {'scope': __doc__, 'sourceConverged': summary.get('converged'),
              'sourceStrictLinearFinal': summary.get('strictLinearFinal'), 'solutionAccepted': False,
              'cells': len(cells), 'faces': len(faces), 'seconds': time.monotonic() - started,
              'inputHashes': {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [args.mesh, summary_path, cell_path, face_path, Path(__file__)]},
              'outputHashes': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs[:2]}}
    outputs[2].write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
