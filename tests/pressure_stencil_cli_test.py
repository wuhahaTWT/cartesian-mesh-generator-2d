#!/usr/bin/env python3
"""Coarse shifted-circle regression: tangential chains need a wider pressure stencil."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools/verification'))
import verify_native_flow as native
import verify_transient_flow as transient


def run(command):
    result = subprocess.run(list(map(str, command)), text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


def main(cli, mesh_cli):
    with tempfile.TemporaryDirectory(prefix='cartmesh-pressure-stencil-') as folder:
        root = Path(folder)
        mesh = root/'mesh.solver.cm2d'
        run([mesh_cli, ROOT/'examples/acceptance/circle.xy', root/'mesh', 8, .25, .15,
             'exterior', root/'foam', 0, 0, '--size-field', '--far-field-spans', 10,
             '--wall-relative-size', .3281250000003282, '--background-relative-size', 1.3125000000013127,
             '--cells-per-level', 22, '--max-safe-wall-level', 11])
        assert len(native.read_cm2d(mesh).cells) == 2504

        def solve(prefix, scheme='face-limited-linear', extra=()):
            run([cli, '--mesh', mesh, '--output', prefix, '--case', 'external', '--nu', .1,
                 '--speed', 1, '--convection', scheme, '--pressure-preconditioner', 'aggregation',
                 '--max-iterations', 1500, '--tolerance', 1e-8, *extra])

        for scheme in ['upwind', 'face-limited-linear']:
            prefix = root/scheme
            solve(prefix, scheme)
            audit = native.verify_case(mesh, prefix, 'external', .1, 1, native.argument_parser().parse_args([]))
            assert audit['valid'], audit['issues']
            summary_path = Path(str(prefix)+'.json')
            raw = summary_path.read_text()
            summary = json.loads(raw)
            assert summary['pressureBoundaryReconstruction'] == 'one-sided-linear-adaptive'
            assert summary['outletBackflow'] == 'reject' and summary['outletBackflowFaces'] == 0
            # Audit must reconstruct the actual wider stencil, not simply
            # trust converged=true or silently apply it to legacy summaries.
            summary['pressureBoundaryReconstruction'] = 'one-sided-linear-2ring'
            summary_path.write_text(json.dumps(summary))
            wrong = native.verify_case(mesh, prefix, 'external', .1, 1, native.argument_parser().parse_args([]))
            assert not wrong['valid'], 'incorrect pressure operator metadata accepted'
            summary_path.write_text(raw)

        whole, first, resumed = (root/name for name in ['whole', 'first', 'resumed'])
        solve(whole, extra=['--time-step', .02, '--steps', 2])
        solve(first, extra=['--time-step', .02, '--steps', 1])
        solve(resumed, extra=['--time-step', .02, '--steps', 1, '--restart', str(first)+'.checkpoint'])
        assert Path(str(whole)+'.checkpoint').read_bytes() == Path(str(resumed)+'.checkpoint').read_bytes()
        for prefix in [whole, resumed]:
            result = transient.verify(mesh, prefix, root/(prefix.name+'-audit.json'))
            assert result['valid'], result
        print(json.dumps({'valid': True, 'cells': 2504, 'steadySchemes': 2,
                          'restartIdentical': True, 'wrongStencilRejected': True}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cli', type=Path, required=True)
    parser.add_argument('--mesh-cli', type=Path, required=True)
    args = parser.parse_args()
    main(args.cli.resolve(), args.mesh_cli.resolve())
