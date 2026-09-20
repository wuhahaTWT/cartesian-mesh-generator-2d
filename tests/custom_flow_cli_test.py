#!/usr/bin/env python3
"""Independently audit general boundary inputs on real rotated Cut-cell meshes."""
import argparse
import csv
import json
import math
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools/verification'))
import verify_native_flow as audit
import verify_transient_flow as transient


def run(command, success=True):
    result = subprocess.run(list(map(str, command)), capture_output=True, text=True, timeout=25)
    assert (result.returncode == 0) == success, result.stdout + result.stderr
    return result


def read_cells(prefix):
    with Path(str(prefix)+'.cells.csv').open() as stream:
        return list(csv.DictReader(stream))


def main(cli, mesh_cli, convection):
    with tempfile.TemporaryDirectory(prefix='cartmesh-custom-') as name:
        root = Path(name)
        mesh_path = root/'nozzle.solver.cm2d'
        run([mesh_cli, ROOT/'examples/complex/nozzle_profile.xy', root/'nozzle',
             4, 1/30, .1, 'interior', root/'foam', 4, 0])
        bc_path = root/'input.boundaries'
        run([cli, '--mesh', mesh_path, '--case', 'duct', '--speed', 1, '--export-boundaries', bc_path])

        def solve(output, extra=(), mesh=mesh_path, bc=bc_path, success=True):
            return run([cli, '--mesh', mesh, '--output', output, '--case', 'custom', '--boundary', bc,
                        '--nu', .1, '--speed', 1, '--tolerance', 1e-9, '--max-iterations', 500,
                        '--pressure-preconditioner', 'aggregation', '--convection', convection, *extra], success)

        def verify(output, mesh=mesh_path):
            result = audit.verify_case(mesh, output, 'custom', .1, 1, audit.argument_parser().parse_args([]))
            assert result['valid'], result['issues']
            assert result['benchmark']['status'] == 'not-qualified'

        base = root/'base'
        solve(base)
        verify(base)
        original = read_cells(base)
        source = bc_path.read_text()
        # Prescribed pressure reference must enter the independent pressure
        # gradients, face forces and momentum equation, not only metadata.
        shifted = root/'shifted.boundaries'
        lines = source.splitlines()
        for i, line in enumerate(lines):
            if ' pressure-outlet ' in line:
                lines[i] = line.rsplit(' ',1)[0]+' 7.25'
        shifted.write_text('\n'.join(lines)+'\n')
        shift = root/'shift'
        solve(shift, bc=shifted)
        verify(shift)
        for a,b in zip(original,read_cells(shift)):
            assert max(abs(float(a[k])-float(b[k])) for k in ('u','v')) < 1e-9
            assert abs(float(b['p'])-float(a['p'])-7.25) < 1e-8

        # Arbitrary orientation: transform actual polygon coordinates and
        # bound face centres, normals and inlet velocities together.
        angle = .63
        c,s = math.cos(angle),math.sin(angle)
        def rotate(x,y): return c*x-s*y, s*x+c*y
        rotated_mesh = root/'rotated.solver.cm2d'
        lines = mesh_path.read_text().splitlines()
        for i in range(2,2+int(lines[1].split()[1])):
            ident,x,y = lines[i].split()
            x,y = rotate(float(x),float(y))
            lines[i] = f'{ident} {x:.17g} {y:.17g}'
        rotated_mesh.write_text('\n'.join(lines)+'\n')
        rotated_bc = root/'rotated.boundaries'
        lines = source.splitlines()
        for i,line in enumerate(lines):
            if not line.startswith('BOUNDARY '): continue
            t = shlex.split(line)
            for j in (3,5,9):
                x,y = rotate(float(t[j]),float(t[j+1]))
                t[j:j+2] = [f'{x:.17g}',f'{y:.17g}']
            t[8] = json.dumps(t[8])
            lines[i] = ' '.join(t)
        rotated_bc.write_text('\n'.join(lines)+'\n')
        rotated = root/'rotated'
        solve(rotated, mesh=rotated_mesh, bc=rotated_bc)
        verify(rotated, rotated_mesh)
        for a,b in zip(original,read_cells(rotated)):
            u,v = rotate(float(a['u']),float(a['v']))
            assert math.hypot(u-float(b['u']),v-float(b['v'])) < 1e-8
            assert abs(float(a['p'])-float(b['p'])) < 1e-7
        rejection = solve(root/'wrong-mesh', mesh=rotated_mesh, success=False)
        assert 'geometry differs' in rejection.stderr, rejection.stderr

        whole, first, resumed = (root/k for k in ('whole','first','resumed'))
        solve(whole, ['--time-step',.02,'--steps',2], bc=shifted)
        solve(first, ['--time-step',.02,'--steps',1], bc=shifted)
        solve(resumed, ['--time-step',.02,'--steps',1,'--restart',str(first)+'.checkpoint'], bc=shifted)
        for output in (whole,resumed):
            assert transient.verify(mesh_path, output, root/(output.name+'-audit.json'))['valid']
        assert Path(str(whole)+'.checkpoint').read_bytes() == Path(str(resumed)+'.checkpoint').read_bytes()
        rejection = solve(root/'wrong-restart', ['--time-step',.02,'--steps',1,
                          '--restart',str(first)+'.checkpoint'], success=False)
        assert 'boundary' in rejection.stderr.lower(), rejection.stderr

        # Both verifier modes reject inconsistent boundary metadata rather
        # than trusting solver-reported success or pressure reference labels.
        for output, is_transient in ((base,False),(whole,True)):
            path = Path(str(output)+'.json')
            original_json = path.read_text()
            for mutate in (
                lambda p: p['boundaryConditions'][0].update(name='wrong patch'),
                lambda p: p['boundaryConditions'].pop(),
                lambda p: p.update(referenceSpeedRole='inlet-speed'),
            ):
                payload = json.loads(original_json)
                mutate(payload)
                path.write_text(json.dumps(payload))
                try:
                    if is_transient:
                        transient.verify(mesh_path,output,root/'bad-audit.json')
                    else:
                        verify(output)
                except audit.VerificationError:
                    pass
                else:
                    raise AssertionError('tampered explicit boundary metadata accepted')
            path.write_text(original_json)
        path = Path(str(base)+'.boundaries')
        path.write_text(source+'TRAILING\n')
        try: verify(base)
        except audit.VerificationError: pass
        else: raise AssertionError('trailing boundary data accepted')

        # A closed moving-wall domain uses a cell pressure gauge and zero
        # throughput, independently audited in both temporal modes.
        square=root/'square.xy'
        square.write_text('0 0\n1 0\n1 1\n0 1\n')
        closed_mesh=root/'closed.solver.cm2d'
        run([mesh_cli,square,root/'closed',4,.1,.1,'interior',root/'closed-foam',4,0])
        closed_bc=root/'closed.boundaries'
        run([cli,'--mesh',closed_mesh,'--case','cavity','--export-boundaries',closed_bc])
        closed=root/'closed-flow'
        solve(closed,mesh=closed_mesh,bc=closed_bc)
        verify(closed,closed_mesh)
        solve(closed,['--time-step',.02,'--steps',2],mesh=closed_mesh,bc=closed_bc)
        assert transient.verify(closed_mesh,closed,root/'closed-audit.json')['valid']
        print(json.dumps({'valid':True,'cells':len(original),'rotation':angle,
                          'pressureOffset':7.25,'restartIdentical':True,'inputBinding':True}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cli', type=Path, required=True)
    parser.add_argument('--mesh-cli', type=Path, required=True)
    parser.add_argument('--convection', choices=['upwind','face-limited-linear'], default='upwind')
    args=parser.parse_args()
    main(args.cli.resolve(),args.mesh_cli.resolve(),args.convection)
