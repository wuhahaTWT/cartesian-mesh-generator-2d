#!/usr/bin/env python3
"""Real channel SST-RANS verifier and fail-closed artifact tamper checks."""
import argparse
import csv
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools" / "verification"))
import verify_native_flow as native  # noqa: E402
import verify_sst_rans as verifier  # noqa: E402


def run(command, timeout=90):
    result = subprocess.run([str(x) for x in command], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=timeout)
    if result.returncode:
        raise RuntimeError(result.stdout)


def make_mesh(mesh_cli, root, output, case="channel"):
    request = native.Request("sst-rans-test", case, 3, 1 / 6, .001, 1.)
    command = native.mesh_command(mesh_cli, output, request, root)
    run(command)
    mesh = output.with_suffix(".solver.cm2d")
    if not mesh.exists():
        raise AssertionError(f"missing generated mesh: {mesh}")
    return mesh


def reject_csv(mesh, prefix, suffix, column, mutate, row_index=0):
    path = Path(str(prefix) + suffix)
    original = path.read_bytes()
    try:
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            records = list(reader)
            fields = reader.fieldnames
        if not fields or column not in fields:
            raise AssertionError(f"missing {column} in {path}")
        mutate(records, row_index, column)
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader(); writer.writerows(records)
        try:
            verifier.audit(mesh, prefix)
        except (ValueError, KeyError, AssertionError):
            return
        raise AssertionError(f"audit accepted tampered {suffix}:{column}")
    finally:
        path.write_bytes(original)


def reject_json(mesh, prefix, mutate):
    path = Path(str(prefix) + ".json")
    original = path.read_bytes()
    try:
        data = json.loads(original)
        mutate(data)
        path.write_text(json.dumps(data) + "\n")
        try:
            verifier.audit(mesh, prefix)
        except (ValueError, KeyError, AssertionError):
            return
        raise AssertionError("audit accepted tampered JSON")
    finally:
        path.write_bytes(original)


def main(args):
    with tempfile.TemporaryDirectory(prefix="cartmesh-sst-rans-") as name:
        root = Path(name)
        channel_mesh = make_mesh(args.mesh_cli, root, root / "channel")
        plate_mesh = make_mesh(args.mesh_cli, root, root / "flatplate", "manufactured")
        for mode in ('channel', 'flatplate', 'flatplate-symmetry', 'channel-sweep', 'flatplate-sweep'):
            is_channel=mode.startswith('channel')
            mesh = channel_mesh if is_channel else plate_mesh
            prefix = root / mode
            run([args.probe, mesh, prefix] + ([] if mode == 'channel' else [mode]))
            accepted = verifier.audit(mesh, prefix)
            assert accepted["valid"] is True and accepted["cells"] > 1

            def add(records, index, key):
                records[index][key] = str(float(records[index][key]) + .1)

            for column in ("u", "p", "speed", "k", "omega", "nuT", "strain", "gradWx"):
                reject_csv(mesh, prefix, ".cells.csv", column, add)

            with (Path(str(prefix) + ".faces.csv")).open(newline="") as stream:
                face_rows = list(csv.DictReader(stream))
            wall_row = next(i for i, row in enumerate(face_rows) if row["wall"] == "1")
            reject_csv(mesh, prefix, ".faces.csv", "omegaBoundary", add, wall_row)
            reject_csv(mesh, prefix, ".faces.csv", "viscosity", add)
            reject_csv(mesh, prefix, ".faces.csv", "kDiffusion", add)
            reject_csv(mesh, prefix, ".faces.csv", "advectionX", add)

            reject_csv(mesh, prefix, ".history.csv", "momentumResidual",
                       lambda records, index, key: records[-1].__setitem__(key, str(float(records[-1][key]) + .1)))
            reject_csv(mesh, prefix, ".history.csv", "kNorm",
                       lambda records, index, key: records[-1].__setitem__(key, '0'))
            for column in ('kCellResidual','omegaCellResidual'):
                reject_csv(mesh,prefix,'.history.csv',column,
                           lambda records,index,key:records[-1].__setitem__(key,'1e-4'))
            for column in ('lossW','sourceK'):
                reject_csv(mesh,prefix,'.cells.csv',column,add)
            reject_csv(mesh, prefix, ".history.csv", "turbulenceIterations",
                       lambda records, index, key: records[index].__setitem__(key, '2'))

            def stale_closure(records, index, key):
                source = max(range(len(records)),
                             key=lambda j: abs(float(records[j][key]) - float(records[index][key])))
                if source == index:
                    raise AssertionError(f"no distinct value for {key}")
                records[index][key] = records[source][key]
            for column in ("sourceW", "nuT", "Dk"):
                reject_csv(mesh, prefix, ".cells.csv", column, stale_closure)

            reject_json(mesh, prefix, lambda data: data.__setitem__("model", "wrong-model"))
            reject_json(mesh, prefix, lambda data: data.__setitem__("pressureConvention", "p/rho+2k/3"))
            reject_json(mesh, prefix, lambda data: data.__setitem__("converged", False))
            reject_json(mesh, prefix, lambda data: data.__setitem__("momentumResidual", 0.0))
            reject_json(mesh, prefix, lambda data: data.__setitem__("scalarRelativeTolerance", 1.0))
            reject_json(mesh, prefix, lambda data: data.__setitem__("globalRelativeImbalance", 1.0))
            for invalid in (0, -1, 501, 1.5, True):
                reject_json(mesh, prefix, lambda data, value=invalid:
                            data.__setitem__("turbulenceUpdatesPerIteration", value))
            for invalid in (-1, 2, 1.5, True):
                reject_json(mesh, prefix, lambda data, value=invalid:
                            data.__setitem__('scalarCorrectionsPerUpdate', value))
            if not is_channel:
                assert abs(accepted['boundarySummary']['wall']['length']-.5) <= 8*math.ulp(.5)
                assert accepted['plateWallSamples']
                if mode in ('flatplate','flatplate-sweep'):
                    assert accepted['boundarySummary']['farfield']['outwardFlux'] > 0
                reject_json(mesh, prefix, lambda data: data.__setitem__('flatPlateLeadingEdge', .51))
                reject_json(mesh, prefix, lambda data: data.__setitem__('flatPlateTop', 'invalid'))
                reject_json(mesh, prefix, lambda data: data.__setitem__('case', 'channel'))
                if mode in ('flatplate','flatplate-sweep'):
                    reject_json(mesh, prefix, lambda data: data.__setitem__('flatPlateTop', 'symmetry'))
        # Similarity transformation: x,y -> 2(x,y)+(-1,2), U -> 2U,
        # nu -> 4nu, k -> 4k, omega unchanged; same Re and dimensionless model.
        similarity=[]
        for label,bounds,options in (
            ('unit',[],[]),
            ('scaled',['-1','2','1','4'],['--speed','2','--nu','.004','--inlet-k','.004','--inlet-omega','2','--leading-edge','0'])):
            mesh_prefix=root/('similarity-'+label);mesh=Path(str(mesh_prefix)+'.cm2d')
            run([args.rect_probe,'8','8','1',mesh_prefix,*bounds])
            prefix=root/('flow-'+label)
            run([args.probe,mesh,prefix,'flatplate-sweep',*options])
            result=verifier.audit(mesh,prefix)
            assert result['valid'] and result['plateReynolds']==500
            with Path(str(prefix)+'.cells.csv').open() as stream:fields=list(csv.DictReader(stream))
            similarity.append((fields,result))
            for key in ('nu','speed','inletK','inletOmega'):
                for invalid in (-1,True,'1',float('nan')):
                    reject_json(mesh,prefix,lambda data,key=key,value=invalid:data.__setitem__(key,value))
                reject_json(mesh,prefix,lambda data,key=key:data.__setitem__(key,data[key]*1.1))
            for options_bad in (['--nu','0'],['--speed','-1'],['--inlet-k','-.1'],
                                ['--inlet-omega','nan'],['--nu','1junk'],['--nu'],
                                ['--nu','.001','--nu','.002'],['--unknown','1'],
                                ['--leading-edge','.12345']):
                result_bad=subprocess.run([str(args.probe),str(mesh),str(root/'invalid'),'flatplate-sweep',*options_bad],
                                          stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
                assert result_bad.returncode!=0, options_bad
        for a,b in zip(similarity[0][0],similarity[1][0]):
            for key,factor in (('u',2),('v',2),('p',4),('k',4),('omega',1)):
                assert abs(float(a[key])-float(b[key])/factor)<1e-5*max(1,abs(float(a[key]))),(key,a[key],b[key])
        for a,b in zip(similarity[0][1]['plateWallSamples'],similarity[1][1]['plateWallSamples']):
            assert abs(a['Cf']-b['Cf'])<1e-5 and abs(a['yPlus']-b['yPlus'])<1e-5
        for bad_bounds in (['0','0','0','1'],['0','1','1','0'],['nan','0','1','1']):
            rejected=subprocess.run([str(args.rect_probe),'8','8','1',str(root/'bad-mesh'),*bad_bounds],
                                    stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
            assert rejected.returncode!=0
        # Anisotropic wall-resolved probes may exceed 256 cells on one axis,
        # but the existing total-cell resource bound must remain enforced.
        run([args.rect_probe,'400','2','0',root/'wide-mesh'])
        for dims in (('400','400'),('4097','2'),('2','4097'),('2.5','2')):
            rejected=subprocess.run([str(args.rect_probe),*dims,'0',str(root/'bad-size')],
                                    stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
            assert rejected.returncode!=0, dims
        ilu=root/'ilu'
        run([args.probe,plate_mesh,ilu,'flatplate-sweep','--scalar-preconditioner','ilu0'])
        assert verifier.audit(plate_mesh,ilu)['scalarPreconditioner']=='ilu0'
        reject_json(plate_mesh,ilu,lambda data:data.__setitem__('scalarPreconditioner','invalid'))
        # Iteration-bounded diagnostics must never masquerade as accepted fields.
        limited=root/'limited'
        stopped=subprocess.run([str(args.probe),str(plate_mesh),str(limited),'flatplate-sweep',
                                '--max-iterations','1'],capture_output=True,timeout=90)
        assert stopped.returncode!=0
        diagnostic=json.loads(Path(str(limited)+'.diagnostics.json').read_text())
        assert diagnostic['converged'] is False and diagnostic['stopReason']=='iteration-limit'
        assert diagnostic['iterations']==diagnostic['maxIterations']==1
        assert diagnostic['performance']['sstUpdates']==1
        assert Path(str(limited)+'.history.csv').exists()
        assert not any(Path(str(limited)+s).exists() for s in ('.json','.cells.csv','.faces.csv'))
        previous=Path(str(ilu)+'.json').read_bytes()
        rerun=subprocess.run([str(args.probe),str(plate_mesh),str(ilu),'flatplate-sweep'],capture_output=True,timeout=90)
        assert rerun.returncode!=0 and b'fresh output prefix' in rerun.stderr
        assert Path(str(ilu)+'.json').read_bytes()==previous
        for options in (['--max-iterations','0'],['--max-iterations','1.5'],['--max-iterations','2001'],
                        ['--scalar-preconditioner','ic0']):
            rejected=subprocess.run([str(args.probe),str(plate_mesh),str(root/'bad-options'),'flatplate-sweep',*options],
                                    capture_output=True,timeout=90)
            assert rejected.returncode!=0
    print("SST-RANS verifier: channel, both flat plate boundaries and bounded corrections audited; all tamper cases rejected.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-cli", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--rect-probe", type=Path, required=True)
    main(parser.parse_args())
