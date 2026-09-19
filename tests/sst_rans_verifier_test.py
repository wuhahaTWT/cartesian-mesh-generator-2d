#!/usr/bin/env python3
"""Real channel SST-RANS verifier and fail-closed artifact tamper checks."""
import argparse
import csv
import json
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


def make_mesh(mesh_cli, root, output):
    request = native.Request("sst-rans-test", "channel", 3, 1 / 6, .001, 1.)
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
        mesh = make_mesh(args.mesh_cli, root, root / "channel")
        prefix = root / "probe"
        run([args.probe, mesh, prefix])
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
    print("SST-RANS verifier: real channel audit passed and all tamper cases rejected.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-cli", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    main(parser.parse_args())
