#!/usr/bin/env python3
"""Real-mesh nonlinear-SST audit and fail-closed closure tamper checks."""
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
import verify_sst_nonlinear as verifier  # noqa: E402


def run(command, timeout=90):
    result = subprocess.run([str(x) for x in command], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=timeout)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {result.stdout}")
    return result


def make_mesh(mesh_cli, root, output, shear=False):
    request = native.Request("sst-test", "manufactured", 3, 1 / 6, .1, 1.)
    command = native.mesh_command(mesh_cli, output, request, root)
    if shear:
        geometry = Path(command[1])
        points = []
        for line in geometry.read_text().splitlines():
            x, y = (float(value) for value in line.split())
            points.append(f"{x + .2 * y:.17g} {y:.17g}")
        geometry.write_text("\n".join(points) + "\n")
        command[2] = str(output)
        command[7] = str(output.parent / (output.name + "-openfoam"))
    run(command)
    mesh = output.with_suffix(".solver.cm2d")
    if not mesh.exists():
        raise AssertionError(f"missing generated solver mesh: {mesh}")
    return mesh


def reject_csv(mesh, prefix, suffix, column, mutate, row_index=0):
    path = Path(str(prefix) + suffix)
    original = path.read_bytes()
    try:
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            rows = list(reader)
            fields = reader.fieldnames
        if not fields or column not in fields or not rows:
            raise AssertionError(f"missing {column} in {path}")
        mutate(rows, row_index, column)
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
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
        raise AssertionError("audit accepted tampered metadata")
    finally:
        path.write_bytes(original)


def main(args):
    with tempfile.TemporaryDirectory(prefix="cartmesh-sst-nonlinear-") as name:
        root = Path(name)
        mesh = make_mesh(args.mesh_cli, root, root / "mesh")
        prefix = root / "probe"
        run([args.probe, mesh, prefix])
        accepted = verifier.audit(mesh, prefix)
        assert accepted["valid"] is True and accepted["cells"] > 1

        # Exercise the independent auditor on a separately generated sheared mesh.
        shear_mesh = make_mesh(args.mesh_cli, root, root / "shear", shear=True)
        shear_prefix = root / "shear-probe"
        run([args.probe, shear_mesh, shear_prefix])
        shear_accepted = verifier.audit(shear_mesh, shear_prefix)
        assert shear_accepted["valid"] is True and shear_accepted["cells"] > 1

        def add(rows, index, key, amount=.1):
            rows[index][key] = str(float(rows[index][key]) + amount)

        # Returned state and independently reconstructed closure quantities.
        for column in ("k", "omega", "gradKx", "gradWx", "nuT", "sourceW", "lossW"):
            reject_csv(mesh, prefix, ".cells.csv", column, add)
        for column in ("kDiffusion", "omegaDiffusion"):
            reject_csv(mesh, prefix, ".faces.csv", column, add)

        with (Path(str(prefix) + ".faces.csv")).open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        wall_row = next(i for i, row in enumerate(rows) if row["wall"] == "1")
        # The face schema names these values kBoundary/omegaBoundary; selecting
        # a wall row makes these the resolved wall boundary checks requested here.
        reject_csv(mesh, prefix, ".faces.csv", "omegaBoundary", add, wall_row)
        reject_csv(mesh, prefix, ".faces.csv", "kBoundary", add, wall_row)

        # A face identity mutation is also rejected before values can be trusted.
        reject_csv(mesh, prefix, ".faces.csv", "face",
                   lambda rows, index, key: rows[index].__setitem__(key, str(index + 1)))

        # Replace final fields with a different cell's closure: this catches stale
        # coefficient columns even if the replacement remains finite and positive.
        def replace_from_other(rows, index, key):
            source = max(range(len(rows)),
                         key=lambda candidate: abs(float(rows[candidate][key]) -
                                                   float(rows[index][key])))
            if source == index:
                raise AssertionError(f"no distinct {key} value for stale-closure test")
            rows[index][key] = rows[source][key]
        for column in ("sourceW", "nuT", "lossW"):
            reject_csv(mesh, prefix, ".cells.csv", column, replace_from_other)

        reject_csv(mesh, prefix, ".history.csv", "kNorm",
                   lambda rows, index, key: rows[-1].__setitem__(key, str(float(rows[-1][key]) + .1)),
                   0)
        reject_csv(mesh, prefix, ".history.csv", "omegaCellResidual",
                   lambda rows, index, key: rows[-1].__setitem__(key, "0.1"))
        reject_json(mesh, prefix,
                    lambda data: data.__setitem__("iterations", int(data["iterations"]) + 1))
        print("SST nonlinear verifier: square/shear audits passed and all tamper cases rejected.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-cli", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    main(parser.parse_args())
