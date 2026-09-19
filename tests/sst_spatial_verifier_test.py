#!/usr/bin/env python3
"""Real-mesh SST spatial audit and fail-closed artifact tamper checks."""
import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools" / "verification"))
import verify_native_flow as native  # noqa: E402
import verify_sst_spatial as verifier  # noqa: E402
import verify_wall_distance as wall_verifier  # noqa: E402


def run(command, timeout=90):
    result = subprocess.run([str(value) for value in command], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=timeout)
    if result.returncode:
        raise RuntimeError(result.stdout)
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


def copy_probe(probe, mesh, prefix):
    run([probe, mesh, prefix], timeout=90)


def tamper_and_reject(mesh, prefix, suffix, column, mutate, audit_fn=verifier.audit):
    path = Path(str(prefix) + suffix)
    original = path.read_bytes()
    try:
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            records = list(reader)
            fields = reader.fieldnames
        assert fields and column in fields, (path, column, fields)
        mutate(records, column)
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(records)
        try:
            audit_fn(mesh, prefix)
        except (ValueError, KeyError, AssertionError):
            return
        raise AssertionError(f"audit accepted tampered {suffix}:{column}")
    finally:
        path.write_bytes(original)


def main(args):
    with tempfile.TemporaryDirectory(prefix="cartmesh-sst-spatial-") as name:
        root = Path(name)
        mesh = make_mesh(args.mesh_cli, root, root / "mesh")
        prefix = root / "probe"
        copy_probe(args.probe, mesh, prefix)
        accepted = verifier.audit(mesh, prefix)
        assert accepted["valid"] is True and accepted["cells"] > 0

        wall_prefix = root / "wall"
        run([args.wall_probe, mesh, wall_prefix], timeout=90)
        wall_accepted = wall_verifier.audit(mesh, wall_prefix)
        assert wall_accepted["valid"] is True
        wall_mesh = native.read_cm2d(mesh)
        internal = next(edge.id for edge in wall_mesh.edges if edge.neighbour >= 0)
        tamper_and_reject(mesh, wall_prefix, ".csv", "distance",
                          lambda rows, key: rows[0].__setitem__(key, str(float(rows[0][key]) + .1)), wall_verifier.audit)
        tamper_and_reject(mesh, wall_prefix, ".csv", "nearestFace",
                          lambda rows, key: rows[0].__setitem__(key, str(internal)), wall_verifier.audit)
        tamper_and_reject(mesh, wall_prefix, ".csv", "cell",
                          lambda rows, key: rows[0].__setitem__(key, "1"), wall_verifier.audit)

        shear_mesh = make_mesh(args.mesh_cli, root, root / "shear", shear=True)
        shear_wall_prefix = root / "shear-wall"
        run([args.wall_probe, shear_mesh, shear_wall_prefix], timeout=90)
        try:
            wall_verifier.audit(shear_mesh, shear_wall_prefix)
        except ValueError:
            pass
        else:
            raise AssertionError("rectangle-only wall verifier accepted shear domain")
        shear_prefix = root / "shear-probe"
        copy_probe(args.probe, shear_mesh, shear_prefix)
        shear_accepted = verifier.audit(shear_mesh, shear_prefix)
        assert shear_accepted["valid"] is True and shear_accepted["cells"] > 0

        # nearestFace must become an internal edge ID, which is known from the
        # independently parsed CM2D rather than from the probe's CSV.
        parsed = native.read_cm2d(mesh)
        internal = next(edge.id for edge in parsed.edges if edge.neighbour >= 0)
        tamper_and_reject(mesh, prefix, ".cells.csv", "distance",
                          lambda rows, key: rows[0].__setitem__(key, str(float(rows[0][key]) + .1)))
        tamper_and_reject(mesh, prefix, ".cells.csv", "nearestFace",
                          lambda rows, key: rows[0].__setitem__(key, str(internal)))
        for column in ("gradKx", "F1", "Dk", "sourceW", "lossK"):
            tamper_and_reject(mesh, prefix, ".cells.csv", column,
                              lambda rows, key: rows[0].__setitem__(key, str(float(rows[0][key]) + .1)))
        tamper_and_reject(mesh, prefix, ".cells.csv", "k",
                          lambda rows, key: rows[0].__setitem__(key, "-1"))
        tamper_and_reject(mesh, prefix, ".faces.csv", "kDiffusion",
                          lambda rows, key: rows[0].__setitem__(key, str(float(rows[0][key]) + .1)))
        tamper_and_reject(mesh, prefix, ".faces.csv", "wall",
                          lambda rows, key: rows[0].__setitem__(key, "0" if rows[0][key] == "1" else "1"))
        tamper_and_reject(mesh, prefix, ".faces.csv", "volumeFlux",
                          lambda rows, key: rows[0].__setitem__(key, str(float(rows[0][key]) + .1)))
    print("SST spatial verifier: real unit/shear meshes and all tamper cases rejected.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-cli", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--wall-probe", type=Path, required=True)
    main(parser.parse_args())
