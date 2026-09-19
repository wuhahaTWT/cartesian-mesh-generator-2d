#!/usr/bin/env python3
"""Small end-to-end contracts for variable-viscosity flow checkpoints."""
import argparse
import csv
import json
import re
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
sys.path.insert(0, str(ROOT / "tools" / "verification"))
from scalar_transport_cli_test import boundary_csv, native_mesh, run  # noqa: E402
from verify_native_flow import argument_parser, read_cm2d, verify_case  # noqa: E402
from verify_scalar_transport import verify as verify_scalar  # noqa: E402
from verify_transient_flow import verify as verify_transient  # noqa: E402


def expect_failure(command, pattern):
    result = run(command)
    assert result.returncode == 1, (result.returncode, result.stdout, result.stderr)
    assert re.search(pattern, result.stderr, re.IGNORECASE), result.stderr


def replace_first_value(contents, value):
    lines = contents.splitlines()
    face = lines[1].split(",", 1)[0]
    lines[1] = f"{face},{value}"
    return "\n".join(lines) + "\n"


def viscosity_csv(mesh_path, flow_prefix, path):
    with flow_prefix.with_suffix(".faces.csv").open(newline="") as stream:
        faces = list(csv.DictReader(stream))
    mesh = read_cm2d(mesh_path)
    with path.open("w", newline="") as stream:
        stream.write("face,viscosity\n")
        for row in faces:
            edge = mesh.edges[int(row["face"])]
            x = (mesh.vertices[edge.v0][0] + mesh.vertices[edge.v1][0]) / 2.
            coefficient = .1 * (1. + x)
            stream.write(f"{row['face']},{coefficient:.17g}\n")
    return faces


def main(args):
    with tempfile.TemporaryDirectory(prefix="cartmesh-variable-viscosity-") as name:
        root = Path(name)
        mesh = native_mesh(args.mesh_cli, root, "unit")
        seed = root / "seed"
        seed_result = run([args.flow_cli, "--mesh", mesh, "--output", seed,
                           "--case", "cavity", "--nu", ".1", "--speed", "1",
                           "--time-step", ".01", "--steps", "1", "--tolerance", "1e-9"])
        assert seed_result.returncode == 0, (seed_result.stdout, seed_result.stderr)
        csv_path = root / "viscosity.csv"
        # The unit mesh has face centres available from the CM2D only; for this
        # contract a positive per-face field is sufficient, while the reader
        # independently binds every ID and rejects missing/duplicate rows.
        faces = viscosity_csv(mesh, seed, csv_path)
        assert len(faces) > 0

        one = root / "one"
        result = run([args.flow_cli, "--mesh", mesh, "--output", one, "--case", "cavity",
                      "--nu", ".1", "--speed", "1", "--face-viscosity", csv_path,
                      "--time-step", ".01", "--steps", "1", "--tolerance", "1e-9"])
        assert result.returncode == 0, (result.stdout, result.stderr)
        two = root / "two"
        result = run([args.flow_cli, "--mesh", mesh, "--output", two, "--case", "cavity",
                      "--nu", ".1", "--speed", "1", "--face-viscosity", csv_path,
                      "--time-step", ".01", "--steps", "2", "--tolerance", "1e-9"])
        assert result.returncode == 0, (result.stdout, result.stderr)
        restart = root / "restart"
        result = run([args.flow_cli, "--mesh", mesh, "--output", restart, "--case", "cavity",
                      "--nu", ".1", "--speed", "1", "--face-viscosity", csv_path,
                      "--restart", one.with_suffix(".checkpoint"), "--time-step", ".01",
                      "--steps", "1", "--tolerance", "1e-9"])
        assert result.returncode == 0, (result.stdout, result.stderr)
        assert two.with_suffix(".checkpoint").read_bytes() == restart.with_suffix(".checkpoint").read_bytes()

        for prefix in (one, two, restart):
            audit = verify_transient(mesh, prefix, root / (prefix.name+"-audit.json"))
            assert audit['valid'], audit

        uniform_file = root / "uniform.csv"
        uniform_file.write_text("face,viscosity\n"+"".join(f"{r['face']},.1\n" for r in faces))
        uniform = root / "uniform"
        result = run([args.flow_cli, "--mesh", mesh, "--output", uniform, "--case", "cavity",
                      "--nu", ".1", "--speed", "1", "--face-viscosity", uniform_file,
                      "--time-step", ".01", "--steps", "1", "--tolerance", "1e-9"])
        assert result.returncode == 0, result.stderr
        for suffix in ('.cells.csv','.fields.json','.residuals.csv','.vtk'):
            assert Path(str(uniform)+suffix).read_bytes() == Path(str(seed)+suffix).read_bytes(), suffix

        changed = root / "changed.csv"
        changed.write_text(replace_first_value(csv_path.read_text(), ".11"))
        expect_failure([args.flow_cli, "--mesh", mesh, "--output", root / "bad-restart",
                        "--case", "cavity", "--nu", ".1", "--speed", "1",
                        "--face-viscosity", changed, "--restart", one.with_suffix(".checkpoint"),
                        "--time-step", ".01", "--steps", "1"], r"face viscosity")

        old = root / "old"
        result = run([args.flow_cli, "--mesh", mesh, "--output", old, "--case", "cavity",
                      "--nu", ".1", "--speed", "1", "--time-step", ".01", "--steps", "1"])
        assert result.returncode == 0, (result.stdout, result.stderr)
        expect_failure([args.flow_cli, "--mesh", mesh, "--output", root / "bad-legacy",
                        "--case", "cavity", "--nu", ".1", "--speed", "1",
                        "--face-viscosity", csv_path, "--restart", old.with_suffix(".checkpoint"),
                        "--time-step", ".01", "--steps", "1"], r"legacy checkpoint")

        boundary = boundary_csv(one, root / "boundary.csv")
        scalar = root / "scalar"
        result = run([args.transport_cli, "--mesh", mesh, "--flow-checkpoint", one.with_suffix(".checkpoint"),
                      "--boundary", root / "boundary.csv", "--output", scalar,
                      "--diffusivity", ".01", "--initial", "0", "--dt", ".01", "--steps", "1"])
        assert result.returncode == 0, (result.stdout, result.stderr)
        verify_scalar(scalar)

        bad_cases = {
            "duplicate": csv_path.read_text() + "0,.1\n",
            "missing": "\n".join(csv_path.read_text().splitlines()[:-1]) + "\n",
            "zero": replace_first_value(csv_path.read_text(), "0"),
            "nan": replace_first_value(csv_path.read_text(), "nan"),
            "extra": replace_first_value(csv_path.read_text(), ".1,extra"),
        }
        for label, contents in bad_cases.items():
            bad = root / (label + ".csv")
            bad.write_text(contents)
            expect_failure([args.flow_cli, "--mesh", mesh, "--output", root / ("bad-" + label),
                            "--case", "cavity", "--nu", ".1", "--speed", "1",
                            "--face-viscosity", bad, "--time-step", ".01", "--steps", "1"],
                           r"viscosity|finite|positive|row|missing|duplicate")

        expect_failure([args.flow_cli, "--mesh", mesh, "--output", root / "bad-file-mms",
                        "--case", "manufactured", "--nu", ".1", "--face-viscosity", csv_path],
                       r"face-viscosity file requires")
        expect_failure([args.flow_cli, "--mesh", mesh, "--output", root / "bad-slope",
                        "--case", "manufactured", "--nu", ".1",
                        "--manufactured-viscosity-slope", "-1"], r"slope > -1")

        mms = root / "mms"
        result = run([args.flow_cli, "--mesh", mesh, "--output", mms, "--case", "manufactured",
                      "--nu", ".1", "--speed", "1", "--tolerance", "1e-9",
                      "--manufactured-viscosity-slope", "1"])
        assert result.returncode == 0, (result.stdout, result.stderr)
        verifier_args = argument_parser().parse_args([])
        verifier_args.manufactured_pressure_slope = 0.
        verified = verify_case(mesh, mms, "manufactured", .1, 1., verifier_args)
        assert verified["valid"] is True, verified
        tampered = root / "tampered"
        for suffix in (".json", ".cells.csv", ".faces.csv", ".residuals.csv", ".fields.json", ".vtk"):
            shutil.copyfile(str(mms) + suffix, str(tampered) + suffix)
        cells = tampered.with_suffix(".cells.csv")
        with cells.open(newline="") as stream:
            reader = csv.DictReader(stream)
            records = list(reader)
            names = reader.fieldnames
        assert names and "sourceX" in names
        records[0]["sourceX"] = str(float(records[0]["sourceX"]) + .123)
        with cells.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=names)
            writer.writeheader(); writer.writerows(records)
        bad_verified = verify_case(mesh, tampered, "manufactured", .1, 1., verifier_args)
        assert bad_verified["valid"] is False
        tampered_face = root / "tampered-face-viscosity"
        for suffix in (".json", ".cells.csv", ".faces.csv", ".residuals.csv", ".fields.json", ".vtk"):
            shutil.copyfile(str(mms) + suffix, str(tampered_face) + suffix)
        faces_path = tampered_face.with_suffix(".faces.csv")
        with faces_path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            face_records = list(reader)
            face_names = reader.fieldnames
        assert face_names and "viscosity" in face_names
        face_records[0]["viscosity"] = str(float(face_records[0]["viscosity"]) * 1.25)
        with faces_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=face_names)
            writer.writeheader(); writer.writerows(face_records)
        bad_face_verified = verify_case(mesh, tampered_face, "manufactured", .1, 1., verifier_args)
        assert bad_face_verified["valid"] is False
    print("Variable-viscosity CLI: v3 restart, frozen transport, CSV rejection, and manufactured verification passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-cli", required=True)
    parser.add_argument("--flow-cli", required=True)
    parser.add_argument("--transport-cli", required=True)
    main(parser.parse_args())
