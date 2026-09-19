#!/usr/bin/env python3
"""Black-box scalar transport CLI checks and input rejection regressions."""
import argparse
import csv
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "verification"))
from verify_scalar_transport import verify as verify_scalar


def run(command, timeout=30):
    return subprocess.run([str(x) for x in command], text=True,
                          capture_output=True, timeout=timeout)


def square(path, width=1.0, height=1.0):
    path.write_text(f"0 0\n{width:.17g} 0\n{width:.17g} {height:.17g}\n0 {height:.17g}\n")


def native_mesh(mesh_cli, root, name, width=1.0, height=1.0):
    boundary = root / (name + ".xy")
    prefix = root / name
    square(boundary, width, height)
    # Interior fluid gives the transport verification an exact rectangular
    # envelope.  The final solver mesh is required; intermediate CM2D is not
    # accepted by the transport CLI.
    solver_case = root / (name + "-openfoam")
    # A case directory requests the mesher's final solver partition and
    # *.solver.cm2d export.  This is the same documented workflow used by the
    # independent native-flow verifier.
    result = run([mesh_cli, boundary, prefix, "4", f"{1.0 / 14.0:.17g}", ".1",
                  "interior", solver_case, "4", "0"])
    assert result.returncode == 0, (result.stdout, result.stderr)
    mesh = prefix.with_suffix(".solver.cm2d")
    assert mesh.exists(), ("native mesher did not produce final solver mesh",
                           result.stdout, result.stderr)
    return mesh


def transport(transport_cli, mesh, output, *args):
    return run([transport_cli, "--mesh", mesh, "--output", output, *args])


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def finite_rows(path):
    result = rows(path)
    assert result
    for row in result:
        for key, value in row.items():
            if not (key in {"exact", "previous"} and value.lower() == "nan"):
                assert math.isfinite(float(value)), (path, key, value)
    return result


def independent_balance(prefix):
    cells = finite_rows(prefix.with_suffix(".cells.csv"))
    faces = finite_rows(prefix.with_suffix(".faces.csv"))
    balance = [float(row["temporalIntegral"]) - float(row["sourceIntegral"])
               for row in cells]
    for row in faces:
        flux = float(row["advectiveFlux"]) + float(row["diffusiveFlux"])
        owner = int(row["owner"])
        balance[owner] += flux
        neighbour = int(row["neighbour"])
        if neighbour >= 0:
            balance[neighbour] -= flux
    return cells, faces, balance


def boundary_csv(face_prefix, path, kind="value", inflow=""):
    face_rows = rows(face_prefix.with_suffix(".faces.csv"))
    boundary = [row for row in face_rows if int(row["neighbour"]) < 0]
    assert boundary
    with path.open("w", newline="") as stream:
        stream.write("face,type,value,inflowValue\n")
        for row in boundary:
            stream.write(f"{row['face']},{kind},1,{inflow}\n")
    return boundary


def diffusivity_csv(face_prefix, path, value=".01"):
    face_rows = rows(face_prefix.with_suffix(".faces.csv"))
    assert face_rows
    with path.open("w", newline="") as stream:
        stream.write("face,diffusivity\n")
        for row in face_rows:
            coefficient = ".02" if int(row["face"]) == int(face_rows[0]["face"]) else value
            stream.write(f"{row['face']},{coefficient}\n")
    return face_rows


def expect_failure(command, pattern, timeout=30):
    result = run(command, timeout)
    assert result.returncode == 1, (result.returncode, result.stdout, result.stderr)
    assert re.search(pattern, result.stderr, re.IGNORECASE), result.stderr
    return result


def copy_prefix(source, destination):
    for suffix in (".json", ".cells.csv", ".faces.csv", ".history.csv", ".vtk"):
        shutil.copyfile(str(source) + suffix, str(destination) + suffix)


def expect_verifier_failure(prefix, pattern):
    try:
        verify_scalar(prefix)
    except (ValueError, KeyError, OSError, ZeroDivisionError) as error:
        assert re.search(pattern, str(error), re.IGNORECASE), str(error)
        return
    raise AssertionError("independent scalar verifier accepted tampered output")


def main(args):
    with tempfile.TemporaryDirectory(prefix="cartmesh-scalar-") as name:
        root = Path(name)
        unit_mesh = native_mesh(args.mesh_cli, root, "unit")

        sine = root / "sine"
        result = transport(args.transport_cli, unit_mesh, sine,
                           "--verification", "sine", "--diffusivity", ".01",
                           "--speed", "1")
        assert result.returncode == 0, (result.stdout, result.stderr)
        sine_summary = json.loads(sine.with_suffix(".json").read_text())
        assert sine_summary["converged"] is True and sine_summary["verification"] == "sine"
        sine_cells, _, sine_balance = independent_balance(sine)
        assert max(abs(v) for v in sine_balance) < 2e-7
        assert len(sine_cells) >= 4
        assert all(math.isfinite(float(row["value"])) for row in sine_cells)
        verify_scalar(sine)

        variable = root / "variable-sine"
        variable_result = transport(args.transport_cli, unit_mesh, variable,
                                    "--verification", "variable-sine", "--diffusivity", ".01",
                                    "--speed", "1")
        assert variable_result.returncode == 0, (variable_result.stdout, variable_result.stderr)
        variable_summary = json.loads(variable.with_suffix(".json").read_text())
        assert variable_summary["verification"] == "variable-sine"
        assert variable_summary["diffusivityModel"] == "face-values"
        assert variable_summary["diffusivityFile"] == ""
        variable_faces = rows(variable.with_suffix(".faces.csv"))
        assert "diffusivity" in variable_faces[0]
        assert {int(row["face"]) for row in variable_faces} == {
            int(row["face"]) for row in rows(sine.with_suffix(".faces.csv"))
        }
        assert all(float(row["diffusivity"]) > 0 for row in variable_faces)
        verify_scalar(variable)

        decay_one = root / "decay_one"
        assert transport(args.transport_cli, unit_mesh, decay_one,
                         "--verification", "decay", "--diffusivity", ".01",
                         "--initial", "0", "--dt", ".01", "--steps", "1").returncode == 0
        decay_two = root / "decay_two"
        assert transport(args.transport_cli, unit_mesh, decay_two,
                         "--verification", "decay", "--diffusivity", ".01",
                         "--initial", "0", "--dt", ".01", "--steps", "2").returncode == 0
        one_cells, _, one_balance = independent_balance(decay_one)
        two_cells, _, two_balance = independent_balance(decay_two)
        assert max(abs(v) for v in one_balance + two_balance) < 2e-7
        assert len(one_cells) == len(two_cells)
        # The second run must export the accepted state from the preceding
        # physical step, rather than repeating the initial state.
        for old, current in zip(one_cells, two_cells):
            assert math.isclose(float(current["previous"]), float(old["value"]),
                                rel_tol=2e-12, abs_tol=2e-12)
            expected = float(current["area"]) * (float(current["value"]) -
                                                   float(current["previous"])) / .01
            assert math.isclose(float(current["temporalIntegral"]), expected,
                                rel_tol=2e-10, abs_tol=2e-10)
        verify_scalar(decay_one)
        verify_scalar(decay_two)

        physical_mesh = native_mesh(args.mesh_cli, root, "channel", 2.0, 1.0)
        flow = root / "flow"
        flow_result = run([args.flow_cli, "--mesh", physical_mesh, "--output", flow,
                           "--case", "channel", "--nu", ".01", "--speed", "1",
                           "--tolerance", "1e-6", "--max-iterations", "100",
                           "--time-step", ".02", "--steps", "1"])
        assert flow_result.returncode == 0, (flow_result.stdout, flow_result.stderr)
        checkpoint = flow.with_suffix(".checkpoint")
        assert checkpoint.exists()
        bc = root / "boundary.csv"
        boundary = boundary_csv(flow, bc)
        physical = root / "physical"
        physical_result = transport(args.transport_cli, physical_mesh, physical,
                                    "--flow-checkpoint", checkpoint,
                                    "--boundary", bc, "--diffusivity", ".01",
                                    "--initial", "0", "--dt", ".01", "--steps", "2")
        assert physical_result.returncode == 0, (physical_result.stdout, physical_result.stderr)
        physical_cells, physical_faces, physical_balance = independent_balance(physical)
        assert len(physical_cells) > 1 and len(physical_faces) >= len(boundary)
        temperatures = [float(row["value"]) for row in physical_cells]
        assert max(temperatures) > .1 and min(temperatures) >= -1e-10 and max(temperatures) <= 1+1e-10
        assert max(abs(v) for v in physical_balance) < 2e-6
        assert physical.with_suffix(".vtk").read_text().count("SCALARS scalar") == 1
        verify_scalar(physical)

        face_csv = root / "face-diffusivity.csv"
        face_rows = diffusivity_csv(physical, face_csv)
        face_values = root / "face-values"
        face_values_result = transport(args.transport_cli, physical_mesh, face_values,
                                       "--flow-checkpoint", checkpoint, "--boundary", bc,
                                       "--diffusivity", ".01", "--face-diffusivity", face_csv,
                                       "--initial", "0", "--dt", ".01", "--steps", "1")
        assert face_values_result.returncode == 0, (face_values_result.stdout, face_values_result.stderr)
        face_summary = json.loads(face_values.with_suffix(".json").read_text())
        assert face_summary["diffusivityModel"] == "face-values"
        assert Path(face_summary["diffusivityFile"]).name == face_csv.name
        assert "diffusivity" in rows(face_values.with_suffix(".faces.csv"))[0]
        verify_scalar(face_values)

        bad_face_header = root / "bad-face-header.csv"
        bad_face_header.write_text("face,value\n" + "\n".join(
            f"{row['face']},.01" for row in face_rows) + "\n")
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", checkpoint,
                        "--boundary", bc, "--output", root / "bad-face-header",
                        "--face-diffusivity", bad_face_header], r"diffusivity CSV header")
        bad_face_duplicate = root / "bad-face-duplicate.csv"
        bad_face_duplicate.write_text(face_csv.read_text() + f"{face_rows[0]['face']},.01\n")
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", checkpoint,
                        "--boundary", bc, "--output", root / "bad-face-duplicate",
                        "--face-diffusivity", bad_face_duplicate], r"duplicate diffusivity face")
        bad_face_value = root / "bad-face-value.csv"
        bad_face_value.write_text(face_csv.read_text().replace(",.02\n", ",0\n", 1))
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", checkpoint,
                        "--boundary", bc, "--output", root / "bad-face-value",
                        "--face-diffusivity", bad_face_value], r"face diffusivity must be positive")
        bad_face_missing = root / "bad-face-missing.csv"
        bad_face_missing.write_text("\n".join(face_csv.read_text().splitlines()[:-1]) + "\n")
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", checkpoint,
                        "--boundary", bc, "--output", root / "bad-face-missing",
                        "--face-diffusivity", bad_face_missing], r"missing diffusivity face")
        bad_face_nan = root / "bad-face-nan.csv"
        bad_face_nan.write_text(face_csv.read_text().replace(",.02\n", ",nan\n", 1))
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", checkpoint,
                        "--boundary", bc, "--output", root / "bad-face-nan",
                        "--face-diffusivity", bad_face_nan], r"expected finite number")
        bad_face_extra = root / "bad-face-extra.csv"
        bad_face_extra.write_text(face_csv.read_text().replace(",.02\n", ",.02,extra\n", 1))
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", checkpoint,
                        "--boundary", bc, "--output", root / "bad-face-extra",
                        "--face-diffusivity", bad_face_extra], r"row requires two fields")
        expect_failure([args.transport_cli, "--mesh", unit_mesh, "--output", root / "bad-variable-dt",
                        "--verification", "variable-sine", "--dt", ".01"],
                       r"sine verification is steady")
        expect_failure([args.transport_cli, "--mesh", unit_mesh, "--output", root / "bad-face-verification",
                        "--verification", "sine", "--face-diffusivity", face_csv],
                       r"face-diffusivity requires frozen carrier without restart or verification")

        tampered_flux = root / "tampered-flux"
        copy_prefix(physical, tampered_flux)
        flux_rows = rows(tampered_flux.with_suffix(".faces.csv"))
        flux_rows[0]["advectiveFlux"] = str(float(flux_rows[0]["advectiveFlux"]) + 1e-4)
        with tampered_flux.with_suffix(".faces.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=flux_rows[0].keys())
            writer.writeheader(); writer.writerows(flux_rows)
        expect_verifier_failure(tampered_flux, r"constitutive flux")
        tampered_previous = root / "tampered-previous"
        copy_prefix(physical, tampered_previous)
        previous_rows = rows(tampered_previous.with_suffix(".cells.csv"))
        previous_rows[0]["previous"] = str(float(previous_rows[0]["previous"]) + 1e-4)
        with tampered_previous.with_suffix(".cells.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=previous_rows[0].keys())
            writer.writeheader(); writer.writerows(previous_rows)
        expect_verifier_failure(tampered_previous, r"time term")
        tampered_source = root / "tampered-source"
        copy_prefix(physical, tampered_source)
        source_rows = rows(tampered_source.with_suffix(".cells.csv"))
        source_rows[0]["sourceIntegral"] = str(float(source_rows[0]["sourceIntegral"]) + 1e-4)
        with tampered_source.with_suffix(".cells.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=source_rows[0].keys())
            writer.writeheader(); writer.writerows(source_rows)
        expect_verifier_failure(tampered_source, r"source integral")
        tampered_diffusivity = root / "tampered-diffusivity"
        copy_prefix(face_values, tampered_diffusivity)
        variable_face_rows = rows(tampered_diffusivity.with_suffix(".faces.csv"))
        variable_face_rows[0]["diffusivity"] = str(float(variable_face_rows[0]["diffusivity"]) * 1.5)
        with tampered_diffusivity.with_suffix(".faces.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=variable_face_rows[0].keys())
            writer.writeheader(); writer.writerows(variable_face_rows)
        expect_verifier_failure(tampered_diffusivity, r"diffus|face|constitutive|metadata")
        tampered_face_source = root / "tampered-face-source"
        copy_prefix(face_values, tampered_face_source)
        face_source_rows = rows(tampered_face_source.with_suffix(".cells.csv"))
        face_source_rows[0]["sourceIntegral"] = str(float(face_source_rows[0]["sourceIntegral"]) + 1e-4)
        with tampered_face_source.with_suffix(".cells.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=face_source_rows[0].keys())
            writer.writeheader(); writer.writerows(face_source_rows)
        expect_verifier_failure(tampered_face_source, r"source integral")

        missing = root / "missing.csv"
        missing.write_text(bc.read_text().rsplit("\n", 2)[0] + "\n")
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", checkpoint,
                        "--boundary", missing, "--output", root / "bad-missing"],
                       r"missing scalar boundary face")
        duplicate = root / "duplicate.csv"
        duplicate.write_text(bc.read_text() + bc.read_text().splitlines()[1] + "\n")
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", checkpoint,
                        "--boundary", duplicate, "--output", root / "bad-duplicate"],
                       r"duplicate or internal boundary face")
        internal = root / "internal.csv"
        lines = bc.read_text().splitlines()
        internal_id = next(row["face"] for row in rows(flow.with_suffix(".faces.csv"))
                           if int(row["neighbour"]) >= 0)
        internal.write_text("\n".join(lines + [f"{internal_id},value,0,"]) + "\n")
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", checkpoint,
                        "--boundary", internal, "--output", root / "bad-internal"],
                       r"duplicate or internal boundary face")
        flux = root / "missing-inflow.csv"
        inflow_face = next(row["face"] for row in rows(flow.with_suffix(".faces.csv"))
                           if int(row["neighbour"]) < 0 and float(row["flux"]) < 0)
        flux.write_text("\n".join([lines[0]] + [
            (line.replace(",value,", ",flux,") if line.startswith(inflow_face + ",") else line)
            for line in lines[1:]]) + "\n")
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", checkpoint,
                        "--boundary", flux, "--output", root / "bad-inflow"],
                       r"inflow requires")

        altered = root / "altered.checkpoint"
        checkpoint_text = checkpoint.read_text()
        altered.write_text(re.sub(r"^(CELL 0 \S+ \S+ )\S+", r"\g<1>0.123",
                                  checkpoint_text, count=1, flags=re.MULTILINE))
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", altered,
                        "--boundary", bc, "--output", root / "bad-checkpoint"],
                       r"mesh mismatch in cell area")
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", checkpoint,
                        "--boundary", bc, "--output", root / "bad-nan", "--diffusivity", "nan"],
                       r"finite number")
        expect_failure([args.transport_cli, "--mesh", physical_mesh, "--flow-checkpoint", checkpoint,
                        "--boundary", bc, "--output", root / "bad-option", "--bogus", "1"],
                       r"unknown option")
        expect_failure([args.transport_cli, "--mesh", unit_mesh, "--output", root / "bad-mix",
                        "--verification", "sine", "--flow-checkpoint", checkpoint],
                       r"explicit carrier.*verification|choose explicit")
        expect_failure([args.transport_cli, "--mesh", unit_mesh, "--output", root / "bad-face-evolve",
                        "--face-diffusivity", face_csv, "--evolve-flow", "external",
                        "--dt", ".01", "--steps", "1"],
                       r"face-diffusivity requires frozen carrier")
        expect_failure([args.transport_cli, "--mesh", unit_mesh, "--output", root / "bad-dt",
                        "--verification", "decay", "--dt", "0", "--steps", "2"],
                       r"decay verification needs dt|invalid time settings")
        expect_failure([args.transport_cli, "--mesh", unit_mesh, "--output", root / "bad-steps",
                        "--verification", "sine", "--steps", "2"],
                       r"invalid time settings|sine verification is steady")

        limited = root / "limited"
        limited_result = transport(args.transport_cli, unit_mesh, limited,
                                   "--verification", "decay", "--diffusivity", ".01",
                                   "--dt", ".01", "--steps", "1", "--max-corrections", "1")
        assert limited_result.returncode == 2, (limited_result.stdout, limited_result.stderr)
        limited_summary = json.loads(limited.with_suffix(".json").read_text())
        assert limited_summary["converged"] is False

    print("Scalar transport CLI: verification, physical checkpoint transport, balance, and rejection paths verified.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--transport-cli", "--cli", dest="transport_cli", required=True)
    parser.add_argument("--mesh-cli", required=True)
    parser.add_argument("--flow-cli", required=True)
    main(parser.parse_args())
