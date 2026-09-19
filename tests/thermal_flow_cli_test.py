#!/usr/bin/env python3
"""Black-box synchronized thermal-flow CLI and joint-checkpoint checks.

The authoritative restart is a single file with this format:
``CARTMESH2D_THERMAL_CHECKPOINT 1`` / ``COUPLING new-time-flux-Euler-v1`` /
``THERMAL_CONFIG`` / ``SOURCES`` / ``BOUNDARIES`` / ``SCALAR`` / ``FLOW``;
the FLOW section embeds the complete versioned ``CARTMESH2D_FLOW_CHECKPOINT`` record,
including geometry, TIME, U/V/P and FLUX.  The file is replaced atomically.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import shutil
import subprocess
import tempfile
import time


sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'/'verification'))
from verify_scalar_transport import verify

def run(command, timeout=30):
    return subprocess.run([str(x) for x in command], text=True,
                          capture_output=True, timeout=timeout)


def square(path):
    path.write_text("0 0\n1 0\n1 1\n0 1\n")


def native_mesh(mesh_cli, root):
    boundary = root / "unit.xy"
    square(boundary)
    prefix = root / "unit"
    case = root / "unit-openfoam"
    result = run([mesh_cli, boundary, prefix, "4", str(1 / 14), ".1",
                  "interior", case, "4", "0"])
    assert result.returncode == 0, (result.stdout, result.stderr)
    mesh = prefix.with_suffix(".solver.cm2d")
    assert mesh.exists() and case.exists()
    return mesh


def invoke(cli, mesh, output, steps, extra=(), restart=None, timeout=30):
    command = [cli, "--mesh", mesh, "--output", output,
               "--verification", "thermal-vortex", "--diffusivity", ".2",
               "--dt", ".01", "--steps", str(steps),
               "--flow-nu", ".1", "--flow-convection", "limited-linear",
               "--convection", "limited-linear"]
    if restart is not None:
        command += ["--restart", restart]
    return run(command + list(extra), timeout)


def csv_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def assert_csv_balance(prefix):
    cells = csv_rows(prefix.with_suffix(".cells.csv"))
    faces = csv_rows(prefix.with_suffix(".faces.csv"))
    assert len(cells) == 196 and faces
    balance = [float(row["temporalIntegral"]) - float(row["sourceIntegral"])
               for row in cells]
    for row in faces:
        flux = float(row["advectiveFlux"]) + float(row["diffusiveFlux"])
        owner = int(row["owner"])
        balance[owner] += flux
        neighbour = int(row["neighbour"])
        if neighbour >= 0:
            balance[neighbour] -= flux
    assert max(abs(value) for value in balance) < 2e-7
    for row in cells:
        for key in ("value", "previous", "temporalIntegral", "sourceIntegral"):
            assert math.isfinite(float(row[key]))
    return cells


def checkpoint_section(text, name):
    match = re.search(r"^" + re.escape(name) + r"\s+(\d+)\s+(.*)$",
                      text, re.MULTILINE)
    assert match, name
    return [float(value) for value in match.group(2).split()]


def expect_failure(command, pattern):
    result = command if isinstance(command, subprocess.CompletedProcess) else run(command)
    assert result.returncode == 1, (result.returncode, result.stdout, result.stderr)
    assert re.search(pattern, result.stderr, re.IGNORECASE), result.stderr
    return result


def main(args):
    with tempfile.TemporaryDirectory(prefix="cartmesh-thermal-flow-") as name:
        root = Path(name)
        mesh = native_mesh(args.mesh_cli, root)

        continuous = root / "continuous"
        result = invoke(args.transport_cli, mesh, continuous, 4)
        assert result.returncode == 0, (result.stdout, result.stderr)
        summary = json.loads(continuous.with_suffix(".json").read_text())
        assert summary["cells"] == 196 and summary["faces"] > 196
        assert math.isclose(summary["time"], .04) and math.isclose(summary["carrierTime"], .04)
        assert math.isclose(summary["acceptedTime"], .04)
        continuous_cells = assert_csv_balance(continuous)
        history = csv_rows(continuous.with_suffix(".thermal-history.csv"))
        assert len(history) == 4 and all(row["accepted"] == "1" for row in history)
        assert [float(row["time"]) for row in history] == [.01, .02, .03, .04]
        for row in history:
            for key in ("heatContent", "scalarGlobalBalance", "flowContinuity"):
                assert math.isfinite(float(row[key]))
        heats = [float(row["heatContent"]) for row in history]
        assert all(a > b for a, b in zip(heats, heats[1:]))
        values = [float(row["value"]) for row in continuous_cells]
        assert max(values) < 1 and min(values) >= -1e-10
        assert float(summary["l2Error"]) < .1
        assert continuous.with_suffix(".thermal.checkpoint").read_text().startswith(
            "CARTMESH2D_THERMAL_CHECKPOINT 1\nCOUPLING new-time-flux-Euler-v1\n")

        first = root / "split-first"
        assert invoke(args.transport_cli, mesh, first, 2).returncode == 0
        split = root / "split-final"
        assert invoke(args.transport_cli, mesh, split, 2,
                      restart=first.with_suffix(".thermal.checkpoint")).returncode == 0
        assert math.isclose(json.loads(split.with_suffix(".json").read_text())["time"], .04)
        for suffix in (".cells.csv", ".faces.csv", ".thermal.checkpoint"):
            assert split.with_suffix(suffix).read_bytes() == continuous.with_suffix(suffix).read_bytes(), suffix
        assert checkpoint_section(first.with_suffix(".thermal.checkpoint").read_text(), "U") != checkpoint_section(
            continuous.with_suffix(".thermal.checkpoint").read_text(), "U")
        assert checkpoint_section(first.with_suffix(".thermal.checkpoint").read_text(), "FLUX") != checkpoint_section(
            continuous.with_suffix(".thermal.checkpoint").read_text(), "FLUX")

        verify(continuous); verify(split)
        tampered=root/'lagged-clock'
        for suffix in ('.json','.cells.csv','.faces.csv','.history.csv','.thermal-history.csv','.thermal.checkpoint'):
            shutil.copyfile(str(split)+suffix,str(tampered)+suffix)
        metadata=json.loads(tampered.with_suffix('.json').read_text());metadata['carrierTime']-=.01
        tampered.with_suffix('.json').write_text(json.dumps(metadata))
        try:
            verify(tampered)
        except ValueError as error:
            assert 'lagged carrier time' in str(error),str(error)
        else:
            raise AssertionError('independent reader accepted lagged flow clock')

        checkpoint = continuous.with_suffix(".thermal.checkpoint")
        before = hashlib.sha256(checkpoint.read_bytes()).digest()
        failed = invoke(args.transport_cli, mesh, continuous, 1,
                        extra=("--flow-max-iterations", "1"), restart=checkpoint)
        assert failed.returncode == 2
        assert json.loads(continuous.with_suffix(".json").read_text())["failedStage"] == "flow"
        assert hashlib.sha256(checkpoint.read_bytes()).digest() == before
        assert json.loads(continuous.with_suffix(".json").read_text())["converged"] is False
        resumed = root / "resumed"
        assert invoke(args.transport_cli, mesh, resumed, 2, restart=checkpoint).returncode == 0
        assert math.isclose(json.loads(resumed.with_suffix(".json").read_text())["time"], .06)

        scalar_failed = invoke(args.transport_cli, mesh, continuous, 1,
                               extra=("--max-corrections", "1"), restart=checkpoint)
        assert scalar_failed.returncode == 2
        scalar_summary = json.loads(continuous.with_suffix(".json").read_text())
        assert scalar_summary["failedStage"] == "scalar"
        assert math.isclose(scalar_summary["acceptedTime"], .04)
        assert hashlib.sha256(checkpoint.read_bytes()).digest() == before
        scalar_resumed = root / "scalar-resumed"
        assert invoke(args.transport_cli, mesh, scalar_resumed, 1,
                      restart=checkpoint).returncode == 0
        assert math.isclose(json.loads(scalar_resumed.with_suffix(".json").read_text())["time"], .05)

        # Setup changes are exact joint-checkpoint compatibility errors.
        expect_failure(invoke(args.transport_cli, mesh, root / "bad-d", 1,
                              extra=("--diffusivity", ".21"), restart=checkpoint),
                       r"thermal diffusivity")
        expect_failure(invoke(args.transport_cli, mesh, root / "bad-nu", 1,
                              extra=("--flow-nu", ".11"), restart=checkpoint),
                       r"incompatible configuration|configuration nu")
        expect_failure(invoke(args.transport_cli, mesh, root / "bad-scheme", 1,
                              extra=("--convection", "upwind"), restart=checkpoint),
                       r"thermal convection")

        text = checkpoint.read_text()
        truncated = root / "truncated.checkpoint"
        truncated.write_text(text[: text.rfind("\n") // 2])
        expect_failure(invoke(args.transport_cli, mesh, root / "bad-truncated", 1,
                              restart=truncated), r"truncated|unexpected end")
        altered = root / "altered.checkpoint"
        altered.write_text(re.sub(r"^(CELL 0 \S+ \S+ )\S+", r"\g<1>0.123",
                                  text, count=1, flags=re.MULTILINE))
        expect_failure(invoke(args.transport_cli, mesh, root / "bad-geometry", 1,
                              restart=altered), r"mesh mismatch in cell area")

        # A bad input fails before the accepted file can be replaced.
        invalid_before = hashlib.sha256(checkpoint.read_bytes()).digest()
        expect_failure(invoke(args.transport_cli, mesh, continuous, 1,
                              extra=("--diffusivity", "nan"), restart=checkpoint),
                       r"expected finite number")
        assert hashlib.sha256(checkpoint.read_bytes()).digest() == invalid_before

        # Observe a live process after it has published a positive-time joint
        # checkpoint, then terminate promptly and restart from that intact file.
        live = root / "cancelled"
        command = [args.transport_cli, "--mesh", mesh, "--output", live,
                   "--verification", "thermal-vortex", "--diffusivity", ".2",
                   "--dt", ".01", "--steps", "1000", "--flow-nu", ".1",
                   "--flow-convection", "limited-linear", "--convection", "limited-linear"]
        proc = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        observed_time = 0.0
        live_checkpoint = live.with_suffix(".thermal.checkpoint")
        try:
            deadline = time.monotonic() + 5
            while time.monotonic() < deadline:
                if live_checkpoint.exists() and proc.poll() is None:
                    match = re.search(r"^TIME\s+([^\s]+)$", live_checkpoint.read_text(), re.MULTILINE)
                    if match and float(match.group(1))>0:
                        observed_time = float(match.group(1))
                        break
                time.sleep(.05)
            assert observed_time>0 and proc.poll() is None, "no live positive-time checkpoint observed"
            proc.terminate()
            proc.communicate(timeout=2)
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.communicate(timeout=2)
        final_text = live_checkpoint.read_text()
        final_match = re.search(r"^TIME\s+([^\s]+)$", final_text, re.MULTILINE)
        assert final_match and float(final_match.group(1)) >= observed_time
        restart_after_cancel = root / "after-cancel"
        assert invoke(args.transport_cli, mesh, restart_after_cancel, 1,
                      restart=live_checkpoint).returncode == 0
        final_summary = json.loads(restart_after_cancel.with_suffix(".json").read_text())
        assert math.isclose(final_summary["acceptedTime"],float(final_match.group(1))+.01,rel_tol=1e-12,abs_tol=1e-14)

    print("Thermal-flow CLI: synchronized evolution, restart identity, CSV balances, rejection, and cancellation verified.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--transport-cli", required=True)
    parser.add_argument("--mesh-cli", required=True)
    main(parser.parse_args())
