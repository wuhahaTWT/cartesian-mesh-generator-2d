#!/usr/bin/env python3
"""Focused negative and sampler tests for the published cavity audit.

These tests deliberately stub the expensive native readback.  They still create
the manifest, producer, mesh and cell artifacts so that provenance checks are
real; no solver or mesh generation process is started.
"""
from __future__ import annotations

import copy
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools" / "verification"))
import audit_cavity_reference as audit  # noqa: E402
import verify_native_flow as native  # noqa: E402


class CavityReferenceAuditTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.binary = self.root / "producer"
        self.binary.write_bytes(b"producer-for-test")
        self.mesh = self.root / "mesh.solver.cm2d"
        self.mesh.write_text("synthetic mesh", encoding="utf-8")
        self.reference = self.root / "reference.json"
        rows_u = [{"coordinate": i / 16, "value": i / 16, "estimatedError": .001}
                  for i in range(1, 16)]
        rows_v = [{"coordinate": i / 16, "value": 0., "estimatedError": .001}
                  for i in range(1, 16)]
        self.reference.write_text(json.dumps({
            "format": "cartmesh2d-published-cavity-reference-v1",
            "reynolds": 100, "lidSpeed": 1, "domain": [0, 0, 1, 1],
            "doi": "10.1590/S1678-58782009000300004", "u": rows_u, "v": rows_v,
        }), encoding="utf-8")

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_compare_reference_reproduces_affine_field(self) -> None:
        cells = [{"x": (i + .5) / 5, "y": (j + .5) / 5,
                  "u": (j + .5) / 5, "v": 0.}
                 for i in range(5) for j in range(5)]
        result = audit.compare_reference(cells, json.loads(self.reference.read_text()))
        self.assertLess(result["centrelineRmse"], 1e-12)
        self.assertLess(result["centrelineMaxError"], 1e-12)

    def test_real_reference_has_center_points_and_known_u_value(self) -> None:
        if not audit.REFERENCE.is_file():
            self.skipTest("bundled published reference has not landed yet")
        reference = audit.load_reference(audit.REFERENCE)
        for component in ("u", "v"):
            centres = [row for row in reference[component] if row["coordinate"] == .5]
            self.assertEqual(len(centres), 1)
        u_low = next(row for row in reference["u"] if row["coordinate"] == .0625)
        self.assertAlmostEqual(u_low["value"], -0.041974991, places=12)

    def test_reference_rejects_missing_nonfinite_and_bad_coordinates(self) -> None:
        valid = json.loads(self.reference.read_text())
        for mutate in (
            lambda ref: ref["u"].pop(),
            lambda ref: ref["u"][2].update(coordinate=0.2),
            lambda ref: ref["v"][4].update(value=float("nan")),
            lambda ref: ref["v"][4].update(estimatedError=-1),
        ):
            candidate = copy.deepcopy(valid)
            mutate(candidate)
            path = self.root / f"bad-{len(list(self.root.glob('bad-*.json')))}.json"
            path.write_text(json.dumps(candidate, allow_nan=True), encoding="utf-8")
            with self.assertRaises(ValueError):
                audit.load_reference(path)

    def _write_cells(self, prefix: Path, offset: float = 0.0) -> None:
        path = Path(str(prefix) + ".cells.csv")
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=("cell", "x", "y", "u", "v", "speed"))
            writer.writeheader()
            for number, (i, j) in enumerate(( (i, j) for i in range(5) for j in range(5) )):
                x, y = (i + .5) / 5, (j + .5) / 5
                writer.writerow({"cell": number, "x": x, "y": y,
                                 "u": y + offset, "v": offset, "speed": 0})

    def _command(self, prefix: Path, *, convection: str = "upwind", method: str = "ic0",
                 tolerance: str = "1e-08", transient: bool = False) -> list[str]:
        command = [str(self.binary), "--mesh", str(self.mesh), "--output", str(prefix),
                   "--case", "cavity", "--nu", "0.01", "--speed", "1",
                   "--convection", convection, "--pressure-preconditioner", method,
                   "--tolerance", tolerance, "--max-iterations", "1000"]
        if transient:
            command += ["--steps", "2"]
        return command

    def _manifest(self, count: int = 3, *, offsets: tuple[float, ...] | None = None) -> Path:
        offsets = offsets or (0.3, 0.2, 0.1)
        runs = []
        for index in range(count):
            prefix = self.root / f"run-{index}"
            self._write_cells(prefix, offsets[index] if index < len(offsets) else 0.0)
            command = self._command(prefix)
            runs.append({"label": f"cavity-l{index + 4}", "returnCode": 0, "timedOut": False,
                         "mesh": str(self.mesh), "prefix": str(prefix),
                         "meshSha256": native.sha256_file(self.mesh),
                         "binarySha256": native.sha256_file(self.binary), "command": command})
        path = self.root / "runs.json"
        path.write_text(json.dumps(runs), encoding="utf-8")
        return path

    def _fake_verify(self, mesh, prefix, case, nu, speed, controls):
        label = Path(prefix).name
        h = {"run-0": 1.0, "run-1": .5, "run-2": .25}[label]
        method = "aggregation" if label == "run-2" and getattr(self, "aggregation_run", False) else "ic0"
        bounds = getattr(self, "fake_bounds", [0., 0., 1., 1.])
        return {"valid": True, "case": "cavity", "issues": [],
                "meshMeasurement": {"bounds": bounds, "characteristicH": h},
                "benchmark": {"centrelineRmse": .2},
                "native": {"case": "cavity", "nu": .01, "speed": 1.,
                           "convection": "upwind", "pressurePreconditioner": method,
                           "tolerance": 1e-8, "viscousStress": "symmetric"}}

    def test_audit_keeps_ghia_failure_when_published_trend_passes(self) -> None:
        # The synthetic cell offsets decrease .3 -> .2 -> .1, so the published
        # trend passes.  The mocked existing Ghia metric is deliberately .2,
        # .2, .2 and therefore remains a visible failure only if changed below.
        manifest = self._manifest()
        def verify(*args):
            result = self._fake_verify(*args)
            result["benchmark"]["centrelineRmse"] = {
                "run-0": .1, "run-1": .2, "run-2": .2}[Path(args[1]).name]
            return result
        with mock.patch.object(native, "verify_case", side_effect=verify):
            report = audit.audit_runs(manifest, self.reference)
        self.assertFalse(report["ghiaSequence"]["valid"])
        self.assertTrue(report["publishedSequence"]["valid"])
        self.assertFalse(report["allChecksPassed"])

    def test_independent_audit_failure_is_not_hidden(self) -> None:
        manifest = self._manifest()
        def verify(*args):
            result = self._fake_verify(*args)
            result["valid"] = False
            result["issues"] = ["controlled independent failure"]
            return result
        with mock.patch.object(native, "verify_case", side_effect=verify):
            report = audit.audit_runs(manifest, self.reference)
        self.assertFalse(report["independentAuditsPassed"])
        self.assertFalse(report["allChecksPassed"])
        self.assertIn("independent geometry/physics audit failed", report["issues"][0])

    def test_rejects_hash_mismatch_and_duplicate_prefix(self) -> None:
        manifest = self._manifest()
        runs = json.loads(manifest.read_text())
        runs[0]["meshSha256"] = "0" * 64
        manifest.write_text(json.dumps(runs), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Mesh hash changed"):
            audit.audit_runs(manifest, self.reference)
        manifest = self._manifest()
        runs = json.loads(manifest.read_text())
        runs[0]["binarySha256"] = "f" * 64
        manifest.write_text(json.dumps(runs), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Producer binary changed"):
            audit.audit_runs(manifest, self.reference)

        manifest = self._manifest()
        runs = json.loads(manifest.read_text())
        runs[1]["prefix"] = runs[0]["prefix"]
        manifest.write_text(json.dumps(runs), encoding="utf-8")
        with mock.patch.object(native, "verify_case", side_effect=self._fake_verify):
            with self.assertRaisesRegex(ValueError, "Duplicate result prefix"):
                audit.audit_runs(manifest, self.reference)

    def test_rejects_transient_and_command_summary_control_mismatch(self) -> None:
        manifest = self._manifest()
        runs = json.loads(manifest.read_text())
        runs[0]["command"] = self._command(Path(runs[0]["prefix"]), transient=True)
        manifest.write_text(json.dumps(runs), encoding="utf-8")
        with self.assertRaises(ValueError):
            audit.audit_runs(manifest, self.reference)

        manifest = self._manifest()
        runs = json.loads(manifest.read_text())
        runs[0]["command"] = self._command(Path(runs[0]["prefix"]), convection="limited-linear")
        manifest.write_text(json.dumps(runs), encoding="utf-8")
        with mock.patch.object(native, "verify_case", side_effect=self._fake_verify):
            with self.assertRaises(ValueError):
                audit.audit_runs(manifest, self.reference)

    def test_rejects_mixed_numerical_controls(self) -> None:
        manifest = self._manifest()
        runs = json.loads(manifest.read_text())
        runs[2]["command"] = self._command(Path(runs[2]["prefix"]), method="aggregation")
        manifest.write_text(json.dumps(runs), encoding="utf-8")
        self.aggregation_run = True
        with mock.patch.object(native, "verify_case", side_effect=self._fake_verify):
            with self.assertRaisesRegex(ValueError, "mixes binary or numerical controls"):
                audit.audit_runs(manifest, self.reference)

    def test_unit_square_rounding_is_accepted_but_real_size_change_is_rejected(self) -> None:
        manifest = self._manifest()
        self.fake_bounds = [0., 0., 1.0000000000000002, 1.]
        with mock.patch.object(native, "verify_case", side_effect=self._fake_verify):
            report = audit.audit_runs(manifest, self.reference)
        self.assertEqual(len(report["cases"]), 3)

        manifest = self._manifest()
        self.fake_bounds = [0., 0., 1.01, 1.]
        with mock.patch.object(native, "verify_case", side_effect=self._fake_verify):
            with self.assertRaisesRegex(ValueError, "unit square"):
                audit.audit_runs(manifest, self.reference)


if __name__ == "__main__":
    unittest.main()
