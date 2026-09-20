#!/usr/bin/env python3
"""Focused runset comparison tests; no solver invocation is required."""
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("compare_transient_steps", ROOT / "tools/verification/compare_transient_steps.py")
COMPARE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(COMPARE)


class TransientStepComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cartmesh-step-compare-")
        self.root = Path(self.temp.name)
        self.mesh = self.root / "mesh.solver.cm2d"
        self.mesh.write_text("mesh fixture\n")
        self.binary = self.root / "flow-cli"
        self.binary.write_bytes(b"binary fixture")
        self.prefixes = []
        for index, (dt, u) in enumerate(((.04, 1.0), (.02, .75), (.01, .625))):
            prefix = self.root / f"run-{index}" / "result"
            prefix.parent.mkdir()
            prefix.with_suffix(".cells.csv").write_text(
                "cell,x,y,area,u,v\n0,0,0,1," + str(u) + ",0\n1,1,0,1," + str(u + .1) + ",0\n")
            self.prefixes.append((prefix, dt))
        self.original_verify = COMPARE.audit.verify
        self.original_sha = COMPARE.audit.native.sha256_file
        def fake_verify(mesh, prefix, report):
            result = {"valid": True, "case": "taylor-green", "time": .2,
                      "dt": next(dt for p, dt in self.prefixes if p == prefix),
                      "controls": {"nu": .1, "speed": 1.0, "tolerance": 1e-9,
                                   "convection": "limited-linear", "viscousStress": "symmetric",
                                   "temporalFaceInterpolation": "old-and-iteration-flux-defect-skew-corrected-v2",
                                   "velocityRelaxation": 1.0}}
            dt = result['dt']
            result['history'] = [{'time': (i+1)*dt} for i in range(round(.2/dt))]
            Path(report).write_text(json.dumps(result))
            return result
        COMPARE.audit.verify = fake_verify

    def tearDown(self):
        COMPARE.audit.verify = self.original_verify
        COMPARE.audit.native.sha256_file = self.original_sha
        self.temp.cleanup()

    def runset(self):
        return {"mesh": str(self.mesh), "meshSha256": hashlib.sha256(self.mesh.read_bytes()).hexdigest(),
                "binarySha256": hashlib.sha256(self.binary.read_bytes()).hexdigest(),
                "runs": [{"command": [str(self.binary), "--mesh", str(self.mesh), "--output", str(prefix),
                                         "--case", "taylor-green", "--nu", ".1", "--speed", "1",
                                         "--convection", "limited-linear", "--time-step", str(dt),
                                         "--steps", str(round(.2 / dt)), "--tolerance", "1e-9"],
                          "dt": dt, "endTime": .2, "returnCode": 0, "valid": True}
                         for prefix, dt in self.prefixes]}

    def write_runset(self, data):
        path = self.root / "runs.json"
        path.write_text(json.dumps(data))
        return path

    def test_valid_mapping_and_diagnostic_order(self):
        report = COMPARE.compare(self.write_runset(self.runset()), self.root / "report.json")
        self.assertTrue(report["valid"])
        self.assertEqual(report["cellIds"], 2)
        self.assertEqual(len(report["adjacentDistances"]), 2)
        self.assertTrue(report["adjacentDistances"][1]["observedOrder"] > 0)

    def test_retained_binary_snapshot_requires_original_hash(self):
        runset = self.write_runset(self.runset())
        retained = self.root / 'saved-cli'
        retained.write_bytes(self.binary.read_bytes())
        self.binary.write_bytes(b'new build')
        with self.assertRaises(COMPARE.audit.native.VerificationError):
            COMPARE.compare(runset, self.root / 'report.json')
        result = COMPARE.compare(runset, self.root / 'report.json', retained)
        self.assertEqual(result['verifiedBinarySnapshot'], str(retained.resolve()))
        retained.write_bytes(b'incorrect snapshot')
        with self.assertRaises(COMPARE.audit.native.VerificationError):
            COMPARE.compare(runset, self.root / 'report.json', retained)

    def test_mutated_mesh_hash_rejected(self):
        data = self.runset()
        self.mesh.write_text("mutated mesh\n")
        with self.assertRaises(COMPARE.audit.native.VerificationError):
            COMPARE.compare(self.write_runset(data), self.root / "report.json")

    def test_config_mismatch_rejected(self):
        data = self.runset()
        calls = [0]
        def mismatched(mesh, prefix, report):
            calls[0] += 1
            result = {"valid": True, "case": "taylor-green", "time": .2,
                      "dt": next(dt for p, dt in self.prefixes if p == prefix),
                      "controls": {"nu": .2 if calls[0] == 2 else .1, "speed": 1.0,
                                   "tolerance": 1e-9, "convection": "limited-linear",
                                   "viscousStress": "symmetric",
                                   "temporalFaceInterpolation": "old-and-iteration-flux-defect-skew-corrected-v2",
                                   "velocityRelaxation": 1.0}}
            dt = result['dt']
            result['history'] = [{'time': (i+1)*dt} for i in range(round(.2/dt))]
            Path(report).write_text(json.dumps(result))
            return result
        COMPARE.audit.verify = mismatched
        with self.assertRaises(COMPARE.audit.native.VerificationError):
            COMPARE.compare(self.write_runset(data), self.root / "report.json")

    def test_wrong_dt_label_and_restart_are_rejected(self):
        data = self.runset()
        data["runs"][0]["dt"] = .02
        with self.assertRaises(COMPARE.audit.native.VerificationError):
            COMPARE.compare(self.write_runset(data), self.root / "report.json")

    def test_command_nu_and_verified_dt_mismatches_are_rejected(self):
        data = self.runset()
        data["runs"][1]["command"][data["runs"][1]["command"].index("--nu") + 1] = ".2"
        with self.assertRaises(COMPARE.audit.native.VerificationError):
            COMPARE.compare(self.write_runset(data), self.root / "report.json")
        data = self.runset()
        original = COMPARE.audit.verify
        def wrong_dt(mesh, prefix, report):
            result = original(mesh, prefix, report)
            if prefix == self.prefixes[1][0]:
                result["dt"] = .03
            dt = result['dt']
            result['history'] = [{'time': (i+1)*dt} for i in range(round(.2/dt))]
            Path(report).write_text(json.dumps(result))
            return result
        COMPARE.audit.verify = wrong_dt
        with self.assertRaises(COMPARE.audit.native.VerificationError):
            COMPARE.compare(self.write_runset(data), self.root / "report.json")
        data = self.runset()
        data["runs"][0]["command"] += ["--restart", str(self.root / "restart.checkpoint")]
        with self.assertRaises(COMPARE.audit.native.VerificationError):
            COMPARE.compare(self.write_runset(data), self.root / "report.json")

    def test_cell_id_mapping_mismatch_rejected(self):
        bad = self.prefixes[2][0].with_suffix(".cells.csv")
        bad.write_text("cell,x,y,area,u,v\n0,0,0,1,.6,0\n2,1,0,1,.7,0\n")
        with self.assertRaises(COMPARE.audit.native.VerificationError):
            COMPARE.compare(self.write_runset(self.runset()), self.root / "report.json")


if __name__ == "__main__":
    unittest.main()
