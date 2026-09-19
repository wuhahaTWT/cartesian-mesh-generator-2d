#!/usr/bin/env python3
"""Synthetic contract tests for fixed-mesh thermal time refinement."""
import importlib.util
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools" / "verification"))
SPEC = importlib.util.spec_from_file_location("verify_thermal_time", ROOT / "tools/verification/verify_thermal_time.py")
TIME = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(TIME)


def case(dt, steps, scalar_error, velocity_error, **updates):
    reference = {
        "meshSha256": "mesh-fixed-sha",
        "cells": 196, "h": 1 / 14, "dt": dt, "steps": steps, "time": dt * steps,
        "nu": .1, "diffusivity": .02, "speed": 1., "velocityRelaxation": .8,
        "flowConvection": "limited-linear", "scalarConvection": "limited-linear",
        "viscousStress": "symmetric", "pressurePreconditioner": "aggregation",
        "coupledFlowTolerance": 1e-9, "coupledFlowToleranceStatus": "explicit-thermal-flow-json",
        "flowTolerance": 1e-9, "flowToleranceStatus": "explicit-flow-json",
        "scalarRelativeTolerance": 1e-9, "scalarAbsoluteTolerance": 1e-12,
        "scalarCellTolerance": 1e-9, "carrierCheckpointByteIdentical": True,
        "scalarContinuousL2": scalar_error, "velocityContinuousL2": velocity_error,
    }
    reference.update(updates)
    return {"valid": True, "reference": reference}


class ThermalTimeVerifierTests(unittest.TestCase):
    def good(self):
        return [case(.1, 1, .30, .20), case(.05, 2, .20, .12), case(.025, 4, .10, .06)]

    def test_monotonic_continuous_reference_errors_pass(self):
        result = TIME.time_series_checks(self.good())
        self.assertTrue(result["valid"])
        self.assertTrue(result["decreases"]["scalarContinuousL2"])
        self.assertTrue(result["decreases"]["velocityContinuousL2"])

    def test_nondecreasing_errors_fail_without_threshold_relaxation(self):
        cases = self.good(); cases[-1]["reference"]["velocityContinuousL2"] = .12
        result = TIME.time_series_checks(cases)
        self.assertFalse(result["valid"])

    def test_rejects_structural_and_audit_failures(self):
        mutations = []
        bad = self.good(); bad[1]["reference"]["meshSha256"] = "other"; mutations.append(bad)
        bad = self.good(); bad[1]["reference"]["steps"] = 3; mutations.append(bad)
        bad = self.good(); bad[1]["reference"].pop("flowConvection"); mutations.append(bad)
        bad = self.good(); bad[1]["reference"]["coupledFlowToleranceStatus"] = "legacy-unknown"; mutations.append(bad)
        bad = self.good(); bad[1]["valid"] = False; mutations.append(bad)
        bad = self.good(); bad[1]["reference"]["scalarContinuousL2"] = math.nan; mutations.append(bad)
        for candidate in mutations:
            with self.assertRaises(ValueError): TIME.time_series_checks(candidate)

    def test_rejects_duplicate_or_reversed_time_and_mixed_controls(self):
        bad = self.good(); bad[1]["reference"]["dt"] = .1; bad[1]["reference"]["steps"] = 1
        with self.assertRaises(ValueError): TIME.time_series_checks(bad)
        bad = self.good()[::-1]
        with self.assertRaises(ValueError): TIME.time_series_checks(bad)
        for key in ("flowConvection", "scalarConvection", "viscousStress", "pressurePreconditioner",
                    "nu", "diffusivity", "speed", "velocityRelaxation", "coupledFlowTolerance",
                    "flowTolerance", "scalarRelativeTolerance"):
            bad = self.good(); bad[1]["reference"][key] = ("upwind" if "Convection" in key else
                                                               "laplacian" if key == "viscousStress" else
                                                               "ic0" if key == "pressurePreconditioner" else
                                                               bad[1]["reference"][key] * 2)
            with self.assertRaises(ValueError): TIME.time_series_checks(bad)

    def test_rejects_coupled_and_standalone_flow_tolerance_mismatch(self):
        bad = self.good(); bad[1]["reference"]["flowTolerance"] = 2e-9
        with self.assertRaises(ValueError): TIME.time_series_checks(bad)
        bad = self.good(); bad[1]["reference"]["coupledFlowTolerance"] = 2e-9
        with self.assertRaises(ValueError): TIME.time_series_checks(bad)
        bad = self.good()
        for item in bad: item["reference"]["flowTolerance"] = 2e-9
        with self.assertRaises(ValueError): TIME.time_series_checks(bad)

    def test_rejects_incomplete_or_unaccepted_trajectory(self):
        with self.assertRaises(ValueError): TIME.time_series_checks(self.good()[:2])
        for key, value in (("dt", 0.), ("dt", math.nan), ("steps", True),
                           ("carrierCheckpointByteIdentical", False)):
            bad = self.good(); bad[1]["reference"][key] = value
            with self.assertRaises(ValueError): TIME.time_series_checks(bad)



if __name__ == "__main__":
    unittest.main()
