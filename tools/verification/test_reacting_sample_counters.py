#!/usr/bin/env python3
"""Exact counter/state semantics for accepted samples and stopped trials."""
import copy
import unittest

from verify_reacting_flame import audit_sample_counters


class SampleCounterTests(unittest.TestCase):
    def setUp(self):
        state = {"time": 2e-6, "steps": 3, "U": [[1, 0, 0, 1, 1]]}
        integration = {"rhsCalls": 20, "rejectedRhsCalls": 0, "jacobianEvaluations": 1,
                       "maximumNonlinearIterations": 3, "relativeTolerance": 1e-7,
                       "sampleResidualEvaluations": 2, "acceptedByBdfOrder": [0, 3],
                       "maximumLocalErrorNorm": .5, "minimumNewtonFraction": 1.,
                       "maximumReflectedSpeciesScaledChange": 0., "lastDampedTrialFailure": "",
                       "canceled": False, "samplingFailure": ""}
        self.sample = {"state": state, "integration": integration}
        self.field = {"complete": True, "failure": "", "final": copy.deepcopy(state),
                      "integration": copy.deepcopy(integration)}

    def stopped(self):
        self.field.update(complete=False, failure="canceled during a trial solve")
        self.field["integration"].update(rhsCalls=24, rejectedRhsCalls=1, jacobianEvaluations=2,
                                         canceled=True, lastDampedTrialFailure="trial diagnostic")

    def test_complete_snapshot_is_exact(self):
        audit_sample_counters(self.field, self.sample, 20, True)
        self.sample["integration"]["rhsCalls"] = 21
        with self.assertRaises(AssertionError):
            audit_sample_counters(self.field, self.sample, 20, True)

    def test_terminal_counters_on_last_accepted_state(self):
        self.stopped()
        self.sample["integration"] = copy.deepcopy(self.field["integration"])
        audit_sample_counters(self.field, self.sample, 20, True)

    def test_regular_sample_before_stopped_trial(self):
        self.stopped()
        audit_sample_counters(self.field, self.sample, 20, True)

    def test_no_accepted_step_before_cancel(self):
        self.sample["state"].update(time=0, steps=0)
        self.field["final"] = copy.deepcopy(self.sample["state"])
        self.sample["integration"].update(rhsCalls=0, jacobianEvaluations=0,
                                           acceptedByBdfOrder=[0, 0], maximumLocalErrorNorm=0.,
                                           sampleResidualEvaluations=1)
        self.field["integration"] = copy.deepcopy(self.sample["integration"])
        self.stopped()
        audit_sample_counters(self.field, self.sample, 0, True)

    def test_intermediate_counter_cannot_use_terminal_exception(self):
        self.stopped()
        self.sample["integration"] = copy.deepcopy(self.field["integration"])
        with self.assertRaises(AssertionError):
            audit_sample_counters(self.field, self.sample, 20, False)

    def test_terminal_exception_requires_same_accepted_state_and_failure(self):
        for corruption in ["state", "failure", "complete"]:
            with self.subTest(corruption=corruption):
                self.setUp(); self.stopped()
                self.sample["integration"] = copy.deepcopy(self.field["integration"])
                if corruption == "state": self.sample["state"]["time"] += 1e-6
                elif corruption == "failure": self.field["failure"] = ""
                else: self.field["complete"] = True
                with self.assertRaises(AssertionError):
                    audit_sample_counters(self.field, self.sample, 20, True)

    def test_controls_and_accepted_statistics_cannot_change(self):
        for key, value in [("relativeTolerance", 1e-4), ("acceptedByBdfOrder", [0, 4]),
                           ("maximumLocalErrorNorm", .9), ("sampleResidualEvaluations", 3)]:
            with self.subTest(key=key):
                self.setUp(); self.stopped(); self.field["integration"][key] = value
                with self.assertRaises(AssertionError):
                    audit_sample_counters(self.field, self.sample, 20, True)

    def test_intermediate_or_decreasing_terminal_statistics_rejected(self):
        self.stopped(); self.sample["integration"]["rhsCalls"] = 22
        with self.assertRaises(AssertionError):
            audit_sample_counters(self.field, self.sample, 20, True)
        self.sample["integration"]["rhsCalls"] = 20
        self.field["integration"]["rhsCalls"] = 19
        with self.assertRaises(AssertionError):
            audit_sample_counters(self.field, self.sample, 20, True)


if __name__ == "__main__":
    unittest.main()
