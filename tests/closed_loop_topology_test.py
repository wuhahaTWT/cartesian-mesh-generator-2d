"""Focused tests of measured feedback and accepted-state semantics."""
from contextlib import redirect_stdout
from dataclasses import dataclass
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools/optimization'))
import closed_loop_topology as loop
from brinkman import Problem
from navier_stokes_brinkman import NavierStokesBrinkman
from optimize_flow import initial_design


class FeedbackTests(unittest.TestCase):
    def test_measured_secants_correct_a_wrong_direction(self):
        x = np.array([.5, .5, 0.])
        samples = [dict(id=1, x=x+[.1, 0, 0], low=.9, high=1.2, topology='A'),
                   dict(id=2, x=x+[0, .1, 0], low=1.1, high=.9, topology='A')]
        correction, info = loop.discrepancy_gradient(x, 1., 1., 'A', samples,
                                                    np.array([True, True, False]))
        np.testing.assert_allclose(correction, [3., -2., 0.], atol=1e-13)
        np.testing.assert_allclose(np.array([-1., 1., 0.])+correction, [2., -1., 0.])
        self.assertEqual(info['rank'], 2)

    def test_failed_and_other_topology_samples_do_not_supply_a_gradient(self):
        x = np.array([.5, .5])
        samples = [dict(id=1, x=x+.1, low=1., high=None, topology='A'),
                   dict(id=2, x=x+.1, low=1., high=999., topology='B')]
        correction, info = loop.discrepancy_gradient(x, 1., 1., 'A', samples,
                                                    np.ones(2, dtype=bool))
        np.testing.assert_array_equal(correction, [0., 0.])
        self.assertEqual(info['rank'], 0)

    def test_stationary_roundoff_probe_does_not_amplify_native_noise(self):
        x = np.array([.5, .5])
        samples = [dict(id=1, x=x+1e-13, low=1., high=1.+1e-8, topology='A')]
        correction, info = loop.discrepancy_gradient(x, 1., 1., 'A', samples,
                                                    np.ones(2, dtype=bool))
        np.testing.assert_array_equal(correction, [0., 0.])
        self.assertEqual(info['rank'], 0)

    def test_topology_probe_preserves_passive_cells_and_volume(self):
        p = Problem(nx=24, ny=16, filter_radius=.08, port_width=.25)
        model = NavierStokesBrinkman(p)
        x = model.feasible_design(initial_design(model, 'geometric'), 6.)
        e = model.evaluate(x, .1, 6., 'total-pressure-power')
        bank = loop.material_probes(model, x, e, 6., .15, 7, 0)
        self.assertIn('open-component-bridge', [name for name, _ in bank])
        for _, candidate in bank:
            np.testing.assert_array_equal(candidate[~model.design], model.fixed_design[~model.design])
            self.assertLessEqual(model.volume(candidate, 6.)[0], p.volume_fraction+1e-10)

    def test_online_feedback_creates_and_accepts_a_new_direction(self):
        @dataclass
        class E:
            objective: float
            gradient: np.ndarray
            volume: float = 0.
            @property
            def volume_gradient(self): return np.zeros(2)
        class Model:
            problem = Problem()
            design = np.ones(2, dtype=bool)
            def enforce_passive(self, x): return x
            def feasible_design(self, x, beta): return x
            def gradient_candidate(self, x, e, beta, move):
                return np.clip(x-move*e.gradient/max(np.max(np.abs(e.gradient)), 1e-30), 0, 1)
            def evaluate(self, x, q, beta, objective):
                return E(2-x[0], np.array([-1., 0.]))
        model = Model()
        class Evaluator:
            def __init__(self, root):
                self.root, self.samples, self.records = root, [], []
            def evaluate(self, x, name, feedback=None, round_index=None):
                i = len(self.records)
                # Low objective decreases in x[0]; high objective increases.
                low, high = 2-x[0], 1+x[0]
                sample = dict(id=i, x=x.copy(), low=low, high=high, topology='A',
                              e=E(low, np.array([-1., 0.])))
                self.samples.append(sample); self.records.append(dict(high=high, feedback=feedback))
                return sample
        def probes(model, x, e, beta, move, seed, iteration):
            return [('porous-gradient', model.gradient_candidate(x, e, beta, move)),
                    ('interface', x.copy())]
        with tempfile.TemporaryDirectory() as tmp, patch.object(loop, 'material_probes', probes), \
             patch.object(loop, 'snapshot'), redirect_stdout(io.StringIO()):
            report = loop.evolve(model, np.array([.5, .5]), .1, 6., Evaluator(Path(tmp)),
                                 rounds=2, move=.1, probes=1)
        self.assertGreater(report['feedbackDrivenAcceptedSteps'], 0)
        self.assertLess(report['finalObjective'], report['initialObjective'])
        self.assertGreater(report['history'][0]['changeFromPorousProposal'], 0)


if __name__ == '__main__':
    unittest.main()
