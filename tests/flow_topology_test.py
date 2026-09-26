#!/usr/bin/env python3
"""Focused research tests; run with tools/optimization/requirements.txt installed.

Checks discretisation/adjoint/constraint contracts, not general CFD accuracy.
The Poiseuille check is one analytic sanity case with a deliberately stated
development tolerance, not a grid-convergence or industrial qualification.
"""
from contextlib import redirect_stdout
import io
import json
import sys
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"tools"/"optimization"))
from brinkman import Problem, StokesBrinkman, port_average
from topology_artifacts import contours, signed_area
from optimize_flow import parser, run, stationarity
from optimize_flow import reference_design
from compare_sharp_designs import assess, match_sharp_area


class FlowTopologyTest(unittest.TestCase):
    def model(self):
        return StokesBrinkman(Problem(nx=18, ny=12, port_width=.25, filter_radius=.12))

    def test_port_flux_is_integrated_and_balanced(self):
        width = .19
        flux = sum(port_average(j/17, (j+1)/17, .27, width)/17 for j in range(17))
        self.assertAlmostEqual(flux, 2*width/3, places=13)
        m = self.model()
        e = m.evaluate(m.uniform_design())
        self.assertLess(e.linear_residual, m.problem.linear_residual_limit)
        self.assertLess(e.continuity, m.problem.linear_residual_limit)

    def test_poiseuille_analytic_sanity(self):
        m = StokesBrinkman(Problem(nx=32, ny=24, case="channel", alpha_max=0))
        e = m.evaluate(np.ones(m.cells))
        u, v = m.cell_velocity(e.velocity)
        ys = (np.arange(m.ny)+.5)*m.dy
        exact = 1-4*(ys-.5)**2
        computed = u.reshape(m.ny, m.nx)[:, m.nx//2]
        velocity_error = np.linalg.norm(computed-exact)/np.linalg.norm(exact)
        pressure = e.pressure.reshape(m.ny, m.nx)
        drop = np.mean(pressure[:, 0])-np.mean(pressure[:, -1])
        exact_drop = 8*m.problem.viscosity*(m.problem.width-m.dx)
        # 2% is a relative analytic-error sanity bound for THIS modest grid.
        # It catches sign, spacing and pressure-scaling mistakes cheaply.
        self.assertLess(velocity_error, .02)
        self.assertLess(abs(drop/exact_drop-1), .02)
        self.assertLess(np.max(np.abs(v)), .02)

    def test_discrete_adjoint_both_objectives_and_projection(self):
        m = self.model()
        rng = np.random.default_rng(7304)
        x = m.enforce_passive(.35+.15*rng.random(m.cells))
        direction = rng.normal(size=m.cells)
        direction[~m.design] = 0
        direction /= np.max(np.abs(direction))
        h = 1e-4
        for objective in ("dissipation", "pressure-power"):
            for beta in (0, 3):
                with self.subTest(objective=objective, beta=beta):
                    e = m.evaluate(x, .1, beta, objective)
                    a = m.evaluate(x+h*direction, .1, beta, objective).objective
                    b = m.evaluate(x-h*direction, .1, beta, objective).objective
                    difference = (a-b)/(2*h)
                    adjoint = e.gradient@direction
                    # Relative directional-gradient agreement; NOT flow accuracy.
                    self.assertLess(abs(difference-adjoint)/max(abs(difference), abs(adjoint), 1e-12), 1e-4)
                    self.assertTrue(np.all(e.gradient[~m.design] == 0))

    def test_volume_derivative_and_continuation_feasibility(self):
        m = self.model()
        rng = np.random.default_rng(12)
        x = m.uniform_design()
        direction = rng.normal(size=m.cells)
        direction[~m.design] = 0
        direction /= np.max(np.abs(direction))
        for beta in (0, 2, 6):
            x = m.feasible_design(x, beta)
            volume, derivative = m.volume(x, beta)
            h = 1e-5
            fd = (m.volume(x+h*direction, beta)[0]-m.volume(x-h*direction, beta)[0])/(2*h)
            self.assertLess(abs(fd-derivative@direction), 1e-8)
            self.assertLessEqual(volume, m.problem.volume_fraction+1e-10)
            self.assertAlmostEqual(volume, m.problem.volume_fraction, places=10)

    def test_constrained_updates_decrease_real_objective(self):
        m = self.model()
        x = m.uniform_design(2)
        first = m.evaluate(x, .1, 2)
        current = first
        for _ in range(5):
            proposal = m.oc_candidate(x, current, 2)
            candidate = m.evaluate(proposal, .1, 2)
            self.assertLessEqual(candidate.objective, current.objective)
            self.assertLessEqual(candidate.volume, m.problem.volume_fraction+1e-10)
            self.assertTrue(np.array_equal(proposal[~m.design], m.fixed_design[~m.design]))
            x, current = proposal, candidate
        self.assertLess(current.objective, first.objective)
        self.assertTrue(0 <= stationarity(m, x, current) <= 1)

    def test_extraction_retains_components_and_holes(self):
        field = np.zeros((12, 16))
        field[1:10, 1:10] = 1
        field[4:7, 4:7] = 0
        field[2:5, 12:15] = 1
        groups = contours(field, 2, 1.5)
        self.assertEqual(len(groups), 2)
        self.assertEqual(sum(len(group)-1 for group in groups), 1)
        self.assertTrue(all(abs(signed_area(loop)) > 0 for group in groups for loop in group))

    def test_invalid_values_and_infeasible_constraints_fail(self):
        for kwargs in (dict(nx=0), dict(ny=3), dict(alpha_max=float("nan")),
                       dict(volume_fraction=0), dict(filter_radius=-1), dict(case="euler")):
            with self.assertRaises(ValueError):
                Problem(**kwargs)
        m = self.model()
        bad = np.full(m.cells, .5); bad[10] = np.nan
        with self.assertRaises(ValueError):
            m.evaluate(bad)
        with self.assertRaises(ValueError):
            m.evaluate(m.uniform_design(), q=0)
        with self.assertRaises(ValueError):
            m.evaluate(m.uniform_design(), beta=float("nan"))
        with self.assertRaises(ValueError):
            StokesBrinkman(Problem(nx=18, ny=12, port_width=.25, volume_fraction=.001)).uniform_design()

    def test_deterministic_analysis_and_update(self):
        m = self.model()
        x = m.uniform_design(2)
        a, b = m.evaluate(x, .1, 2), m.evaluate(x, .1, 2)
        self.assertTrue(np.array_equal(a.velocity, b.velocity))
        self.assertTrue(np.array_equal(a.gradient, b.gradient))
        self.assertTrue(np.array_equal(m.oc_candidate(x, a, 2), m.oc_candidate(x, b, 2)))

    def test_iteration_budget_does_not_claim_convergence(self):
        with tempfile.TemporaryDirectory() as temporary, redirect_stdout(io.StringIO()):
            args = parser().parse_args(["--output", str(Path(temporary)/"run"),
                                        "--nx", "18", "--ny", "12", "--port-width", ".25",
                                        "--filter-radius", ".12", "--iterations", "1", "--no-plot"])
            report = run(args)
            self.assertEqual(report["status"], "iteration-limit")
            self.assertFalse(report["optimizationConverged"])
            self.assertEqual(sum(stage["accepted"] for stage in report["stages"]), 3)
            self.assertEqual(report["finalParameters"], {"q": .1, "beta": 6})
            with np.load(args.output/"accepted-design.npz") as accepted, np.load(args.output/"final.npz") as final:
                for name in ("design", "rho", "u", "v", "p", "q", "beta", "objective"):
                    np.testing.assert_array_equal(accepted[name], final[name])
            with self.assertRaises(ValueError):
                run(args)  # A second invocation must not overwrite the evidence.

    def test_failed_candidate_or_stage_preserves_last_accepted_state(self):
        original = StokesBrinkman.evaluate
        # Call three is either a trial in stage zero or the next stage's first
        # analysis. Neither failure may publish an unevaluated design or beta.
        for budget in (2, 1):
            with self.subTest(budget=budget), tempfile.TemporaryDirectory() as temporary:
                count = 0
                saved = {}
                def fail_third(model, x, q=.1, beta=0, objective="dissipation"):
                    nonlocal count
                    count += 1
                    if count == 3:
                        raise ArithmeticError("injected failed analysis")
                    result = original(model, x, q, beta, objective)
                    saved.update(x=x.copy(), rho=result.rho.copy(), objective=result.objective, beta=beta)
                    return result
                args = parser().parse_args(["--output", str(Path(temporary)/"run"),
                                            "--nx", "18", "--ny", "12", "--port-width", ".25",
                                            "--filter-radius", ".12", "--iterations", str(budget), "--no-plot"])
                with patch.object(StokesBrinkman, "evaluate", fail_third), redirect_stdout(io.StringIO()):
                    report = run(args)
                self.assertEqual(report["status"], "analysis-failed")
                self.assertFalse(report["optimizationConverged"])
                self.assertEqual(report["finalParameters"]["beta"], saved["beta"])
                with np.load(args.output/"final.npz") as final:
                    np.testing.assert_array_equal(final["design"].ravel(), saved["x"])
                    np.testing.assert_array_equal(final["rho"].ravel(), saved["rho"])
                    self.assertEqual(float(final["objective"]), saved["objective"])
                self.assertEqual(json.loads((args.output/"summary.json").read_text())["status"], "analysis-failed")

    def test_sharp_area_matching_preserves_passive_ports_and_source(self):
        m = self.model()
        rho = m.physical(reference_design(m, 6), 6)[0].reshape(m.ny, m.nx)
        original = rho.copy()
        target = m.problem.volume_fraction*m.problem.width*m.problem.height
        matched, report = match_sharp_area(rho, m, target)
        self.assertLessEqual(abs(report["matchedArea"]-target), report["areaTolerance"])
        np.testing.assert_array_equal(rho, original)
        passive = ~m.design.reshape(m.ny, m.nx)
        np.testing.assert_array_equal(matched[passive], rho[passive])
        with self.assertRaises(ValueError):
            match_sharp_area(rho, m, 100)

    def test_grid_change_cannot_be_hidden_by_one_improved_result(self):
        def row(level, baseline, candidate):
            return dict(level=level, baseline=dict(inletFlux=1, fluxWeightedPressureDrop=baseline),
                        candidate=dict(inletFlux=1, fluxWeightedPressureDrop=candidate))
        self.assertFalse(assess([row(6, 10, 8)])["meshRobustImprovementObserved"])
        unstable = assess([row(6, 10, 8), row(7, 11, 10.5)])
        self.assertFalse(unstable["meshRobustImprovementObserved"])
        stable = assess([row(6, 10, 8), row(7, 9.9, 7.9)])
        self.assertTrue(stable["meshRobustImprovementObserved"])
        self.assertFalse(stable["physicalAccuracyQualified"])
        incomplete = assess([row(6, 10, 8), row(7, 9.9, 7.9), dict(level=8, baseline=None)])
        self.assertFalse(incomplete["meshRobustImprovementObserved"])
        self.assertEqual(incomplete["unpairedLevels"], [8])
        unequal = row(6, 10, 8)
        unequal["candidate"]["inletFlux"] = .5
        with self.assertRaises(ValueError):
            assess([unequal])


if __name__ == "__main__":
    unittest.main(verbosity=2)
