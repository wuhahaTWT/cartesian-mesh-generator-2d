#!/usr/bin/env python3
"""Focused research tests; run with tools/optimization/requirements.txt installed.

Checks discretisation/adjoint/constraint contracts, not general CFD accuracy.
The Poiseuille check is one analytic sanity case with a deliberately stated
development tolerance, not a grid-convergence or industrial qualification.
"""
from contextlib import redirect_stdout
import io
import json
from dataclasses import asdict
import shlex
import sys
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"tools"/"optimization"))
from brinkman import Problem, StokesBrinkman, port_average
from topology_artifacts import contours, port_connectivity, signed_area
from optimize_flow import parser, run, stationarity
from optimize_flow import initial_design, reference_design
from compare_sharp_designs import assess, match_sharp_area, reuse_baseline
from verify_extracted_flow import parabolic_boundaries
from render_connectivity_study import grouped_sensitivity


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

    def test_real_oc_stall_has_feasible_gradient_descent(self):
        # Preserve the observed failure's design only, not a cached flow field.
        fixture = Path(__file__).with_name("fixtures")/"topology_oc_stall.npz"
        with np.load(fixture) as saved:
            m = StokesBrinkman(Problem(**json.loads(str(saved["problem"]))))
            x = saved["design"].ravel()
            q, beta = float(saved["q"]), float(saved["beta"])
        e = m.evaluate(x, q, beta)
        oc = m.oc_candidate(x, e, beta)
        self.assertGreater(e.gradient@(oc-x), 0)
        self.assertGreater(m.evaluate(oc, q, beta).objective, e.objective)
        candidate = m.gradient_candidate(x, e, beta)
        updated = m.evaluate(candidate, q, beta)
        self.assertLess(e.gradient@(candidate-x), 0)
        self.assertLess(updated.objective, e.objective)
        self.assertLessEqual(updated.volume, m.problem.volume_fraction+1e-10)
        self.assertLessEqual(np.max(np.abs(candidate-x)), .15+1e-12)
        np.testing.assert_array_equal(candidate[~m.design], x[~m.design])

    def test_initialization_records_explicit_connected_seed(self):
        m = StokesBrinkman(Problem(nx=48, ny=24, width=2))
        for name, count in (("geometric", 2), ("merged", 1)):
            x = initial_design(m, name)
            rho = m.physical(x, 0)[0].reshape(m.ny, m.nx)
            groups = contours(rho, 2, 1)
            self.assertEqual(len(groups), count)
            reachability = port_connectivity(groups, asdict(m.problem))
            self.assertTrue(reachability["allPortsCovered"])
            self.assertEqual(reachability["inletToOutletReachability"],
                             [[True, False], [False, True]] if name == "geometric" else [[True, True], [True, True]])
            self.assertLessEqual(m.volume(x, 0)[0], m.problem.volume_fraction+1e-10)
            np.testing.assert_array_equal(x[~m.design], m.fixed_design[~m.design])
        with self.assertRaises(ValueError):
            initial_design(StokesBrinkman(Problem(case="bend")), "merged")

    def test_driver_recovers_from_non_descent_oc_proposal(self):
        with tempfile.TemporaryDirectory() as temporary, redirect_stdout(io.StringIO()):
            args = parser().parse_args(["--output", str(Path(temporary)/"run"),
                                        "--nx", "18", "--ny", "12", "--port-width", ".25",
                                        "--filter-radius", ".12", "--iterations", "1", "--no-plot"])
            def stalled(model, x, evaluation, beta, move):
                return x.copy()
            with patch.object(StokesBrinkman, "oc_candidate", stalled):
                report = run(args)
            self.assertEqual(sum(s["gradientSteps"] for s in report["stages"]), 3)
            for stage in report["stages"]:
                self.assertLess(stage["final"]["objective"], stage["initial"]["objective"])
                self.assertLessEqual(stage["final"]["volumeFraction"], args.volume+1e-10)

    def test_internal_component_cannot_receive_artificial_ports(self):
        with tempfile.TemporaryDirectory() as temporary:
            source, target = Path(temporary)/"template", Path(temporary)/"boundaries"
            source.write_text('CARTMESH2D_FLOW_BOUNDARIES 1\nCOUNTS 1 5 5\n'
                'BOUNDARY 0 0 0.5 0.25 -0.1 0 velocity-inlet "inlet" .02 0 0\n'
                'BOUNDARY 1 0 1.5 0.25 0.1 0 pressure-outlet "outlet" 0 0 0\n'
                'BOUNDARY 2 0 0 0.25 -0.1 0 velocity-inlet "inlet" .02 0 0\n'
                'BOUNDARY 3 0 2 0.75 0.1 0 pressure-outlet "outlet" 0 0 0\n'
                'BOUNDARY 4 0 2 0.95 0.02 0 pressure-outlet "outlet" 0 0 0\n')
            counts = parabolic_boundaries(source, target, asdict(Problem(width=2)), .02)
            self.assertEqual(counts, dict(inletFaces=1, outletFaces=1, closedArtificialOpenings=3))
            rows = {int(row[1]):row for line in target.read_text().splitlines()
                    if line.startswith('BOUNDARY ') for row in [shlex.split(line)]}
            for index in (0, 1, 4):
                self.assertEqual(rows[index][7], "wall")
                self.assertEqual(rows[index][9:12], ["0", "0", "0"])

    def test_cached_baseline_cannot_cross_physical_controls_or_geometry(self):
        # A cache from another experiment must fail before reading any field.
        controls = dict(speed=.02, nu=1., tolerance=1e-8, iterations=4000, smallAlpha=.25)
        design = dict(problem=asdict(Problem(width=2)), extraction=dict(boundarySha256="exact-boundary"))
        old = dict(controls=controls, designs=dict(baseline=design), problem=design["problem"], rows=[])
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory/"summary.json").write_text(json.dumps(old))
            for key, value in (("speed", .01), ("nu", .5), ("tolerance", 1e-6), ("smallAlpha", .1)):
                with self.subTest(key=key), self.assertRaisesRegex(ValueError, "different "+key):
                    reuse_baseline(directory, design, dict(controls, **{key:value}), 6, None)
            other = dict(design, extraction=dict(boundarySha256="other-boundary"))
            with self.assertRaisesRegex(ValueError, "different boundary"):
                reuse_baseline(directory, other, controls, 6, None)
            other = dict(design, problem=dict(design["problem"], port_width=.2))
            with self.assertRaisesRegex(ValueError, "different port parameter"):
                reuse_baseline(directory, other, controls, 6, None)
            with self.assertRaisesRegex(ValueError, "completed baseline"):
                reuse_baseline(directory, design, controls, 6, None)

    def test_sensitivity_keeps_rejections_and_separates_mesh_controls(self):
        def report(level, alpha, candidate):
            return dict(controls=dict(levels=[level], smallAlpha=alpha), problem=asdict(Problem(width=2)),
                        source=dict(path=f"level-{level}-{alpha}"),
                        designs={key:dict(extraction=dict(boundarySha256=key)) for key in ("baseline", "candidate")},
                        rows=[dict(level=level, baseline=dict(inletFlux=1, fluxWeightedPressureDrop=10),
                                   candidate=dict(inletFlux=1, fluxWeightedPressureDrop=candidate) if candidate else None)])
        data = grouped_sensitivity([report(5, .25, None), report(6, .25, 6), report(7, .25, 5.9), report(7, .4, 5.8)])
        self.assertEqual(len(data), 2)
        first = data[0]
        self.assertFalse(first["allAttemptedGrids"]["meshRobustImprovementObserved"])
        self.assertEqual(first["allAttemptedGrids"]["unpairedLevels"], [5])
        self.assertEqual(first["acceptedPairsOnly"]["levels"], [6, 7])
        self.assertTrue(first["acceptedPairsOnly"]["assessment"]["meshRobustImprovementObserved"])
        self.assertFalse(data[1]["acceptedPairsOnly"]["assessment"]["meshRobustImprovementObserved"])
        with self.assertRaisesRegex(ValueError, "duplicate grid"):
            grouped_sensitivity([report(6, .25, 6), report(6, .25, 5)])


if __name__ == "__main__":
    unittest.main(verbosity=2)
