#!/usr/bin/env python3
"""Focused nonlinear operator/adjoint contracts; no physical qualification."""
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"tools/optimization"))
from brinkman import Problem, StokesBrinkman
from navier_stokes_brinkman import NavierStokesBrinkman
from topology_cases import CASES,problem_parameters
from engineering_baselines import family
from optimize_flow import reference_design


class NavierStokesTopologyTest(unittest.TestCase):
    def model(self, reynolds):
        return NavierStokesBrinkman(Problem(nx=18, ny=12, port_width=.25,
                                           filter_radius=.12, reynolds=reynolds))

    def test_re_zero_is_exact_legacy_path(self):
        m = self.model(0)
        x = m.uniform_design(3)
        a, b = StokesBrinkman(m.problem).evaluate(x, beta=3), m.evaluate(x, beta=3)
        for name in ("velocity", "pressure", "gradient", "objective"):
            np.testing.assert_array_equal(getattr(a,name), getattr(b,name))

    def test_quadratic_convection_jacobian(self):
        m = self.model(50)
        rng = np.random.default_rng(714)
        u, d = rng.normal(size=(2,m.velocities))
        value, jac = m.convection(u)
        np.testing.assert_allclose(jac@u, 2*value, atol=1e-13)
        h = 1e-4
        fd = (m.convection(u+h*d, False)-m.convection(u-h*d, False))/(2*h)
        np.testing.assert_allclose(jac@d, fd, atol=1e-10, rtol=1e-9)

    def test_nonlinear_adjoint_and_continuity(self):
        rng = np.random.default_rng(316)
        for reynolds in (10,50,100,200):
            m = self.model(reynolds)
            x = m.enforce_passive(.35+.2*rng.random(m.cells))
            d = rng.normal(size=m.cells); d[~m.design]=0; d/=np.max(abs(d))
            for objective in ("dissipation", "pressure-power"):
                with self.subTest(reynolds=reynolds,objective=objective):
                    e = m.evaluate(x, beta=3, objective=objective)
                    h=1e-4
                    fd=(m.evaluate(x+h*d,beta=3,objective=objective).objective-
                        m.evaluate(x-h*d,beta=3,objective=objective).objective)/(2*h)
                    self.assertLess(abs(fd-e.gradient@d)/max(abs(fd),abs(e.gradient@d),1e-12),1e-5)
                    self.assertLessEqual(e.linear_residual,m.problem.newton_residual_limit)
                    self.assertLessEqual(e.continuity,m.problem.linear_residual_limit)

    def test_channel_has_zero_inertia_effect(self):
        # Fully developed u(y), v=0: convective acceleration vanishes.
        for reynolds in (0,100):
            m=NavierStokesBrinkman(Problem(nx=24,ny=16,case="channel",alpha_max=0,reynolds=reynolds))
            e=m.evaluate(np.ones(m.cells))
            self.assertLess(e.linear_residual,m.problem.linear_residual_limit)
            self.assertTrue(np.isfinite(e.velocity).all())

    def test_invalid_controls_rejected(self):
        for controls in (dict(reynolds=-1),dict(reynolds=float("nan")),dict(newton_iterations=0),dict(newton_residual_limit=0)):
            with self.assertRaises(ValueError):
                Problem(**controls)

    def test_benchmark_ports_volume_and_total_power_adjoint(self):
        rng=np.random.default_rng(901)
        for case in CASES:
            with self.subTest(case=case):
                m=NavierStokesBrinkman(Problem(case=case,**problem_parameters(case),reynolds=50))
                self.assertAlmostEqual(float(np.sum(m.pressure_weights)),0,places=13)
                for label,x in family(m,3):
                    self.assertAlmostEqual(m.volume(x,3)[0],m.problem.volume_fraction,places=10)
                x=m.enforce_passive(.35+.2*rng.random(m.cells))
                d=rng.normal(size=m.cells);d[~m.design]=0;d/=np.max(abs(d))
                e=m.evaluate(x,beta=3,objective="total-pressure-power")
                h=1e-4
                fd=(m.evaluate(x+h*d,beta=3,objective="total-pressure-power").objective-
                    m.evaluate(x-h*d,beta=3,objective="total-pressure-power").objective)/(2*h)
                self.assertLess(abs(fd-e.gradient@d)/max(abs(fd),abs(e.gradient@d),1e-12),1e-5)

    def test_high_drag_pressure_adjoint_residual_correction(self):
        # T01-B real failure: direct LU transpose residual was 4.56e-8 for
        # this geometric reference. Retain the same 1e-8 gate, correct LU.
        m=NavierStokesBrinkman(Problem(nx=36,ny=24,width=1.5,alpha_max=360000,filter_radius=.055))
        e=m.evaluate(reference_design(m,6),.1,6,"total-pressure-power")
        self.assertLessEqual(e.adjoint_residual,m.problem.linear_residual_limit)


if __name__ == "__main__":
    unittest.main(verbosity=2)
