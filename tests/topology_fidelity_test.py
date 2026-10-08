#!/usr/bin/env python3
"""Matched controls, true ports, GCI and non-cherry-picked ranking contracts."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tools/optimization"))
from brinkman import Problem
from navier_stokes_brinkman import NavierStokesBrinkman
from fidelity_metrics import leakage,grid_uncertainty,ranking,paired_grid_orderings
from native_flow import matched_controls,prescribed_boundaries,neutral_boundary_template,native
from optimize_flow import initial_design
from dataclasses import asdict
from compare_sharp_designs import match_sharp_area,area
from topology_artifacts import contours
import native_flow as bridge


class FidelityTest(unittest.TestCase):
    def test_area_search_rejects_plateau_intermediate_without_losing_components(self):
        fixture=Path(__file__).parent/'fixtures/topology_area_plateau.npz'
        with np.load(fixture) as f:
            specs=json.loads(str(f['specifications']))
            for case,p in specs.items():
                with self.subTest(case=case):
                    model=NavierStokesBrinkman(Problem(**p))
                    rho=f[case].copy()
                    original=rho.copy()
                    plateau=rho.copy()
                    mask=model.design.reshape(rho.shape)
                    plateau[mask]=np.clip(plateau[mask]+.5,0,1)
                    with self.assertRaisesRegex(ValueError,'zero-area contour'):
                        contours(plateau,p['width'],p['height'])
                    target=p['width']*p['height']*p['volume_fraction']
                    field,record=match_sharp_area(rho,model,target)
                    self.assertTrue(record['rejectedSearchSamples'])
                    np.testing.assert_array_equal(rho,original)
                    np.testing.assert_array_equal(field.ravel()[~model.design],model.fixed_design[~model.design])
                    groups=contours(field,p['width'],p['height'])
                    self.assertLessEqual(abs(area(groups)-target),record['areaTolerance'])
                    self.assertGreaterEqual(len(groups),1)
    def test_seeded_random_design_is_reproducible_and_feasible(self):
        m=NavierStokesBrinkman(Problem(nx=18,ny=12,port_width=.25,filter_radius=.12))
        a=initial_design(m,"random",41)
        np.testing.assert_array_equal(a,initial_design(m,"random",41))
        self.assertGreater(np.max(abs(a-initial_design(m,"random",42))),.01)
        self.assertAlmostEqual(m.volume(a,0)[0],m.problem.volume_fraction,places=10)
        np.testing.assert_array_equal(a[~m.design],m.fixed_design[~m.design])

    def test_mean_re_mapping_and_exact_stokes(self):
        p=asdict(Problem())
        for re in (0,10,50,100):
            c=matched_controls(p,re,.02)
            self.assertAlmostEqual(c['momentumInertia']*(2/3)*c['speed']*p['port_width']/c['nu'],re,places=12)
        self.assertEqual(matched_controls(p,0,.02)['momentumInertia'],0)
        with self.assertRaises(ValueError):matched_controls(p,-1,.02)

    def test_outlet_on_bottom_and_no_false_internal_ports(self):
        with tempfile.TemporaryDirectory() as d:
            src,dst=Path(d)/'in',Path(d)/'out'
            src.write_text('CARTMESH2D_FLOW_BOUNDARIES 1\nCOUNTS 1 3 3\n'
                'BOUNDARY 0 0 0 0.8 -0.2 0 wall "wall" 0 0 0\n'
                'BOUNDARY 1 0 0.8 0 0 -0.2 wall "wall" 0 0 0\n'
                'BOUNDARY 2 0 0.4 0.8 -0.2 0 velocity-inlet "inlet" .02 0 0\n')
            c=prescribed_boundaries(src,dst,asdict(Problem(case='elbow',width=1,port_width=.2)),.02)
            self.assertEqual(c,dict(inletFaces=1,outletFaces=1,closedArtificialOpenings=1))
            text=dst.read_text()
            self.assertIn('velocity-outlet "outlet_0"',text)
            self.assertIn('BOUNDARY 2 0 0.4 0.8 -0.2 0 wall',text)

    def test_neutral_template_does_not_require_duct_openings(self):
        # Minimal two-port corner fixture: its right side is a wall, not an
        # outlet. Building a template must not instantiate a duct solver.
        vertices=((0.,0.),(1.,0.),(1.,1.),(0.,1.))
        edges=tuple(native.Edge(i,i,(i+1)%4,0,-1,0) for i in range(4))
        mesh=native.Mesh(Path('corner.cm2d'),vertices,edges,
                         (native.Cell(0,1.,(0,1,2,3),(0,1,2,3)),),(0,)*8)
        measured=native.measure(mesh)
        with tempfile.TemporaryDirectory() as directory:
            template,target=Path(directory)/'template',Path(directory)/'assigned'
            neutral_boundary_template(mesh,measured,template)
            result=prescribed_boundaries(template,target,asdict(Problem(case='elbow',width=1,port_width=.2)),.02)
            self.assertEqual(result['inletFaces'],1)
            self.assertEqual(result['outletFaces'],1)
            self.assertIn('BOUNDARY 1 0 1 0.5 1 0 wall',target.read_text())
            self.assertTrue(target.read_text().endswith('END\n'))

    def row(self,level,value,cells):
        return dict(level=level,metrics=dict(dimensionlessTotalPower=value,cells=cells,area=1.))

    def test_gci_reference_and_two_grid_boundary(self):
        rows=[self.row(4,.96178,1),self.row(5,.96854,4),self.row(6,.97050,16)]
        result=grid_uncertainty(rows)
        self.assertEqual(result['method'],'three-grid-GCI')
        # Rounded published input table, not the extra digits in its Fortran
        # sample output: independent equal-ratio analytic expression.
        self.assertAlmostEqual(result['observedOrder'],np.log((.96854-.96178)/(.97050-.96854))/np.log(2),places=10)
        self.assertAlmostEqual(result['gciRelative'],.0010308,places=7)
        two=grid_uncertainty(rows[1:]);self.assertIsNone(two['uncertainty'])
        self.assertEqual(two['method'],'two-grid-change')
        oscillatory=grid_uncertainty([self.row(4,1,1),self.row(5,2,4),self.row(6,1.5,16)])
        self.assertIsNone(oscillatory['uncertainty'])
        missing=grid_uncertainty(rows+[dict(level=7,metrics=None,status='failed')])
        self.assertFalse(missing['allAttemptedGridsAccepted'])
        rows[-1]['meshPaddingFraction']=1/37
        mixed=grid_uncertainty(rows)
        self.assertIsNone(mixed['uncertainty'])
        self.assertIn('padding families',mixed['gciUnavailableReason'])

    def test_ranking_does_not_promote_two_grid_change_to_uncertainty(self):
        a=dict(id='a',porousObjective=1.,grid=dict(J_T=3.,uncertainty=.1,twoGridChange=.2))
        b=dict(id='b',porousObjective=2.,grid=dict(J_T=2.,uncertainty=.1,twoGridChange=.2))
        failed=dict(id='failed',porousObjective=0,grid={})
        r=ranking([a,b,failed]);self.assertEqual(r['kendallTau'],-1)
        self.assertEqual(r['attempted'],3);self.assertEqual(len(r['conditionalGciReversals']),1)
        self.assertEqual(len(r['resolvedReversals']),0)
        for row in (a,b):
            row['porousEvaluationContract']={'alphaMax':100.,'Re':50}
            row['grid']['asymptoticRangeIndependentlyEstablished']=True
        self.assertEqual(len(ranking([a,b])['resolvedReversals']),1)
        b['porousEvaluationContract']['alphaMax']=1000.
        self.assertEqual(len(ranking([a,b])['resolvedReversals']),0)
        self.assertEqual(len(ranking([a,b])['conditionalGciReversals']),1)
        b['grid']['uncertainty']=None
        r=ranking([a,b]);self.assertEqual(len(r['resolvedReversals']),0)
        self.assertEqual(len(r['twoGridSensitivityReversals']),1)

    def test_rank_comparisons_use_shared_levels_and_keep_order_changes(self):
        a=dict(porousObjective=1.,native=[self.row(3,4,20),self.row(4,2,80),self.row(7,1,1000)])
        b=dict(porousObjective=2.,native=[self.row(3,3,20),self.row(4,3,80),self.row(6,5,500)])
        r=paired_grid_orderings(a,b)
        self.assertEqual(r['commonLevels'],[3,4])
        self.assertTrue(r['nativeOrderChangesAcrossCommonGrids'])
        self.assertFalse(r['reversalObservedOnAtLeastTwoCommonGrids'])
        # Individual finest levels have opposite ordering to the shared fine
        # level. They must not replace the actual paired comparison.
        self.assertEqual(r['comparisons'][-1]['difference'],-1.)
        b['native'][1]['meshPaddingFraction']=1/37
        self.assertEqual(paired_grid_orderings(a,b)['commonLevels'],[3])

    def test_leakage_energy_partition(self):
        m=NavierStokesBrinkman(Problem(nx=18,ny=12,port_width=.25,filter_radius=.12,reynolds=50))
        e=m.evaluate(m.uniform_design(3),beta=3)
        values=leakage(m,e,.1)
        self.assertLess(values['dissipationReconstructionRelativeError'],1e-12)
        for name in ('solidFluxFraction','solidDissipationFraction'):
            self.assertTrue(0<=values[name]<=1)


if __name__=='__main__':unittest.main(verbosity=2)
