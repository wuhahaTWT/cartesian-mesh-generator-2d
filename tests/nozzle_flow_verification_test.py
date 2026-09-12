#!/usr/bin/env python3
"""Check CFD evidence parsing and timeout isolation without running a solver."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch, Mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools/verification'))
import openfoam_nozzle_flow as nozzle
import run_nozzle_flow as driver


class NozzleEvidenceTests(unittest.TestCase):
    def test_common_profile_integral(self):
        def linear(ys, offset=0):
            return [{'y': y, 'normalized_ux': 2*y+offset} for y in ys]
        profiles = [linear([-1,0,1]), linear([-.9,-.2,.3,.9]), linear([-.8,0,.8])]
        report = nozzle.compare_profiles(profiles)
        self.assertEqual(report['common_y_interval'], [-.8,.8])
        self.assertEqual(len(report['sample_y']), 64)
        self.assertTrue(all(value < 1e-14 for value in report['adjacent_l2']))
        shifted = nozzle.compare_profiles([linear([-1,0,1]),linear([-1,0,1],1)])
        self.assertAlmostEqual(shifted['adjacent_l2'][0],1)
        with self.assertRaises(ValueError):
            nozzle.compare_profiles([linear([0,0,1]),linear([0,1])])

    def test_last_iteration_is_authoritative(self):
        block = ('Time = 279\nSolving for Ux, Initial residual = 1e-9\n'
                 'Solving for Uy, Initial residual = 1e-9\n'
                 'Solving for p, Initial residual = 1e-9\n')
        self.assertTrue(nozzle.solver_convergence(block+'SIMPLE solution converged')['converged_at_requested_residual'])
        bad = block.replace('for p, Initial residual = 1e-9','for p, Initial residual = 1e-4')
        self.assertFalse(nozzle.solver_convergence(bad+'SIMPLE solution converged')['converged_at_requested_residual'])
        self.assertFalse(nozzle.solver_convergence(block+'SIMPLE solution converged\n'+block)['converged_at_requested_residual'])

    def test_field_cardinality(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'U'
            p.write_text(nozzle.header('volVectorField','0','U')+'internalField nonuniform List<vector> 2 ((1 0 0)(2 0 0));\n')
            self.assertEqual(nozzle.vector_field(p),[(1,0,0),(2,0,0)])
            p.write_text(p.read_text().replace('List<vector> 2','List<vector> 3'))
            with self.assertRaises(ValueError): nozzle.vector_field(p)
            p.write_text(nozzle.header('volScalarField','0','p')+'internalField uniform nan;\n')
            with self.assertRaises(ValueError): nozzle.scalar_field(p)

    def test_timeout_cleans_only_named_container(self):
        with tempfile.TemporaryDirectory() as folder:
            process = Mock(pid=123456)
            process.wait.side_effect=[subprocess.TimeoutExpired(['fixture'],1),0]
            with patch.object(driver.subprocess,'Popen',return_value=process), \
                 patch.object(driver.subprocess,'run',return_value=SimpleNamespace(returncode=0,stdout='',stderr='')) as cleanup, \
                 patch.object(driver.os,'killpg',create=True) as kill:
                result=driver.run_stage(Path(folder)/'stage',['fixture'],1,'owned-fixture')
                self.assertEqual(result['status'],'timeout')
                self.assertEqual(cleanup.call_args.args[0],['docker','rm','-f','owned-fixture'])
                if driver.os.name == 'posix':
                    kill.assert_called_once()
                    process.kill.assert_not_called()
                else:
                    kill.assert_not_called()
                    process.kill.assert_called_once()
                self.assertEqual(json.loads((Path(folder)/'stage/stage.json').read_text())['returncode'],124)

    def test_existing_run_is_preserved(self):
        with tempfile.TemporaryDirectory() as folder:
            sentinel=Path(folder)/'evidence';sentinel.write_text('keep')
            result=subprocess.run([sys.executable,str(Path(driver.__file__)),
                '--output-root',folder,'--cli','/usr/bin/true'],capture_output=True)
            self.assertNotEqual(result.returncode,0)
            self.assertEqual(sentinel.read_text(),'keep')
            self.assertEqual(len(list(Path(folder).iterdir())),1)


if __name__=='__main__':
    unittest.main()
