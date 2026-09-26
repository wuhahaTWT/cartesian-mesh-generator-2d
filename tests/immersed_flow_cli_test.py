#!/usr/bin/env python3
import importlib.util
import json
import os
import signal
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

CLI=Path(sys.argv.pop(1)).resolve()
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('immersed_verify',ROOT/'tools/verification/verify_immersed_flow.py')
verify=importlib.util.module_from_spec(spec)
spec.loader.exec_module(verify)

class ImmersedFlowTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
    def run_case(self,name,*args,expected=0):
        path=self.root/name
        result=subprocess.run([str(CLI),'--output',str(path),*map(str,args)],text=True,capture_output=True,timeout=60)
        self.assertEqual(result.returncode,expected,result.stdout+result.stderr)
        return path,result
    def test_poiseuille_from_rest(self):
        path,_=self.run_case('channel','--nx',16,'--ny',16,'--nu',.1,'--drive',1.2,'--max-steps',10000)
        audit=verify.verify(path,True)
        # One inexpensive development grid, not a mesh-convergence criterion.
        # L2 velocity error relative to the analytic pressure-driven parabola.
        self.assertLess(audit['channel_relative_l2'],.01)
        self.assertEqual(audit['continuity'],0)
    def test_cylinder_deflects_and_slows_flow(self):
        path,_=self.run_case('cylinder','--case','cylinder','--nx',48,'--ny',16,'--nu',.1,'--drive',1.2,'--max-steps',10000)
        audit=verify.verify(path,True)
        result=json.loads((path/'summary.json').read_text())
        self.assertGreater(result['grid']['classification_counts'][1],0)
        self.assertGreater(result['grid']['classification_counts'][2],0)
        self.assertLess(result['metrics']['mean_velocity'],.9)
        self.assertGreater(result['metrics']['penalty_drag_per_density'],0)
        self.assertLess(abs(result['metrics']['penalty_lift_per_density']),1e-5)
        self.assertLess(audit['wall_speed_max'],.2) # coarse-grid smoke only; actual error is reported
        self.assertGreater(max(abs(row['v']) for row in verify.table(path/'cells.csv')),.05)
    def test_budget_is_not_convergence_and_reproducible(self):
        a,_=self.run_case('a','--case','cylinder','--nx',48,'--ny',16,'--max-steps',3,expected=2)
        b,_=self.run_case('b','--case','cylinder','--nx',48,'--ny',16,'--max-steps',3,expected=2)
        audit=verify.verify(a)
        self.assertEqual(audit['stop_reason'],'iteration-limit')
        self.assertEqual(audit['steps'],3)
        for name in ['u.csv','v.csv','cells.csv','walls.csv','history.csv','field.vtk']:
            self.assertEqual((a/name).read_bytes(),(b/name).read_bytes())
    def test_failed_candidate_preserves_initial_state(self):
        path,_=self.run_case('failure','--case','cylinder','--nx',48,'--ny',16,'--max-steps',2,'--continuity-tolerance',1e-30,expected=2)
        report=json.loads((path/'summary.json').read_text())
        self.assertEqual(report['stop_reason'],'candidate-failed')
        self.assertEqual(report['state_status'],'initial-only')
        self.assertEqual(report['steps'],0)
        verify.verify(path)
    def test_invalid_geometry_is_explicit(self):
        for name,text in [('cross','1 .3\n2 .7\n1 .7\n2 .3\n'),('zero','1 .3\n1.5 .3\n2 .3\n'),('duplicate','1 .3\n2 .3\n2 .3\n2 .7\n1 .7\n')]:
            xy=self.root/(name+'.xy');xy.write_text(text)
            path,result=self.run_case(name,'--case','custom','--boundary',xy,'--nx',48,'--ny',16,expected=1)
            self.assertIn('invalid solid geometry',result.stderr)
            self.assertFalse(path.exists())
    def test_custom_polygon_preserves_coordinates(self):
        xy=self.root/'solid.xy';xy.write_text('1.21 .31\n1.63 .29\n1.69 .65\n1.28 .69\n')
        path,_=self.run_case('custom','--case','custom','--boundary',xy,'--nx',64,'--ny',24,'--max-steps',2,expected=2)
        verify.verify(path)
        points=lambda p:[tuple(map(float,line.split())) for line in p.read_text().splitlines() if line.strip()]
        self.assertEqual(points(xy),points(path/'boundary.xy'))
    def test_nested_loops_keep_fluid_hole(self):
        xy=self.root/'hole.xy'
        xy.write_text('1 .25\n2 .25\n2 .75\n1 .75\n\n1.3 .4\n1.7 .4\n1.7 .6\n1.3 .6\n')
        path,_=self.run_case('hole','--case','custom','--boundary',xy,'--nx',64,'--ny',32,'--max-steps',2,expected=2)
        verify.verify(path)
        report=json.loads((path/'summary.json').read_text())
        self.assertAlmostEqual(report['grid']['solid_area_geometry'],.42)
    def test_input_and_output_guards(self):
        for index,args in enumerate([('--nx','-2'),('--nu','nan'),('--penalty-time','0'),('--linear-solver','unknown'),('--case','unknown'),('--max-steps','0'),('--nx','3')]):
            self.run_case(str(index),*args,expected=1)
        path,_=self.run_case('existing','--nx',8,'--ny',8,'--max-steps',1,expected=2)
        before=(path/'summary.json').read_bytes()
        _,result=self.run_case('existing',expected=1)
        self.assertIn('already exists',result.stderr)
        self.assertEqual(before,(path/'summary.json').read_bytes())
    @unittest.skipUnless(sys.platform == 'darwin', 'system Cholesky is a macOS option')
    def test_system_cholesky_matches_portable_pressure_solve(self):
        a,_=self.run_case('ic0','--case','cylinder','--nx',48,'--ny',16,'--max-steps',3,expected=2)
        b,_=self.run_case('cholesky','--case','cylinder','--nx',48,'--ny',16,'--max-steps',3,'--linear-solver','cholesky',expected=2)
        for directory in (a,b):
            audit=verify.verify(directory)
            self.assertEqual(audit['stop_reason'],'iteration-limit')
            self.assertEqual(audit['steps'],3)
        for name,field in [('u.csv','u'),('v.csv','v'),('cells.csv','pressure_fluctuation')]:
            difference=max(abs(x[field]-y[field]) for x,y in zip(verify.table(a/name),verify.table(b/name)))
            self.assertLess(difference,1e-7) # same equations, allowance above requested 1e-8 linear tolerance
    @unittest.skipIf(os.name == 'nt', 'POSIX signal check; Windows termination needs separate validation')
    def test_cancellation_retains_accepted_state(self):
        path=self.root/'cancelled'
        process=subprocess.Popen([str(CLI),'--output',str(path),'--case','cylinder','--nx','128','--ny','32','--max-steps','1000000'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        try:
            self.assertIn('step=1 ',process.stdout.readline())
            process.send_signal(signal.SIGTERM)
            stdout,stderr=process.communicate(timeout=30)
            self.assertEqual(process.returncode,130,stdout+stderr)
            self.assertEqual(verify.verify(path)['stop_reason'],'cancelled')
        finally:
            if process.poll() is None:
                process.kill()
                process.communicate()
    def test_surface_coupling_reduces_dense_wall_error(self):
        xy=self.root/'skew.xy';xy.write_text('1.21 .31\n1.63 .29\n1.69 .65\n1.28 .69\n')
        common=['--case','custom','--boundary',xy,'--nx',32,'--ny',16,'--nu',.1,'--drive',1.2,'--max-steps',10000]
        base,_=self.run_case('base',*common)
        wall,_=self.run_case('wall',*common,'--wall-method','surface-penalty','--wall-penalty-time',1e-6)
        a,b=verify.verify(base,True),verify.verify(wall,True)
        # Development guard: dense original-wall samples must improve by >=4x,
        # not just quadrature markers. Not an engineering wall-accuracy gate.
        self.assertLess(b['wall_speed_max'],a['wall_speed_max']/4)
        self.assertLess(b['wall_normal_speed_max'],a['wall_normal_speed_max']/4)
        # On a quadratic interval the three Gauss interpolation basis functions
        # have sum(abs(L_i)) <= 7/3. An omitted interval/endpoint breaks this
        # bound (the old sparse-marker skew-wall case exceeded it by >50x).
        self.assertLessEqual(b['wall_speed_max'],(7/3)*b['marker_speed_max']+1e-12)
        self.assertLess(b['surface_power_per_density'],0)
        self.assertLess(b['force_balance'],1e-4) # same normalized steady balance budget
        report=json.loads((wall/'summary.json').read_text())
        self.assertEqual(report['controls']['linear_solver'],'jacobi')
        self.assertGreater(report['metrics']['mean_velocity'],.05) # not a stagnant fake solution
        self.assertEqual((base/'boundary.xy').read_bytes(),(wall/'boundary.xy').read_bytes())
    def test_surface_quadrature_ignores_collinear_vertex_subdivision(self):
        # The wall-energy integral should depend on geometry, not how a straight
        # input edge was divided. Reversing the same loop also preserves physics.
        points=[(1.125,.375),(1.625,.375),(1.625,.625),(1.125,.625)]
        subdivided=[]
        for a,b in zip(points,points[1:]+points[:1]):
            subdivided.extend([a,((2*a[0]+b[0])/3,(2*a[1]+b[1])/3),((a[0]+2*b[0])/3,(a[1]+2*b[1])/3)])
        outputs=[]
        for name,vertices in [('plain',points),('split',subdivided),('reversed',list(reversed(points)))]:
            xy=self.root/(name+'.xy');xy.write_text(''.join(f'{x:.17g} {y:.17g}\n' for x,y in vertices))
            path,_=self.run_case(name,'--case','custom','--boundary',xy,'--nx',48,'--ny',16,'--wall-method','surface-penalty','--wall-penalty-time',1e-6,'--max-steps',3,expected=2)
            self.assertEqual(verify.verify(path)['steps'],3)
            outputs.append(path)
        for other in outputs[1:]:
            for name,field in [('u.csv','u'),('v.csv','v'),('cells.csv','pressure_fluctuation')]:
                delta=max(abs(a[field]-b[field]) for a,b in zip(verify.table(outputs[0]/name),verify.table(other/name)))
                self.assertLess(delta,1e-6) # allowance for independent PCG paths at relative 1e-8
    def test_coupled_backend_guards_and_rejected_state(self):
        _,result=self.run_case('ic0','--wall-method','surface-penalty','--linear-solver','ic0',expected=1)
        self.assertIn('requires jacobi or cholesky',result.stderr)
        for name,args in [('method',['--wall-method','unknown']),('time',['--wall-penalty-time',0])]:
            self.run_case(name,*args,expected=1)
        path,_=self.run_case('fail','--case','cylinder','--nx',48,'--ny',16,'--wall-method','surface-penalty','--max-steps',2,'--continuity-tolerance',1e-30,expected=2)
        report=json.loads((path/'summary.json').read_text())
        self.assertEqual(report['stop_reason'],'candidate-failed')
        self.assertEqual(report['steps'],0)
        verify.verify(path)
        self.assertTrue(all(row['wall_force']==0 for row in verify.table(path/'u.csv')))
    @unittest.skipUnless(sys.platform == 'darwin', 'system Cholesky is a macOS option')
    def test_coupled_cholesky_matches_portable_jacobi(self):
        outputs=[]
        for solver in ('jacobi','cholesky'):
            path,_=self.run_case(solver,'--case','cylinder','--nx',48,'--ny',16,'--wall-method','surface-penalty','--wall-penalty-time',1e-6,'--max-steps',3,'--linear-solver',solver,expected=2)
            self.assertEqual(verify.verify(path)['steps'],3)
            outputs.append(path)
        for name,field in [('u.csv','u'),('v.csv','v'),('cells.csv','pressure_fluctuation')]:
            delta=max(abs(a[field]-b[field]) for a,b in zip(verify.table(outputs[0]/name),verify.table(outputs[1]/name)))
            self.assertLess(delta,1e-6)
    def test_surface_mode_without_solids_preserves_channel(self):
        outputs=[]
        for method in ('brinkman','surface-penalty'):
            path,_=self.run_case(method,'--nx',16,'--ny',16,'--wall-method',method,'--max-steps',3,expected=2)
            self.assertEqual(verify.verify(path)['steps'],3);outputs.append(path)
        for name in ('cells.csv','u.csv','v.csv','walls.csv','field.vtk'):
            self.assertEqual((outputs[0]/name).read_bytes(),(outputs[1]/name).read_bytes())
    def test_surface_nested_loop_orientation(self):
        xy=self.root/'nested.xy'
        xy.write_text('1 .25\n2 .25\n2 .75\n1 .75\n\n1.3 .4\n1.3 .6\n1.7 .6\n1.7 .4\n')
        path,_=self.run_case('nested','--case','custom','--boundary',xy,'--nx',64,'--ny',32,'--wall-method','surface-penalty','--max-steps',2,expected=2)
        self.assertEqual(verify.verify(path)['steps'],2)
        self.assertAlmostEqual(json.loads((path/'summary.json').read_text())['grid']['solid_area_geometry'],.42)
    def test_surface_reader_detects_force_corruption(self):
        path,_=self.run_case('force','--case','cylinder','--nx',48,'--ny',16,'--wall-method','surface-penalty','--max-steps',2,expected=2)
        verify.verify(path)
        p=path/'u.csv';original=p.read_text();lines=original.splitlines();values=lines[1].split(',');values[-1]='1';lines[1]=','.join(values);p.write_text('\n'.join(lines)+'\n')
        with self.assertRaisesRegex(verify.VerificationError,'spread wall force mismatch'):
            verify.verify(path)
        p.write_text(original)
        p=path/'walls.csv';lines=p.read_text().splitlines();values=lines[1].split(',')
        values[2]=str(-float(values[2]));values[3]=str(-float(values[3]));lines[1]=','.join(values);p.write_text('\n'.join(lines)+'\n')
        with self.assertRaisesRegex(verify.VerificationError,'wall normal orientation mismatch|summary mismatch: wall_normal_flux_net'):
            verify.verify(path)
    def test_reader_detects_corruption(self):
        path,_=self.run_case('corrupt','--nx',8,'--ny',8,'--max-steps',1,expected=2)
        text=(path/'field.vtk').read_text().replace('DIMENSIONS 9 9 1','DIMENSIONS 8 9 1')
        (path/'field.vtk').write_text(text)
        with self.assertRaisesRegex(verify.VerificationError,'VTK dimensions'):
            verify.verify(path)

if __name__=='__main__':
    unittest.main()
