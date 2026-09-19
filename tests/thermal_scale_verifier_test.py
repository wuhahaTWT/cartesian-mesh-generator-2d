#!/usr/bin/env python3
"""Bounded contract tests for the thermal scale verifier."""
import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
VERIFY_PATH = ROOT / "tools" / "verification" / "verify_thermal_scale.py"
sys.path.insert(0, str(ROOT / "tools" / "verification"))
SPEC = importlib.util.spec_from_file_location("verify_thermal_scale", VERIFY_PATH)
VERIFY = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(VERIFY)


class ThermalScaleMathTests(unittest.TestCase):
    def test_process_group_timeout_kills_spawned_descendant(self):
        with tempfile.TemporaryDirectory(prefix="cartmesh-thermal-timeout-") as name:
            root=Path(name); marker=root/'delayed-marker'; log=root/'run'
            child_path=root/'child.py'
            child_path.write_text("import pathlib,time; print('CHILD_STARTED',flush=True); time.sleep(.8); "
                                 "pathlib.Path("+repr(str(marker))+").write_text('late')")
            parent=("import subprocess,sys,time; "
                    "print('PRETIMEOUT',flush=True); "
                    "p=subprocess.Popen([sys.executable,"+repr(str(child_path))+"]); "
                    "time.sleep(2.0); p.wait()")
            result=VERIFY.run_process_group([sys.executable,'-c',parent],log,.5)
            self.assertTrue(result['timedOut'])
            output=Path(result['stdout']).read_text()
            self.assertEqual(output.count('PRETIMEOUT'),1)
            self.assertEqual(output.count('CHILD_STARTED'),1)
            time.sleep(.85)
            self.assertFalse(marker.exists())

    def test_backward_euler_amplitude_and_continuous_limit(self):
        rate, dt, steps = 2.0, 0.01, 7
        expected = (1.0 + rate * dt) ** (-steps)
        self.assertAlmostEqual(VERIFY.amplitude(rate, dt, steps), expected, places=14)
        self.assertLess(
            abs(VERIFY.amplitude(rate, 1.0e-5, 100000) - math.exp(-rate)),
            1.0e-5,
        )

    def test_amplitude_rejects_invalid_controls(self):
        for rate, dt, steps in ((0.0, .1, 1), (-1., .1, 1), (1., 0., 1),
                                (1., .1, 0), (1., .1, True), (math.nan, .1, 1)):
            with self.assertRaises(ValueError):
                VERIFY.amplitude(rate, dt, steps)

    def test_series_checks_rejects_mixed_controls_bad_grids_and_nondecreasing_errors(self):
        def case(h, scalar_error, velocity_error, **overrides):
            reference = dict(h=h, dt=.01, steps=2, time=.02, nu=.1,
                             diffusivity=.02, speed=1., velocityRelaxation=.6,
                             coupledFlowTolerance=1e-8, flowTolerance=1e-8,
                             coupledFlowToleranceStatus='explicit-thermal-flow-json',
                             flowToleranceStatus='explicit-flow-json',
                             scalarRelativeTolerance=1e-9, scalarAbsoluteTolerance=1e-12, scalarCellTolerance=1e-9,
                             scalarBackwardEulerL2=scalar_error,
                             velocityBackwardEulerL2=velocity_error)
            reference.update(overrides)
            return {"valid": True, "reference": reference}

        good = VERIFY.series_checks([case(.2, .3, .2), case(.1, .2, .1),
                                     case(.05, .1, .05)])
        self.assertTrue(good["valid"])
        for bad in (
            [case(.2, .3, .2), case(.1, .2, .1, dt=.02), case(.05, .1, .05)],
            [case(.2, .3, .2), case(.2, .2, .1), case(.05, .1, .05)],
        ):
            with self.assertRaises(ValueError):
                VERIFY.series_checks(bad)
        for key,value in [('steps',3),('time',.03),('nu',.2),('diffusivity',.04),('speed',2.),('flowTolerance',2e-8),('coupledFlowTolerance',2e-9),('scalarRelativeTolerance',2e-9)]:
            with self.assertRaises(ValueError):
                VERIFY.series_checks([case(.2,.3,.2),case(.1,.2,.1,**{key:value}),case(.05,.1,.05)])
        nondecreasing = VERIFY.series_checks(
            [case(.2, .3, .2), case(.1, .2, .1), case(.05, .1, .1)])
        self.assertFalse(nondecreasing["valid"])

    def test_series_checks_requires_valid_audits_and_finite_errors(self):
        def item(h, valid=True, scalar_error=.3, velocity_error=.2):
            return {"valid": valid, "reference": {
                "h": h, "dt": .01, "steps": 2, "time": .02, "nu": .1,
                "diffusivity": .02, "speed": 1., "velocityRelaxation": .6,
                "coupledFlowTolerance": 1e-8, "flowTolerance": 1e-8,
                "coupledFlowToleranceStatus": "explicit-thermal-flow-json",
                "flowToleranceStatus": "explicit-flow-json",
                "scalarRelativeTolerance": 1e-9, "scalarAbsoluteTolerance": 1e-12, "scalarCellTolerance": 1e-9,
                "scalarBackwardEulerL2": scalar_error,
                "velocityBackwardEulerL2": velocity_error,
            }}
        with self.assertRaises(ValueError):
            VERIFY.series_checks([item(.2), item(.1), item(.05, valid=False)])
        with self.assertRaises(ValueError):
            VERIFY.series_checks([item(.2), item(.1, scalar_error=math.nan), item(.05)])

    def test_command_line_rejects_invalid_study_arguments(self):
        with tempfile.TemporaryDirectory(prefix="cartmesh-thermal-scale-args-") as name:
            for extra in (("--cells-across", "4", "5"), ("--dt", "0"),
                          ("--steps", "0"), ("--timeout", "0"), ("--timeout", "601"),
                          ("--timeout", "nan")):
                result = subprocess.run(
                    [sys.executable, str(VERIFY_PATH), "--output", str(Path(name) / "study"),
                     *extra], text=True, capture_output=True, timeout=10)
                self.assertEqual(result.returncode, 2, result.stderr)
            result = subprocess.run([sys.executable, str(VERIFY_PATH), "--output", str(Path(name) / "bad-flow-tol"), "--flow-tolerance", "nan"], text=True, capture_output=True, timeout=10)
            self.assertEqual(result.returncode, 2, result.stderr)


class ThermalScaleArtifactTests(unittest.TestCase):
    """Use the real small CLI fixture when CMake supplies its three binaries."""

    @classmethod
    def setUpClass(cls):
        if not os.environ.get("CARTMESH_MESH_CLI"):
            raise unittest.SkipTest("CLI paths were not supplied")
        def option(name, env, default):
            return Path(os.environ.get(env, default)).resolve()
        cls.mesh_cli = option("--mesh-cli", "CARTMESH_MESH_CLI", "build/cartmesh2d_cli")
        cls.flow_cli = option("--flow-cli", "CARTMESH_FLOW_CLI", "build/cartmesh2d_flow_cli")
        cls.transport_cli = option("--transport-cli", "CARTMESH_TRANSPORT_CLI", "build/cartmesh2d_transport_cli")
        cls.tmp = tempfile.TemporaryDirectory(prefix="cartmesh-thermal-scale-")
        cls.root = Path(cls.tmp.name)
        boundary = cls.root / "unit.xy"
        boundary.write_text("0 0\n1 0\n1 1\n0 1\n")
        prefix = cls.root / "mesh"
        result = subprocess.run(
            [str(cls.mesh_cli), boundary, prefix, "4", str(1 / 14), ".1",
             "interior", cls.root / "mesh-openfoam", "4", "0"],
            text=True, capture_output=True, timeout=30)
        if result.returncode:
            raise RuntimeError(result.stderr or result.stdout)
        cls.mesh = prefix.with_suffix(".solver.cm2d")
        cls.scalar = cls.root / "thermal"
        cls.flow = cls.root / "flow"
        common = ["--mesh", cls.mesh, "--dt", ".01", "--steps", "2"]
        result = subprocess.run(
            [str(cls.transport_cli), *common, "--output", cls.scalar,
             "--verification", "thermal-vortex", "--diffusivity", ".02",
             "--flow-nu", ".1", "--flow-convection", "limited-linear",
             "--convection", "limited-linear", "--pressure-preconditioner", "aggregation"],
            text=True, capture_output=True, timeout=30)
        if result.returncode:
            raise RuntimeError(result.stderr or result.stdout)
        result = subprocess.run(
            [str(cls.flow_cli), "--mesh", cls.mesh, "--output", cls.flow,
             "--case", "taylor-green", "--time-step", ".01", "--steps", "2",
             "--nu", ".1", "--speed", "1", "--tolerance", "1e-8",
             "--convection", "limited-linear", "--pressure-preconditioner", "aggregation"],
            text=True, capture_output=True, timeout=30)
        if result.returncode:
            raise RuntimeError(result.stderr or result.stdout)

    def test_relaxation_flag_validation_and_steady_frozen_rejection(self):
        for value in ("0", "-1", "1.1", "nan"):
            flow = subprocess.run(
                [str(self.flow_cli), "--mesh", self.mesh, "--output", self.root / ("bad-flow-" + value),
                 "--case", "manufactured", "--velocity-relaxation", value],
                text=True, capture_output=True, timeout=10)
            self.assertNotEqual(flow.returncode, 0, value)
            transport = subprocess.run(
                [str(self.transport_cli), "--mesh", self.mesh, "--output", self.root / ("bad-transport-" + value),
                 "--verification", "decay", "--dt", ".01", "--steps", "1",
                 "--flow-velocity-relaxation", value],
                text=True, capture_output=True, timeout=10)
            self.assertNotEqual(transport.returncode, 0, value)
        frozen = subprocess.run(
            [str(self.transport_cli), "--mesh", self.mesh, "--output", self.root / "frozen-relaxation",
             "--verification", "decay", "--dt", ".01", "--steps", "1",
             "--flow-velocity-relaxation", ".8"],
            text=True, capture_output=True, timeout=10)
        self.assertNotEqual(frozen.returncode, 0)
        self.assertIn('requires evolving flow',frozen.stderr)
        steady=subprocess.run([
            str(self.flow_cli),'--mesh',self.mesh,'--output',self.root/'steady-relaxation',
            '--case','manufactured','--velocity-relaxation','.8'],
            text=True,capture_output=True,timeout=10)
        self.assertNotEqual(steady.returncode,0)
        self.assertIn('requires transient flow',steady.stderr)

    def test_valid_relaxation_evolution_matches_standalone_carrier(self):
        thermal = self.root / "thermal-alpha08"
        flow = self.root / "flow-alpha08"
        common = ["--mesh", self.mesh, "--dt", ".01", "--steps", "2"]
        result = subprocess.run(
            [str(self.transport_cli), *common, "--output", thermal,
             "--verification", "thermal-vortex", "--diffusivity", ".02",
             "--flow-nu", ".1", "--flow-velocity-relaxation", ".8",
             "--flow-convection", "limited-linear", "--convection", "limited-linear",
             "--pressure-preconditioner", "aggregation"],
            text=True, capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
        result = subprocess.run(
            [str(self.flow_cli), "--mesh", self.mesh, "--output", flow,
             "--case", "taylor-green", "--time-step", ".01", "--steps", "2",
             "--nu", ".1", "--speed", "1", "--velocity-relaxation", ".8",
             "--tolerance", "1e-8", "--convection", "limited-linear",
             "--pressure-preconditioner", "aggregation"],
            text=True, capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
        thermal_info = json.loads(Path(str(thermal) + ".json").read_text())
        flow_info = json.loads(Path(str(flow) + ".json").read_text())
        self.assertAlmostEqual(thermal_info["flowVelocityRelaxation"], .8)
        self.assertAlmostEqual(flow_info["velocityRelaxation"], .8)
        self.assertEqual(Path(str(thermal) + ".carrier.checkpoint").read_bytes(),
                         Path(str(flow) + ".checkpoint").read_bytes())
        audit = VERIFY.audit(thermal, flow, self.root / "alpha08-audit.json")
        self.assertTrue(audit["valid"], audit["issues"])

    def test_explicit_default_preserves_carrier_checkpoint(self):
        prefix=self.root/'explicit-default'
        result=subprocess.run([
            str(self.flow_cli),'--mesh',self.mesh,'--output',prefix,'--case','taylor-green',
            '--time-step','.01','--steps','2','--nu','.1','--speed','1','--tolerance','1e-8',
            '--convection','limited-linear','--pressure-preconditioner','aggregation',
            '--velocity-relaxation','.6'],text=True,capture_output=True,timeout=30)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(Path(str(prefix)+'.checkpoint').read_bytes(),
                         Path(str(self.flow)+'.checkpoint').read_bytes())

    def test_generate_small_nondefault_flow_tolerance_is_propagated_and_series_checked(self):
        study_dir = self.root / 'generated-tolerance-1e-9'
        command = [sys.executable, str(VERIFY_PATH), '--output', study_dir,
                   '--cells-across', '6', '10', '14', '--flow-tolerance', '1e-9',
                   '--timeout', '30', '--mesh-cli', self.mesh_cli,
                   '--flow-cli', self.flow_cli, '--transport-cli', self.transport_cli]
        result = subprocess.run(command, text=True, capture_output=True, timeout=180)
        self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
        report = json.loads((study_dir / 'study.json').read_text())
        self.assertEqual(report['flowTolerance'], 1e-9)
        self.assertTrue(report['spatial']['valid'])
        for case in report['cases']:
            reference = case['reference']
            self.assertEqual(reference['coupledFlowTolerance'], 1e-9)
            self.assertEqual(reference['flowTolerance'], 1e-9)
            self.assertEqual(reference['coupledFlowToleranceStatus'], 'explicit-thermal-flow-json')
            self.assertEqual(reference['flowToleranceStatus'], 'explicit-flow-json')
        flow_tolerances=[]; solver_tolerances=[]
        for run in report['runs']:
            command=[str(x) for x in run['command']]
            for i,token in enumerate(command[:-1]):
                if token == '--flow-tolerance': flow_tolerances.append(float(command[i+1]))
                if token == '--tolerance': solver_tolerances.append(float(command[i+1]))
        self.assertGreaterEqual(sum(math.isclose(x,1e-9,rel_tol=0,abs_tol=1e-15) for x in flow_tolerances), 3)
        self.assertGreaterEqual(sum(math.isclose(x,1e-9,rel_tol=0,abs_tol=1e-15) for x in solver_tolerances), 3)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "tmp"):
            cls.tmp.cleanup()

    def copy_case(self, name):
        dst = self.root / name
        for suffix in (".json", ".cells.csv", ".faces.csv", ".history.csv",
                       ".thermal-history.csv", ".thermal.checkpoint", ".carrier.checkpoint"):
            shutil.copyfile(str(self.scalar) + suffix, str(dst) + suffix)
        flow_dst = self.root / (name + "-flow")
        for suffix in (".json", ".cells.csv", ".faces.csv", ".residuals.csv",
                       ".time-history.csv", ".checkpoint"):
            shutil.copyfile(str(self.flow) + suffix, str(flow_dst) + suffix)
        return dst, flow_dst

    def test_real_artifacts_pass_reference_and_binding_tampering_is_rejected(self):
        scalar_prefix, flow_prefix = self.copy_case("pass")
        result = VERIFY.reference_errors(scalar_prefix, flow_prefix)
        self.assertEqual(result["cells"], 196)
        self.assertGreater(result["scalarContinuousL2"], result["scalarBackwardEulerL2"])
        self.assertGreater(result["velocityContinuousL2"], result["velocityBackwardEulerL2"])
        self.assertEqual(result["coupledFlowTolerance"], 1e-8)
        self.assertEqual(result["flowTolerance"], 1e-8)

        legacy, legacy_flow = self.copy_case("legacy-missing-flow-tolerance")
        legacy_metadata = json.loads(Path(str(legacy) + ".json").read_text())
        legacy_metadata.pop("flowTolerance", None)
        Path(str(legacy) + ".json").write_text(json.dumps(legacy_metadata))
        legacy_result = VERIFY.reference_errors(legacy, legacy_flow)
        self.assertEqual(legacy_result["coupledFlowToleranceStatus"], "legacy-unknown")
        explicit_null, explicit_null_flow = self.copy_case("explicit-null-flow-tolerance")
        null_metadata = json.loads(Path(str(explicit_null) + ".json").read_text())
        null_metadata["flowTolerance"] = None
        Path(str(explicit_null) + ".json").write_text(json.dumps(null_metadata))
        with self.assertRaises(ValueError):
            VERIFY.reference_errors(explicit_null, explicit_null_flow)

        for suffix, mutate in (
            (".cells.csv", self._swap_first_two_cells),
            (".json", self._tamper_scalar_speed),
        ):
            bad, bad_flow = self.copy_case("bad-" + suffix[1:-4])
            mutate(bad, bad_flow)
            with self.assertRaises(ValueError):
                VERIFY.reference_errors(bad, bad_flow)

        bad, bad_flow = self.copy_case("bad-checkpoint")
        checkpoint = Path(str(bad) + ".carrier.checkpoint")
        checkpoint.write_text(checkpoint.read_text() + "\n")
        with self.assertRaises(ValueError):
            VERIFY.reference_errors(bad, bad_flow)

        bad, bad_flow = self.copy_case("bad-flow-tolerance")
        flow_metadata = json.loads(Path(str(bad_flow) + ".json").read_text())
        flow_metadata["tolerance"] = 2e-8
        Path(str(bad_flow) + ".json").write_text(json.dumps(flow_metadata))
        with self.assertRaises(ValueError):
            VERIFY.reference_errors(bad, bad_flow, expected_flow_tolerance=1e-8)

        bad, bad_flow = self.copy_case("bad-relaxation")
        metadata = json.loads(Path(str(bad) + ".json").read_text())
        metadata["flowVelocityRelaxation"] = .8
        Path(str(bad) + ".json").write_text(json.dumps(metadata))
        with self.assertRaises(ValueError):
            VERIFY.reference_errors(bad, bad_flow)

    def test_independent_audit_and_iteration_history_contract(self):
        result=VERIFY.audit(self.scalar,self.flow,self.root/'independent.json')
        self.assertTrue(result['valid'],result['issues'])
        counts=result['iterations']
        self.assertEqual(counts['steps'],2)
        self.assertGreater(counts['flowNonlinearIterations'],0)
        self.assertGreater(counts['scalarCorrections'],0)
        self.assertGreater(counts['scalarLinearIterations'],0)
        for label,kind in [('missing-step','coverage'),('negative-linear','negative'),('wrong-clock','clock')]:
            bad,flow=self.copy_case(label)
            path=Path(str(bad)+'.history.csv')
            with path.open() as f:reader=csv.DictReader(f);fields=reader.fieldnames;rows=list(reader)
            if kind=='coverage':rows=[r for r in rows if r['step']=='1']
            elif kind=='negative':rows[0]['linearIterations']='-1'
            else:rows[0]['time']='0.125'
            with path.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
            with self.assertRaises(ValueError):VERIFY.iteration_counts(bad,flow)

    def test_audit_converts_missing_and_bad_metadata_to_invalid_report(self):
        missing=self.root/'missing-metadata'; flow=self.root/'missing-metadata-flow'
        shutil.copyfile(str(self.scalar)+'.json',str(missing)+'.json')
        shutil.copyfile(str(self.flow)+'.json',str(flow)+'.json')
        Path(str(missing)+'.json').unlink()
        report=VERIFY.audit(missing,flow,self.root/'missing-audit.json')
        self.assertFalse(report['valid'])
        self.assertTrue(report['issues'])
        bad,flow=self.copy_case('bad-metadata')
        metadata=json.loads(Path(str(bad)+'.json').read_text());metadata['steps']=[]
        Path(str(bad)+'.json').write_text(json.dumps(metadata))
        report=VERIFY.audit(bad,flow,self.root/'bad-metadata-audit.json')
        self.assertFalse(report['valid'])
        self.assertTrue(report['issues'])

    def _swap_first_two_cells(self, prefix, _flow):
        path = Path(str(prefix) + ".cells.csv")
        lines = path.read_text().splitlines()
        lines[1], lines[2] = lines[2], lines[1]
        path.write_text("\n".join(lines) + "\n")

    def _tamper_scalar_speed(self, prefix, _flow):
        path = Path(str(prefix) + ".json")
        data = json.loads(path.read_text())
        data["flowSpeed"] = 2.
        path.write_text(json.dumps(data))


if __name__ == "__main__":
    for option, env in (("--mesh-cli", "CARTMESH_MESH_CLI"),
                        ("--flow-cli", "CARTMESH_FLOW_CLI"),
                        ("--transport-cli", "CARTMESH_TRANSPORT_CLI")):
        if option in sys.argv:
            index = sys.argv.index(option)
            os.environ[env] = sys.argv[index + 1]
            del sys.argv[index:index + 2]
    unittest.main()
