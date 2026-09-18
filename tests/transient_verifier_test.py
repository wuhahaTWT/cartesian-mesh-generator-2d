#!/usr/bin/env python3
"""Portable independent transient-verifier tests using a real CLI result."""
import csv, importlib.util, json, os, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
from transient_flow_cli_test import rectangle
SPEC = importlib.util.spec_from_file_location("verify_transient_flow", ROOT / "tools/verification/verify_transient_flow.py")
VERIFIER = importlib.util.module_from_spec(SPEC); assert SPEC.loader is not None; SPEC.loader.exec_module(VERIFIER)

class TransientVerifierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cli=sys.argv[sys.argv.index("--cli")+1] if "--cli" in sys.argv else os.environ.get("CARTMESH_FLOW_CLI","build/cartmesh2d_flow_cli")
        cls.cli=Path(cli).resolve(); cls.temp=tempfile.TemporaryDirectory(prefix="cartmesh-transient-verifier-"); cls.root=Path(cls.temp.name)
        cls.mesh=cls.root/"taylor.solver.cm2d"; rectangle(cls.mesh); cls.source=cls.root/"source"
        cmd=[str(cls.cli),"--mesh",str(cls.mesh),"--output",str(cls.source),"--case","taylor-green","--nu",".01","--speed","1","--tolerance","1e-9","--max-iterations","150","--time-step",".02","--steps","2"]
        result=subprocess.run(cmd,text=True,capture_output=True,timeout=30)
        if result.returncode: raise RuntimeError(result.stderr or result.stdout)
    @classmethod
    def tearDownClass(cls): cls.temp.cleanup()
    def copy_case(self,name):
        p=self.root/name
        for s in (".cells.csv",".faces.csv",".json",".residuals.csv",".time-history.csv"): shutil.copyfile(str(self.source)+s,str(p)+s)
        return p
    def assert_rejected(self,mutate,name):
        p=self.copy_case(name)
        self.assertTrue(VERIFIER.verify(self.mesh,p,self.root/(name+".before.json"))["valid"])
        mutate(p)
        with self.assertRaises(VERIFIER.native.VerificationError): VERIFIER.verify(self.mesh,p,self.root/(name+".audit.json"))
    def test_failed_cli_verification_overwrites_stale_pass(self):
        prefix = self.copy_case("cli_failure")
        output = self.root / "cli-failure.audit.json"
        VERIFIER.verify(self.mesh, prefix, output)
        summary_path = Path(str(prefix)+".json")
        summary = json.loads(summary_path.read_text())
        summary["dt"] = .03
        summary_path.write_text(json.dumps(summary))
        completed = subprocess.run([sys.executable, str(ROOT / "tools/verification/verify_transient_flow.py"),
            "--mesh", str(self.mesh), "--prefix", str(prefix), "--output", str(output)], capture_output=True, text=True, timeout=20)
        self.assertEqual(completed.returncode, 1)
        self.assertIs(json.loads(output.read_text())["valid"], False)

    def test_real_result_passes(self): self.assertTrue(VERIFIER.verify(self.mesh,self.source,self.root/"pass.audit.json")["valid"])
    def test_temporal_csv_and_missing_schema_rejected(self):
        def temporal(p):
            f=Path(str(p)+".cells.csv")
            with f.open() as stream: rows=list(csv.DictReader(stream))
            rows[0]["temporalX"]=str(float(rows[0]["temporalX"])+.1)
            with f.open("w",newline="") as s: w=csv.DictWriter(s,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
        self.assert_rejected(temporal,"bad_temporal")
        def missing(p):
            f=Path(str(p)+".cells.csv"); lines=f.read_text().splitlines(); lines[0]=lines[0].replace(",temporalY",""); f.write_text("\n".join(lines)+"\n")
        self.assert_rejected(missing,"missing_temporal")
        def face(p):
            f=Path(str(p)+".faces.csv"); rows=list(csv.DictReader(f.read_text().splitlines())); rows[0]["advectionX"]=str(float(rows[0]["advectionX"])+.1)
            with f.open("w",newline="") as s: w=csv.DictWriter(s,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
        self.assert_rejected(face,"bad_face_momentum")
    def test_summary_cfl_force_dt_nan_and_status_rejected(self):
        for key,val,name in (("maxCourant",9.,"bad_cfl"),("forceX",1.,"bad_force"),("dt",.03,"bad_dt")):
            def mutate(p,key=key,val=val):
                f=Path(str(p)+".json"); d=json.loads(f.read_text()); d[key]=val; f.write_text(json.dumps(d))
            self.assert_rejected(mutate,name)
        def nan(p):
            f=Path(str(p)+".time-history.csv"); f.write_text(f.read_text().replace(",0.02,",",nan,",1))
        self.assert_rejected(nan,"bad_nan")
        def status(p):
            f=Path(str(p)+".json"); d=json.loads(f.read_text()); d["status"]="time_step_not_converged"; d["converged"]=False; f.write_text(json.dumps(d))
        self.assert_rejected(status,"bad_status")
    def test_negative_first_time_and_false_failed_last_row_rejected(self):
        def negative(p):
            f=Path(str(p)+".time-history.csv"); ls=f.read_text().splitlines(); x=ls[1].split(","); x[1]="-0.02"; ls[1]=",".join(x); f.write_text("\n".join(ls)+"\n")
        self.assert_rejected(negative,"bad_time")
        def failed(p):
            f=Path(str(p)+".time-history.csv"); ls=f.read_text().splitlines(); x=ls[-1].split(","); x[3]="0"; ls[-1]=",".join(x); f.write_text("\n".join(ls)+"\n")
        self.assert_rejected(failed,"bad_failed")
        def stale(p):
            failed(p)
            f=Path(str(p)+".json"); d=json.loads(f.read_text()); d["status"]="converged"; d["converged"]=True; f.write_text(json.dumps(d))
        self.assert_rejected(stale,"stale_valid_summary")

if __name__ == "__main__":
    # Keep the runner portable while accepting the repository test convention.
    if "--cli" in sys.argv:
        i = sys.argv.index("--cli"); os.environ["CARTMESH_FLOW_CLI"] = sys.argv[i + 1]; del sys.argv[i:i + 2]
    unittest.main()
