#!/usr/bin/env python3
"""Focused negative tests for the independent native-flow momentum audit."""

import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools/verification"))
import verify_native_flow as verifier


def _flow_cli() -> Path:
    if "--cli" in sys.argv:
        index = sys.argv.index("--cli")
        try:
            value = sys.argv[index + 1]
        except IndexError as exc:
            raise SystemExit("--cli requires a path") from exc
        del sys.argv[index:index + 2]
        return Path(value).resolve()
    return Path(os.environ.get("CARTMESH_FLOW_CLI", "build/cartmesh2d_flow_cli")).resolve()


FLOW_CLI = _flow_cli()


def tiny_mesh() -> verifier.Mesh:
    return verifier.Mesh(
        Path("tiny.solver.cm2d"),
        ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
        (verifier.Edge(0, 0, 1, 0, -1, 1),),
        (verifier.Cell(0, 1.0, (0, 1, 2), (0,)),),
        (0, 0, 0, 0, 0, 0, 0),
    )


class FlowVerifierSchemaTests(unittest.TestCase):
    def test_legacy_faces_are_explicitly_unavailable(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "faces.csv"
            path.write_text("face,owner,neighbour,flux\n0,0,-1,0\n", encoding="utf-8")
            rows, available = verifier.read_faces(path, tiny_mesh())
            self.assertFalse(available)
            self.assertEqual(rows[0]["flux"], 0.0)

    def test_partial_face_schema_fails(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "faces.csv"
            path.write_text(
                "face,owner,neighbour,flux,pressure,advectionX,diffusionX,diffusionY\n"
                "0,0,-1,0,0,0,0,0\n", encoding="utf-8")
            with self.assertRaises(verifier.VerificationError):
                verifier.read_faces(path, tiny_mesh())


class FlowVerifierTamperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(tempfile.mkdtemp(prefix="cartmesh-flow-verifier-"))
        cls.mesh = cls.root / "channel.solver.cm2d"
        cls.prefix = cls.root / "channel"
        nx, ny, width = 32, 8, 4.0
        vertices = [(i * width / nx, j / ny)
                    for j in range(ny + 1) for i in range(nx + 1)]
        edges, cells, by_pair = [], [], {}
        for j in range(ny):
            for i in range(nx):
                ids = [j * (nx + 1) + i, j * (nx + 1) + i + 1,
                       (j + 1) * (nx + 1) + i + 1, (j + 1) * (nx + 1) + i]
                incident = []
                for a, b in zip(ids, ids[1:] + ids[:1]):
                    pair = tuple(sorted((a, b)))
                    if pair not in by_pair:
                        by_pair[pair] = len(edges)
                        edges.append([len(edges), a, b, len(cells), -1, 1])
                    else:
                        edges[by_pair[pair]][4:] = [len(cells), 0]
                    incident.append(by_pair[pair])
                cells.append([len(cells), 0, 0, width / nx / ny, 4, *ids, 4, *incident])
        lines = ["CM2D 1", f"VERTICES {len(vertices)}"]
        lines.extend(f"{i} {x:.17g} {y:.17g}" for i, (x, y) in enumerate(vertices))
        lines.extend([f"EDGES {len(edges)}", *(" ".join(map(str, edge)) for edge in edges),
                      f"CELLS {len(cells)}", *(" ".join(map(str, cell)) for cell in cells),
                      "AUDIT 0 0 0 0 0 0 0", "END"])
        cls.mesh.write_text("\n".join(lines) + "\n", encoding="utf-8")
        result = subprocess.run([str(FLOW_CLI), "--mesh", str(cls.mesh), "--output", str(cls.prefix),
                                 "--case", "channel", "--nu", ".01", "--speed", "1",
                                 "--max-iterations", "700"], capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            raise RuntimeError(f"flow CLI fixture failed: {result.stdout[-1000:]} {result.stderr[-1000:]}")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.root, ignore_errors=True)

    def _fixture(self):
        return self.mesh, self.prefix

    @staticmethod
    def _args():
        return SimpleNamespace(
            geometry_absolute_tolerance=1e-11, geometry_relative_tolerance=1e-9,
            continuity_absolute_tolerance=1e-10, continuity_relative_tolerance=1e-7,
            max_reported_continuity=1e-8, max_iterations=1500,
            channel_velocity_l2=.08, channel_transverse_l2=.05,
            channel_pressure_l2=.12, channel_gradient_error=.12,
            channel_flux_error=.08, cavity_rmse=.12, cavity_max_error=.30,
            max_speed_ratio=4.0, external_lift_drag_ratio=.30,
        )

    def _copy_fixture(self, mesh, prefix, folder):
        copied_mesh = Path(folder) / "channel.solver.cm2d"
        shutil.copyfile(mesh, copied_mesh)
        copied_prefix = Path(folder) / "channel"
        for suffix in (".cells.csv", ".faces.csv", ".residuals.csv", ".json"):
            shutil.copyfile(Path(str(prefix) + suffix), Path(str(copied_prefix) + suffix))
        return copied_mesh, copied_prefix

    def test_untampered_and_legacy_status(self):
        mesh, prefix = self._fixture()
        result = verifier.verify_case(mesh, prefix, "channel", .01, 1.0, self._args())
        self.assertTrue(result["valid"], result["issues"])
        self.assertTrue(result["momentumAudit"]["valid"])
        with tempfile.TemporaryDirectory() as folder:
            copied_mesh, copied_prefix = self._copy_fixture(mesh, prefix, folder)
            faces = Path(str(copied_prefix) + ".faces.csv")
            with faces.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            columns = ["face", "owner", "neighbour", "flux"]
            with faces.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
                writer.writeheader(); writer.writerows(rows)
            summary = Path(str(copied_prefix) + ".json")
            data = json.loads(summary.read_text())
            for key in ("pressureForceX", "pressureForceY", "discreteForceX", "discreteForceY",
                        "pressureDiscretization", "convection"):
                del data[key]
            summary.write_text(json.dumps(data))
            legacy = verifier.verify_case(copied_mesh, copied_prefix, "channel", .01, 1.0, self._args())
            self.assertTrue(legacy["valid"], legacy["issues"])
            self.assertEqual(legacy["momentumAudit"]["status"], "unavailable")
            self.assertFalse(legacy["momentumAudit"]["valid"])

    def test_face_terms_tamper_are_rejected(self):
        mesh, prefix = self._fixture()
        for field in ("pressure", "advectionX", "diffusionY"):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as folder:
                copied_mesh, copied_prefix = self._copy_fixture(mesh, prefix, folder)
                faces = Path(str(copied_prefix) + ".faces.csv")
                with faces.open(newline="") as stream:
                    rows = list(csv.DictReader(stream))
                rows[0][field] = str(float(rows[0][field]) + .01)
                with faces.open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
                    writer.writeheader(); writer.writerows(rows)
                result = verifier.verify_case(copied_mesh, copied_prefix, "channel", .01, 1.0, self._args())
                self.assertFalse(result["valid"])
                self.assertFalse(result["momentumAudit"]["valid"])

    def test_cell_field_tamper_is_rejected(self):
        mesh, prefix = self._fixture()
        with tempfile.TemporaryDirectory() as folder:
            copied_mesh, copied_prefix = self._copy_fixture(mesh, prefix, folder)
            cells = Path(str(copied_prefix) + ".cells.csv")
            with cells.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["p"] = str(float(rows[0]["p"]) + .01)
            with cells.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
                writer.writeheader(); writer.writerows(rows)
            result = verifier.verify_case(copied_mesh, copied_prefix, "channel", .01, 1.0, self._args())
            self.assertFalse(result["valid"])
            self.assertFalse(result["momentumAudit"]["valid"])

    def test_deleting_new_face_column_fails(self):
        mesh, prefix = self._fixture()
        with tempfile.TemporaryDirectory() as folder:
            copied_mesh, copied_prefix = self._copy_fixture(mesh, prefix, folder)
            faces = Path(str(copied_prefix) + ".faces.csv")
            with faces.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            rows = [{key: value for key, value in row.items() if key != "diffusionY"} for row in rows]
            with faces.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
                writer.writeheader(); writer.writerows(rows)
            with self.assertRaises(verifier.VerificationError):
                verifier.verify_case(copied_mesh, copied_prefix, "channel", .01, 1.0, self._args())


if __name__ == "__main__":
    unittest.main()
