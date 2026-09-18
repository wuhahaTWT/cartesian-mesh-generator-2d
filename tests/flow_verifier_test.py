#!/usr/bin/env python3
"""Focused negative tests for the independent native-flow momentum audit."""

import csv
import json
import math
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


class CavitySamplingTests(unittest.TestCase):
    @staticmethod
    def samples(n, function, scale=1.0):
        return [dict(x=scale * (i + .5) / n, y=scale * (j + .5) / n,
                     u=function((i + .5) / n, (j + .5) / n))
                for j in range(n) for i in range(n)]

    def test_affine_field_at_benchmark_points_across_scales(self):
        f = lambda x, y: .7 + 1.3 * x - .8 * y
        for scale in (1e-6, 1., 1e6):
            for n in (14, 30, 62):
                rows = self.samples(n, f, scale)
                for coordinate, _ in verifier.GHIA_U[1:-1]:
                    x, y = .5, coordinate
                    walls = [dict(x=x*scale, y=endpoint*scale, u=f(x, endpoint))
                             for endpoint in (0., 1.)]
                    with self.subTest(scale=scale, n=n, y=y):
                        self.assertAlmostEqual(verifier.affine_sample(
                            rows, x*scale, y*scale, "u", boundary=walls), f(x, y), places=13)
        # Keep the minimal evidence for why the former average was replaced.
        rows = self.samples(14, f)
        walls = [dict(x=.5, y=e, u=f(.5, e)) for e in (0., 1.)]
        legacy_errors = [abs(verifier.idw(rows, .5, y, "u", boundary=walls)-f(.5, y))
                         for y, _ in verifier.GHIA_U[1:-1]]
        self.assertGreater(max(legacy_errors), .01)

    def test_skewed_samples_and_exact_points(self):
        f = lambda x, y: 2. + 3. * x - 4. * y
        rows = [dict(x=x+.3*y, y=y, u=f(x+.3*y, y))
                for x, y in ((0., 0.), (1., 0.), (0., 1.), (1., 1.), (.4, .8))]
        self.assertAlmostEqual(verifier.affine_sample(rows, .4, .3, "u"), f(.4, .3), places=14)
        self.assertEqual(verifier.affine_sample(rows, 0., 0., "u"), 2.)
        walls = [dict(x=.5, y=1., u=1.)]
        self.assertEqual(verifier.affine_sample(rows, .5, 1., "u", boundary=walls), 1.)
        self.assertEqual(verifier.affine_sample(list(reversed(rows)), .4, .3, "u"),
                         verifier.affine_sample(rows, .4, .3, "u"))

    def test_smooth_quadratic_sampling_refines(self):
        f = lambda x, y: x*x + y*y
        errors = []
        for n in (8, 16, 32):
            rows = self.samples(n, f)
            walls = [dict(x=.5, y=e, u=f(.5, e)) for e in (0., 1.)]
            differences = [verifier.affine_sample(rows, .5, y, "u", boundary=walls)-f(.5, y)
                           for y, _ in verifier.GHIA_U[1:-1]]
            errors.append(math.sqrt(sum(e*e for e in differences)/len(differences)))
        self.assertGreater(errors[0]/errors[1], 3.)
        self.assertGreater(errors[1]/errors[2], 3.)

    def test_invalid_or_rank_deficient_samples_fail(self):
        line = [dict(x=float(i), y=0., u=1.) for i in range(8)]
        with self.assertRaisesRegex(verifier.VerificationError, "rank-deficient"):
            verifier.affine_sample(line, .5, .2, "u")
        with self.assertRaises(verifier.VerificationError):
            verifier.affine_sample([], .5, .2, "u")
        with self.assertRaises(verifier.VerificationError):
            verifier.affine_sample(line, math.nan, .2, "u")
        with self.assertRaises(verifier.VerificationError):
            verifier.affine_sample(line, .5, .2, "u", count=2)
        line[0]["u"] = math.nan
        with self.assertRaises(verifier.VerificationError):
            verifier.affine_sample(line, .5, .2, "u")
        duplicate = [dict(x=0., y=0., u=float(i)) for i in range(3)]
        with self.assertRaisesRegex(verifier.VerificationError, "conflicting"):
            verifier.affine_sample(duplicate, 0., 0., "u")

    def test_cavity_sequence_rejects_mixed_sampling_methods(self):
        cases = [dict(case="cavity", label=str(i), meshMeasurement=dict(characteristicH=h),
                      benchmark=dict(centrelineRmse=e, samplingMethod=method))
                 for i, h, e, method in ((0, .1, .02, "legacy"), (1, .05, .01, "affine-v2"))]
        result = verifier.sequence_checks(cases)
        self.assertFalse(result["valid"])
        self.assertTrue(any("sampling" in issue for issue in result["issues"]))


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
                        "pressureDiscretization", "convection", "viscousStress",
                        "forceDefinition", "reconstructedForceX", "reconstructedForceY",
                        "wallForceX", "wallForceY", "wallViscousForceX", "wallViscousForceY",
                        "wallForceDefinition"):
                data.pop(key, None)
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

    def test_selected_force_tamper_is_rejected(self):
        mesh, prefix = self._fixture()
        with tempfile.TemporaryDirectory() as folder:
            copied_mesh, copied_prefix = self._copy_fixture(mesh, prefix, folder)
            summary = Path(str(copied_prefix) + ".json")
            payload = json.loads(summary.read_text())
            if "viscousStress" not in payload:
                self.skipTest("fixture is from legacy flow CLI without viscousStress metadata")
            payload["forceX"] = float(payload["forceX"]) + .01
            summary.write_text(json.dumps(payload), encoding="utf-8")
            result = verifier.verify_case(copied_mesh, copied_prefix, "channel", .01, 1.0, self._args())
            self.assertFalse(result["valid"])
            self.assertFalse(result["momentumAudit"]["valid"])

    def test_unknown_viscous_stress_is_rejected(self):
        mesh, prefix = self._fixture()
        with tempfile.TemporaryDirectory() as folder:
            copied_mesh, copied_prefix = self._copy_fixture(mesh, prefix, folder)
            summary = Path(str(copied_prefix) + ".json")
            payload = json.loads(summary.read_text())
            payload["viscousStress"] = "bogus"
            summary.write_text(json.dumps(payload), encoding="utf-8")
            result = verifier.verify_case(copied_mesh, copied_prefix, "channel", .01, 1.0, self._args())
            self.assertFalse(result["valid"])
            self.assertTrue(any("viscousStress" in issue for issue in result["issues"]))

    def test_toggled_viscous_stress_metadata_is_rejected(self):
        mesh, prefix = self._fixture()
        with tempfile.TemporaryDirectory() as folder:
            copied_mesh, copied_prefix = self._copy_fixture(mesh, prefix, folder)
            summary = Path(str(copied_prefix) + ".json")
            payload = json.loads(summary.read_text())
            if "viscousStress" not in payload:
                self.skipTest("fixture is from legacy flow CLI without viscousStress metadata")
            payload["viscousStress"] = "laplacian" if payload["viscousStress"] == "symmetric" else "symmetric"
            summary.write_text(json.dumps(payload), encoding="utf-8")
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

class FlowVerifierManufacturedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(tempfile.mkdtemp(prefix="cartmesh-manufactured-verifier-"))
        cls.mesh = cls.root / "manufactured.solver.cm2d"
        cls.prefix = cls.root / "manufactured"
        cls.slope_prefix = cls.root / "manufactured-slope1"
        nx = ny = 16
        vertices = [(i / nx, j / ny) for j in range(ny + 1) for i in range(nx + 1)]
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
                cells.append([len(cells), 0, 0, 1.0 / nx / ny, 4, *ids, 4, *incident])
        lines = ["CM2D 1", f"VERTICES {len(vertices)}"]
        lines.extend(f"{i} {x:.17g} {y:.17g}" for i, (x, y) in enumerate(vertices))
        lines.extend([f"EDGES {len(edges)}", *(
            " ".join(map(str, edge)) for edge in edges),
            f"CELLS {len(cells)}", *(
                " ".join(map(str, cell)) for cell in cells),
            "AUDIT 0 0 0 0 0 0 0", "END"])
        cls.mesh.write_text("\n".join(lines) + "\n", encoding="utf-8")
        result = subprocess.run(
            [str(FLOW_CLI), "--mesh", str(cls.mesh), "--output", str(cls.prefix),
             "--case", "manufactured", "--nu", ".1", "--speed", "1",
             "--convection", "limited-linear", "--max-iterations", "1200"],
            capture_output=True, text=True, timeout=90)
        if result.returncode != 0:
            raise RuntimeError(f"manufactured CLI fixture failed: {result.stdout[-1000:]} {result.stderr[-1000:]}")
        slope_result = subprocess.run(
            [str(FLOW_CLI), "--mesh", str(cls.mesh), "--output", str(cls.slope_prefix),
             "--case", "manufactured", "--nu", ".1", "--speed", "1",
             "--manufactured-pressure-slope", "1", "--convection", "limited-linear",
             "--max-iterations", "1200"],
            capture_output=True, text=True, timeout=90)
        if slope_result.returncode != 0:
            raise RuntimeError(f"manufactured slope fixture failed: {slope_result.stdout[-1000:]} {slope_result.stderr[-1000:]}")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.root, ignore_errors=True)

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

    def _copy(self, folder):
        mesh = Path(folder) / "manufactured.solver.cm2d"
        prefix = Path(folder) / "manufactured"
        shutil.copyfile(self.mesh, mesh)
        for suffix in (".cells.csv", ".faces.csv", ".residuals.csv", ".json"):
            shutil.copyfile(Path(str(self.prefix) + suffix), Path(str(prefix) + suffix))
        return mesh, prefix

    def _copy_slope(self, folder):
        mesh = Path(folder) / "manufactured.solver.cm2d"
        prefix = Path(folder) / "manufactured-slope1"
        shutil.copyfile(self.mesh, mesh)
        for suffix in (".cells.csv", ".faces.csv", ".residuals.csv", ".json"):
            shutil.copyfile(Path(str(self.slope_prefix) + suffix), Path(str(prefix) + suffix))
        return mesh, prefix

    def test_independent_mms_and_fields(self):
        result = verifier.verify_case(self.mesh, self.prefix, "manufactured", .1, 1.0, self._args())
        self.assertTrue(result["valid"], result["issues"])
        self.assertTrue(result["benchmark"]["boundaryAudit"]["allStationaryNoSlipWalls"])
        self.assertLess(result["benchmark"]["exactColumnMaxAbsolute"]["exactU"], 1e-10)
        self.assertGreater(result["benchmark"]["velocityL2Relative"], 0.0)
        self.assertEqual(result["benchmark"]["pressureSlope"], 0.0)
        self.assertIn(result["pressureBoundaryReconstruction"], ("zero-normal", "one-sided-linear"))
        self.assertTrue(result["momentumAudit"]["valid"])

    def test_nonzero_pressure_slope_fixture(self):
        mesh, prefix = self.mesh, self.slope_prefix
        args = self._args()
        args.manufactured_pressure_slope = 1.0
        result = verifier.verify_case(mesh, prefix, "manufactured", .1, 1.0, args)
        self.assertTrue(result["valid"], result["issues"])
        self.assertEqual(result["benchmark"]["pressureSlope"], 1.0)

    def test_pressure_slope_payload_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            mesh = Path(folder) / "manufactured.solver.cm2d"
            prefix = Path(folder) / "manufactured-slope1"
            shutil.copyfile(self.mesh, mesh)
            for suffix in (".cells.csv", ".faces.csv", ".residuals.csv", ".json"):
                shutil.copyfile(Path(str(self.slope_prefix) + suffix), Path(str(prefix) + suffix))
            summary = Path(str(prefix) + ".json")
            payload = json.loads(summary.read_text())
            payload["manufacturedPressureSlope"] = 2.0
            summary.write_text(json.dumps(payload), encoding="utf-8")
            args = self._args()
            args.manufactured_pressure_slope = 1.0
            result = verifier.verify_case(mesh, prefix, "manufactured", .1, 1.0, args)
            self.assertFalse(result["valid"])
            self.assertTrue(any("manufacturedPressureSlope" in issue for issue in result["issues"]))

    def test_pressure_boundary_reconstruction_metadata_is_checked(self):
        with tempfile.TemporaryDirectory() as folder:
            mesh = Path(folder) / "manufactured.solver.cm2d"
            prefix = Path(folder) / "manufactured-slope1"
            shutil.copyfile(self.mesh, mesh)
            for suffix in (".cells.csv", ".faces.csv", ".residuals.csv", ".json"):
                shutil.copyfile(Path(str(self.slope_prefix) + suffix), Path(str(prefix) + suffix))
            summary = Path(str(prefix) + ".json")
            payload = json.loads(summary.read_text())
            self.assertIn(payload.get("pressureBoundaryReconstruction"), ("zero-normal", "one-sided-linear"))
            payload["pressureBoundaryReconstruction"] = "unknown-mode"
            summary.write_text(json.dumps(payload), encoding="utf-8")
            args = self._args()
            args.manufactured_pressure_slope = 1.0
            result = verifier.verify_case(mesh, prefix, "manufactured", .1, 1.0, args)
            self.assertFalse(result["valid"])
            self.assertTrue(any("pressureBoundaryReconstruction" in issue for issue in result["issues"]))

    def test_pressure_gradient_two_ring_recovers_linear_tip_field(self):
        polygons = [
            ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
            ((1.0, 0.0), (1.0, 1.0), (0.0, 1.0)),
            ((1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0)),
            ((0.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)),
        ]
        vertex_ids = {}
        vertices = []
        edges = []
        edge_ids = {}
        cell_edges = []
        for cell_id, polygon in enumerate(polygons):
            ids = []
            for point in polygon:
                if point not in vertex_ids:
                    vertex_ids[point] = len(vertices)
                    vertices.append(point)
                ids.append(vertex_ids[point])
            incident = []
            for a, b in zip(ids, ids[1:] + ids[:1]):
                pair = tuple(sorted((a, b)))
                if pair not in edge_ids:
                    edge_ids[pair] = len(edges)
                    edges.append(verifier.Edge(len(edges), a, b, cell_id, -1, 0))
                else:
                    edge_id = edge_ids[pair]
                    old = edges[edge_id]
                    edges[edge_id] = verifier.Edge(old.id, old.v0, old.v1, old.owner, cell_id, 0)
                incident.append(edge_ids[pair])
            area, _ = verifier.polygon(list(polygon))
            cell_edges.append(verifier.Cell(cell_id, area, tuple(ids), tuple(incident)))
        mesh = verifier.Mesh(Path("tip.solver.cm2d"), tuple(vertices), tuple(edges),
                             tuple(cell_edges), (0, 0, 0, 0, 0, 0, 0))
        measured = verifier.measure(mesh, 1e-11, 1e-9)
        geometries = verifier.face_geometry(mesh, measured)
        values = [2.0 * x - 3.0 * y for x, y in measured.centroids]
        gradients = verifier.reconstruct_gradient(
            mesh, measured, geometries, values, [0.0] * len(edges), [False] * len(edges),
            skip_unknown_boundary=True)
        self.assertAlmostEqual(gradients[0][0], 2.0, places=12)
        self.assertAlmostEqual(gradients[0][1], -3.0, places=12)

    def test_nonzero_pressure_slope_is_rejected_for_ordinary_case(self):
        args = self._args()
        args.manufactured_pressure_slope = 1.0
        with self.assertRaises(verifier.VerificationError):
            verifier.verify_case(self.mesh, self.prefix, "channel", .1, 1.0, args)

    def test_source_and_exact_tampering_are_rejected(self):
        for field in ("sourceX", "exactU"):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as folder:
                mesh, prefix = self._copy(folder)
                cells = Path(str(prefix) + ".cells.csv")
                with cells.open(newline="") as stream:
                    rows = list(csv.DictReader(stream))
                rows[0][field] = str(float(rows[0][field]) + 1e-4)
                with cells.open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
                    writer.writeheader()
                    writer.writerows(rows)
                result = verifier.verify_case(mesh, prefix, "manufactured", .1, 1.0, self._args())
                self.assertFalse(result["valid"])
                self.assertIn("manufactured benchmark checks failed", result["issues"])

    def test_ordinary_case_rejects_manufactured_columns(self):
        with tempfile.TemporaryDirectory() as folder:
            mesh, prefix = self._copy(folder)
            cells = Path(str(prefix) + ".cells.csv")
            with cells.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["sourceX"] = "1e-4"
            with cells.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            parsed_mesh = verifier.read_cm2d(mesh)
            with self.assertRaises(verifier.VerificationError):
                verifier.read_cells(cells, parsed_mesh,
                                    verifier.measure(parsed_mesh, 1e-11, 1e-9), "channel")

    def test_manufactured_pressure_gauge_is_checked_numerically(self):
        mesh = verifier.read_cm2d(self.mesh)
        measured = verifier.measure(mesh, 1e-11, 1e-9)
        rows = verifier.read_cells(Path(str(self.prefix) + ".cells.csv"), mesh,
                                   measured, "manufactured")
        for row in rows:
            row["p"] += .01  # Gradients cannot detect an arbitrary pressure offset.
        result = verifier.manufactured_checks(mesh, measured, rows, .1, 1.0)
        self.assertFalse(result["valid"])
        self.assertFalse(result["pressureGauge"]["valid"])

    def test_manufactured_requires_all_exact_and_source_columns(self):
        with tempfile.TemporaryDirectory() as folder:
            mesh, prefix = self._copy(folder)
            cells = Path(str(prefix) + ".cells.csv")
            with cells.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            rows = [{key: value for key, value in row.items() if key != "sourceY"} for row in rows]
            with cells.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaises(verifier.VerificationError):
                verifier.verify_case(mesh, prefix, "manufactured", .1, 1.0, self._args())

    def test_cli_rejects_non_unit_manufactured_domain(self):
        parsed = verifier.read_cm2d(self.mesh)
        bad_mesh = self.root / "bad-domain.solver.cm2d"
        lines = ["CM2D 1", f"VERTICES {len(parsed.vertices)}"]
        lines.extend(f"{i} {4.0 * x:.17g} {y:.17g}" for i, (x, y) in enumerate(parsed.vertices))
        lines.extend([f"EDGES {len(parsed.edges)}", *(
            " ".join(map(str, (edge.id, edge.v0, edge.v1, edge.owner, edge.neighbour, edge.patch)))
            for edge in parsed.edges),
            f"CELLS {len(parsed.cells)}"])
        lines.extend(" ".join(map(str, (cell.id, 0, 0, 4.0 * cell.stored_area, len(cell.vertices),
                                         *cell.vertices, len(cell.edges), *cell.edges)))
                     for cell in parsed.cells)
        lines.extend(["AUDIT 0 0 0 0 0 0 0", "END"])
        bad_mesh.write_text("\n".join(lines) + "\n", encoding="utf-8")
        result = subprocess.run(
            [str(FLOW_CLI), "--mesh", str(bad_mesh), "--output", str(self.root / "bad"),
             "--case", "manufactured", "--nu", ".1", "--speed", "1", "--max-iterations", "10"],
            capture_output=True, text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
