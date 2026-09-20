#!/usr/bin/env python3
"""Focused negative tests for the independent native-flow momentum audit."""

import csv
import decimal
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


class OutletBackflowVerifierTests(unittest.TestCase):
    def setUp(self):
        self.mesh = verifier.Mesh(
            Path("square.solver.cm2d"), ((0., 0.), (1., 0.), (1., 1.), (0., 1.)),
            tuple(verifier.Edge(i, i, (i+1) % 4, 0, -1, 2) for i in range(4)),
            (verifier.Cell(0, 1., (0, 1, 2, 3), (0, 1, 2, 3)),), (0,) * 7)
        self.measured = verifier.measure(self.mesh, 1e-11, 1e-9)
        self.geo = verifier.face_geometry(self.mesh, self.measured)

    def test_final_flux_sign_controls_incoming_outlet_tangent(self):
        for q in (-.1, -1e-20, 0., .1):
            bc = verifier.flow_boundaries(self.mesh, self.measured, "counterflow", 1.,
                                          "normal-inlet", [0., q, 0., -q])
            self.assertFalse(bc["fixedU"][1])
            self.assertTrue(bc["fixedP"][1])
            self.assertEqual(bc["fixedV"][1], q < 0)
            self.assertEqual(bc["constantV"][1], q < 0)
            self.assertEqual(bc["v"][1], 0.)
            self.assertEqual(bc["u"][3], -1.)  # signed left profile at y=.5
            self.assertTrue(bc["fixedV"][0])

    def test_reverse_normal_advection_and_tangential_stress(self):
        bc = verifier.flow_boundaries(self.mesh, self.measured, "external", 1.,
                                      "normal-inlet", [0., -.1, 0., .1])
        edge, geo = self.mesh.edges[1], self.geo[1]
        gu, gv = [(9., 2.)], [(7., 8.)]
        self.assertEqual(verifier._viscous_face_gradient(self.mesh, self.measured, edge, geo,
            [-.2], gu, bc["u"], bc["fixedU"], bc["constantU"]), (0., 2.))
        self.assertEqual(verifier._viscous_face_gradient(self.mesh, self.measured, edge, geo,
            [.3], gv, bc["v"], bc["fixedV"], bc["constantV"]), (-.6, 0.))
        for limiter in (None, [1.]):
            self.assertEqual(verifier._advective_value(self.mesh, self.measured, edge, geo,
                -.1, [-.2], gu, limiter, bc["fixedU"], bc["u"], normal_inlet=True), -.2)
            self.assertEqual(verifier._advective_value(self.mesh, self.measured, edge, geo,
                -.1, [.3], gv, limiter, bc["fixedV"], bc["v"], normal_inlet=True), 0.)

    def test_audit_uses_actual_flux_and_rejects_default_or_unknown_reverse_policy(self):
        cells = [{"u": -.2, "v": .3, "p": 0.}]
        records = [dict(flux=q, wall=0., pressure=0., advectionX=0.,
                        advectionY=0., diffusionX=0., diffusionY=0.)
                   for q in (0., -.1, 0., .1)]
        payload = dict(outletBackflow="normal-inlet", convection="limited-linear",
                       viscousStress="symmetric", momentumResidual=0.,
                       pressureForceX=0., pressureForceY=0., discreteForceX=0.,
                       discreteForceY=0., forceX=0., forceY=0.)
        def audit():
            return verifier.reconstruct_momentum_audit(self.mesh, self.measured, cells,
                records, .1, 1., "external", payload)
        result = audit()
        self.assertEqual(result["outletBackflowAudit"],
                         dict(faceCount=1, inwardVolumeFlux=.1, faceIds=[1]))
        self.assertEqual(result["maxFaceDeviation"]["advectionY"], 0.)
        payload.update(outletBackflowFaces=1, outletInflow=.1)
        self.assertEqual(audit()["outletBackflowAudit"], result["outletBackflowAudit"])
        for key, bad in (("outletBackflowFaces", 0), ("outletInflow", .2)):
            old = payload[key]
            payload[key] = bad
            with self.assertRaisesRegex(verifier.VerificationError, key):
                audit()
            payload[key] = old
        records[1]["advectionY"] = -.03
        self.assertEqual(audit()["maxFaceDeviation"]["advectionY"], .03)
        for mode in ("reject", "unknown"):
            payload["outletBackflow"] = mode
            with self.assertRaisesRegex(verifier.VerificationError, "backflow|outletBackflow"):
                audit()
        self.assertEqual(verifier.outlet_backflow_mode({}), "reject")

    def test_counterflow_exact_source_formula(self):
        # Values at extrema independently establish sign and amplitude.
        nu, speed = .2, 3.
        for y, factor in ((0., 1.), (.5, -1.), (1., 1.)):
            sample = verifier.counterflow_sample(y, speed, nu)
            self.assertEqual(sample["u"], speed*(1+2*factor))
            self.assertEqual(sample["v"], 0.)
            self.assertEqual(sample["p"], 0.)
            self.assertAlmostEqual(sample["sourceX"], factor*8*math.pi**2*nu*speed)
            self.assertEqual(sample["sourceY"], 0.)

    def test_normal_inlet_is_numerically_inert_for_closed_cavity(self):
        flux = [0.]*4
        default = verifier.flow_boundaries(self.mesh, self.measured, "cavity", 1.)
        enabled = verifier.flow_boundaries(self.mesh, self.measured, "cavity", 1., "normal-inlet", flux)
        self.assertEqual(default, enabled)
        cells = [{"u": .2, "v": .3, "p": 0.}]
        records = [dict(flux=0., wall=1., pressure=0., advectionX=0.,
                        advectionY=0., diffusionX=0., diffusionY=0.) for _ in flux]
        payload = dict(convection="limited-linear", viscousStress="symmetric", momentumResidual=0.,
                       pressureForceX=0., pressureForceY=0., discreteForceX=0.,
                       discreteForceY=0., forceX=0., forceY=0.)
        results = []
        for mode in ("reject", "normal-inlet"):
            payload["outletBackflow"] = mode
            result = verifier.reconstruct_momentum_audit(self.mesh, self.measured, cells,
                records, .1, 1., "cavity", payload)
            result.pop("outletBackflow")
            results.append(result)
        self.assertEqual(*results)

    def test_re20_reference_is_not_compared_to_other_reynolds_numbers(self):
        args = SimpleNamespace(external_lift_drag_ratio=.3, max_speed_ratio=4.)
        cells = [dict(x=1.1, y=0., u=1., speed=1.)]
        payload = dict(forceX=1.0225, forceY=0.)
        for nu, matches in ((.05, True), (.025, False)):
            result = verifier.external_checks(tiny_mesh(), self.measured, cells, payload, nu, 1., args)
            reference = result["openCylinderReference"]
            self.assertEqual(reference["reynoldsMatches"], matches)
            self.assertFalse(reference["acceptanceGate"])
            if matches:
                self.assertEqual(reference["relativeDifference"], 0.)
            else:
                self.assertIsNone(reference["relativeDifference"])


class PressureBoundaryStencilTests(unittest.TestCase):
    def setUp(self):
        # Gradient-only fixture: wall-adjacent cell 0 has two almost parallel
        # interior neighbour directions. Their shared second ring spans y.
        centres = ((0., 0.), (-1., .01), (-2., -.01), (-1., 1.), (-1., -1.))
        pairs = ((0, -1), (0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 4))
        edges = tuple(verifier.Edge(i, 0, 1, owner, other, int(other < 0))
                      for i, (owner, other) in enumerate(pairs))
        cells = tuple(verifier.Cell(i, 1., (), tuple(e.id for e in edges
                       if e.owner == i or e.neighbour == i)) for i in range(5))
        self.mesh = verifier.Mesh(Path('pressure-stencil.cm2d'), (), edges, cells, (0,) * 7)
        self.measured = verifier.Measurement((1.,)*5, centres, (-2., -1., 0., 1.), 5., 1., ())
        self.geo = [verifier.FaceGeometry((0., .5), (0., 1.), (0., 0.), 1., .5)
                    for _ in edges]
        self.fixed = [False] * len(edges)
        self.bc = [0.] * len(edges)
        self.values = [2. + 3.*x - 4.*y for x, y in centres]

    def gradients(self, values=None, skip=True, fixed=None, second_ring=True):
        return verifier.reconstruct_gradient(self.mesh, self.measured, self.geo,
            self.values if values is None else values, self.bc,
            self.fixed if fixed is None else fixed, skip_unknown_boundary=skip,
            boundary_second_ring=second_ring)

    def test_affine_pressure_survives_nearly_collinear_boundary_stencil(self):
        for gx, gy in self.gradients():
            self.assertAlmostEqual(gx, 3., places=11)
            self.assertAlmostEqual(gy, -4., places=11)

    def test_second_ring_reduces_boundary_noise_amplification(self):
        values = self.values.copy()
        values[1] += 1e-4
        values[2] -= 1e-4
        gx, gy = self.gradients(values)[0]
        self.assertLess(math.hypot(gx - 3., gy + 4.), 1e-4)
        # Two direct samples alone solve -gx+.01gy=delta1 and
        # -2gx-.01gy=delta2, amplifying this noise into a .01 y error.
        d1, d2 = values[1] - values[0], values[2] - values[0]
        direct_gy = (2.*d1 - d2)/.03
        direct_gx = .01*direct_gy - d1
        direct_error = math.hypot(direct_gx - 3., direct_gy + 4.)
        self.assertGreater(direct_error, .009)
        self.assertLess(math.hypot(gx - 3., gy + 4.), direct_error/100.)
        old_gx, old_gy = self.gradients(values, second_ring=False)[0]
        self.assertAlmostEqual(old_gx, direct_gx, places=10)
        self.assertAlmostEqual(old_gy, direct_gy, places=10)

    def test_second_ring_is_used_only_when_pressure_boundary_is_omitted(self):
        changed = self.values.copy()
        changed[3] += 1.
        self.assertNotEqual(self.gradients(changed)[0], self.gradients()[0])
        # Velocity / legacy zero-normal reconstruction keeps the direct ring.
        self.assertEqual(self.gradients(changed, skip=False)[0],
                         self.gradients(skip=False)[0])
        # Prescribed pressure retains its boundary sample and direct ring.
        fixed = self.fixed.copy()
        fixed[0] = True
        self.assertEqual(self.gradients(changed, fixed=fixed)[0],
                         self.gradients(fixed=fixed)[0])
        # Interior full-rank cell 1 must not sample second-ring cell 2.
        changed = self.values.copy()
        changed[2] += 1.
        self.assertEqual(self.gradients(changed)[1], self.gradients()[1])


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


class PolygonMeasurementTests(unittest.TestCase):
    @staticmethod
    def high_precision(points):
        # Reference uses the global-coordinate formula in 80-digit arithmetic
        # on the exact binary input values, not native/exported cell metadata.
        with decimal.localcontext() as context:
            context.prec = 80
            points = [tuple(decimal.Decimal.from_float(x) for x in p) for p in points]
            pairs = list(zip(points, points[1:] + points[:1]))
            cross = [x*v-u*y for ((x, y), (u, v)) in pairs]
            twice = sum(cross)
            centre = tuple(float(sum((a[k]+b[k])*c for (a, b), c in zip(pairs, cross))
                                 / (3*twice)) for k in (0, 1))
            return float(twice/2), centre

    def assert_precise(self, points):
        area, centre = verifier.polygon(points)
        expected_area, expected_centre = self.high_precision(points)
        self.assertLessEqual(abs(area-expected_area), 4*math.ulp(expected_area))
        for actual, expected in zip(centre, expected_centre):
            self.assertLessEqual(abs(actual-expected), 4*math.ulp(expected))

    def test_translated_small_polygons_match_high_precision(self):
        for shape in (((0., 0.), (1., 0.), (.75, 1.), (0., .8)),
                      ((0., 0.), (1., 0.), (1., .2), (.2, .2), (.2, 1.), (0., 1.))):
            for size, offset in ((1., (0., 0.)), (1e-6, (.6, -.8)),
                                 (1e-4, (1e6, -1e6)), (1e3, (-1e9, 2e9))):
                with self.subTest(shape=shape, size=size, offset=offset):
                    self.assert_precise([(offset[0]+size*x, offset[1]+size*y) for x, y in shape])

    def test_actual_small_cutcell_rejects_global_origin_cancellation(self):
        points = [(0.6325000000000001, -0.81375),
                  (0.6325000000000001, -0.7957824995171311),
                  (0.62557023302, -0.801469612303)]
        self.assert_precise(points)
        pairs = list(zip(points, points[1:]+points[:1]))
        cross = [a[0]*b[1]-b[0]*a[1] for a, b in pairs]
        old_area = .5*math.fsum(cross)
        exact_area, _ = self.high_precision(points)
        self.assertGreater(abs(old_area-exact_area), 100*math.ulp(exact_area))

    def test_thin_centres_do_not_shift_tangentially(self):
        height=9.531188624900715e-5
        for x in (.3125,.9375,1024.3125):
            points=[(x,0.),(x+.0625,0.),(x+.0625,height),(x,height)]
            for start in range(4):
                _,centre=verifier.polygon(points[start:]+points[:start])
                self.assertEqual(centre,(x+.03125,height/2))
        # The independent high-precision path also handles an asymmetric thin
        # polygon; its centre must not be replaced by a box midpoint.
        self.assert_precise([(.3125,0.),(.375,0.),(.32,height)])

    def test_actual_centroid_halfway_below_old_precision_trigger(self):
        # Actual cell 4001, aspect ~29.36: old aspect>32 trigger missed it.
        # The exact x moment is halfway between binary64 values. Expected
        # tie-to-even result is from exact input fractions, not exported metadata.
        x0=float.fromhex('-0x1.f99999999999ap-2');x1=float.fromhex('-0x1.f333333333333p-2')
        y0=float.fromhex('0x1.07a798a899018p-11');y1=float.fromhex('0x1.7767a6cfbae7fp-11')
        points=[(x0,y0),(x1,y0),(x1,y1),(x0,y1)]
        self.assertLess((x1-x0)/(y1-y0),32)
        for start in range(4):
            _,centre=verifier.polygon(points[start:]+points[:start])
            self.assertEqual(centre[0].hex(),'-0x1.f666666666666p-2')
            self.assertEqual(centre[1].hex(),'0x1.3f879fbc29f4cp-11')

    def test_invalid_polygons_still_fail(self):
        for points in ([], [(1., 1.)], [(0., 0.), (1., 0.), (2., 0.)],
                       [(0., 0.), (0., 1.), (1., 0.)]):
            with self.assertRaises(verifier.VerificationError):
                verifier.polygon(points)


class Cm2dReaderTests(unittest.TestCase):
    @staticmethod
    def _mesh_text(cell="0 0 0 0.5 3 0 1 2 3 0 1 2"):
        return "\n".join((
            "CM2D 1 VERTICES 3",
            "0 0 0",
            "1 1 0",
            "2 0 1",
            "EDGES 3",
            "0 0 1 0 -1 1",
            "1 1 2 0 -1 1",
            "2 2 0 0 -1 1",
            "CELLS 1",
            cell,
            "AUDIT 0 0 0 0 0 0 0 END",
            ""))

    def _write(self, folder, text):
        path = Path(folder) / "mesh.solver.cm2d"
        path.write_text(text, encoding="utf-8-sig")
        return path

    def test_stream_reader_accepts_wrapped_whitespace_and_bom(self):
        with tempfile.TemporaryDirectory() as folder:
            path = self._write(folder, self._mesh_text().replace(" ", "\t", 4))
            mesh = verifier.read_cm2d(path)
            self.assertEqual(len(mesh.vertices), 3)
            self.assertEqual(len(mesh.cells), 1)

    def test_stream_reader_rejects_truncation(self):
        with tempfile.TemporaryDirectory() as folder:
            path = self._write(folder, self._mesh_text().rsplit(" END", 1)[0])
            with self.assertRaisesRegex(verifier.VerificationError, "truncated CM2D record"):
                verifier.read_cm2d(path)

    def test_stream_reader_rejects_trailing_tokens(self):
        with tempfile.TemporaryDirectory() as folder:
            path = self._write(folder, self._mesh_text() + "extra\n")
            with self.assertRaisesRegex(verifier.VerificationError, "trailing or missing"):
                verifier.read_cm2d(path)

    def test_stream_reader_preserves_id_and_arity_checks(self):
        with tempfile.TemporaryDirectory() as folder:
            bad_id = self._mesh_text().replace("1 1 0\n", "7 1 0\n")
            with self.assertRaisesRegex(verifier.VerificationError, "non-contiguous vertex ids"):
                verifier.read_cm2d(self._write(folder, bad_id))
            bad_arity = self._mesh_text(cell="0 0 0 0.5 3 0 1 2 2 0 1")
            with self.assertRaisesRegex(verifier.VerificationError, "invalid cell 0 loop"):
                verifier.read_cm2d(self._write(folder, bad_arity))


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
        self.assertIn(result["pressureBoundaryReconstruction"], ("zero-normal", "one-sided-linear", "one-sided-linear-2ring", "one-sided-linear-adaptive"))
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
            self.assertIn(payload.get("pressureBoundaryReconstruction"), ("zero-normal", "one-sided-linear", "one-sided-linear-2ring", "one-sided-linear-adaptive"))
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

    def test_counterflow_native_fixture_and_tampering(self):
        prefix = self.root / "counterflow"
        completed = subprocess.run(
            [str(FLOW_CLI), "--mesh", str(self.mesh), "--output", str(prefix),
             "--case", "counterflow", "--nu", ".1", "--speed", "1",
             "--convection", "limited-linear", "--outlet-backflow", "normal-inlet",
             "--max-iterations", "1500"], capture_output=True, text=True, timeout=90)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        result = verifier.verify_case(self.mesh, prefix, "counterflow", .1, 1., self._args())
        self.assertTrue(result["valid"], result["issues"])
        self.assertGreater(result["benchmark"]["outlet"]["incomingFaces"], 0)
        self.assertGreater(result["benchmark"]["outlet"]["outgoingFaces"], 0)
        for suffix, column in ((".cells.csv", "sourceX"), (".cells.csv", "exactU"),
                               (".faces.csv", "advectionY"), (".faces.csv", "flux")):
            with self.subTest(column=column), tempfile.TemporaryDirectory() as folder:
                copied = Path(folder) / "counterflow"
                for ending in (".cells.csv", ".faces.csv", ".json", ".residuals.csv"):
                    shutil.copyfile(Path(str(prefix)+ending), Path(str(copied)+ending))
                path = Path(str(copied)+suffix)
                with path.open(newline="") as stream:
                    rows = list(csv.DictReader(stream))
                idx = (result["momentumAudit"]["outletBackflowAudit"]["faceIds"][0]
                       if suffix == ".faces.csv" else 0)
                rows[idx][column] = str(float(rows[idx][column]) + .01)
                with path.open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
                    writer.writeheader(); writer.writerows(rows)
                tampered = verifier.verify_case(self.mesh, copied, "counterflow", .1, 1., self._args())
                self.assertFalse(tampered["valid"])
        path = Path(str(prefix)+".json")
        payload = json.loads(path.read_text())
        for mode in ("reject", None, "unknown"):
            with self.subTest(mode=mode):
                if mode is None:
                    payload.pop("outletBackflow", None)
                else:
                    payload["outletBackflow"] = mode
                path.write_text(json.dumps(payload))
                if mode == "unknown":
                    with self.assertRaisesRegex(verifier.VerificationError, "outletBackflow"):
                        verifier.verify_case(self.mesh, prefix, "counterflow", .1, 1., self._args())
                else:
                    self.assertFalse(verifier.verify_case(self.mesh, prefix, "counterflow", .1, 1., self._args())["valid"])


if __name__ == "__main__":
    unittest.main()
