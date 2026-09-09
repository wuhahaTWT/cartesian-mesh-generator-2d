"""Analytic and actual-failure tests for the independent determinant reader."""
import math
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "verification"))
from check_directional_connectivity import determinant, measure, verify_native_report


class DirectionalConnectivityTest(unittest.TestCase):
    def test_serialized_four_cell_mesh(self):
        # Four unit squares, each with two orthogonal internal side faces.
        # Physical outer edges are excluded from the directional tensor.
        mesh = """CM2D 1
VERTICES 9
0 0 0
1 1 0
2 2 0
3 0 1
4 1 1
5 2 1
6 0 2
7 1 2
8 2 2
EDGES 12
0 0 1 0 -1 2
1 1 4 0 1 0
2 4 3 0 2 0
3 3 0 0 -1 2
4 1 2 1 -1 2
5 2 5 1 -1 2
6 5 4 1 3 0
7 4 7 2 3 0
8 7 6 2 -1 2
9 6 3 2 -1 2
10 5 8 3 -1 2
11 8 7 3 -1 2
CELLS 4
0 0 0 1 4 0 1 4 3 4 0 1 2 3
1 1 1 1 4 1 2 5 4 4 4 5 6 1
2 2 2 1 4 3 4 7 6 4 2 7 8 9
3 3 3 1 4 4 5 8 7 4 6 10 11 7
AUDIT 0 0 0 0 0 0 0
END
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "squares.cm2d"
            path.write_text(mesh)
            result = measure(path)
            self.assertTrue(result["valid"])
            self.assertEqual(result["cell_count"], 4)
            self.assertEqual(result["minimum_measured"], 0.125)
            self.assertEqual(result["failed_cells"], [])
            for corrupt in (mesh.replace("1 1 4 0 1 0", "1 1 4 0 0 0"),
                            mesh.replace("1 1 4 0 1 0", "1 1 4 0 3 0")):
                path.write_text(corrupt)
                with self.assertRaises(ValueError):
                    measure(path)

    def test_rank_and_square(self):
        self.assertEqual(determinant([]), 0)
        self.assertEqual(determinant([(1, 0)]), 0)
        self.assertEqual(determinant([(1, 0), (-1, 0)]), 0)
        self.assertEqual(determinant([(1, 0), (0, 1)]), 0.125)
        self.assertEqual(determinant([(1, 0), (0, 1), (-1, 0), (0, -1)]), 0.5)

    def test_actual_naca_failure_and_finite_transforms(self):
        # Actual internal side faces of OpenFOAM cell 108446, NACA r03.
        a = (1.000083814, 0.0012572093)
        b = (0.9998148997045184, 0.0008471631535156218)
        c = (1.0020160673496092, 0.0028003700031249967)
        vectors = [(c[0] - b[0], c[1] - b[1]), (c[0] - a[0], c[1] - a[1])]
        expected = 0.0003310395831  # Actual OpenFOAM 2606 checkMesh output.
        for scale in (0.001, 1, 1000):
            for angle in (0, math.radians(17)):
                transformed = [(scale * (math.cos(angle) * x - math.sin(angle) * y),
                                scale * (math.sin(angle) * x + math.cos(angle) * y))
                               for x, y in vectors]
                self.assertAlmostEqual(determinant(transformed), expected, delta=1e-12)

    def test_invalid_vectors_are_not_a_pass(self):
        for vectors in [[(0, 0)], [(math.inf, 1)], [(1, math.nan)]]:
            with self.assertRaises(ValueError):
                determinant(vectors)

    def test_native_failure_report_cannot_hide_bad_cells(self):
        measured = {"valid": False, "minimum_measured": 0.000331,
                    "failed_cells": [{"cell_id": 108446}]}
        native = {"scope": "uncoupled_planar_uniform_extrusion_empty_front_back",
                  "valid": False, "threshold": 0.001, "minimum": 0.000331,
                  "failed_cell_ids": [108446], "input_issue_count": 0}
        verify_native_report(measured, {"directional_connectivity": native})
        for key, wrong in (("valid", True), ("threshold", 0.0001),
                           ("minimum", None), ("minimum", 0.001),
                           ("failed_cell_ids", []), ("input_issue_count", 1)):
            with self.assertRaises(ValueError):
                verify_native_report(measured, {"directional_connectivity": {**native, key: wrong}})


if __name__ == "__main__":
    unittest.main()
