import os
import unittest
import numpy
from pyKVFinder.grid import (
    get_vertices,
    _get_dimensions,
    _get_sincos,
    _get_cavity_name,
    _get_cavity_label,
    _select_cavities,
)


class TestGetCavityName(unittest.TestCase):
    def test_indexes(self):
        indexes = list(range(0, 100, 10))
        cavity_names = [_get_cavity_name(index) for index in indexes]
        self.assertListEqual(
            cavity_names,
            ["KAA", "KAK", "KAU", "KBE", "KBO", "KBY", "KCI", "KCS", "KDC", "KDM"],
        )


class TestGetCavityLabel(unittest.TestCase):
    def test_cavity_names(self):
        cavity_names = [
            "KAA",
            "KAK",
            "KAU",
            "KBE",
            "KBO",
            "KBY",
            "KCI",
            "KCS",
            "KDC",
            "KDM",
        ]
        labels = [_get_cavity_label(name) for name in cavity_names]
        self.assertListEqual(labels, list(range(2, 100, 10)))


class TestSelectCavities(unittest.TestCase):
    def test_selection(self):
        # Create a dummy grid
        grid = numpy.arange(0, 27).reshape(3, 3, 3)
        # Selection: list of cavity labels
        selected = _select_cavities(grid, selection=[2, 3, 4])
        self.assertListEqual(
            selected.tolist(),
            [
                [[0, 1, 2], [3, 4, 1], [1, 1, 1]],
                [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
                [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
            ],
        )


class TestGetVertices(unittest.TestCase):
    def setUp(self):
        # Define valid xyzr, padding and step
        self.xyzr = [[0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0]]
        self.probe_out = 4.0
        self.step = 0.6

    def test_xyzr(self):
        # Test a xyzrc value without padding and step
        result = get_vertices(self.xyzr, 0.0, self.step).tolist()
        self.assertListEqual(
            result,
            [
                [-0.6, -0.6, -0.6],
                [1.6, -0.6, -0.6],
                [-0.6, 1.6, -0.6],
                [-0.6, -0.6, 1.6],
            ],
        )

    def test_probe_out(self):
        # Test a padding addition
        result = get_vertices(self.xyzr, self.probe_out, self.step).tolist()
        self.assertListEqual(
            result,
            [
                [-4.6, -4.6, -4.6],
                [5.6, -4.6, -4.6],
                [-4.6, 5.6, -4.6],
                [-4.6, -4.6, 5.6],
            ],
        )

    def test_step(self):
        # Test a step addition
        result = get_vertices(self.xyzr, self.probe_out, 0.5).tolist()
        self.assertListEqual(
            result,
            [
                [-4.5, -4.5, -4.5],
                [5.5, -4.5, -4.5],
                [-4.5, 5.5, -4.5],
                [-4.5, -4.5, 5.5],
            ],
        )

    def test_bad_type(self):
        # bad xyzr
        self.assertRaises(TypeError, get_vertices, "string", self.probe_out, self.step)
        self.assertRaises(TypeError, get_vertices, True, self.probe_out, self.step)
        self.assertRaises(TypeError, get_vertices, 1, self.probe_out, self.step)
        self.assertRaises(TypeError, get_vertices, 1.0, self.probe_out, self.step)
        # bad probe_out
        self.assertRaises(TypeError, get_vertices, self.xyzr, "string", self.step)
        self.assertRaises(TypeError, get_vertices, self.xyzr, True, self.step)
        self.assertRaises(TypeError, get_vertices, self.xyzr, [4.0, 4.0], self.step)
        # bad step
        self.assertRaises(TypeError, get_vertices, self.xyzr, self.probe_out, "string")
        self.assertRaises(TypeError, get_vertices, self.xyzr, self.probe_out, True)
        self.assertRaises(
            TypeError, get_vertices, self.xyzr, self.probe_out, [0.6, 0.6]
        )

    def test_bad_values(self):
        # vertices
        # shape (4,)
        self.assertRaises(
            ValueError, get_vertices, [1.0, 1.0, 1.0, 1.0], self.probe_out, self.step
        )
        # shape (1, 5)
        self.assertRaises(
            ValueError,
            get_vertices,
            [[1.0, 1.0, 1.0, 1.0, 1.0]],
            self.probe_out,
            self.step,
        )
        # shape (1, 1, 4)
        self.assertRaises(
            ValueError,
            get_vertices,
            [[[1.0, 1.0, 1.0, 1.0]]],
            self.probe_out,
            self.step,
        )
        # probe_out
        self.assertRaises(ValueError, get_vertices, self.xyzr, -1, self.step)
        # step
        self.assertRaises(ValueError, get_vertices, self.xyzr, self.probe_out, 0.0)
        self.assertRaises(ValueError, get_vertices, self.xyzr, self.probe_out, -1.0)


class TestGetDimensions(unittest.TestCase):
    def setUp(self):
        self.vertices = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]

    def test_bad_type(self):
        # bad vertices
        self.assertRaises(TypeError, _get_dimensions, "string")
        self.assertRaises(TypeError, _get_dimensions, True)
        self.assertRaises(TypeError, _get_dimensions, 1)
        self.assertRaises(TypeError, _get_dimensions, 1.0)
        # bad step
        self.assertRaises(TypeError, _get_dimensions, self.vertices, "string")
        self.assertRaises(TypeError, _get_dimensions, self.vertices, True)
        self.assertRaises(TypeError, _get_dimensions, self.vertices, [1, 1.0])

    def test_bad_values(self):
        # bad vertices
        # shape (3,)
        self.assertRaises(ValueError, _get_dimensions, [1.0, 1.0, 1.0])
        # shape (1, 3)
        self.assertRaises(ValueError, _get_dimensions, [[1.0, 1.0, 1.0]])
        # shape (1, 1, 3)
        self.assertRaises(ValueError, _get_dimensions, [[[1.0, 1.0, 1.0]]])
        # bad step
        self.assertRaises(ValueError, _get_dimensions, self.vertices, 0.0)
        self.assertRaises(ValueError, _get_dimensions, self.vertices, -1)

    def test_vertices(self):
        # Test vertices input
        result = _get_dimensions(self.vertices)
        self.assertEqual(result, (2, 2, 2))


class TestGetSincos(unittest.TestCase):
    def test_bad_type(self):
        # bad vertices
        self.assertRaises(TypeError, _get_sincos, "string")
        self.assertRaises(TypeError, _get_sincos, True)
        self.assertRaises(TypeError, _get_sincos, 1)
        self.assertRaises(TypeError, _get_sincos, 1.0)

    def test_bad_shape(self):
        # shape (3,)
        self.assertRaises(ValueError, _get_sincos, [1.0, 1.0, 1.0])
        # shape (1, 3)
        self.assertRaises(ValueError, _get_sincos, [[1.0, 1.0, 1.0]])
        # shape (1, 1, 3)
        self.assertRaises(ValueError, _get_sincos, [[[1.0, 1.0, 1.0]]])

    def test_vertices(self):
        # Aligned vertices
        data = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        result = _get_sincos(data).tolist()
        self.assertListAlmostEqual(result, [0.0, 1.0, 0.0, 1.0], 3)
        # a = 45o rotation
        data = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.707, -0.707],
            [0.0, 0.707, 0.707],
        ]
        result = _get_sincos(data).tolist()
        self.assertListAlmostEqual(result, [0.707, 0.707, 0.0, 1.0], 3)
        # b = 45o rotation
        data = [
            [0.0, 0.0, 0.0],
            [0.707, 0.0, -0.707],
            [0.0, 1.0, 0.0],
            [0.707, 0.0, 0.707],
        ]
        result = _get_sincos(data).tolist()
        self.assertListAlmostEqual(result, [0.0, 1.0, -0.707, 0.707], 3)

    def assertListAlmostEqual(self, list1, list2, tol):
        self.assertEqual(len(list1), len(list2))
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(a, b, tol)
