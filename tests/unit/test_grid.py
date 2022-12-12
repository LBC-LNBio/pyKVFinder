from typing import Union
import unittest
import os
import toml
import numpy
import pyKVFinder
from pyKVFinder.grid import (
    get_vertices,
    _get_vertices_from_box,
    _get_vertices_from_residues,
    get_vertices_from_file,
    _get_dimensions,
    _get_sincos,
    _get_cavity_name,
    _get_cavity_label,
    _select_cavities,
    _process_spatial,
    _process_depth,
    _process_residues,
    _process_hydropathy,
    detect,
    spatial,
    depth,
    constitutional,
    hydropathy,
)
from pyKVFinder.utils import read_pdb, read_cavity

PYKVFINDER_TESTS_DIR = os.path.join(
    os.path.dirname(pyKVFinder.__file__), "data", "tests"
)
UNIT_TESTS_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


class TestGetVerticesFromBox(unittest.TestCase):
    def test_custom_box(self):
        # Expected
        expected = [
            [-0.89, 3.34, -2.41],
            [15.51, 3.34, -2.41],
            [-0.89, 14.74, -2.41],
            [-0.89, 3.34, 10.19],
        ]
        # Prepare data
        box = toml.load(os.path.join(PYKVFINDER_TESTS_DIR, "custom-box.toml"))["box"]
        # Result
        result = _get_vertices_from_box(box).tolist()
        self.assertListEqual(expected, result)


class TestGetVerticesFromResidues(unittest.TestCase):
    def test_residues_box(self):
        # Expected
        expected = [
            [-2.634, 2.843, -3.68],
            [24.422, 2.843, -3.68],
            [-2.634, 22.698, -3.68],
            [-2.634, 2.843, 16.204],
        ]
        # Prepare data
        atomic = read_pdb(os.path.join(PYKVFINDER_TESTS_DIR, "1FMO.pdb"))
        atominfo = numpy.asarray(
            ([[f"{atom[0]}_{atom[1]}_{atom[2]}", atom[3]] for atom in atomic[:, :4]])
        )
        xyzr = atomic[:, 4:].astype(numpy.float64)
        box = toml.load(os.path.join(PYKVFINDER_TESTS_DIR, "residues-box.toml"))["box"]
        # Result
        result = _get_vertices_from_residues(box, atominfo, xyzr).tolist()
        self.assertListEqual(expected, result)


class TestGetVerticesFromFile(unittest.TestCase):
    def test_custom_box(self):
        # Expected vertices
        expected = [
            [-0.89, 3.34, -2.41],
            [15.51, 3.34, -2.41],
            [-0.89, 14.74, -2.41],
            [-0.89, 3.34, 10.19],
        ]
        # Dummy atomic
        atomic = [
            ["13", "E", "GLU", "C", "3.5", "8.0", "4.0", "1.908"],
            ["13", "E", "GLU", "C", "-6.73", "-14.62", "-15.897", "1.908"],
        ]
        # Get vertices from file
        vertices, selected = get_vertices_from_file(
            os.path.join(PYKVFINDER_TESTS_DIR, "custom-box.toml"), atomic
        )
        # Atom selection
        self.assertListEqual(selected.tolist(), [atomic[0]])
        # Vertices
        self.assertListEqual(vertices.tolist(), expected)

    def test_residues_box(self):
        # Expected vertices
        expected = [
            [-4.0, 0.5, -3.5],
            [11.0, 0.5, -3.5],
            [-4.0, 15.5, -3.5],
            [-4.0, 0.5, 11.5],
        ]
        # Dummy atomic
        atomic = [
            ["49", "E", "LEU", "C", "3.5", "8.0", "4.0", "1.908"],
            ["50", "E", "GLY", "C", "3.5", "8.0", "4.0", "1.908"],
            ["51", "E", "THR", "C", "3.5", "8.0", "4.0", "1.908"],
            ["52", "E", "GLU", "C", "-6.73", "-14.62", "-15.897", "1.908"],
        ]
        # Get vertices from file
        vertices, selected = get_vertices_from_file(
            os.path.join(PYKVFINDER_TESTS_DIR, "residues-box.toml"), atomic
        )
        # Atom selection
        self.assertListEqual(selected.tolist(), atomic[0:3])
        # Vertices
        self.assertListEqual(vertices.tolist(), expected)


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
        # Define valid atomic, padding and step
        self.atomic = [
            ["1", "A", "GLU", "C", "0.0", "0.0", "0.0", "0.0"],
            ["1", "A", "GLU", "C", "1.0", "1.0", "1.0", "1.0"],
        ]
        self.probe_out = 4.0
        self.step = 0.6

    def test_atomic(self):
        # Test a atomic value without padding and step
        result = get_vertices(self.atomic, 0.0, self.step).tolist()
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
        result = get_vertices(self.atomic, self.probe_out, self.step).tolist()
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
        result = get_vertices(self.atomic, self.probe_out, 0.5).tolist()
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
        # bad atomic
        self.assertRaises(TypeError, get_vertices, "string", self.probe_out, self.step)
        self.assertRaises(TypeError, get_vertices, True, self.probe_out, self.step)
        self.assertRaises(TypeError, get_vertices, 1, self.probe_out, self.step)
        self.assertRaises(TypeError, get_vertices, 1.0, self.probe_out, self.step)
        # bad probe_out
        self.assertRaises(TypeError, get_vertices, self.atomic, "string", self.step)
        self.assertRaises(TypeError, get_vertices, self.atomic, True, self.step)
        self.assertRaises(TypeError, get_vertices, self.atomic, [4.0, 4.0], self.step)
        # bad step
        self.assertRaises(
            TypeError, get_vertices, self.atomic, self.probe_out, "string"
        )
        self.assertRaises(TypeError, get_vertices, self.atomic, self.probe_out, True)
        self.assertRaises(
            TypeError, get_vertices, self.atomic, self.probe_out, [0.6, 0.6]
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
        self.assertRaises(ValueError, get_vertices, self.atomic, -1, self.step)
        # step
        self.assertRaises(ValueError, get_vertices, self.atomic, self.probe_out, 0.0)
        self.assertRaises(ValueError, get_vertices, self.atomic, self.probe_out, -1.0)


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


class TestProcessSpatial(unittest.TestCase):
    def test_raw_data(self):
        # Dummy volume and area
        volume = numpy.array([1.0, 2.0, 3.0])
        area = numpy.array([1.0, 2.0, 3.0])
        # Process volume and area
        volume, area = _process_spatial(volume, area, len(volume), None)
        # Assert results
        self.assertDictEqual(volume, {"KAA": 1.0, "KAB": 2.0, "KAC": 3.0})
        self.assertDictEqual(area, {"KAA": 1.0, "KAB": 2.0, "KAC": 3.0})

    def test_selection(self):
        # Dummy volume and area
        volume = numpy.array([1.0, 2.0, 3.0])
        area = numpy.array([1.0, 2.0, 3.0])
        # Process volume and area
        volume, area = _process_spatial(volume, area, len(volume), selection=[2])
        # Assert results
        self.assertDictEqual(volume, {"KAA": 1.0})
        self.assertDictEqual(area, {"KAA": 1.0})


class TestProcessSpatial(unittest.TestCase):
    def test_raw_data(self):
        # Dummy volume and area
        volume = numpy.array([1.0, 2.0, 3.0])
        area = numpy.array([1.0, 2.0, 3.0])
        # Process volume and area
        volume, area = _process_spatial(volume, area, len(volume), None)
        # Assert results
        self.assertDictEqual(volume, {"KAA": 1.0, "KAB": 2.0, "KAC": 3.0})
        self.assertDictEqual(area, {"KAA": 1.0, "KAB": 2.0, "KAC": 3.0})

    def test_selection(self):
        # Dummy volume and area
        volume = numpy.array([1.0, 2.0, 3.0])
        area = numpy.array([1.0, 2.0, 3.0])
        # Process volume and area
        volume, area = _process_spatial(volume, area, len(volume), selection=[2])
        # Assert results
        self.assertDictEqual(volume, {"KAA": 1.0})
        self.assertDictEqual(area, {"KAA": 1.0})


class TestProcessDepth(unittest.TestCase):
    def test_raw_data(self):
        # Dummy max_depth and avg_depth
        max_depth = numpy.array([1.0, 2.0, 3.0])
        avg_depth = numpy.array([1.0, 2.0, 3.0])
        # Process max_depth and avg_depth
        max_depth, avg_depth = _process_depth(
            max_depth, avg_depth, len(max_depth), None
        )
        # Assert results
        self.assertDictEqual(max_depth, {"KAA": 1.0, "KAB": 2.0, "KAC": 3.0})
        self.assertDictEqual(avg_depth, {"KAA": 1.0, "KAB": 2.0, "KAC": 3.0})

    def test_selection(self):
        # Dummy max_depth and avg_depth
        max_depth = numpy.array([1.0, 2.0, 3.0])
        avg_depth = numpy.array([1.0, 2.0, 3.0])
        # Process max_depth and avg_depth
        max_depth, avg_depth = _process_depth(
            max_depth, avg_depth, len(max_depth), selection=[2]
        )
        # Assert results
        self.assertDictEqual(max_depth, {"KAA": 1.0})
        self.assertDictEqual(avg_depth, {"KAA": 1.0})


class TestProcessResidues(unittest.TestCase):
    def test_raw_data(self):
        # Dummy residues
        residues = [
            "14_A_SER",
            "-1",
            "15_A_GLU",
            "-1",
            "16_A_ALA",
            "-1",
        ]
        # Process residues
        residues = _process_residues(residues, 3, None)
        # Assert results
        self.assertDictEqual(
            residues,
            {
                "KAA": [["14", "A", "SER"]],
                "KAB": [["15", "A", "GLU"]],
                "KAC": [["16", "A", "ALA"]],
            },
        )

    def test_selection(self):
        # Dummy residues
        residues = [
            "14_A_SER",
            "-1",
            "-1",
            "-1",
        ]
        # Process residues
        residues = _process_residues(residues, 3, selection=[2])
        # Assert results
        self.assertDictEqual(residues, {"KAA": [["14", "A", "SER"]]})


class TestProcessHydropathy(unittest.TestCase):
    def test_raw_data(self):
        # Dummy avg_hydropathy
        avg_hydropathy = numpy.array([1.0, 2.0, 3.0])
        # Process avg_hydropathy
        avg_hydropathy = _process_hydropathy(avg_hydropathy, len(avg_hydropathy), None)
        # Assert results
        self.assertDictEqual(avg_hydropathy, {"KAA": 1.0, "KAB": 2.0, "KAC": 3.0})

    def test_selection(self):
        # Dummy max_depth and avg_depth
        # Dummy avg_hydropathy
        avg_hydropathy = numpy.array([1.0, 2.0, 3.0])
        # Process avg_hydropathy
        avg_hydropathy = _process_hydropathy(
            avg_hydropathy, len(avg_hydropathy), selection=[2]
        )
        # Assert results
        self.assertDictEqual(avg_hydropathy, {"KAA": 1.0})


class TestDetect(unittest.TestCase):
    def test_detect(self):
        # Prepare data
        atomic = read_pdb(os.path.join(PYKVFINDER_TESTS_DIR, "1FMO.pdb"))
        vertices = get_vertices(atomic)
        # Detect cavities
        ncav, cavities = detect(atomic, vertices)
        cavities[cavities == 1] = 0
        # Expected grid
        expected = read_cavity(
            os.path.join(PYKVFINDER_TESTS_DIR, "1FMO.KVFinder.output.pdb"),
            os.path.join(PYKVFINDER_TESTS_DIR, "1FMO.pdb"),
        )
        # Assert number of cavities
        self.assertEqual(ncav, int(expected.max() - 1))
        # Grid similarity
        tol = 100
        self.assertEqual((cavities - expected).sum() < tol, True)


class TestSpatial(unittest.TestCase):
    def test_spatial(self):
        # Prepare data
        cavities = numpy.full((5, 5, 5), -1)
        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    cavities[i, j, k] = 0

        for i in range(2, 4):
            for j in range(2, 4):
                for k in range(2, 4):
                    cavities[i, j, k] = 2

        # Surface
        surface, volume, area = spatial(cavities)
        # Assert
        self.assertDictEqual(volume, {"KAA": 1.73})
        self.assertDictEqual(area, {"KAA": 3.13})
        cavities[3, 3, 3] = -1
        self.assertListEqual(cavities.tolist(), surface.tolist())


class TestDepth(unittest.TestCase):
    def test_depth(self):
        # Expected
        expected = numpy.zeros((5, 5, 5))
        expected[2, 2, 2] = 0.6
        # Prepare data
        cavities = numpy.full((5, 5, 5), -1)
        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    cavities[i, j, k] = 0

        cavities[2, 2, 2] = 2
        cavities[2, 2, 1] = 2
        # Surface
        depths, max_depth, avg_depth = depth(cavities)
        # Assert
        self.assertDictEqual(max_depth, {"KAA": 0.6})
        self.assertDictEqual(avg_depth, {"KAA": 0.3})
        self.assertListEqual(depths.tolist(), expected.tolist())


class TestConstitutional(unittest.TestCase):
    def test_constitutional(self):
        # Cavities
        cavities = numpy.full((5, 5, 5), -1)
        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    cavities[i, j, k] = 0

        cavities[2, 2, 1] = 2
        cavities[2, 2, 4] = 3
        # Dummy atomic
        atomic = [
            ["13", "E", "GLU", "C", "1.0", "1.0", "1.0", "1.908"],
        ]
        # Dummy vertices
        vertices = numpy.array(
            [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
        )
        # Surface
        residues = constitutional(cavities, atomic, vertices)
        # Assert
        self.assertDictEqual(
            residues, {"KAA": [["13", "E", "GLU"]], "KAB": [["13", "E", "GLU"]]}
        )


class TestConstitutional(unittest.TestCase):
    def test_constitutional(self):
        # Expected
        expected = numpy.zeros((5, 5, 5))
        expected[2, 2, 1] = 0.76
        expected[2, 2, 4] = 0.76
        # Surface
        surface = numpy.full((5, 5, 5), -1)
        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    surface[i, j, k] = 0

        surface[2, 2, 1] = 2
        surface[2, 2, 4] = 3
        # Dummy atomic
        atomic = [
            ["13", "E", "GLU", "C", "1.0", "1.0", "1.0", "1.908"],
        ]
        # Dummy vertices
        vertices = numpy.array(
            [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
        )
        # Surface
        scales, avg_hydropathy = hydropathy(surface, atomic, vertices)
        # Assert
        self.assertListEqual(scales.tolist(), expected.tolist())
        self.assertDictEqual(
            avg_hydropathy, {"KAA": 0.38, "KAB": 0.38, "EisenbergWeiss": [-1.42, 2.6]}
        )
