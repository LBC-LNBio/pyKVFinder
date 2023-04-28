import os
import unittest
from unittest import mock

import numpy
import toml

import pyKVFinder
from pyKVFinder.grid import (
    _get_cavity_label,
    _get_cavity_name,
    _get_dimensions,
    _get_opening_name,
    _get_opening_label,
    _get_sincos,
    _get_vertices_from_box,
    _get_vertices_from_residues,
    _process_depth,
    _process_hydropathy,
    _process_residues,
    _process_spatial,
    _select_cavities,
    constitutional,
    depth,
    detect,
    export,
    get_vertices,
    get_vertices_from_file,
    hydropathy,
    spatial,
)
from pyKVFinder.utils import read_cavity, read_pdb

DATADIR = os.path.join(os.path.dirname(pyKVFinder.__file__), "data")
FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")


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
        box = toml.load(os.path.join(DATADIR, "tests", "custom-box.toml"))["box"]
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
        atomic = read_pdb(os.path.join(DATADIR, "tests", "1FMO.pdb"))
        atominfo = numpy.asarray(
            ([[f"{atom[0]}_{atom[1]}_{atom[2]}", atom[3]] for atom in atomic[:, :4]])
        )
        xyzr = atomic[:, 4:].astype(numpy.float64)
        box = toml.load(os.path.join(DATADIR, "tests", "residues-box.toml"))["box"]
        # Result
        result = _get_vertices_from_residues(box, atominfo, xyzr).tolist()
        self.assertListEqual(expected, result)


class TestGetVerticesFromFile(unittest.TestCase):
    def setUp(self):
        self.atomic1 = [
            ["13", "E", "GLU", "C", "3.5", "8.0", "4.0", "1.908"],
            ["13", "E", "GLU", "C", "-6.73", "-14.62", "-15.897", "1.908"],
        ]
        self.atomic2 = [
            ["49", "E", "LEU", "C", "3.5", "8.0", "4.0", "1.908"],
            ["50", "E", "GLY", "C", "3.5", "8.0", "4.0", "1.908"],
            ["51", "E", "THR", "C", "3.5", "8.0", "4.0", "1.908"],
            ["52", "E", "GLU", "C", "-6.73", "-14.62", "-15.897", "1.908"],
        ]

    def test_custom_box(self):
        # Expected vertices
        expected = [
            [-0.89, 3.34, -2.41],
            [15.51, 3.34, -2.41],
            [-0.89, 14.74, -2.41],
            [-0.89, 3.34, 10.19],
        ]
        # Get vertices from file
        vertices, selected = get_vertices_from_file(
            os.path.join(DATADIR, "tests", "custom-box.toml"), self.atomic1
        )
        # Atom selection
        self.assertListEqual(selected.tolist(), [self.atomic1[0]])
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
        # Get vertices from file
        vertices, selected = get_vertices_from_file(
            os.path.join(DATADIR, "tests", "residues-box.toml"), self.atomic2
        )
        # Atom selection
        self.assertListEqual(selected.tolist(), self.atomic2[0:3])
        # Vertices
        self.assertListEqual(vertices.tolist(), expected)

    def test_wrong_fn_format(self):
        # Check wrong fn formats
        for fn in [1, 1.0, [1], {"vdw": 1}, numpy.ones(1)]:
            self.assertRaises(TypeError, get_vertices_from_file, fn, self.atomic1)

    def test_wrong_atomic_format(self):
        # Check wrong atomic format
        for atomic in [True, 4, 4.0, {"atomic": []}, "4.0"]:
            self.assertRaises(
                TypeError,
                get_vertices_from_file,
                os.path.join(DATADIR, "tests", "custom-box.toml"),
                atomic,
            )

    def test_invalid_atomic(self):
        # Check invalid atomic
        for atomic in [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],  # shape (8,)
            [[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],  # shape (1, 9)
            [[[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]],  # shape (1, 1, 8)
        ]:
            self.assertRaises(
                ValueError,
                get_vertices_from_file,
                os.path.join(DATADIR, "tests", "custom-box.toml"),
                atomic,
            )

    def test_wrong_step_format(self):
        # Check wrong step format
        for step in [True, [0.6], {"step": 0.6}, "0.6", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                get_vertices_from_file,
                os.path.join(DATADIR, "tests", "custom-box.toml"),
                self.atomic1,
                step=step,
            )

    def test_invalid_step(self):
        # Check invalid step
        for step in [-1.0, 0.0]:
            self.assertRaises(
                ValueError,
                get_vertices_from_file,
                os.path.join(DATADIR, "tests", "custom-box.toml"),
                self.atomic1,
                step=step,
            )

    def test_wrong_probe_in_format(self):
        # Check wrong probe in format
        for probe_in in [True, [1.4], {"probe_in": 1.4}, "1.4", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                get_vertices_from_file,
                os.path.join(DATADIR, "tests", "custom-box.toml"),
                self.atomic1,
                probe_in=probe_in,
            )

    def test_invalid_probe_in(self):
        # Check invalid probe in
        self.assertRaises(
            ValueError,
            get_vertices_from_file,
            os.path.join(DATADIR, "tests", "custom-box.toml"),
            self.atomic1,
            probe_in=-1.0,
        )

    def test_wrong_probe_out_format(self):
        # Check wrong probe out format
        for probe_out in [True, [4.0], {"probe_out": 4.0}, "4.0", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                get_vertices_from_file,
                os.path.join(DATADIR, "tests", "custom-box.toml"),
                self.atomic1,
                probe_out=probe_out,
            )

    def test_invalid_probe_out(self):
        # Check invalid probe out
        self.assertRaises(
            ValueError,
            get_vertices_from_file,
            os.path.join(DATADIR, "tests", "custom-box.toml"),
            self.atomic1,
            probe_out=-1.0,
        )

    def test_invalid_probe_pair(self):
        # Check probe in greater than probe out
        self.assertRaises(
            ValueError,
            get_vertices_from_file,
            os.path.join(DATADIR, "tests", "custom-box.toml"),
            self.atomic1,
            probe_in=4,
            probe_out=1.4,
        )

    def test_wrong_nthreads_format(self):
        # Check wrong nthreads format
        for nthreads in [1.0, [1], {"nthreads": 1}, numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                get_vertices_from_file,
                os.path.join(DATADIR, "tests", "custom-box.toml"),
                self.atomic1,
                nthreads=nthreads,
            )

    def test_invalid_nthreads(self):
        # Check invalid nthreads
        for nthreads in [-1, 0]:
            self.assertRaises(
                ValueError,
                get_vertices_from_file,
                os.path.join(DATADIR, "tests", "custom-box.toml"),
                self.atomic1,
                nthreads=nthreads,
            )

    @mock.patch(
        "toml.load",
        return_value={
            "SETTINGS": {
                "visiblebox": {
                    "p1": {"x": 3.11, "y": 7.34, "z": 1.59},
                    "p2": {"x": 11.51, "y": 7.34, "z": 1.59},
                    "p3": {"x": 3.11, "y": 10.74, "z": 1.59},
                    "p4": {"x": 3.11, "y": 7.34, "z": 6.19},
                }
            }
        },
    )
    def test_parKVFinder_box_file(self, _):
        # Expected vertices
        expected = [
            [-0.89, 3.34, -2.41],
            [15.51, 3.34, -2.41],
            [-0.89, 14.74, -2.41],
            [-0.89, 3.34, 10.19],
        ]
        # Get vertices from file
        vertices, _ = get_vertices_from_file("mocked", self.atomic1)
        self.assertListEqual(vertices.tolist(), expected)

    @mock.patch(
        "toml.load",
        return_value={
            "residues": [
                [
                    "49",
                    "E",
                    "LEU",
                ],
            ],
        },
    )
    def test_missing_padding(self, _):
        # Expected vertices
        expected = [
            [-4.0, 0.5, -3.5],
            [11.0, 0.5, -3.5],
            [-4.0, 15.5, -3.5],
            [-4.0, 0.5, 11.5],
        ]
        # Get vertices from file
        vertices, _ = get_vertices_from_file(
            os.path.join(DATADIR, "tests", "residues-box.toml"), self.atomic2
        )
        self.assertListEqual(vertices.tolist(), expected)

    @mock.patch(
        "toml.load",
        return_value={
            "box": {
                "p1": {"x": 3.11, "y": 7.34, "z": 1.59},
                "p2": {"x": 11.51, "y": 7.34, "z": 1.59},
                "p3": {"x": 3.11, "y": 10.74, "z": 1.59},
                "p4": {"x": 3.11, "y": 7.34, "z": 6.19},
                "padding": 3.5,
            }
        },
    )
    def test_invalid_box_keys(self, _):
        self.assertRaises(ValueError, get_vertices_from_file, "mocked", self.atomic1)

    @mock.patch(
        "toml.load",
        return_value={
            "residues": [
                [
                    "49",
                    "E",
                    "LEU",
                ],
            ],
            "padding": 3.5,
            "p1": {"x": 3.11, "y": 7.34, "z": 1.59},
        },
    )
    def test_invalid_residues_keys(self, _):
        self.assertRaises(ValueError, get_vertices_from_file, "mocked", self.atomic1)

    @mock.patch(
        "toml.load",
        return_value={"wrong_key": []},
    )
    def test_invalid_box(self, _):
        self.assertRaises(ValueError, get_vertices_from_file, "mocked", self.atomic1)


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

    def test_invalid_cavity_name(self):
        for cavity_name in ["AAA", "BAA", "CAA", "kaa", "1"]:
            self.assertRaises(ValueError, _get_cavity_label, cavity_name)


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

    def test_wrong_atomic_format(self):
        # Check wrong atomic format
        for atomic in [True, 4, 4.0, {"step": 4.0}, "4.0"]:
            self.assertRaises(
                TypeError, get_vertices, atomic, self.probe_out, self.step
            )

    def test_invalid_atomic(self):
        # Check invalid atomic
        for atomic in [
            [1.0, 1.0, 1.0, 1.0],  # shape (4,)
            [[1.0, 1.0, 1.0, 10, 1.0]],  # shape (1, 5)
            [[[1.0, 1.0, 1.0, 1.0]]],  # shape (1, 1, 4)
        ]:
            self.assertRaises(
                ValueError, get_vertices, atomic, self.probe_out, self.step
            )

    def test_wrong_probe_out_format(self):
        # Check wrong probe out format
        for probe_out in [True, [4.0], {"step": 4.0}, "4.0", numpy.ones(1)]:
            self.assertRaises(
                TypeError, get_vertices, self.atomic, probe_out, self.step
            )

    def test_invalid_probe_out(self):
        # Check invalid probe_out
        self.assertRaises(ValueError, get_vertices, self.atomic, -1, self.step)

    def test_wrong_step_format(self):
        # Check wrong step format
        for step in [True, [0.6], {"step": 0.6}, "0.6", numpy.ones(1)]:
            self.assertRaises(
                TypeError, get_vertices, self.atomic, self.probe_out, step
            )

    def test_invalid_step(self):
        # Check invalid step
        for step in [-1.0, 0.0]:
            self.assertRaises(
                ValueError, get_vertices, self.atomic, self.probe_out, step
            )


class TestGetDimensions(unittest.TestCase):
    def setUp(self):
        self.vertices = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]

    def test_vertices(self):
        # Test vertices input
        result = _get_dimensions(self.vertices)
        self.assertEqual(result, (2, 2, 2))

    def test_step_as_integer(self):
        result = _get_dimensions(self.vertices, step=1)
        self.assertEqual(result, (1, 1, 1))

    def test_wrong_vertices_format(self):
        for vertices in ["vertices", True, 1, 1.0, {"vertices": []}]:
            self.assertRaises(TypeError, _get_dimensions, vertices)

    def test_invalid_vertices(self):
        for vertices in [
            [1.0, 1.0, 1.0],  # shape (3,)
            [[1.0, 1.0, 1.0]],  # shape (1, 3)
            [[[1.0, 1.0, 1.0]]],  # shape (1, 1, 3)
        ]:
            self.assertRaises(ValueError, _get_dimensions, vertices)

    def test_wrong_step_format(self):
        for step in ["step", True, {"step": [0.6]}, [0.6, 0.6]]:
            self.assertRaises(TypeError, _get_dimensions, self.vertices, step=step)

    def test_invalid_step(self):
        for step in [-1.0, 0.0]:
            self.assertRaises(ValueError, _get_dimensions, self.vertices, step=step)


class TestGetSincos(unittest.TestCase):
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

    def test_wrong_vertices_format(self):
        for vertices in ["vertices", True, 1, 1.0, {"vertices": []}]:
            self.assertRaises(TypeError, _get_sincos, vertices)

    def test_invalid_vertices(self):
        for vertices in [
            [1.0, 1.0, 1.0],  # shape (3,)
            [[1.0, 1.0, 1.0]],  # shape (1, 3)
            [[[1.0, 1.0, 1.0]]],  # shape (1, 1, 3)
        ]:
            self.assertRaises(ValueError, _get_sincos, vertices)


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
    def setUp(self):
        self.atomic = read_pdb(os.path.join(DATADIR, "tests", "1FMO.pdb"))
        self.vertices = get_vertices(self.atomic)

    def test_detect(self):
        # Detect cavities
        ncav, cavities = detect(self.atomic, self.vertices)
        cavities[cavities == 1] = 0
        # Expected grid
        expected = read_cavity(
            os.path.join(DATADIR, "tests", "1FMO.KVFinder.output.pdb"),
            os.path.join(DATADIR, "tests", "1FMO.pdb"),
        )
        # Assert number of cavities
        self.assertEqual(ncav, int(expected.max() - 1))
        # Grid similarity
        tol = 100
        self.assertEqual((cavities - expected).sum() < tol, True)

    def test_detect_w_SAS(self):
        # Detect cavities
        ncav, _ = detect(self.atomic, self.vertices, surface="SAS")
        # Assert number of cavities
        self.assertEqual(ncav, 4)

    def test_atomic_as_list(self):
        # Detect cavities
        ncav, cavities = detect(self.atomic.tolist(), self.vertices)
        cavities[cavities == 1] = 0
        # Expected grid
        expected = read_cavity(
            os.path.join(DATADIR, "tests", "1FMO.KVFinder.output.pdb"),
            os.path.join(DATADIR, "tests", "1FMO.pdb"),
        )
        # Assert number of cavities
        self.assertEqual(ncav, int(expected.max() - 1))

    def test_vertices_as_list(self):
        # Detect cavities
        ncav, cavities = detect(self.atomic, self.vertices.tolist())
        cavities[cavities == 1] = 0
        # Expected grid
        expected = read_cavity(
            os.path.join(DATADIR, "tests", "1FMO.KVFinder.output.pdb"),
            os.path.join(DATADIR, "tests", "1FMO.pdb"),
        )
        # Assert number of cavities
        self.assertEqual(ncav, int(expected.max() - 1))

    def test_step_as_integer(self):
        # Detect cavities
        ncav, _ = detect(self.atomic, self.vertices, step=1)
        # Assert number of cavities
        self.assertEqual(ncav, 9)

    def test_probe_in_as_integer(self):
        # Detect cavities
        ncav, _ = detect(self.atomic, self.vertices, probe_in=1)
        # Assert number of cavities
        self.assertEqual(ncav, 30)

    def test_probe_out_as_integer(self):
        # Detect cavities
        ncav, cavities = detect(self.atomic, self.vertices, probe_out=4)
        cavities[cavities == 1] = 0
        # Expected grid
        expected = read_cavity(
            os.path.join(DATADIR, "tests", "1FMO.KVFinder.output.pdb"),
            os.path.join(DATADIR, "tests", "1FMO.pdb"),
        )
        # Assert number of cavities
        self.assertEqual(ncav, int(expected.max() - 1))

    def test_removal_distance_as_integer(self):
        # Detect cavities
        ncav, cavities = detect(self.atomic, self.vertices, removal_distance=2)
        cavities[cavities == 1] = 0
        # Expected grid
        expected = read_cavity(
            os.path.join(DATADIR, "tests", "1FMO.KVFinder.output.pdb"),
            os.path.join(DATADIR, "tests", "1FMO.pdb"),
        )
        # Assert number of cavities
        self.assertEqual(ncav, int(expected.max() - 1))

    def test_volume_cutoff_as_integer(self):
        # Detect cavities
        ncav, cavities = detect(self.atomic, self.vertices, volume_cutoff=5)
        cavities[cavities == 1] = 0
        # Expected grid
        expected = read_cavity(
            os.path.join(DATADIR, "tests", "1FMO.KVFinder.output.pdb"),
            os.path.join(DATADIR, "tests", "1FMO.pdb"),
        )
        # Assert number of cavities
        self.assertEqual(ncav, int(expected.max() - 1))

    def test_ligand_cutoff_as_integer(self):
        # Detect cavities
        ncav, cavities = detect(self.atomic, self.vertices, ligand_cutoff=5)
        cavities[cavities == 1] = 0
        # Expected grid
        expected = read_cavity(
            os.path.join(DATADIR, "tests", "1FMO.KVFinder.output.pdb"),
            os.path.join(DATADIR, "tests", "1FMO.pdb"),
        )
        # Assert number of cavities
        self.assertEqual(ncav, int(expected.max() - 1))

    def test_latomic_as_list(self):
        # Detect cavities
        ncav, cavities = detect(
            self.atomic, self.vertices, latomic=self.atomic.tolist()
        )
        cavities[cavities == 1] = 0
        # Expected grid
        expected = read_cavity(
            os.path.join(DATADIR, "tests", "1FMO.KVFinder.output.pdb"),
            os.path.join(DATADIR, "tests", "1FMO.pdb"),
        )
        # Assert number of cavities
        self.assertEqual(ncav, int(expected.max() - 1))

    def test_wrong_atomic_format(self):
        # Check wrong atomic format
        for atomic in [True, 4, 4.0, {"step": 4.0}, "4.0"]:
            self.assertRaises(TypeError, detect, atomic, self.vertices)

    def test_invalid_atomic(self):
        # Check invalid atomic
        for atomic in [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],  # shape (8,)
            [[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],  # shape (1, 9)
            [[[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]],  # shape (1, 1, 8)
        ]:
            self.assertRaises(ValueError, detect, atomic, self.vertices)

    def test_wrong_vertices_format(self):
        # Check wrong vertices format
        for vertices in ["vertices", True, 1, 1.0, {"vertices": []}]:
            self.assertRaises(TypeError, detect, self.atomic, vertices)

    def test_invalid_vertices(self):
        for vertices in [
            [1.0, 1.0, 1.0],  # shape (3,)
            [[1.0, 1.0, 1.0]],  # shape (1, 3)
            [[[1.0, 1.0, 1.0]]],  # shape (1, 1, 3)
        ]:
            self.assertRaises(ValueError, detect, self.atomic, vertices)

    def test_wrong_step_format(self):
        # Check wrong step format
        for step in [True, [0.6], {"step": 0.6}, "0.6", numpy.ones(1)]:
            self.assertRaises(TypeError, detect, self.atomic, self.vertices, step=step)

    def test_invalid_step(self):
        # Check invalid step
        for step in [-1.0, 0.0]:
            self.assertRaises(
                ValueError,
                detect,
                self.atomic,
                self.vertices,
                step=step,
            )

    def test_wrong_probe_in_format(self):
        # Check wrong probe out format
        for probe_in in [True, [1.4], {"step": 1.4}, "1.4", numpy.ones(1)]:
            self.assertRaises(
                TypeError, detect, self.atomic, self.vertices, probe_in=probe_in
            )

    def test_invalid_probe_in(self):
        # Check invalid probe in
        self.assertRaises(
            ValueError,
            detect,
            self.atomic,
            self.vertices,
            probe_in=-1,
        )

    def test_wrong_probe_out_format(self):
        # Check wrong probe out format
        for probe_out in [True, [4.0], {"step": 4.0}, "4.0", numpy.ones(1)]:
            self.assertRaises(
                TypeError, detect, self.atomic, self.vertices, probe_out=probe_out
            )

    def test_invalid_probe_out(self):
        # Check invalid probe out
        self.assertRaises(
            ValueError,
            detect,
            self.atomic,
            self.vertices,
            probe_out=-1,
        )

    def test_invalid_probe_pair(self):
        # Check probe in greater than probe out
        self.assertRaises(
            ValueError,
            detect,
            self.atomic,
            self.vertices,
            probe_in=4,
            probe_out=1.4,
        )

    def test_wrong_removal_distance_format(self):
        # Check wrong removal distance format
        for removal_distance in [
            True,
            [2.4],
            {"removal_distance": 2.4},
            "2.4",
            numpy.ones(1),
        ]:
            self.assertRaises(
                TypeError,
                detect,
                self.atomic,
                self.vertices,
                removal_distance=removal_distance,
            )

    def test_invalid_removal_distance(self):
        # Check invalid removal distance
        self.assertRaises(
            ValueError,
            detect,
            self.atomic,
            self.vertices,
            removal_distance=-1.0,
        )

    def test_wrong_volume_cutoff_format(self):
        # Check wrong volume cutoff format
        for volume_cutoff in [
            True,
            [5.0],
            {"volume_cutoff": 5.0},
            "5.0",
            numpy.ones(1),
        ]:
            self.assertRaises(
                TypeError,
                detect,
                self.atomic,
                self.vertices,
                volume_cutoff=volume_cutoff,
            )

    def test_invalid_volume_cutoff(self):
        # Check invalid volume cutoff
        self.assertRaises(
            ValueError,
            detect,
            self.atomic,
            self.vertices,
            volume_cutoff=-1,
        )

    def test_wrong_latomic_format(self):
        # Check wrong latomic format
        for latomic in [True, 4, 4.0, {"latomic": []}, "4.0"]:
            self.assertRaises(
                TypeError, detect, self.atomic, self.vertices, latomic=latomic
            )

    def test_invalid_latomic(self):
        # Check invalid latomic
        for latomic in [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],  # shape (8,)
            [[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],  # shape (1, 9)
            [[[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]],  # shape (1, 1, 8)
        ]:
            self.assertRaises(
                ValueError, detect, self.atomic, self.vertices, latomic=latomic
            )

    def test_wrong_ligand_cutoff_format(self):
        # Check wrong ligand cutoff format
        for ligand_cutoff in [
            True,
            [5.0],
            {"ligand_cutoff": 5.0},
            "5.0",
            numpy.ones(1),
        ]:
            self.assertRaises(
                TypeError,
                detect,
                self.atomic,
                self.vertices,
                ligand_cutoff=ligand_cutoff,
            )

    def test_invalid_ligand_cutoff(self):
        # Check invalid ligand cutoff
        self.assertRaises(
            ValueError,
            detect,
            self.atomic,
            self.vertices,
            ligand_cutoff=-1,
        )

    def test_wrong_box_adjustment_format(self):
        # Check wrong box adjustment format
        for box_adjustment in [
            1,
            1.0,
            [1.0],
            {"box_adjustment": 1.0},
            "1.0",
            numpy.ones(1),
        ]:
            self.assertRaises(
                TypeError,
                detect,
                self.atomic,
                self.vertices,
                box_adjustment=box_adjustment,
            )

    def test_wrong_surface_format(self):
        # Check wrong surface format
        for surface in [True, [4.0], {"surface": "SES"}, numpy.ones(1)]:
            self.assertRaises(
                TypeError, detect, self.atomic, self.vertices, surface=surface
            )

    def test_invalid_surface(self):
        # Check probe in greater than probe out
        self.assertRaises(
            ValueError,
            detect,
            self.atomic,
            self.vertices,
            surface="invalid",
        )

    def test_wrong_nthreads_format(self):
        # Check wrong nthreads format
        for nthreads in [1.0, [1.0], {"nthreads": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError, detect, self.atomic, self.vertices, nthreads=nthreads
            )

    def test_invalid_nthreads(self):
        # Check invalid nthreads
        for nthreads in [-1, 0]:
            self.assertRaises(
                ValueError,
                detect,
                self.atomic,
                self.vertices,
                nthreads=nthreads,
            )

    def test_wrong_verbose_format(self):
        # Check wrong verbose format
        for verbose in [1, 1.0, [4.0], {"verbose": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError, detect, self.atomic, self.vertices, verbose=verbose
            )


class TestSpatial(unittest.TestCase):
    def setUp(self):
        # Prepare data
        self.cavities = numpy.full((5, 5, 5), -1)
        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    self.cavities[i, j, k] = 0

        for i in range(2, 4):
            for j in range(2, 4):
                for k in range(2, 4):
                    self.cavities[i, j, k] = 2

    def test_spatial(self):
        # Surface
        surface, volume, area = spatial(self.cavities)
        # Assert
        self.assertDictEqual(volume, {"KAA": 1.73})
        self.assertDictEqual(area, {"KAA": 2.99})
        self.cavities[3, 3, 3] = -1
        self.assertListEqual(self.cavities.tolist(), surface.tolist())

    def test_step_as_integer(self):
        # Surface
        surface, volume, area = spatial(self.cavities, step=1)
        # Assert
        self.assertDictEqual(volume, {"KAA": 8.0})
        self.assertDictEqual(area, {"KAA": 8.29})
        self.cavities[3, 3, 3] = -1
        self.assertListEqual(self.cavities.tolist(), surface.tolist())

    def test_selection(self):
        cavities = numpy.random.randint(-1, 3, 27).reshape(3, 3, 3)
        for selection in [[2], ["KAA"]]:
            surface, _, _ = spatial(cavities, selection=selection)
            self.assertListEqual(numpy.unique(surface).tolist(), [-1, 0, 2])

    def test_wrong_cavities_format(self):
        for cavities in [1, 1.0, "1.0", [1.0], {"cavities": [1.0]}]:
            self.assertRaises(TypeError, spatial, cavities)

    def test_invalid_cavities(self):
        for cavities in [
            numpy.zeros((1)),  # shape (1,)
            numpy.zeros((1, 1)),  # shape (1, 1)
            numpy.zeros((1, 1, 1, 1)),  # shape (1, 1, 1, 1)
        ]:
            self.assertRaises(ValueError, spatial, cavities)

    def test_wrong_step_format(self):
        # Check wrong step format
        for step in [True, [0.6], {"step": 0.6}, "0.6", numpy.ones(1)]:
            self.assertRaises(TypeError, spatial, self.cavities, step=step)

    def test_invalid_step(self):
        # Check invalid step
        for step in [-1.0, 0.0]:
            self.assertRaises(
                ValueError,
                spatial,
                self.cavities,
                step=step,
            )

    def test_wrong_selection_format(self):
        self.assertRaises(TypeError, spatial, self.cavities, selection=["KAA", 2])

    def test_invalid_selection(self):
        for selection in [[1], [0], [-1]]:
            self.assertRaises(ValueError, spatial, self.cavities, selection=selection)

    def test_wrong_nthreads_format(self):
        # Check wrong nthreads format
        for nthreads in [1.0, [1.0], {"nthreads": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                spatial,
                self.cavities,
                nthreads=nthreads,
            )

    def test_invalid_nthreads(self):
        # Check invalid nthreads
        for nthreads in [-1, 0]:
            self.assertRaises(
                ValueError,
                spatial,
                self.cavities,
                nthreads=nthreads,
            )

    def test_wrong_verbose_format(self):
        # Check wrong verbose format
        for verbose in [1, 1.0, [4.0], {"verbose": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                spatial,
                self.cavities,
                verbose=verbose,
            )


class TestDepth(unittest.TestCase):
    def setUp(self):
        self.cavities = numpy.full((5, 5, 5), -1)
        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    self.cavities[i, j, k] = 0
        self.cavities[2, 2, 2] = 2
        self.cavities[2, 2, 1] = 2

    def test_depth(self):
        # Expected
        expected = numpy.zeros((5, 5, 5))
        expected[2, 2, 2] = 0.6
        # Depth
        depths, max_depth, avg_depth = depth(self.cavities)
        # Assert
        self.assertDictEqual(max_depth, {"KAA": 0.6})
        self.assertDictEqual(avg_depth, {"KAA": 0.3})
        self.assertListEqual(depths.tolist(), expected.tolist())

    def test_step_as_integer(self):
        # Expected
        expected = numpy.zeros((5, 5, 5))
        expected[2, 2, 2] = 1.0
        # Depth
        depths, max_depth, avg_depth = depth(self.cavities, step=1)
        # Assert
        self.assertDictEqual(max_depth, {"KAA": 1.0})
        self.assertDictEqual(avg_depth, {"KAA": 0.5})
        self.assertListEqual(depths.tolist(), expected.tolist())

    def test_selection(self):
        for selection in [[2], ["KAA"]]:
            depths, _, _ = depth(self.cavities, selection=selection)
            self.assertListEqual(
                numpy.argwhere(depths != 0).tolist(),
                [[2, 2, 2]],
            )

    def test_wrong_cavities_format(self):
        for cavities in [1, 1.0, "1.0", [1.0], {"cavities": [1.0]}]:
            self.assertRaises(TypeError, depth, cavities)

    def test_invalid_cavities(self):
        for cavities in [
            numpy.zeros((1)),  # shape (1,)
            numpy.zeros((1, 1)),  # shape (1, 1)
            numpy.zeros((1, 1, 1, 1)),  # shape (1, 1, 1, 1)
        ]:
            self.assertRaises(ValueError, depth, cavities)

    def test_wrong_step_format(self):
        # Check wrong step format
        for step in [True, [0.6], {"step": 0.6}, "0.6", numpy.ones(1)]:
            self.assertRaises(TypeError, depth, self.cavities, step=step)

    def test_invalid_step(self):
        # Check invalid step
        for step in [-1.0, 0.0]:
            self.assertRaises(
                ValueError,
                depth,
                self.cavities,
                step=step,
            )

    def test_wrong_selection_format(self):
        self.assertRaises(TypeError, depth, self.cavities, selection=["KAA", 2])

    def test_invalid_selection(self):
        for selection in [[1], [0], [-1]]:
            self.assertRaises(ValueError, depth, self.cavities, selection=selection)

    def test_wrong_nthreads_format(self):
        # Check wrong nthreads format
        for nthreads in [1.0, [1.0], {"nthreads": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                depth,
                self.cavities,
                nthreads=nthreads,
            )

    def test_invalid_nthreads(self):
        # Check invalid nthreads
        for nthreads in [-1, 0]:
            self.assertRaises(
                ValueError,
                depth,
                self.cavities,
                nthreads=nthreads,
            )

    def test_wrong_verbose_format(self):
        # Check wrong verbose format
        for verbose in [1, 1.0, [4.0], {"verbose": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                depth,
                self.cavities,
                verbose=verbose,
            )


class TestConstitutional(unittest.TestCase):
    def setUp(self):
        # Dummy cavities
        self.cavities = numpy.full((5, 5, 5), -1)
        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    self.cavities[i, j, k] = 0

        self.cavities[2, 2, 1] = 2
        self.cavities[2, 2, 4] = 3

        # Dummy atomic
        self.atomic = numpy.array(
            [
                ["13", "E", "GLU", "C", "1.0", "1.0", "1.0", "1.908"],
            ]
        )

        # Dummy vertices
        self.vertices = numpy.array(
            [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
        )

    def test_constitutional(self):
        # Constitutional
        residues = constitutional(self.cavities, self.atomic, self.vertices)
        # Assert
        self.assertDictEqual(
            residues, {"KAA": [["13", "E", "GLU"]], "KAB": [["13", "E", "GLU"]]}
        )

    def test_atomic_as_list(self):
        # Constitutional
        residues = constitutional(self.cavities, self.atomic.tolist(), self.vertices)
        # Assert
        self.assertDictEqual(
            residues, {"KAA": [["13", "E", "GLU"]], "KAB": [["13", "E", "GLU"]]}
        )

    def test_vertices_as_list(self):
        # Constitutional
        residues = constitutional(self.cavities, self.atomic, self.vertices.tolist())
        # Assert
        self.assertDictEqual(
            residues, {"KAA": [["13", "E", "GLU"]], "KAB": [["13", "E", "GLU"]]}
        )

    def test_step_as_integer(self):
        # Constitutional
        residues = constitutional(self.cavities, self.atomic, self.vertices, step=1)
        # Assert
        self.assertDictEqual(residues, {"KAA": [["13", "E", "GLU"]]})

    def test_probe_in_as_integer(self):
        # Constitutional
        residues = constitutional(self.cavities, self.atomic, self.vertices, probe_in=1)
        # Assert
        self.assertDictEqual(
            residues, {"KAA": [["13", "E", "GLU"]], "KAB": [["13", "E", "GLU"]]}
        )

    def test_selection(self):
        for selection in [[2], ["KAA"]]:
            residues = constitutional(
                self.cavities, self.atomic, self.vertices, selection=selection
            )
            self.assertDictEqual(residues, {"KAA": [["13", "E", "GLU"]]})

    def test_wrong_cavities_format(self):
        for cavities in [1, 1.0, "1.0", [1.0], {"cavities": [1.0]}]:
            self.assertRaises(
                TypeError, constitutional, cavities, self.atomic, self.vertices
            )

    def test_invalid_cavities(self):
        for cavities in [
            numpy.zeros((1)),  # shape (1,)
            numpy.zeros((1, 1)),  # shape (1, 1)
            numpy.zeros((1, 1, 1, 1)),  # shape (1, 1, 1, 1)
        ]:
            self.assertRaises(
                ValueError, constitutional, cavities, self.atomic, self.vertices
            )

    def test_wrong_atomic_format(self):
        # Check wrong atomic format
        for atomic in [True, 4, 4.0, {"step": 4.0}, "4.0"]:
            self.assertRaises(
                TypeError, constitutional, self.cavities, atomic, self.vertices
            )

    def test_invalid_atomic(self):
        # Check invalid atomic
        for atomic in [
            [1.0, 1.0, 1.0, 1.0],  # shape (4,)
            [[1.0, 1.0, 1.0, 10, 1.0]],  # shape (1, 5)
            [[[1.0, 1.0, 1.0, 1.0]]],  # shape (1, 1, 4)
        ]:
            self.assertRaises(
                ValueError, constitutional, self.cavities, atomic, self.vertices
            )

    def test_wrong_vertices_format(self):
        for vertices in ["vertices", True, 1, 1.0, {"vertices": []}]:
            self.assertRaises(
                TypeError, constitutional, self.cavities, self.atomic, vertices
            )

    def test_invalid_vertices(self):
        for vertices in [
            [1.0, 1.0, 1.0],  # shape (3,)
            [[1.0, 1.0, 1.0]],  # shape (1, 3)
            [[[1.0, 1.0, 1.0]]],  # shape (1, 1, 3)
        ]:
            self.assertRaises(
                ValueError, constitutional, self.cavities, self.atomic, vertices
            )

    def test_wrong_step_format(self):
        # Check wrong step format
        for step in [True, [0.6], {"step": 0.6}, "0.6", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                constitutional,
                self.cavities,
                self.atomic,
                self.vertices,
                step=step,
            )

    def test_invalid_step(self):
        # Check invalid step
        for step in [-1.0, 0.0]:
            self.assertRaises(
                ValueError,
                constitutional,
                self.cavities,
                self.atomic,
                self.vertices,
                step=step,
            )

    def test_wrong_probe_in_format(self):
        # Check wrong probe out format
        for probe_in in [True, [1.4], {"probe_in": 1.4}, "1.4", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                constitutional,
                self.cavities,
                self.atomic,
                self.vertices,
                probe_in=probe_in,
            )

    def test_invalid_probe_in(self):
        # Check invalid probe in
        self.assertRaises(
            ValueError,
            constitutional,
            self.cavities,
            self.atomic,
            self.vertices,
            probe_in=-1,
        )

    def test_wrong_ignore_backbone_format(self):
        for ignore_backbone in [1, 1.0, "1", [True], {"ignore_backbone": True}]:
            self.assertRaises(
                TypeError,
                constitutional,
                self.cavities,
                self.atomic,
                self.vertices,
                ignore_backbone=ignore_backbone,
            )

    def test_wrong_selection_format(self):
        self.assertRaises(
            TypeError,
            constitutional,
            self.cavities,
            self.atomic,
            self.vertices,
            selection=["KAA", 2],
        )

    def test_invalid_selection(self):
        for selection in [[1], [0], [-1]]:
            self.assertRaises(
                ValueError,
                constitutional,
                self.cavities,
                self.atomic,
                self.vertices,
                selection=selection,
            )

    def test_wrong_nthreads_format(self):
        # Check wrong nthreads format
        for nthreads in [1.0, [1.0], {"nthreads": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                constitutional,
                self.cavities,
                self.atomic,
                self.vertices,
                nthreads=nthreads,
            )

    def test_invalid_nthreads(self):
        # Check invalid nthreads
        for nthreads in [-1, 0]:
            self.assertRaises(
                ValueError,
                constitutional,
                self.cavities,
                self.atomic,
                self.vertices,
                nthreads=nthreads,
            )

    def test_wrong_verbose_format(self):
        # Check wrong verbose format
        for verbose in [1, 1.0, [4.0], {"verbose": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                constitutional,
                self.cavities,
                self.atomic,
                self.vertices,
                verbose=verbose,
            )


class TestHydropathy(unittest.TestCase):
    def setUp(self):
        # Surface
        self.surface = numpy.full((5, 5, 5), -1)
        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    self.surface[i, j, k] = 0

        self.surface[2, 2, 1] = 2
        self.surface[2, 2, 4] = 3
        # Dummy atomic
        self.atomic = numpy.array(
            [
                ["13", "E", "GLU", "C", "1.0", "1.0", "1.0", "1.908"],
            ]
        )
        # Dummy vertices
        self.vertices = numpy.array(
            [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
        )
        # Expected
        self.expected = numpy.zeros((5, 5, 5))
        self.expected[2, 2, 1] = 0.76
        self.expected[2, 2, 4] = 0.76

    def test_hydropathy(self):
        # Surface
        scales, avg_hydropathy = hydropathy(self.surface, self.atomic, self.vertices)
        # Assert
        self.assertListEqual(scales.tolist(), self.expected.tolist())
        self.assertDictEqual(
            avg_hydropathy, {"KAA": 0.76, "KAB": 0.76, "EisenbergWeiss": [-1.42, 2.6]}
        )

    def test_atomic_as_list(self):
        # Surface
        scales, avg_hydropathy = hydropathy(
            self.surface, self.atomic.tolist(), self.vertices
        )
        # Assert
        self.assertListEqual(scales.tolist(), self.expected.tolist())
        self.assertDictEqual(
            avg_hydropathy, {"KAA": 0.76, "KAB": 0.76, "EisenbergWeiss": [-1.42, 2.6]}
        )

    def test_vertices_as_list(self):
        # Surface
        scales, avg_hydropathy = hydropathy(
            self.surface, self.atomic, self.vertices.tolist()
        )
        # Assert
        self.assertListEqual(scales.tolist(), self.expected.tolist())
        self.assertDictEqual(
            avg_hydropathy, {"KAA": 0.76, "KAB": 0.76, "EisenbergWeiss": [-1.42, 2.6]}
        )

    def test_step_as_integer(self):
        # Surface
        scales, avg_hydropathy = hydropathy(
            self.surface, self.atomic, self.vertices, step=1
        )
        # Assert
        self.assertListEqual(scales.tolist(), self.expected.tolist())
        self.assertDictEqual(
            avg_hydropathy, {"KAA": 0.76, "KAB": 0.76, "EisenbergWeiss": [-1.42, 2.6]}
        )

    def test_probe_in_as_integer(self):
        # Surface
        scales, avg_hydropathy = hydropathy(
            self.surface, self.atomic, self.vertices, probe_in=1
        )
        # Assert
        self.assertListEqual(scales.tolist(), self.expected.tolist())
        self.assertDictEqual(
            avg_hydropathy, {"KAA": 0.76, "KAB": 0.76, "EisenbergWeiss": [-1.42, 2.6]}
        )

    def test_selection(self):
        for selection in [[2], ["KAA"]]:
            scales, avg_hydropathy = hydropathy(
                self.surface, self.atomic, self.vertices, selection=selection
            )
            self.assertListEqual(
                scales.tolist(),
                numpy.where(self.surface == 2, self.expected, 0).tolist(),
            )
            self.assertDictEqual(
                avg_hydropathy, {"KAA": 0.76, "EisenbergWeiss": [-1.42, 2.6]}
            )

    def test_ignore_backbone(self):
        # Surface
        scales, avg_hydropathy = hydropathy(
            self.surface, self.atomic, self.vertices, ignore_backbone=True
        )
        # Assert
        self.assertListEqual(scales.tolist(), numpy.zeros((5, 5, 5)).tolist())
        self.assertDictEqual(
            avg_hydropathy, {"KAA": 0.0, "KAB": 0.0, "EisenbergWeiss": [-1.42, 2.6]}
        )

    def test_wrong_surface_format(self):
        for surface in [1, 1.0, "1.0", [1.0], {"cavities": [1.0]}]:
            self.assertRaises(
                TypeError, hydropathy, surface, self.atomic, self.vertices
            )

    def test_invalid_surface(self):
        for surface in [
            numpy.zeros((1)),  # shape (1,)
            numpy.zeros((1, 1)),  # shape (1, 1)
            numpy.zeros((1, 1, 1, 1)),  # shape (1, 1, 1, 1)
        ]:
            self.assertRaises(
                ValueError, hydropathy, surface, self.atomic, self.vertices
            )

    def test_wrong_atomic_format(self):
        # Check wrong atomic format
        for atomic in [True, 4, 4.0, {"step": 4.0}, "4.0"]:
            self.assertRaises(
                TypeError, hydropathy, self.surface, atomic, self.vertices
            )

    def test_invalid_atomic(self):
        # Check invalid atomic
        for atomic in [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],  # shape (8,)
            [[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],  # shape (1, 9)
            [[[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]],  # shape (1, 1, 8)
        ]:
            self.assertRaises(
                ValueError, hydropathy, self.surface, atomic, self.vertices
            )

    def test_wrong_vertices_format(self):
        # Check wrong vertices format
        for vertices in ["vertices", True, 1, 1.0, {"vertices": []}]:
            self.assertRaises(
                TypeError, hydropathy, self.surface, self.atomic, vertices
            )

    def test_invalid_vertices(self):
        for vertices in [
            [1.0, 1.0, 1.0],  # shape (3,)
            [[1.0, 1.0, 1.0]],  # shape (1, 3)
            [[[1.0, 1.0, 1.0]]],  # shape (1, 1, 3)
        ]:
            self.assertRaises(
                ValueError, hydropathy, self.surface, self.atomic, vertices
            )

    def test_wrong_step_format(self):
        # Check wrong step format
        for step in [True, [0.6], {"step": 0.6}, "0.6", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                hydropathy,
                self.surface,
                self.atomic,
                self.vertices,
                step=step,
            )

    def test_invalid_step(self):
        # Check invalid step
        for step in [-1.0, 0.0]:
            self.assertRaises(
                ValueError,
                hydropathy,
                self.surface,
                self.atomic,
                self.vertices,
                step=step,
            )

    def test_wrong_probe_in_format(self):
        # Check wrong probe out format
        for probe_in in [True, [1.4], {"probe_in": 1.4}, "1.4", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                hydropathy,
                self.surface,
                self.atomic,
                self.vertices,
                probe_in=probe_in,
            )

    def test_invalid_probe_in(self):
        # Check invalid probe in
        self.assertRaises(
            ValueError,
            hydropathy,
            self.surface,
            self.atomic,
            self.vertices,
            probe_in=-1,
        )

    def test_wrong_hydrophobicity_scale_format(self):
        for hydrophobicity_scale in [
            True,
            1,
            1.0,
            [""],
            {"hydrophobicity_scale": True},
        ]:
            self.assertRaises(
                TypeError,
                hydropathy,
                self.surface,
                self.atomic,
                self.vertices,
                hydrophobicity_scale=hydrophobicity_scale,
            )

    def test_wrong_ignore_backbone_format(self):
        for ignore_backbone in [1, 1.0, "1", [True], {"ignore_backbone": True}]:
            self.assertRaises(
                TypeError,
                hydropathy,
                self.surface,
                self.atomic,
                self.vertices,
                ignore_backbone=ignore_backbone,
            )

    def test_wrong_selection_format(self):
        self.assertRaises(
            TypeError,
            hydropathy,
            self.surface,
            self.atomic,
            self.vertices,
            selection=["KAA", 2],
        )

    def test_invalid_selection(self):
        for selection in [[1], [0], [-1]]:
            self.assertRaises(
                ValueError,
                hydropathy,
                self.surface,
                self.atomic,
                self.vertices,
                selection=selection,
            )

    def test_wrong_nthreads_format(self):
        # Check wrong nthreads format
        for nthreads in [1.0, [1.0], {"nthreads": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                hydropathy,
                self.surface,
                self.atomic,
                self.vertices,
                nthreads=nthreads,
            )

    def test_invalid_nthreads(self):
        # Check invalid nthreads
        for nthreads in [-1, 0]:
            self.assertRaises(
                ValueError,
                hydropathy,
                self.surface,
                self.atomic,
                self.vertices,
                nthreads=nthreads,
            )

    def test_wrong_verbose_format(self):
        # Check wrong verbose format
        for verbose in [1, 1.0, [4.0], {"verbose": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                hydropathy,
                self.surface,
                self.atomic,
                self.vertices,
                verbose=verbose,
            )


class TestExport(unittest.TestCase):
    def setUp(self):
        self.cavities = 2 * numpy.ones((1, 1, 1), dtype=numpy.int32)
        self.surface = 2 * numpy.ones((1, 1, 1), dtype=numpy.int32)
        self.depths = numpy.random.randint(0, 10, (1, 1, 1)) * numpy.random.rand(1)
        self.scales = ((-1.42) + numpy.random.rand(1) * (2.6 - (-1.42))).reshape(
            1, 1, 1
        )
        self.vertices = numpy.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )

    def test_cavity_wo_surface(self):
        # Cavity without surface
        export(
            "tests/cavities_wo_surface.pdb",
            self.cavities,
            None,
            self.vertices,
            step=1.0,
            B=None,
            Q=None,
        )
        with open("tests/cavities_wo_surface.pdb", "r") as f:
            self.assertEqual(
                f.read(),
                "ATOM      1  H   KAA   259       0.000   0.000   0.000  1.00  0.00\n",
            )
        os.remove("tests/cavities_wo_surface.pdb")

    def test_cavity_w_surface(self):
        # Cavity with surface
        export(
            "tests/cavities.pdb",
            self.cavities,
            self.surface,
            self.vertices,
            step=1.0,
            B=None,
            Q=None,
        )
        with open("tests/cavities.pdb", "r") as f:
            self.assertEqual(
                f.read(),
                "ATOM      1  HA  KAA   259       0.000   0.000   0.000  1.00  0.00\n",
            )
        os.remove("tests/cavities.pdb")

    def test_cavity_with_depth(self):
        # Cavity with depth
        export(
            "tests/cavities_with_depth.pdb",
            self.cavities,
            self.surface,
            self.vertices,
            step=1.0,
            B=self.depths,
            Q=None,
        )
        with open("tests/cavities_with_depth.pdb", "r") as f:
            self.assertEqual(
                f.read(),
                f"ATOM      1  HA  KAA   259       0.000   0.000   0.000  1.00{self.depths.item():6.2f}\n",
            )
        os.remove("tests/cavities_with_depth.pdb")

    def test_cavity_and_hydropathy(self):
        # Cavity and hydropathy
        export(
            "tests/cavities.pdb",
            self.cavities,
            self.surface,
            self.vertices,
            step=1.0,
            B=self.depths,
            Q=self.scales
        )
        # Assert cavity
        with open("tests/cavities.pdb", "r") as f:
            self.assertEqual(
                f.read(),
                f"ATOM      1  HA  KAA   259       0.000   0.000   0.000{self.scales.item():6.2f}{self.depths.item():6.2f}\n",
            )
        os.remove("tests/cavities.pdb")

    def test_hydropathy(self):
        # Hydropathy
        export(
            "tests/hydropathy.pdb",
            self.surface,
            self.surface,
            self.vertices,
            step=1.0,
            Q=self.scales,
        )
        # Assert hydropathy
        with open("tests/hydropathy.pdb", "r") as f:
            self.assertEqual(
                f.read(),
                f"ATOM      1  HA  KAA   259       0.000   0.000   0.000{self.scales.item():6.2f}  0.00\n",
            )
        os.remove("tests/hydropathy.pdb")

    def test_selection(self):
        # Selection
        for selection in [[2], ["KAA"]]:
            export(
                "tests/cavities.pdb",
                self.cavities,
                self.surface,
                self.vertices,
                step=1.0,
                selection=selection,
            )
            # Assert selection
            with open("tests/cavities.pdb", "r") as f:
                self.assertEqual(
                    f.read(),
                    f"ATOM      1  HA  KAA   259       0.000   0.000   0.000  1.00  0.00\n",
                )
            os.remove("tests/cavities.pdb")

        for selection in [[3], ["KAB"]]:
            export(
                "tests/cavities.pdb",
                self.cavities,
                self.surface,
                self.vertices,
                step=1.0,
                selection=selection,
            )
            # Assert selection
            with open("tests/cavities.pdb", "r") as f:
                self.assertEqual(
                    f.read(),
                    f"",
                )
            os.remove("tests/cavities.pdb")

    def test_vertices_as_list(self):
        # Vertices as list
        export(
            "tests/cavities.pdb",
            self.cavities,
            self.surface,
            self.vertices.tolist(),
            step=1.0,
        )
        with open("tests/cavities.pdb", "r") as f:
            self.assertEqual(
                f.read(),
                "ATOM      1  HA  KAA   259       0.000   0.000   0.000  1.00  0.00\n",
            )
        os.remove("tests/cavities.pdb")

    def test_step_as_integer(self):
        # Step as integer
        export(
            "tests/cavities.pdb",
            self.cavities,
            self.surface,
            self.vertices.tolist(),
            step=1,
        )
        with open("tests/cavities.pdb", "r") as f:
            self.assertEqual(
                f.read(),
                "ATOM      1  HA  KAA   259       0.000   0.000   0.000  1.00  0.00\n",
            )
        os.remove("tests/cavities.pdb")

    def test_wrong_fn_format(self):
        # Check wrong fn format
        for fn in [1.0, [1], {"nthreads": 1}, numpy.ones(1)]:
            self.assertRaises(
                TypeError, export, fn, self.cavities, self.surface, self.vertices
            )

    def test_wrong_cavities_format(self):
        for cavities in [1, 1.0, "1.0", [1.0], {"cavities": [1.0]}]:
            self.assertRaises(
                TypeError, export, "any", cavities, self.surface, self.vertices
            )

    def test_invalid_cavities(self):
        for cavities in [
            numpy.zeros((1)),  # shape (1,)
            numpy.zeros((1, 1)),  # shape (1, 1)
            numpy.zeros((1, 1, 1, 1)),  # shape (1, 1, 1, 1)
        ]:
            self.assertRaises(
                ValueError, export, "any", cavities, self.surface, self.vertices
            )

    def test_wrong_surface_format(self):
        for surface in [1, 1.0, "1.0", [1.0], {"cavities": [1.0]}]:
            self.assertRaises(
                TypeError, export, "any", self.cavities, surface, self.vertices
            )

    def test_invalid_surface(self):
        for surface in [
            numpy.zeros((1)),  # shape (1,)
            numpy.zeros((1, 1)),  # shape (1, 1)
            numpy.zeros((1, 1, 1, 1)),  # shape (1, 1, 1, 1)
        ]:
            self.assertRaises(
                ValueError, export, "any", self.cavities, surface, self.vertices
            )

    def test_wrong_vertices_format(self):
        for vertices in ["vertices", True, 1, 1.0, {"vertices": []}]:
            self.assertRaises(
                TypeError, export, "any", self.cavities, self.surface, vertices
            )

    def test_invalid_vertices(self):
        for vertices in [
            [1.0, 1.0, 1.0],  # shape (3,)
            [[1.0, 1.0, 1.0]],  # shape (1, 3)
            [[[1.0, 1.0, 1.0]]],  # shape (1, 1, 3)
        ]:
            self.assertRaises(
                ValueError, export, "any", self.cavities, self.surface, vertices
            )

    def test_wrong_step_format(self):
        # Check wrong step format
        for step in [True, [0.6], {"step": 0.6}, "0.6", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                step=step,
            )

    def test_invalid_step(self):
        # Check invalid step
        for step in [-1.0, 0.0]:
            self.assertRaises(
                ValueError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                step=step,
            )

    def test_wrong_output_hydropathy_format(self):
        # Check wrong fn format
        for output_hydropathy in [1.0, [1], {"nthreads": 1}, numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                output_hydropathy=output_hydropathy,
            )

    def test_wrong_B_format(self):
        for B in [1, 1.0, "1.0", [1.0], {"cavities": [1.0]}]:
            self.assertRaises(
                TypeError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                B=B,
            )

    def test_invalid_B(self):
        for B in [
            numpy.zeros((1)),  # shape (1,)
            numpy.zeros((1, 1)),  # shape (1, 1)
            numpy.zeros((1, 1, 1, 1)),  # shape (1, 1, 1, 1)
        ]:
            self.assertRaises(
                ValueError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                B=B,
            )

    def test_wrong_Q_format(self):
        for scales in [1, 1.0, "1.0", [1.0], {"cavities": [1.0]}]:
            self.assertRaises(
                TypeError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                Q=scales,
            )

    def test_invalid_scales(self):
        for scales in [
            numpy.zeros((1)),  # shape (1,)
            numpy.zeros((1, 1)),  # shape (1, 1)
            numpy.zeros((1, 1, 1, 1)),  # shape (1, 1, 1, 1)
        ]:
            self.assertRaises(
                ValueError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                Q=scales,
            )

    def test_wrong_selection_format(self):
        self.assertRaises(
            TypeError,
            export,
            "any",
            self.cavities,
            self.surface,
            self.vertices,
            selection=["KAA", 2],
        )

    def test_invalid_selection(self):
        for selection in [[1], [0], [-1]]:
            self.assertRaises(
                ValueError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                selection=selection,
            )

    def test_wrong_nthreads_format(self):
        # Check wrong nthreads format
        for nthreads in [1.0, [1.0], {"nthreads": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                nthreads=nthreads,
            )

    def test_invalid_nthreads(self):
        # Check invalid nthreads
        for nthreads in [-1, 0]:
            self.assertRaises(
                ValueError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                nthreads=nthreads,
            )

    def test_wrong_verbose_format(self):
        # Check wrong verbose format
        for verbose in [1, 1.0, [4.0], {"verbose": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                verbose=verbose,
            )

    def test_wrong_append_format(self):
        # Check wrong append format
        for append in [1, 1.0, [4.0], {"append": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                append=append,
            )

    def test_wrong_model_format(self):
        # Check wrong model format
        for model in [1.0, [1.0], {"model": 1}, "1", numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                export,
                "any",
                self.cavities,
                self.surface,
                self.vertices,
                model=model,
            )

    def test_cavity_and_surface_not_defined(self):
        self.assertRaises(
            TypeError,
            export,
            "any",
            None,
            None,
            self.vertices,
        )


class TestMolecule(unittest.TestCase):
    def setUp(self):
        self.C = pyKVFinder.Molecule(
            os.path.join(FIXTURES, "C.pdb"), radii={"GEN": {"C": 1.66}}
        )
        self.H = pyKVFinder.Molecule(os.path.join(FIXTURES, "H.pdb"))
        self.N = pyKVFinder.Molecule(os.path.join(FIXTURES, "N.pdb"))

    def test_atomic(self):
        self.assertListEqual(
            self.C.atomic.tolist(),
            [["1", "A", "UNK", "C", "0.0", "0.0", "0.0", "1.66"]],
        )

    def test_radii(self):
        self.assertDictEqual(self.C.radii, {"GEN": {"C": 1.66}})

    def test_get_padding(self):
        self.assertEqual(self.C._get_padding(), 1.8)
        self.assertEqual(self.H._get_padding(), 1.0)
        self.assertEqual(self.N._get_padding(), 2.2)

    def test_vdw(self):
        expected = [
            [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]],
            [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]],
            [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 0, 1], [1, 1, 1, 1]],
            [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]],
        ]
        self.H.vdw(1.0)
        self.assertListEqual(self.H.grid.tolist(), expected)

    def test_vdw_volume(self):
        self.H.vdw(1.0)
        self.assertEqual(self.H.volume(), 1.0)

    def test_ses(self):
        expected = [
            [
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
            ],
            [
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
            ],
            [
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 1, 1, 1, 1],
            ],
            [
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 1, 1, 1, 1],
            ],
            [
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 1, 1, 1, 1],
            ],
            [
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
            ],
        ]
        self.H.surface(1.0, 1.0, "SES", None)
        self.assertListEqual(self.H.grid.tolist(), expected)

    def test_ses_volume(self):
        self.H.surface(1.0, 1.0, "SES", None)
        self.assertEqual(self.H.volume(), 27.0)

    def test_sas(self):
        expected = [
            [
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
            ],
            [
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
            ],
            [
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 1, 1, 1, 1],
            ],
            [
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 1, 1, 1, 1],
            ],
            [
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 0, 0, 0, 1],
                [1, 1, 1, 1, 1, 1],
            ],
            [
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
            ],
        ]
        self.H.surface(1.0, 1.0, "SAS", None)
        self.assertListEqual(self.H.grid.tolist(), expected)

    def test_sas_volume(self):
        self.H.surface(1.0, 1.0, "SAS", None)
        self.assertEqual(self.H.volume(), 27.0)


class TestGetOpeningName(unittest.TestCase):
    def test_indexes(self):
        indexes = list(range(0, 100, 10))
        opening_names = [_get_opening_name(index) for index in indexes]
        self.assertListEqual(
            opening_names,
            ["OAA", "OAK", "OAU", "OBE", "OBO", "OBY", "OCI", "OCS", "ODC", "ODM"],
        )


class TestGetOpeningLabel(unittest.TestCase):
    def test_opening_names(self):
        opening_names = [
            "OAA",
            "OAK",
            "OAU",
            "OBE",
            "OBO",
            "OBY",
            "OCI",
            "OCS",
            "ODC",
            "ODM",
        ]
        labels = [_get_opening_label(name) for name in opening_names]
        self.assertListEqual(labels, list(range(2, 100, 10)))

    def test_invalid_opening_name(self):
        self.assertRaises(ValueError, _get_opening_label, "AAA")


class TestProcessOpenings(unittest.TestCase):
    pass


class TestOpenings(unittest.TestCase):
    pass


class TestExportOpenings(unittest.TestCase):
    pass
