import io
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
        self.assertRaises(ValueError, _get_cavity_label, "AAA")


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

    def test_wrong_probe_out_format(self):
        # Check wrong probe out format
        for probe_out in [True, [4.0], {"step": 4.0}, "4.0", numpy.ones(1)]:
            self.assertRaises(
                TypeError, get_vertices, self.atomic, probe_out, self.step
            )

    def test_wrong_step_format(self):
        # Check wrong step format
        for step in [True, [0.6], {"step": 0.6}, "0.6", numpy.ones(1)]:
            self.assertRaises(
                TypeError, get_vertices, self.atomic, self.probe_out, step
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

    def test_invalid_probe_out(self):
        # Check invalid probe_out
        self.assertRaises(ValueError, get_vertices, self.atomic, -1, self.step)

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

    def test_wrong_step_format(self):
        for step in ["step", True, {"step": [0.6]}, [0.6, 0.6]]:
            self.assertRaises(TypeError, _get_dimensions, self.vertices, step=step)

    def test_invalid_vertices(self):
        for vertices in [
            [1.0, 1.0, 1.0],  # shape (3,)
            [[1.0, 1.0, 1.0]],  # shape (1, 3)
            [[[1.0, 1.0, 1.0]]],  # shape (1, 1, 3)
        ]:
            self.assertRaises(ValueError, _get_dimensions, vertices)

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

    def test_verbose_prints(self):
        # Check verbose prints to stdout
        ncav, cavities = detect(self.atomic, self.vertices)
        # surface='SES'
        with mock.patch("sys.stdout", new_callable=io.StringIO) as stdout:
            with mock.patch("_pyKVFinder._detect", return_value=(ncav, cavities)) as _:
                # Detect cavities
                _, _ = detect(self.atomic, self.vertices, verbose=True)
        expected = f"> Surface representation: Solvent Excluded Surface (SES)\n"
        self.assertEqual(stdout.getvalue(), expected)
        # surface='SAS'
        with mock.patch("sys.stdout", new_callable=io.StringIO) as stdout:
            with mock.patch("_pyKVFinder._detect", return_value=(ncav, cavities)) as _:
                # Detect cavities
                _, _ = detect(self.atomic, self.vertices, surface="SAS", verbose=True)
        expected = f"> Surface representation: Solvent Accessible Surface (SAS)\n"
        self.assertEqual(stdout.getvalue(), expected)


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
