import unittest
import os
import numpy
import argparse
from pyKVFinder.utils import (
    read_vdw,
    _process_pdb_line,
    read_pdb,
    read_xyz,
    _read_cavity,
    read_cavity,
    _process_box,
    calculate_frequencies,
)

UNIT_TESTS_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


class TestReadVdw(unittest.TestCase):
    def test_valid_vdw(self):
        # Check valid files formats
        expected = {
            "ALA": {"N": 1.824, "C": 1.908, "H": 0.6, "O": 1.6612},
            "GEN": {"N": 1.824, "C": 1.908, "H": 0.6, "O": 1.6612},
        }
        result = read_vdw(os.path.join(UNIT_TESTS_DIR, "vdw1.dat"))
        self.assertEqual(result, expected)

    def test_invalid_vdw(self):
        # Check invalid files formats
        self.assertRaises(
            ValueError, read_vdw, os.path.join(UNIT_TESTS_DIR, "vdw2.dat")
        )
        self.assertRaises(
            ValueError, read_vdw, os.path.join(UNIT_TESTS_DIR, "vdw3.dat")
        )


class TestProcessPdbLine(unittest.TestCase):
    def setUp(self):
        self.lines = [
            "ATOM      1   C  UNK A   1       0.000   0.000   0.000  1.00  0.00      A    C  ",
            "ATOM      1   N  UNK A   1       0.000   0.000   0.000  1.00  0.00      A    N  ",
            "ATOM      1   O  UNK A   1       0.000   0.000   0.000  1.00  0.00      A    O  ",
            "ATOM      1   H  UNK A   1       0.000   0.000   0.000  1.00  0.00      A    H  ",
            "HETATM    1  N9  UNK A   1       0.000   0.000   0.000  1.00  0.00      C    N  ",
            "HETATM    1  C5' UNK A   1       0.000   0.000   0.000  1.00  0.00      C    C  ",
        ]
        self.vdw = {
            "GEN": {
                "N": 1.824,
                "C": 1.66,
                "O": 1.69,
                "H": 0.91,
            }
        }
        self.expected = [
            [1, 'A', 'UNK', 'C', 0.000, 0.000, 0.000, 1.66],
            [1, 'A', 'UNK', 'N', 0.000, 0.000, 0.000, 1.824],
            [1, 'A', 'UNK', 'O', 0.000, 0.000, 0.000, 1.69],
            [1, 'A', 'UNK', 'H', 0.000, 0.000, 0.000, 0.91],
            [1, 'A', 'UNK', 'N9', 0.000, 0.000, 0.000, 1.824],
            [1, 'A', 'UNK', 'C5\'', 0.000, 0.000, 0.000, 1.66],
        ]

    def test_line(self):
        # Read different lines of a pdb file
        for data, expected in zip(self.lines, self.expected):
            result = _process_pdb_line(data, self.vdw)
            self.assertListEqual(result, expected)

class TestReadPdb(unittest.TestCase):
    def test_atom(self):
        expected = [
            ["13", "E", "GLU", "N", "-6.693", "-15.642", "-14.858", "1.824"],
            ["13", "E", "GLU", "CA", "-6.73", "-14.62", "-15.897", "1.908"],
            ["13", "E", "GLU", "C", "-7.49", "-13.357", "-15.508", "1.908"],
            ["13", "E", "GLU", "O", "-7.388", "-12.338", "-16.197", "1.6612"],
            ["13", "E", "GLU", "CB", "-7.166", "-15.114", "-17.306", "1.908"],
            ["13", "E", "GLU", "CG", "-7.485", "-16.623", "-17.424", "1.908"],
            ["13", "E", "GLU", "CD", "-6.267", "-17.503", "-17.618", "1.908"],
            ["13", "E", "GLU", "OE1", "-5.439", "-17.684", "-16.719", "1.6612"],
            ["13", "E", "GLU", "OE2", "-6.178", "-18.047", "-18.817", "1.6612"],
        ]
        result = read_pdb(os.path.join(UNIT_TESTS_DIR, "atom.pdb")).tolist()
        self.assertListEqual(result, expected)

    def test_hetatm(self):
        # Test HETATM entries
        expected = [
            ["351", "E", "ADN", "C5'", "11.087", "9.79", "2.052", "1.66"],
            ["351", "E", "ADN", "O5'", "11.545", "8.52", "1.545", "1.69"],
            ["351", "E", "ADN", "C4'", "10.688", "9.68", "3.523", "1.66"],
            ["351", "E", "ADN", "O4'", "9.714", "10.725", "3.81", "1.69"],
            ["351", "E", "ADN", "C3'", "9.973", "8.374", "3.903", "1.66"],
            ["351", "E", "ADN", "O3'", "10.879", "7.361", "4.304", "1.69"],
            ["351", "E", "ADN", "C2'", "9.115", "8.82", "5.059", "1.66"],
            ["351", "E", "ADN", "O2'", "9.887", "9.034", "6.232", "1.69"],
            ["351", "E", "ADN", "C1'", "8.625", "10.16", "4.5", "1.66"],
            ["351", "E", "ADN", "N1", "3.499", "10.104", "4.402", "1.97"],
            ["351", "E", "ADN", "C2", "4.376", "10.259", "5.387", "1.66"],
            ["351", "E", "ADN", "N3", "5.705", "10.249", "5.351", "1.97"],
            ["351", "E", "ADN", "C4", "6.136", "10.087", "4.094", "1.66"],
            ["351", "E", "ADN", "C5", "5.353", "9.952", "2.974", "1.66"],
            ["351", "E", "ADN", "C6", "3.957", "9.957", "3.146", "1.66"],
            ["351", "E", "ADN", "N6", "3.083", "9.826", "2.142", "1.97"],
            ["351", "E", "ADN", "N7", "6.146", "9.791", "1.843", "1.97"],
            ["351", "E", "ADN", "C8", "7.374", "9.872", "2.291", "1.66"],
            ["351", "E", "ADN", "N9", "7.444", "10.056", "3.646", "1.97"],
        ]
        result = read_pdb(os.path.join(UNIT_TESTS_DIR, "hetatm.pdb")).tolist()
        self.assertListEqual(result, expected)

    def test_skippable_entries(self):
        # Test skippable entries
        expected = [
            ["13", "A", "GLY", "N", "12.681", "37.302", "-25.211", "1.824"],
            ["13", "A", "GLY", "CA", "11.982", "37.996", "-26.241", "1.908"],
        ]
        result = read_pdb(os.path.join(UNIT_TESTS_DIR, "skippable.pdb")).tolist()
        self.assertListEqual(result, expected)


class TestReadXyz(unittest.TestCase):
    def test_xyz(self):
        expected = [
            ["1", "A", "UNK", "N", "-6.693", "-15.642", "-14.858", "1.97"],
            ["2", "A", "UNK", "C", "-6.73", "-14.62", "-15.897", "1.66"],
            ["3", "A", "UNK", "C", "-7.49", "-13.357", "-15.508", "1.66"],
            ["4", "A", "UNK", "O", "-7.388", "-12.338", "-16.197", "1.69"],
            ["5", "A", "UNK", "C", "-7.166", "-15.114", "-17.306", "1.66"],
            ["6", "A", "UNK", "C", "-7.485", "-16.623", "-17.424", "1.66"],
            ["7", "A", "UNK", "C", "-6.267", "-17.503", "-17.618", "1.66"],
            ["8", "A", "UNK", "O", "-5.439", "-17.684", "-16.719", "1.69"],
            ["9", "A", "UNK", "O", "-6.178", "-18.047", "-18.817", "1.69"],
        ]
        result = read_xyz(os.path.join(UNIT_TESTS_DIR, "xyz.xyz")).tolist()
        self.assertListEqual(result, expected)


class TestReadCavity(unittest.TestCase):
    def setUp(self):
        self.cavity = None

    def test__read_cavity(self):
        expected = [
            [0.0, 0.0, 0.0, 2.0],
            [0.0, 0.0, 0.6, 2.0],
            [0.0, 0.6, 0.0, 2.0],
            [0.6, 0.0, 0.0, 2.0],
            [0.0, 0.6, 0.6, 2.0],
            [0.6, 0.6, 0.0, 2.0],
            [0.6, 0.0, 0.6, 2.0],
            [0.6, 0.6, 0.6, 2.0],
        ]
        result = _read_cavity(os.path.join(UNIT_TESTS_DIR, "cavity.pdb")).tolist()
        self.assertListEqual(result, expected)

    def test_read_cavity(self):
        expected = [[[2, 2], [2, 2]], [[2, 2], [2, 2]]]
        result = read_cavity(
            os.path.join(UNIT_TESTS_DIR, "cavity.pdb"),
            os.path.join(UNIT_TESTS_DIR, "receptor.pdb"),
        )[:2, :2, :2].tolist()
        self.assertListEqual(result, expected)


class TestProcessBox(unittest.TestCase):
    def setUp(self):
        self.args = argparse.Namespace(
            box=False,
            probe_out=4.0,
            sincos=numpy.array([0.0, 1.0, 0.0, 1.0]),
            vertices=numpy.array(
                [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
            ),
        )

    def test_grid(self):
        expected = {
            "p1": [0.0, 0.0, 0.0],
            "p2": [10.0, 0.0, 0.0],
            "p3": [0.0, 10.0, 0.0],
            "p4": [0.0, 0.0, 10.0],
        }
        result = _process_box(self.args)
        self.assertDictEqual(result, expected)

    def test_custom_box(self):
        self.args.box = True
        expected = {
            "p1": [4.0, 4.0, 4.0],
            "p2": [6.0, 4.0, 4.0],
            "p3": [4.0, 6.0, 4.0],
            "p4": [4.0, 4.0, 6.0],
        }
        result = _process_box(self.args)
        self.assertDictEqual(result, expected)


class TestCalculateFrequencies(unittest.TestCase):
    def test_residues(self):
        residues = {
            "KAA": [
                ["49", "E", "LEU"],
                ["50", "E", "GLY"],
                ["51", "E", "THR"],
                ["52", "E", "GLY"],
                ["53", "E", "SER"],
            ]
        }
        expected = {
            "KAA": {
                "RESIDUES": {"GLY": 2, "LEU": 1, "SER": 1, "THR": 1},
                "CLASS": {"R1": 3, "R2": 0, "R3": 2, "R4": 0, "R5": 0, "RX": 0},
            }
        }
        result = calculate_frequencies(residues)
        self.assertDictEqual(expected, result)


if __name__ == "__main__":
    unittest.main()
