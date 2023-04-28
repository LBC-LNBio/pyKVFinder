import argparse
import io
import os
import unittest
from unittest import mock

import numpy

from pyKVFinder.utils import (
    _process_box,
    _process_pdb_line,
    _read_cavity,
    calculate_frequencies,
    plot_frequencies,
    read_cavity,
    read_pdb,
    read_vdw,
    read_xyz,
    write_results,
)

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")


class TestReadVdw(unittest.TestCase):
    def test_valid_vdw(self):
        # Check valid files formats
        expected = {
            "ALA": {"N": 1.824, "C": 1.908, "H": 0.6, "O": 1.6612},
            "GEN": {"N": 1.824, "C": 1.908, "H": 0.6, "O": 1.6612},
        }
        result = read_vdw(os.path.join(FIXTURES, "vdw1.dat"))
        self.assertEqual(result, expected)

    def test_invalid_vdw(self):
        # Check invalid files formats
        self.assertRaises(ValueError, read_vdw, os.path.join(FIXTURES, "vdw2.dat"))
        self.assertRaises(ValueError, read_vdw, os.path.join(FIXTURES, "vdw3.dat"))

    def test_wrong_fn_format(self):
        # Check wrong fn formats
        for fn in [1, 1.0, [1], {"vdw": 1}, numpy.ones(1)]:
            self.assertRaises(TypeError, read_vdw, fn)

    def test_default_input(self):
        self.assertIsInstance(read_vdw(None), dict)


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
            [1, "A", "UNK", "C", 0.000, 0.000, 0.000, 1.66],
            [1, "A", "UNK", "N", 0.000, 0.000, 0.000, 1.824],
            [1, "A", "UNK", "O", 0.000, 0.000, 0.000, 1.69],
            [1, "A", "UNK", "H", 0.000, 0.000, 0.000, 0.91],
            [1, "A", "UNK", "N9", 0.000, 0.000, 0.000, 1.824],
            [1, "A", "UNK", "C5'", 0.000, 0.000, 0.000, 1.66],
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
        result = read_pdb(os.path.join(FIXTURES, "atom.pdb")).tolist()
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
        result = read_pdb(os.path.join(FIXTURES, "hetatm.pdb")).tolist()
        self.assertListEqual(result, expected)

    def test_skippable_entries(self):
        # Test skippable entries
        expected = [
            ["13", "A", "GLY", "N", "12.681", "37.302", "-25.211", "1.824"],
            ["13", "A", "GLY", "CA", "11.982", "37.996", "-26.241", "1.908"],
        ]
        result = read_pdb(os.path.join(FIXTURES, "skippable.pdb")).tolist()
        self.assertListEqual(result, expected)

    def test_nmr_models(self):
        # Test MODEL entries
        expected = [
            ["13", "E", "GLU", "N", "-1.111", "-15.642", "-14.858", "1.824"],
            ["13", "E", "GLU", "N", "-2.222", "-15.642", "-14.858", "1.824"],
            ["13", "E", "GLU", "N", "-3.333", "-15.642", "-14.858", "1.824"],
            ["13", "E", "GLU", "N", "-4.444", "-15.642", "-14.858", "1.824"],
            ["13", "E", "GLU", "N", "-5.555", "-15.642", "-14.858", "1.824"],
        ]
        result = read_pdb(os.path.join(FIXTURES, "nmr.pdb")).tolist()
        self.assertListEqual(result, expected)

    def test_one_nmr_model(self):
        # Model number
        model = numpy.random.randint(1, 6)
        # Test one model in NMR entry
        expected = [
            ["13", "E", "GLU", "N", "-1.111", "-15.642", "-14.858", "1.824"],
            ["13", "E", "GLU", "N", "-2.222", "-15.642", "-14.858", "1.824"],
            ["13", "E", "GLU", "N", "-3.333", "-15.642", "-14.858", "1.824"],
            ["13", "E", "GLU", "N", "-4.444", "-15.642", "-14.858", "1.824"],
            ["13", "E", "GLU", "N", "-5.555", "-15.642", "-14.858", "1.824"],
        ]
        result = read_pdb(os.path.join(FIXTURES, "nmr.pdb"), model=model).tolist()
        self.assertListEqual(result, [expected[model - 1]])

    def test_wrong_fn_format(self):
        # Check wrong fn formats
        for fn in [1, 1.0, [1], {"vdw": 1}, numpy.ones(1)]:
            self.assertRaises(TypeError, read_pdb, fn, None, None)

    def test_wrong_model_format(self):
        # Check wrong model formats
        for model in [1.0, [1], {"model": 1}, numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError, read_pdb, os.path.join(FIXTURES, "nmr.pdb"), None, model
            )


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
        result = read_xyz(os.path.join(FIXTURES, "xyz.xyz")).tolist()
        self.assertListEqual(result, expected)

    def test_wrong_fn_format(self):
        # Check wrong fn formats
        for fn in [1, 1.0, [1], {"xyz": 1}, numpy.ones(1)]:
            self.assertRaises(TypeError, read_xyz, fn, None)


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
        result = _read_cavity(os.path.join(FIXTURES, "cavity.pdb")).tolist()
        self.assertListEqual(result, expected)

    def test_read_cavity(self):
        expected = [[[2, 2], [2, 2]], [[2, 2], [2, 2]]]
        result = read_cavity(
            os.path.join(FIXTURES, "cavity.pdb"),
            os.path.join(FIXTURES, "receptor.pdb"),
        )[:2, :2, :2].tolist()
        self.assertListEqual(result, expected)

    def test_xyz_receptor(self):
        expected = [[[2, 2], [2, 2]], [[2, 2], [2, 2]]]
        result = read_cavity(
            os.path.join(FIXTURES, "cavity.pdb"),
            os.path.join(FIXTURES, "receptor.xyz"),
        )[:2, :2, :2].tolist()
        self.assertListEqual(result, expected)

    def test_parameters_as_integers(self):
        for kwargs in [{"step": 1, "probe_in": 1, "probe_out": 4}]:
            expected = [[[2, 2], [2, 2]], [[2, 2], [2, 2]]]
            result = read_cavity(
                os.path.join(FIXTURES, "cavity.pdb"),
                os.path.join(FIXTURES, "receptor.pdb"),
                **kwargs,
            )[:2, :2, :2].tolist()
            self.assertListEqual(result, expected)

    def test_wrong_cavity_format(self):
        # Check wrong cavity formats
        for cavity in [1, 1.0, [1], {"cavity": 1}, numpy.ones(1)]:
            self.assertRaises(
                TypeError, read_cavity, cavity, os.path.join(FIXTURES, "receptor.pdb")
            )

    def test_wrong_receptor_format(self):
        # Check wrong receptor formats
        for receptor in [1, 1.0, [1], {"receptor": 1}, numpy.ones(1)]:
            self.assertRaises(
                TypeError, read_cavity, os.path.join(FIXTURES, "cavity.pdb"), receptor
            )

    def test_wrong_receptor_extension(self):
        # Check wrong receptor extension
        for receptor in ["receptor.mol", "receptor.cif"]:
            self.assertRaises(
                TypeError, read_cavity, os.path.join(FIXTURES, "cavity.pdb"), receptor
            )

    def test_wrong_step_format(self):
        # Check wrong step format
        for step in [[1], {"step": 1}, numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                read_cavity,
                os.path.join(FIXTURES, "cavity.pdb"),
                os.path.join(FIXTURES, "receptor.pdb"),
                step=step,
            )

    def test_invalid_step(self):
        # Check invalid step
        for step in [-1.0, 0.0]:
            self.assertRaises(
                ValueError,
                read_cavity,
                os.path.join(FIXTURES, "cavity.pdb"),
                os.path.join(FIXTURES, "receptor.pdb"),
                step=step,
            )

    def test_wrong_probe_in_format(self):
        # Check wrong probe in format
        for probe_in in [[1], {"step": 1}, numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                read_cavity,
                os.path.join(FIXTURES, "cavity.pdb"),
                os.path.join(FIXTURES, "receptor.pdb"),
                probe_in=probe_in,
            )

    def test_invalid_probe_in(self):
        # Check invalid probe in
        self.assertRaises(
            ValueError,
            read_cavity,
            os.path.join(FIXTURES, "cavity.pdb"),
            os.path.join(FIXTURES, "receptor.pdb"),
            probe_in=-1,
        )

    def test_wrong_probe_out_format(self):
        # Check wrong probe out format
        for probe_out in [[1], {"step": 1}, numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                read_cavity,
                os.path.join(FIXTURES, "cavity.pdb"),
                os.path.join(FIXTURES, "receptor.pdb"),
                probe_out=probe_out,
            )

    def test_invalid_probe_out(self):
        # Check invalid probe out
        self.assertRaises(
            ValueError,
            read_cavity,
            os.path.join(FIXTURES, "cavity.pdb"),
            os.path.join(FIXTURES, "receptor.pdb"),
            probe_out=-1,
        )

    def test_wrong_nthreads_format(self):
        # Check wrong nthreads format
        for nthreads in [1.0, [1], {"nthreads": 1}, numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                read_cavity,
                os.path.join(FIXTURES, "cavity.pdb"),
                os.path.join(FIXTURES, "receptor.pdb"),
                nthreads=nthreads,
            )

    def test_invalid_nthreads(self):
        # Check invalid nthreads
        for nthreads in [-1, 0]:
            self.assertRaises(
                ValueError,
                read_cavity,
                os.path.join(FIXTURES, "cavity.pdb"),
                os.path.join(FIXTURES, "receptor.pdb"),
                nthreads=nthreads,
            )

    def test_invalid_probe_pair(self):
        # Check probe in greater than probe out
        self.assertRaises(
            ValueError,
            read_cavity,
            os.path.join(FIXTURES, "cavity.pdb"),
            os.path.join(FIXTURES, "receptor.pdb"),
            probe_in=4,
            probe_out=1.4,
        )

    def test_wrong_surface_format(self):
        # Check wrong surface format
        for surface in [1.0, [1], {"nthreads": 1}, numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                read_cavity,
                os.path.join(FIXTURES, "cavity.pdb"),
                os.path.join(FIXTURES, "receptor.pdb"),
                surface=surface,
            )

    def test_invalid_surface(self):
        # Check probe in greater than probe out
        self.assertRaises(
            ValueError,
            read_cavity,
            os.path.join(FIXTURES, "cavity.pdb"),
            os.path.join(FIXTURES, "receptor.pdb"),
            surface="invalid",
        )

    def test_wrong_verbose_format(self):
        # Check wrong verbose format
        for verbose in [1.0, [1], {"nthreads": 1}, numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                read_cavity,
                os.path.join(FIXTURES, "cavity.pdb"),
                os.path.join(FIXTURES, "receptor.pdb"),
                verbose=verbose,
            )


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


class TestFrequencies(unittest.TestCase):
    def setUp(self):
        self.residues = {
            "KAA": [
                ["49", "E", "LEU"],
                ["50", "E", "GLY"],
                ["51", "E", "THR"],
                ["52", "E", "GLY"],
                ["53", "E", "SER"],
                ["54", "E", "HYP"],
            ]
        }

    def test_calculate_frequencies(self):
        frequencies = calculate_frequencies(self.residues)
        expected = {
            "KAA": {
                "RESIDUES": {"GLY": 2, "LEU": 1, "SER": 1, "THR": 1, "HYP": 1},
                "CLASS": {"R1": 3, "R2": 0, "R3": 2, "R4": 0, "R5": 0, "RX": 1},
            }
        }
        self.assertDictEqual(expected, frequencies)

    def test_plot_frequencies(self):
        frequencies = calculate_frequencies(self.residues)
        plot_frequencies(frequencies, "tests/frequencies.png")
        self.assertEqual(os.path.exists("tests/frequencies.png"), True)
        os.remove("tests/frequencies.png")

    def test_plot_frequencies_wrong_fn_format(self):
        frequencies = calculate_frequencies(self.residues)
        for fn in [1, 1.0, [1], {"fn": 1}, numpy.ones(1)]:
            self.assertRaises(TypeError, plot_frequencies, frequencies, fn)


class TestWriteResults(unittest.TestCase):
    def test_positional_arguments(self):
        write_results(
            "tests/results.toml", "input.pdb", "ligand.pdb", "output.pdb", step=0.1
        )
        expected = f'# pyKVFinder results\n\n[FILES]\nINPUT = "{os.path.abspath("input.pdb")}"\nLIGAND = "{os.path.abspath("ligand.pdb")}"\nOUTPUT = "{os.path.abspath("output.pdb")}"\n\n[PARAMETERS]\nSTEP = 0.1\n'
        with open("tests/results.toml", "r") as f:
            self.assertEqual(f.read(), expected)
        os.remove("tests/results.toml")

    def test_volume(self):
        write_results("tests/volume.toml", "input.pdb", None, None, volume={"KAA": 100})
        expected = f'# pyKVFinder results\n\n[FILES]\nINPUT = "{os.path.abspath("input.pdb")}"\n\n[PARAMETERS]\nSTEP = 0.6\n\n[RESULTS.VOLUME]\nKAA = 100\n'
        with open("tests/volume.toml", "r") as f:
            self.assertEqual(f.read(), expected)
        os.remove("tests/volume.toml")

    def test_area(self):
        write_results("tests/area.toml", "input.pdb", None, None, area={"KAA": 100})
        expected = f'# pyKVFinder results\n\n[FILES]\nINPUT = "{os.path.abspath("input.pdb")}"\n\n[PARAMETERS]\nSTEP = 0.6\n\n[RESULTS.AREA]\nKAA = 100\n'
        with open("tests/area.toml", "r") as f:
            self.assertEqual(f.read(), expected)
        os.remove("tests/area.toml")

    def test_max_depth(self):
        write_results(
            "tests/max_depth.toml", "input.pdb", None, None, max_depth={"KAA": 3.00}
        )
        expected = f'# pyKVFinder results\n\n[FILES]\nINPUT = "{os.path.abspath("input.pdb")}"\n\n[PARAMETERS]\nSTEP = 0.6\n\n[RESULTS.MAX_DEPTH]\nKAA = 3.0\n'
        with open("tests/max_depth.toml", "r") as f:
            self.assertEqual(f.read(), expected)
        os.remove("tests/max_depth.toml")

    def test_avg_depth(self):
        write_results(
            "tests/avg_depth.toml", "input.pdb", None, None, avg_depth={"KAA": 2.87}
        )
        expected = f'# pyKVFinder results\n\n[FILES]\nINPUT = "{os.path.abspath("input.pdb")}"\n\n[PARAMETERS]\nSTEP = 0.6\n\n[RESULTS.AVG_DEPTH]\nKAA = 2.87\n'
        with open("tests/avg_depth.toml", "r") as f:
            self.assertEqual(f.read(), expected)
        os.remove("tests/avg_depth.toml")

    def test_avg_hydropathy(self):
        write_results(
            "tests/avg_hydropathy.toml",
            "input.pdb",
            None,
            None,
            avg_hydropathy={"KAA": -0.81},
        )
        expected = f'# pyKVFinder results\n\n[FILES]\nINPUT = "{os.path.abspath("input.pdb")}"\n\n[PARAMETERS]\nSTEP = 0.6\n\n[RESULTS.AVG_HYDROPATHY]\nKAA = -0.81\n'
        with open("tests/avg_hydropathy.toml", "r") as f:
            self.assertEqual(f.read(), expected)
        os.remove("tests/avg_hydropathy.toml")

    def test_residues(self):
        write_results(
            "tests/residues.toml",
            "input.pdb",
            None,
            None,
            residues={"KAA": [["49", "E", "LEU"]]},
        )
        expected = f'# pyKVFinder results\n\n[FILES]\nINPUT = "{os.path.abspath("input.pdb")}"\n\n[PARAMETERS]\nSTEP = 0.6\n\n[RESULTS.RESIDUES]\nKAA = [ [ "49", "E", "LEU",],]\n'
        with open("tests/residues.toml", "r") as f:
            self.assertEqual(f.read(), expected)
        os.remove("tests/residues.toml")

    def test_frequencies(self):
        write_results(
            "tests/frequencies.toml",
            "input.pdb",
            None,
            None,
            frequencies={
                "KAA": {
                    "RESIDUES": {"LEU": 1},
                    "CLASS": {"R1": 1, "R2": 0, "R3": 0, "R4": 0, "R5": 0, "RX": 0},
                }
            },
        )
        expected = f'# pyKVFinder results\n\n[FILES]\nINPUT = "{os.path.abspath("input.pdb")}"\n\n[PARAMETERS]\nSTEP = 0.6\n\n[RESULTS.FREQUENCY.KAA.RESIDUES]\nLEU = 1\n\n[RESULTS.FREQUENCY.KAA.CLASS]\nR1 = 1\nR2 = 0\nR3 = 0\nR4 = 0\nR5 = 0\nRX = 0\n'
        with open("tests/frequencies.toml", "r") as f:
            self.assertEqual(f.read(), expected)
        os.remove("tests/frequencies.toml")

    def test_wrong_fn_format(self):
        # Check wrong fn format
        for fn in [1.0, [1], {"nthreads": 1}, numpy.ones(1)]:
            self.assertRaises(TypeError, write_results, fn, "input.pdb", None, None)

    def test_wrong_input_format(self):
        # Check wrong input format
        for input in [1.0, [1], {"nthreads": 1}, numpy.ones(1)]:
            self.assertRaises(
                TypeError, write_results, "results.toml", input, None, None
            )

    def test_wrong_ligand_format(self):
        # Check wrong ligand format
        for ligand in [1.0, [1], {"nthreads": 1}, numpy.ones(1)]:
            self.assertRaises(
                TypeError, write_results, "results.toml", "input.pdb", ligand, None
            )

    def test_wrong_output_format(self):
        # Check wrong output format
        for output in [1.0, [1], {"nthreads": 1}, numpy.ones(1)]:
            self.assertRaises(
                TypeError, write_results, "results.toml", "input.pdb", None, output
            )

    def test_wrong_output_hydropathy_format(self):
        # Check wrong output_hydropathy format
        for output_hydropathy in [1.0, [1], {"nthreads": 1}, numpy.ones(1)]:
            self.assertRaises(
                TypeError,
                write_results,
                "results.toml",
                "input.pdb",
                None,
                None,
                output_hydropathy=output_hydropathy,
            )

    def test_wrong_volume_format(self):
        # Check wrong volume format
        for volume in [1.0, [1], numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                write_results,
                "results.toml",
                "input.pdb",
                None,
                None,
                volume=volume,
            )

    def test_wrong_area_format(self):
        # Check wrong area format
        for area in [1.0, [1], numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                write_results,
                "results.toml",
                "input.pdb",
                None,
                None,
                area=area,
            )

    def test_wrong_max_depth_format(self):
        # Check wrong max_depth format
        for max_depth in [1.0, [1], numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                write_results,
                "results.toml",
                "input.pdb",
                None,
                None,
                max_depth=max_depth,
            )

    def test_wrong_avg_depth_format(self):
        # Check wrong avg_depth format
        for avg_depth in [1.0, [1], numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                write_results,
                "results.toml",
                "input.pdb",
                None,
                None,
                avg_depth=avg_depth,
            )

    def test_wrong_avg_hydropathy_format(self):
        # Check wrong avg_hydropathy format
        for avg_hydropathy in [1.0, [1], numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                write_results,
                "results.toml",
                "input.pdb",
                None,
                None,
                avg_hydropathy=avg_hydropathy,
            )

    def test_wrong_residues_format(self):
        # Check wrong residues format
        for residues in [1.0, [1], numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                write_results,
                "results.toml",
                "input.pdb",
                None,
                None,
                residues=residues,
            )

    def test_wrong_frequencies_format(self):
        # Check wrong frequencies format
        for frequencies in [1.0, [1], numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                write_results,
                "results.toml",
                "input.pdb",
                None,
                None,
                frequencies=frequencies,
            )

    def test_wrong_step_format(self):
        # Check wrong step format
        for step in [{"step": 1}, [1], numpy.ones(1), "1"]:
            self.assertRaises(
                TypeError,
                write_results,
                "results.toml",
                "input.pdb",
                None,
                None,
                step=step,
            )

    def test_invalid_step(self):
        # Check invalid step
        for step in [-1.0, 0.0]:
            self.assertRaises(
                ValueError,
                write_results,
                "results.toml",
                "input.pdb",
                None,
                None,
                step=step,
            )

    def test_step_as_integer(self):
        # Check step as integer
        write_results("tests/step.toml", "input.pdb", None, None, step=1)
        expected = f'# pyKVFinder results\n\n[FILES]\nINPUT = "{os.path.abspath("input.pdb")}"\n\n[PARAMETERS]\nSTEP = 1.0\n'
        with open("tests/step.toml", "r") as f:
            self.assertEqual(f.read(), expected)
        os.remove("tests/step.toml")


if __name__ == "__main__":
    unittest.main()
