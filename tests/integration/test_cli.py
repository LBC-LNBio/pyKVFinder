import argparse
import io
import os
import unittest
from unittest import mock

import numpy
import toml

import pyKVFinder

DATADIR = os.path.join(os.path.dirname(pyKVFinder.__file__), "data")


class TestCLI(unittest.TestCase):
    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy=False,
            box=None,
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_standard_mode(self, _):
        # Run pyKVFinder CLI with standard mode
        # $ pyKVFinder <.pdb>
        self.assertEqual(pyKVFinder.main.cli(), 0)

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.xyz"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy=False,
            box=None,
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_standard_mode_with_xyz(self, _):
        # Run pyKVFinder CLI with standard mode
        # $ pyKVFinder <.pdb>
        self.assertEqual(pyKVFinder.main.cli(), 0)

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy=False,
            box=os.path.join(DATADIR, "tests", "custom-box.toml"),
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_custom_box_mode(self, _):
        # Run pyKVFinder CLI with box mode
        # $ pyKVFinder <.pdb> -B <.toml>
        # - custom_box.toml
        # >>> [box]
        # >>> p1 = [x1, y1, z1]
        # >>> p2 = [x2, y2, z2]
        # >>> p3 = [x3, y3, z3]
        # >>> p4 = [x4, y4, z4]
        self.assertEqual(pyKVFinder.main.cli(), 0)

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy=False,
            box=os.path.join(DATADIR, "tests", "residues-box.toml"),
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_residues_box_mode(self, _):
        # - residues_box.toml
        # >>> [box]
        # >>> residues = [ ["resnum", "chain", "resname",], ["resnum", "chain",
        # ... "resname"], ]
        # >>> padding =  3.5
        self.assertEqual(pyKVFinder.main.cli(), 0)

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy=False,
            box=None,
            ligand=os.path.join(DATADIR, "tests", "ADN.pdb"),
            ligand_cutoff=5.0,
        ),
    )
    def test_ligand_mode(self, _):
        # Run pyKVFinder CLI with ligand mode
        # $ pyKVFinder <.pdb> -L <.pdb>
        self.assertEqual(pyKVFinder.main.cli(), 0)

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy=False,
            box=None,
            ligand=os.path.join(DATADIR, "tests", "ADN.xyz"),
            ligand_cutoff=5.0,
        ),
    )
    def test_ligand_mode_with_xyz(self, _):
        # Run pyKVFinder CLI with ligand mode
        # $ pyKVFinder <.pdb> -L <.pdb>
        self.assertEqual(pyKVFinder.main.cli(), 0)

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy="EisenbergWeiss",
            box=None,
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_hydropathy_mode(self, _):
        # Run pyKVFinder CLI with hydropathy mode
        # $ pyKVFinder <.pdb> --hydropathy
        self.assertEqual(pyKVFinder.main.cli(), 0)

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=True,
            plot_frequencies=True,
            hydropathy=False,
            box=None,
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_plot_frequencies(self, _):
        # Run pyKVFinder CLI with depth mode
        # $ pyKVFinder <.pdb> --depth
        self.assertEqual(pyKVFinder.main.cli(), 0)

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=1000000,  # Force ncav = 0
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=True,
            plot_frequencies=True,
            hydropathy=False,
            box=None,
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_no_cavities_detected(self, _):
        # Run pyKVFinder CLI with depth mode
        # $ pyKVFinder <.pdb> --depth
        with mock.patch("sys.stdout", new_callable=io.StringIO) as stdout:
            self.assertEqual(pyKVFinder.main.cli(), 0)
        self.assertEqual(stdout.getvalue().split("\n")[1], "> No cavities detected!")

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=True,
            plot_frequencies=False,
            hydropathy=False,
            box=None,
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_depth_mode(self, _):
        # Run pyKVFinder CLI with depth mode
        # $ pyKVFinder <.pdb> --depth
        self.assertEqual(pyKVFinder.main.cli(), 0)

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input="non-existing.pdb",
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy=False,
            box=None,
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_non_existing_receptor(self, _):
        # non-existing <.pdb>
        # $ pyKVFinder non-existing.pdb
        self.assertRaises(FileNotFoundError, pyKVFinder.main.cli)  # FileNotFoundError

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out="string",
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy=False,
            box=None,
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_invalid_float(self, _):
        # bad float (step, probe_in, probe_out, removal_distance, volume_cutoff, ligand_cutoff)
        # $ pyKVFinder -o string
        self.assertRaises(TypeError, pyKVFinder.main.cli)  # TypeError

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="A",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy=False,
            box=None,
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_invalid_surface(self, _):
        # bad surface
        # $ pyKVFinder -S A
        # self.assertEqual(pyKVFinder.main.cli(), 2)  # argparse.ArgumentTypeError
        self.assertRaises(ValueError, pyKVFinder.main.cli)  # ValueError

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy=False,
            box=None,
            ligand="non-existing.pdb",
            ligand_cutoff=5.0,
        ),
    )
    def test_non_existing_ligand(self, _):
        # bad surface
        # $ pyKVFinder -L non-existing.pdb
        self.assertRaises(FileNotFoundError, pyKVFinder.main.cli)  # FileNotFoundError

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy="non-existing-dictionary.dat",
            box=None,
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_non_existing_hydropathy_file(self, _):
        # bad surface
        # $ pyKVFinder --hydopathy non-existing-dictionary.dat
        self.assertRaises(FileNotFoundError, pyKVFinder.main.cli)  # FileNotFoundError

    @mock.patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            input=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            verbose=False,
            base_name=None,
            output_directory=os.path.join(DATADIR, "tests", "output"),
            model=None,
            nthreads=11,
            dictionary=os.path.join(DATADIR, "vdw.dat"),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            removal_distance=2.4,
            surface="SES",
            ignore_backbone=False,
            depth=False,
            plot_frequencies=False,
            hydropathy=os.path.join(DATADIR, "tests", "1FMO.pdb"),
            box=None,
            ligand=None,
            ligand_cutoff=5.0,
        ),
    )
    def test_invalid_hydropathy_file(self, _):
        # bad surface
        # $ pyKVFinder --hydopathy <.pdb>
        self.assertRaises(
            toml.decoder.TomlDecodeError, pyKVFinder.main.cli
        )  # toml.decoder.TomlDecodeError
