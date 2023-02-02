import os
import subprocess
import unittest

import pyKVFinder


class TestCLI(unittest.TestCase):
    def test_standard_mode(self):
        # Run pyKVFinder CLI with standard mode
        # $ pyKVFinder <.pdb>
        result = subprocess.run(
            [
                "pyKVFinder",
                f"{os.path.dirname(pyKVFinder.__file__)}/data/tests/1FMO.pdb",
                "-O",
                f"{os.path.dirname(pyKVFinder.__file__)}/data/tests/output",
            ],
            capture_output=True,
        )
        print(result)
        self.assertEqual(result.returncode, 0)

    def test_box_mode(self):
        # Run pyKVFinder CLI with box mode
        # $ pyKVFinder <.pdb> -B <.toml>
        # - custom_box.toml
        # >>> [box]
        # >>> p1 = [x1, y1, z1]
        # >>> p2 = [x2, y2, z2]
        # >>> p3 = [x3, y3, z3]
        # >>> p4 = [x4, y4, z4]
        result = subprocess.run(
            [
                "pyKVFinder",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/1FMO.pdb",
                "-O",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/output",
                "-B",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/custom-box.toml",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0)
        # - residues_box.toml
        # >>> [box]
        # >>> residues = [ ["resnum", "chain", "resname",], ["resnum", "chain",
        # ... "resname"], ]
        # >>> padding =  3.5
        result = subprocess.run(
            [
                "pyKVFinder",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/1FMO.pdb",
                "-O",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/output",
                "-B",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/residues-box.toml",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0)

    def test_lig_mode(self):
        # Run pyKVFinder CLI with ligand mode
        # $ pyKVFinder <.pdb> -L <.pdb>
        result = subprocess.run(
            [
                "pyKVFinder",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/1FMO.pdb",
                "-O",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/output",
                "-L",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/ADN.pdb",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0)

    def test_hydropathy_mode(self):
        # Run pyKVFinder CLI with hydropathy mode
        # $ pyKVFinder <.pdb> --hydropathy
        result = subprocess.run(
            [
                "pyKVFinder",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/1FMO.pdb",
                "-O",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/output",
                "--hydropathy",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0)

    def test_depth_mode(self):
        # Run pyKVFinder CLI with depth mode
        # $ pyKVFinder <.pdb> --depth
        result = subprocess.run(
            [
                "pyKVFinder",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/1FMO.pdb",
                "-O",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/output",
                "--depth",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0)

    def test_bad_args(self):
        # non-existing <.pdb>
        # $ pyKVFinder non-existing.pdb
        result = subprocess.run(
            [
                "pyKVFinder",
                "non-existing.pdb",
                "-O",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/output",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 1)  # FileNotFoundError
        # bad float (step, probe_in, probe_out, removal_distance, volume_cutoff, ligand_cutoff)
        result = subprocess.run(
            ["pyKVFinder", f"{os.path.dirname(pyKVFinder.__file__)}/tests/1FMO.pdb", "-o", "string"],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 2)  # argparse.ArgumentTypeError
        # bad surface
        result = subprocess.run(
            [
                "pyKVFinder",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/1FMO.pdb",
                "-S",
                "A",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 2)  # argparse.ArgumentTypeError
        # non-existing ligand <.pdb>
        result = subprocess.run(
            [
                "pyKVFinder",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/1FMO.pdb",
                "-O",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/output",
                "-L",
                "non-existing.pdb",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 1)  # FileNotFoundError
        # bad hydropathy file (not a TOML file)
        result = subprocess.run(
            [
                "pyKVFinder",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/1FMO.pdb",
                "--hydropathy",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/1FMO.pdb",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 2)  # argparse.ArgumentTypeError
        # non-existing dictionary file
        result = subprocess.run(
            [
                "pyKVFinder",
                f"{os.path.dirname(pyKVFinder.__file__)}/tests/1FMO.pdb",
                "--hydropathy",
                "non-existing-dictionary.dat",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 2)  # argparse.ArgumentTypeError
