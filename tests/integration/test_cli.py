import os
import unittest
import subprocess


class TestCLI(unittest.TestCase):
    def test_standard_mode(self):
        # Run pyKVFinder CLI with standard mode
        result = subprocess.run(
            [
                "pyKVFinder",
                "pyKVFinder/data/tests/1FMO.pdb",
                "-O",
                "pyKVFinder/data/tests/output",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0)

    def test_box_mode(self):
        # Run pyKVFinder CLI with box mode
        # - custom_box.toml
        result = subprocess.run(
            [
                "pyKVFinder",
                "pyKVFinder/data/tests/1FMO.pdb",
                "-O",
                "pyKVFinder/data/tests/output",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0)
        # - residues_box.toml
        result = subprocess.run(
            [
                "pyKVFinder",
                "pyKVFinder/data/tests/1FMO.pdb",
                "-O",
                "pyKVFinder/data/tests/output",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0)

    def test_lig_mode(self):
        pass

    def test_hydropathy_mode(self):
        pass

    def test_depth_mode(self):
        pass

    def test_bad_args(self):
        pass
