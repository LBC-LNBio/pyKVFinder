import unittest
import os
import pyKVFinder


def repeat(times):
    def repeatHelper(f):
        def callHelper(*args):
            for _ in range(0, times):
                f(*args)

        return callHelper

    return repeatHelper


class TestPackage(unittest.TestCase):
    def setUp(self):
        self.vdw = os.path.join(os.path.dirname(pyKVFinder.__file__), "data", "vdw.dat")
        self.pdb = os.path.join(
            os.path.dirname(pyKVFinder.__file__), "data", "tests", "1FMO.pdb"
        )
        self.xyz = os.path.join(
            os.path.dirname(pyKVFinder.__file__), "data", "tests", "1FMO.xyz"
        )
        self.ligand = os.path.join(
            os.path.dirname(pyKVFinder.__file__), "data", "tests", "ADN.pdb"
        )
        self.residues_box = os.path.join(
            os.path.dirname(pyKVFinder.__file__), "data", "tests", "residues-box.toml"
        )
        self.custom_box = os.path.join(
            os.path.dirname(pyKVFinder.__file__), "data", "tests", "custom-box.toml"
        )
        self.cavity = os.path.join(
            os.path.dirname(pyKVFinder.__file__),
            "data",
            "tests",
            "1FMO.KVFinder.output.pdb",
        )
        self.atomic = pyKVFinder.read_pdb(self.pdb)
        self.vertices = pyKVFinder.get_vertices(self.atomic)
        _, self.cavities = pyKVFinder.detect(self.atomic, self.vertices)
        self.surface, _, _ = pyKVFinder.spatial(self.cavities)

    def test_read_vdw(self):
        vdw = pyKVFinder.read_vdw(self.vdw)
        self.assertEqual(len(vdw), 36)
        self.assertListEqual(
            list(vdw.keys()),
            [
                "ALA",
                "ARG",
                "ASH",
                "ASN",
                "ASP",
                "CYM",
                "CYS",
                "CYX",
                "GLH",
                "GLN",
                "GLU",
                "GLY",
                "HID",
                "HIE",
                "HIP",
                "ILE",
                "LEU",
                "LYN",
                "LYS",
                "MET",
                "PHE",
                "PRO",
                "SER",
                "THR",
                "TRP",
                "TYR",
                "VAL",
                "HIS",
                "PTR",
                "SEP",
                "TPO",
                "H2D",
                "Y1P",
                "T1P",
                "S1P",
                "GEN",
            ],
        )

    def test_read_pdb(self):
        self.assertEqual(len(self.atomic), 2792)

    def test_read_xyz(self):
        atomic = pyKVFinder.read_xyz(self.xyz)
        self.assertEqual(len(atomic), 2792)

    def test_read_cavity(self):
        grid = pyKVFinder.read_cavity(self.cavity, self.pdb)
        self.assertEqual(grid.max() - 1 > 0, True)

    def test_get_vertices(self):
        # Expected
        expected = [
            [-19.911, -32.125, -30.806],
            [40.188, -32.125, -30.806],
            [-19.911, 43.446, -30.806],
            [-19.911, -32.125, 27.352],
        ]
        # Assert
        self.assertListEqual(self.vertices.tolist(), expected)

    def test_get_grid_from_file(self):
        # Residues box
        expected = [
            [-2.634, 2.843, -3.68],
            [24.422, 2.843, -3.68],
            [-2.634, 22.698, -3.68],
            [-2.634, 2.843, 16.204],
        ]
        vertices, atomic = pyKVFinder.get_vertices_from_file(
            self.residues_box, self.atomic
        )
        self.assertListEqual(vertices.tolist(), expected)
        self.assertEqual(len(atomic), 723)
        # Custom box
        expected = [
            [-0.89, 3.34, -2.41],
            [15.51, 3.34, -2.41],
            [-0.89, 14.74, -2.41],
            [-0.89, 3.34, 10.19],
        ]
        vertices, atomic = pyKVFinder.get_vertices_from_file(
            self.custom_box, self.atomic
        )
        self.assertListEqual(vertices.tolist(), expected)
        self.assertEqual(len(atomic), 388)

    def test_detect(self):
        self.assertEqual(self.cavities.max() - 1 > 0, True)

    @repeat(100)
    def test_cavity_tags_within_bounds(self):
        # Check if cavity tags are within expeted bounds [-1, ncav+1]
        # NOTE: A segmentation fault error has occurred in ~1/70 times.
        # So we run the test 100 times to make sure it's not happening.
        # Detection
        ncav, cavities = pyKVFinder.detect(self.atomic, self.vertices)
        self.assertEqual(((cavities >= -1) & (cavities <= ncav + 1)).all(), True)
        # Depth characterization
        depths, max_depth, avg_depth = pyKVFinder.depth(cavities)
        self.assertEqual(((cavities >= -1) & (cavities <= ncav + 1)).all(), True)
        # Constitutional characterization
        residues = pyKVFinder.constitutional(cavities, self.atomic, self.vertices)
        self.assertEqual(((cavities >= -1) & (cavities <= ncav + 1)).all(), True)
        # Surface characterization
        surface, volume, area = pyKVFinder.spatial(cavities)
        self.assertEqual(((cavities >= -1) & (cavities <= ncav + 1)).all(), True)
        # Hydropathy characterization
        pyKVFinder.hydropathy(surface, self.atomic, self.vertices)
        self.assertEqual(((cavities >= -1) & (cavities <= ncav + 1)).all(), True)

    def test_surface(self):
        self.assertEqual(self.surface.max() - 1 > 0, True)

    def test_depth(self):
        self.depths, max_depth, avg_depth = pyKVFinder.depth(self.cavities)
        self.assertEqual(self.depths.sum() > 0.0, True)
        # Apply selection
        depths, max_depth, avg_depth = pyKVFinder.depth(
            self.cavities, selection=["KAA"]
        )
        self.assertListEqual(list(max_depth.keys()), ["KAA"])

    def test_constitutional(self):
        # Full constitutional (interface residues, frequencies, bar charts)
        residues = pyKVFinder.constitutional(self.cavities, self.atomic, self.vertices)
        cavity_names = list(residues.keys())
        self.assertEqual(len(cavity_names) > 0, True)
        frequencies = pyKVFinder.calculate_frequencies(residues)
        cavity_names = list(residues.keys())
        self.assertEqual(len(cavity_names) > 0, True)
        pyKVFinder.plot_frequencies(
            frequencies,
            os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "barplots.pdf",
            ),
        )
        # Ignore backbone
        interface_residues = sum([len(r) for r in residues.values()])
        residues = pyKVFinder.constitutional(
            self.cavities, self.atomic, self.vertices, ignore_backbone=True
        )
        self.assertEqual(
            interface_residues > sum([len(r) for r in residues.values()]), True
        )

    def test_hydropathy(self):
        scales, avg_hydropathy = pyKVFinder.hydropathy(
            self.surface, self.atomic, self.vertices
        )
        cavity_names = list(avg_hydropathy.keys())
        self.assertEqual(len(cavity_names) > 0, True)
        self.assertEqual(scales.sum() != 0.0, True)

    def test_export(self):
        pyKVFinder.export(
            os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "cavities-test.pdb",
            ),
            self.cavities,
            self.surface,
            self.vertices,
        )

    def test_write_results(self):
        pyKVFinder.write_results(
            os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "results.toml",
            ),
            input=self.pdb,
            ligand=self.ligand,
            output=os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "cavities-test.pdb",
            ),
            output_hydropathy=None,
            volume=None,
            area=None,
            max_depth=None,
            avg_depth=None,
            avg_hydropathy=None,
            residues=None,
            frequencies=None,
        )


class TestPackageWorkflow(unittest.TestCase):
    def setUp(self):
        # Pdb
        self.pdb = os.path.join(
            os.path.dirname(pyKVFinder.__file__), "data", "tests", "1FMO.pdb"
        )
        # Full workflow
        self.results = pyKVFinder.run_workflow(
            self.pdb,
            include_depth=True,
            include_hydropathy=True,
            hydrophobicity_scale="EisenbergWeiss",
        )

    def test_standard_workflow(self):
        results = pyKVFinder.run_workflow(self.pdb)
        self.assertEqual(results.ncav > 0, True)

    def test_full_workflow(self):
        results = pyKVFinder.run_workflow(
            self.pdb,
            include_depth=True,
            include_hydropathy=True,
            hydrophobicity_scale="EisenbergWeiss",
        )
        self.assertEqual(results.ncav > 0, True)

    def test_ligand_mode(self):
        results = pyKVFinder.run_workflow(
            self.pdb,
            os.path.join(
                os.path.dirname(pyKVFinder.__file__), "data", "tests", "ADN.pdb"
            ),
        )
        self.assertEqual(results.ncav > 0, True)

    def test_residues_box(self):
        results = pyKVFinder.run_workflow(
            self.pdb,
            box=os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "residues-box.toml",
            ),
        )
        self.assertEqual(results.ncav > 0, True)

    def test_custom_box(self):
        results = pyKVFinder.run_workflow(
            self.pdb,
            box=os.path.join(
                os.path.dirname(pyKVFinder.__file__), "data", "tests", "custom-box.toml"
            ),
        )
        self.assertEqual(results.ncav, 1)

    def test_pyKVFinderResults_methods(self):
        # export
        self.results.export(
            output=os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "cavities.pdb",
            ),
            output_hydropathy=os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "hydropathy.pdb",
            ),
        )
        # write
        self.results.write(
            fn=os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "results.toml",
            ),
            output=os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "cavities.pdb",
            ),
            output_hydropathy=os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "hydropathy.pdb",
            ),
        )
        # plot_frequencies
        self.results.plot_frequencies(
            pdf=os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "barplots.pdf",
            )
        )
        # export_all
        self.results.export_all(
            fn=os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "results.toml",
            ),
            output=os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "cavities.pdb",
            ),
            output_hydropathy=os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "hydropathy.pdb",
            ),
            include_frequencies_pdf=True,
            pdf=os.path.join(
                os.path.dirname(pyKVFinder.__file__),
                "data",
                "tests",
                "output",
                "barplots.pdf",
            ),
        )
