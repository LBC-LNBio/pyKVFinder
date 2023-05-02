import os
import unittest

import pyKVFinder


def repeat(times):
    def repeatHelper(f):
        def callHelper(*args):
            for _ in range(0, times):
                f(*args)

        return callHelper

    return repeatHelper


DATADIR = os.path.join(os.path.dirname(pyKVFinder.__file__), "data")


class TestPackage(unittest.TestCase):
    def setUp(self):
        self.vdw = os.path.join(DATADIR, "vdw.dat")
        self.pdb = os.path.join(DATADIR, "tests", "1FMO.pdb")
        self.xyz = os.path.join(DATADIR, "tests", "1FMO.xyz")
        self.ligand = os.path.join(DATADIR, "tests", "ADN.pdb")
        self.residues_box = os.path.join(DATADIR, "tests", "residues-box.toml")
        self.custom_box = os.path.join(DATADIR, "tests", "custom-box.toml")
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

    @repeat(1)
    def test_cavity_tags_within_bounds(self):
        # Check if cavity tags are within expeted bounds [-1, ncav+1]
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
                DATADIR,
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
                DATADIR,
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
                DATADIR,
                "tests",
                "output",
                "results.toml",
            ),
            input=self.pdb,
            ligand=self.ligand,
            output=os.path.join(
                DATADIR,
                "tests",
                "output",
                "cavities-test.pdb",
            ),
            volume=None,
            area=None,
            max_depth=None,
            avg_depth=None,
            avg_hydropathy=None,
            residues=None,
            frequencies=None,
        )


class TestMolecule(unittest.TestCase):
    def setUp(self):
        self.molecule = pyKVFinder.Molecule(os.path.join(DATADIR, "tests", "ClO4.pdb"))

    def test_vdw(self):
        self.molecule.vdw(0.1)
        self.assertEqual((self.molecule.grid == 0).sum() > 0, True)

    def test_vdw_volume(self):
        self.molecule.vdw(0.1)
        self.assertAlmostEqual(self.molecule.volume(), 83.64, 2)

    def test_ses(self):
        self.molecule.surface(0.1, 1.4, "SES")
        self.assertEqual((self.molecule.grid == 0).sum() > 0, True)

    def test_ses_volume(self):
        self.molecule.surface(0.1, 1.4, "SES")
        self.assertAlmostEqual(self.molecule.volume(), 90.80, 2)

    def test_sas(self):
        self.molecule.surface(0.1, 1.4, "SAS")
        self.assertEqual((self.molecule.grid == 0).sum() > 0, True)

    def test_sas_volume(self):
        self.molecule.surface(0.1, 1.4, "SAS")
        self.assertAlmostEqual(self.molecule.volume(), 340.28, 2)

    def test_export(self):
        self.molecule.vdw(0.1)
        self.molecule.export(os.path.join(DATADIR, "tests", "output", "molecule.pdb"))


class TestPackageWorkflow(unittest.TestCase):
    def setUp(self):
        # PDB
        self.pdb = os.path.join(DATADIR, "tests", "1FMO.pdb")
        # XYZ
        self.xyz = os.path.join(DATADIR, "tests", "1FMO.xyz")
        # Full workflow
        self.results = pyKVFinder.run_workflow(
            self.pdb,
            include_depth=True,
            include_hydropathy=True,
            hydrophobicity_scale="EisenbergWeiss",
        )

    def test_custom_vdw(self):
        results = pyKVFinder.run_workflow(
            self.pdb, vdw=os.path.join(DATADIR, "vdw.dat")
        )
        self.assertEqual(results.ncav, self.results.ncav)

    def test_standard_workflow(self):
        results = pyKVFinder.run_workflow(self.pdb)
        self.assertEqual(results.ncav, self.results.ncav)

    def test_standard_workflow_with_xyz_input(self):
        results = pyKVFinder.run_workflow(self.xyz)
        self.assertEqual(results.ncav > 0, True)

    def test_full_workflow(self):
        results = pyKVFinder.run_workflow(
            self.pdb,
            include_depth=True,
            include_hydropathy=True,
            hydrophobicity_scale="EisenbergWeiss",
        )
        self.assertEqual(results.ncav, self.results.ncav)

    def test_ligand_mode(self):
        results = pyKVFinder.run_workflow(
            self.pdb,
            os.path.join(DATADIR, "tests", "ADN.pdb"),
        )
        self.assertEqual(results.ncav > 0, True)

    def test_ligand_mode_with_xyz_input(self):
        results = pyKVFinder.run_workflow(
            self.pdb,
            os.path.join(DATADIR, "tests", "ADN.xyz"),
        )
        self.assertEqual(results.ncav > 0, True)

    def test_residues_box(self):
        results = pyKVFinder.run_workflow(
            self.pdb,
            box=os.path.join(
                DATADIR,
                "tests",
                "residues-box.toml",
            ),
        )
        self.assertEqual(results.ncav > 0, True)

    def test_custom_box(self):
        results = pyKVFinder.run_workflow(
            self.pdb,
            box=os.path.join(DATADIR, "tests", "custom-box.toml"),
        )
        self.assertEqual(results.ncav, 1)

    def test_no_cavities_detected(self):
        results = pyKVFinder.run_workflow(
            self.pdb,
            volume_cutoff=1000000,  # Force ncav = 0
        )
        self.assertEqual(results, None)

    def test_invalid_input_extension(self):
        self.assertRaises(TypeError, pyKVFinder.run_workflow, "any.mol")

    def test_invalid_ligand_extension(self):
        self.assertRaises(
            TypeError, pyKVFinder.run_workflow, self.pdb, ligand="any.mol"
        )


class TestPyKVFinderResults(unittest.TestCase):
    def setUp(self):
        # PDB
        self.pdb = os.path.join(DATADIR, "tests", "1FMO.pdb")
        # Full workflow
        self.results = pyKVFinder.run_workflow(
            self.pdb,
            include_depth=True,
            include_hydropathy=True,
            hydrophobicity_scale="EisenbergWeiss",
        )

    def test_repr(self):
        self.assertEqual(self.results.__repr__(), "<pyKVFinderResults object>")

    def test_export(self):
        output = os.path.join(
            DATADIR,
            "tests",
            "output",
            "cavities.pdb",
        )
        # export
        self.results.export(output=output)
        self.assertEqual(os.path.exists(output), True)
        os.remove(output)

    def test_write(self):
        results = os.path.join(
            DATADIR,
            "tests",
            "output",
            "results.toml",
        )
        output = os.path.join(
            DATADIR,
            "tests",
            "output",
            "cavities.pdb",
        )
        # write
        self.results.write(fn=results, output=output)
        self.assertEqual(os.path.exists(results), True)
        os.remove(results)

    def test_plot_frequencies(self):
        pdf = os.path.join(
            DATADIR,
            "tests",
            "output",
            "barplots.pdf",
        )
        # plot_frequencies
        self.results.plot_frequencies(pdf=pdf)
        self.assertEqual(os.path.exists(pdf), True)
        os.remove(pdf)

    def test_export_all(self):
        results = os.path.join(
            DATADIR,
            "tests",
            "output",
            "results.toml",
        )
        output = os.path.join(
            DATADIR,
            "tests",
            "output",
            "cavities.pdb",
        )
        pdf = os.path.join(
            DATADIR,
            "tests",
            "output",
            "barplots.pdf",
        )
        # export_all
        self.results.export_all(
            fn=results,
            output=output,
            include_frequencies_pdf=True,
            pdf=pdf,
        )
        self.assertEqual(os.path.exists(results), True)
        os.remove(results)
        self.assertEqual(os.path.exists(output), True)
        os.remove(output)
        self.assertEqual(os.path.exists(pdf), True)
        os.remove(pdf)
