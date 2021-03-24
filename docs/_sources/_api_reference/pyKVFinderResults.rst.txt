pyKVFinder.pyKVFinderResults
============================

A class containing pyKVFinder detection and characterization results.

.. code-block:: python

    pyKVFinder.pyKVFinderResults(cavities, surface, depths, scales, volume, area, max_depth, avg_depth, avg_hydropathy, residues, _vertices, _step, _ncav, _pdb = None, _ligand = None)

:Attributes:

    ``cavities`` : *numpy.ndarray*
        A numpy array with integer labels in each position, that are:
            * -1: bulk point.
            * 0: biomolecule point.
            * 1: empty space.
            * >=2: cavity point.
    ``surface`` : *numpy.ndarray*
        A numpy array with integer labels in each position, that are:
            * -1: bulk point.
            * 0: biomolecule or empty space point.
            * >=2: cavity point.
    ``depths`` : *numpy.ndarray*
        A numpy array with depth of cavity points (depth[nx][ny][nz])
    ``scales``: *numpy.ndarray*
        Hydrophobicity scale values mapped at surface points (scales[nx][ny][nz])
    ``volume`` : *dict*
        A dictionary with volume of each detected cavity
    ``area`` : *dict*
        A dictionary with area of each detected cavity
    ``max_depth`` : *dict*
        A dictionary with maximum depth of each detected cavity
    ``avg_depth`` : *dict*
        A dictionary with average depth of each detected cavity
    ``avg_hydropathy`` : *dict*
        A dictionary with average hydropathy of each detected cavity and range of the hydrophobicity scale mapped
    ``residues`` : *dict*
        A dictionary with interface residues of each detected cavity
    ``frequency`` : *dict*
        A dictionary with frequency of residues and class of residues of each detected cavity
    ``_vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)
    ``_step`` : *float*
        Grid spacing (A)
    ``_ncav`` : *int*
        Number of cavities in ``cavities`` and ``surface`` numpy arrays
    ``_pdb`` : *str, default None*
        A path to input PDB file
    ``_ligand`` : *str, default None*
        A path to ligand PDB file

.. raw:: html

    <h3>Methods</h3>

.. toctree::
    :hidden:
    :maxdepth: 1

    pyKVFinderResults.export <pyKVFinderResults.export>
    pyKVFinderResults.write <pyKVFinderResults.write>
    pyKVFinderResults.plot_frequencies <pyKVFinderResults.plot_frequencies>
    pyKVFinderResults.export_all <pyKVFinderResults.export_all>

* `pyKVFinder.pyKVFinderResults.export <pyKVFinderResults.export.html>`_: Exports cavities and characterizations to PDB-formatted files.
* `pyKVFinder.pyKVFinderResults.write <pyKVFinderResults.write.html>`_: Writes file paths and characterizations to a TOML-formatted results file.
* `pyKVFinder.pyKVFinderResults.plot_frequencies <pyKVFinderResults.plot_frequencies.html>`_: Plot histograms of frequencies (residues and classes of residues) in a PDF file.
* `pyKVFinder.pyKVFinderResults.export_all <pyKVFinderResults.export_all.html>`_: Exports cavities and characterizations to PDB-formatted files, writes file paths and characterizations to a TOML-formatted results file, and optionally plot histograms of frequencies in a PDF file.
