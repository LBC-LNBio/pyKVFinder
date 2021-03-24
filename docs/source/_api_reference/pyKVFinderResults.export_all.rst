pyKVFinder.pyKVFinderResults.export_all
---------------------------------------

Exports cavities and characterization to PDB-formatted files, writes file paths and characterization to a TOML-formatted file, and optionally plot histograms of frequencies (residues and classes of residues) in a PDF file.

.. code-block:: python

    pyKVFinder.pyKVFinderResults.export_all(fn = 'results.toml', output = 'cavity.pdb', output_hydropathy = 'hydropathy.pdb', include_frequencies_pdf = False, pdf = 'histogtrams.pdf', nthreads = os.cpu_count() - 1)
  
:Args:

    ``fn`` : *str, default 'results.toml'*
        A path to TOML-formatted file for writing file paths and cavity characterization (volume, area and interface residues) per cavity detected
    ``output`` : *str, default 'cavity.pdb'*
        A path to PDB file for writing cavities
    ``output_hydropathy`` : *str, default 'hydropathy.pdb'*
        A path to PDB file for writing hydropathy at surface points
    ``include_frequencies_pdf`` : *bool, default False*
        Whether to plot frequencies (residues and classes of residues) to PDF file
    ``pdf`` : *str, default 'histograms.pdf'*
        A path to a PDF file
    ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads

:Returns:
    
    File with TOML-formatted data corresponding to file paths and cavity characterization per detected cavity.
    
    File with PDB-formatted data corresponding to cavity points (H), surface points (HA) and a target variable (B) in B-factor column.
    
    (Optional) File with PDB-formatted data corresponding to hydropathy mapped in B-factor column at surface points (HA).
    
    (Optional) PDF file with histograms of calculated frequencies (residues and classes of residues) of each detected cavity.

.. note::

  The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on. 

  The classes of residues are:

    * ``R1`` : Alipathic apolar
        Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine
    * ``R2`` : Aromatic
        Phenylalanine, Tryptophan, Tyrosine
    * ``R3`` : Polar Uncharged
        Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine
    * ``R4`` : Negatively charged
        Aspartate, Glutamate
    * ``R5`` : Positively charged
        Arginine, Histidine, Lysine
    * ``RX`` : Non-standard
        Non-standard residues

.. seealso::

    * `pyKVFinder.pyKVFinderResults <pyKVFinderResults.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

.. code-block:: python

  >>> from pyKVFinder import pyKVFinder
  >>> import os
  >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
  >>> results = pyKVFinder(pdb)
  >>> results.export_all()

Yet, we can set a ``include_frequencies_pdf`` flag to True to plot the histograms of the frequencies in a PDF file.

.. code-block:: python

  >>> results.export_all(include_frequencies_pdf=True)
