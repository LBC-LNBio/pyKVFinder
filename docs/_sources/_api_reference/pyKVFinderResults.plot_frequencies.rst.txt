pyKVFinder.pyKVFinderResults.plot_frequencies
---------------------------------------------

Plot histograms of frequencies (residues and classes of residues) in a PDF file.

.. code-block:: python
  
  pyKVFinder.pyKVFinderResults.plot_frequencies(pdf = 'histogram.pdf')

:Args:
  
  ``pdf`` : *str, default 'histograms.pdf'*
    A path to a PDF file

:Returns:
  
  A PDF file with histograms of calculated frequencies (residues and classes of residues) of each detected cavity.

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
  * `pyKVFinder.pyKVFinderResults.export_all <pyKVFinderResults.export_all.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

.. code-block:: python

  >>> from pyKVFinder import pyKVFinder
  >>> import os
  >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
  >>> results = pyKVFinder(pdb)
  >>> results.plot_frequencies()
