pyKVFinder.pyKVFinderResults.write
----------------------------------

Writes file paths and cavity characterization to TOML-formatted file.

.. code-block:: python
  
  pyKVFinder.pyKVFinderResults.write(fn = 'results.toml, output = None, output_hydropathy = None)

:Args:

  ``fn`` : *str, default 'results.toml'*
    A path to TOML-formatted file for writing file paths and cavity characterization (volume, area, depth and interface residues) per cavity detected
  ``output`` : *str, default None*
    A path to a cavity PDB file
  ``output_hydropathy`` : *str, default None*
    A path to PDB file for writing hydropathy at surface points

:Returns:

  File with TOML-formatted data corresponding to file paths and cavity characterization per detected cavity.

.. note::

  The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on. 

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
  >>> results.write()
