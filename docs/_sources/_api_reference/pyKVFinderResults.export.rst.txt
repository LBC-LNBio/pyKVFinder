pyKVFinder.pyKVFinderResults.export
-----------------------------------

Exports cavitiy (H) and surface (HA) points to PDB-formatted file with a variable (B; optional) in B-factor column, and hydropathy to PDB-formatted file in B-factor column at surface points (HA).

.. code-block:: python
  
  pyKVFinder.pyKVFinderResults.export(output = 'cavity.pdb', output_hydropathy = 'hydropathy.pdb', nthreads = os.cpu_count() - 1)

:Args:

  ``output`` : *str, default 'cavity.pdb'*
    A path to PDB file for writing cavities
  ``output_hydropathy`` : *str, default 'hydropathy.pdb'*
    A path to PDB file for writing hydropathy at surface points
  ``nthreads`` : *int, default 'number of cpus - 1'*
    Number of threads

:Returns:

  File with PDB-formatted data corresponding to cavity points (H), surface points (HA) and a target variable (B) in B-factor column.
  
  (Optional) File with PDB-formatted data corresponding to hydropathy mapped in B-factor column at surface points (HA).

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
  >>> results.export()
