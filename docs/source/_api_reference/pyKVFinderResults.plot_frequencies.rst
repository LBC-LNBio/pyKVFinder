pyKVFinder.pyKVFinderResults.plot_frequencies
---------------------------------------------

.. autofunction:: pyKVFinder.pyKVFinderResults.plot_frequencies

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
