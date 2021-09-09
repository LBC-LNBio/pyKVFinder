pyKVFinder.pyKVFinderResults.export_all
---------------------------------------

.. autofunction:: pyKVFinder.pyKVFinderResults.export_all

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

Yet, we can set a ``include_frequencies_pdf`` flag to True to plot the bar charts of the frequencies in a PDF file.

.. code-block:: python

  >>> results.export_all(include_frequencies_pdf=True)
