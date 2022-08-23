pyKVFinder.export_openings
==========================

.. autofunction:: pyKVFinder.export_openings

.. seealso::

    * `pyKVFinder.export <export.html>`_
    * `pyKVFinder.detect <detect.html>`_
    * `pyKVFinder.depths <depth.html>`_
    * `pyKVFinder.openings <openings.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

With the opening points identified with ``pyKVFinder.openings``, we can export them to a PDB-formatted file:

.. code-block:: python

    >>> from pyKVFinder import export_openings
    >>> export('openings.pdb', openings, vertices)
