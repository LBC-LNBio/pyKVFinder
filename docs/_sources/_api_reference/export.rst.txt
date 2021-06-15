pyKVFinder.export
=================

.. autofunction:: pyKVFinder.export

.. seealso::

    * `pyKVFinder.detect <detect.html>`_
    * `pyKVFinder.spatial <spatial.html>`_
    * `pyKVFinder.depths <depth.html>`_
    * `pyKVFinder.hydropathy <hydropathy.html>`_
    * `pyKVFinder.write_results <write_results.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

With the cavity and surface points identified and depth and hydrophobicity scale mapped in the 3D grid, we can:

* Export cavity points

.. code-block:: python

    >>> from pyKVFinder import export
    >>> export('cavity_wo_surface.pdb', cavities, None, vertices)

* Export cavity and surface points

.. code-block:: python

    >>> export('cavities.pdb', cavities, surface, vertices)

* Export cavity and surface points with depth mapped on them

.. code-block:: python

    >>> export('cavities_with_depth.pdb', cavities, surface, vertices, B=depths)

* Export surface points with hydrophobicity_scale mapped on them

.. code-block:: python

    >>> export(None, None, surface, vertices, output_hydropathy='hydropathy.pdb', scales=scales)

* Export all

.. code-block:: python

    >>> export('cavities.pdb', cavities, surface, vertices, B=depths, output_hydropathy='hydropathy.pdb', scales=scales)
