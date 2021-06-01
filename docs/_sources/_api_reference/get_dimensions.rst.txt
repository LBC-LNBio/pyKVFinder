pyKVFinder.get_dimensions
=========================

.. autofunction:: pyKVFinder.get_dimensions

.. seealso::

    * `pyKVFinder.get_vertices <get_vertices.html>`_
    * `pyKVFinder.detect <detect.html>`_

.. raw:: html

    <h4><u>Example</u></h4>

With the coordinates of the 3D grid vertices calculated with ``pyKVFinder.get_vertices``, we can get dimensions of the 3D grid:

.. code-block:: python

    >>> from pyKVFinder import get_dimensions
    >>> nx, ny, nz = get_dimensions(vertices)
    >>> nx, ny, nz
    (101, 126, 97)
