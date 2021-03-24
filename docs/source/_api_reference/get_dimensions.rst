pyKVFinder.get_dimensions
=========================

Gets dimensions of 3D grid from vertices.

.. code-block:: python

    pyKVFinder.get_dimensions(vertices, step = 0.6)

:Args:

    ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)
    ``step`` : *float*
        Grid spacing (A)

:Returns:

    ``nx`` : *int*
        x grid units
    ``ny`` : *int*
        y grid units 
    ``nz`` : *int*
        z grid units 

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
