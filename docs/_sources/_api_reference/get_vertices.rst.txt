pyKVFinder.get_vertices
=======================

Gets 3D grid vertices.

.. code-block:: python

    pyKVFinder.get_vertices(xyzr, probe_out = 4.0, step = 0.6)

:Args:

    ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz atomic coordinates and radii values (x, y, z, radius) 
    ``probe_out`` : *float, default 4.0*
        Probe Out size (A)
    ``step`` : *float, default 0.6*
        Grid spacing (A)

:Returns:

    ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)

.. seealso::

    * `pyKVFinder.read_pdb <read_pdb.html>`_
    * `pyKVFinder.get_sincos <get_sincos.html>`_
    * `pyKVFinder.get_dimensions <get_dimensions.html>`_
    * `pyKVFinder.detect <detect.html>`_
    * `pyKVFinder.constitutional <constitutional.html>`_
    * `pyKVFinder.hydropathy <hydropathy.html>`_
    * `pyKVFinder.export <export.html>`_

.. raw:: html

    <h4><u>Example</u></h4>

With the atomic coordinates read with ``pyKVFinder.read_pdb``, we can get the coordinates of 3D grid vertices (origin, X-axis, Y-axis, Z-axis):

.. code-block:: python

    >>> from pyKVFinder import get_vertices
    >>> vertices = get_vertices(xyzr)
    >>> vertices
    array([[-19.911, -32.125, -30.806],
       [ 40.188, -32.125, -30.806],
       [-19.911,  43.446, -30.806],
       [-19.911, -32.125,  27.352]])
