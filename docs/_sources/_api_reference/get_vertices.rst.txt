pyKVFinder.get_vertices
=======================

.. autofunction:: pyKVFinder.get_vertices

.. seealso::

    * `pyKVFinder.read_pdb <read_pdb.html>`_
    * `pyKVFinder.read_xyz <read_xyz.html>`_
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
