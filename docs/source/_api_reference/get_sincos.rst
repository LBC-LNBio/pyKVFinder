pyKVFinder.get_sincos
=====================

Gets sine and cossine of the grid rotation angles.

.. code-block:: python
    
    pyKVFinder.get_sincos(vertices)
  
:Args:

    ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)

:Returns:

    ``sincos`` : *numpy.ndarray*
        A numpy array with sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb)

.. seealso::

    * `pyKVFinder.get_vertices <get_vertices.html>`_
    * `pyKVFinder.detect <detect.html>`_
    * `pyKVFinder.constitutional <constitutional.html>`_
    * `pyKVFinder.hydropathy <hydropathy.html>`_
    * `pyKVFinder.export <export.html>`_

.. raw:: html

    <h4><u>Example</u></h4>

With the coordinates of the 3D grid vertices calculated with ``pyKVFinder.get_vertices``, we can get sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb):

.. code-block:: python

    >>> from pyKVFinder import get_sincos
    >>> sincos = get_sincos(vertices)
    >>> sincos
    array([0., 1., 0., 1.])
