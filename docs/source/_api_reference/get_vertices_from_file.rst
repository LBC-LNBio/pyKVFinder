pyKVFinder.get_vertices_from_file
=================================

.. autofunction:: pyKVFinder.get_vertices_from_file

.. seealso::

    * `pyKVFinder.read_pdb <read_pdb.html>`_
    * `pyKVFinder.read_xyz <read_xyz.html>`_
    * `pyKVFinder.detect <detect.html>`_
    * `pyKVFinder.constitutional <constitutional.html>`_
    * `pyKVFinder.hydropathy <hydropathy.html>`_
    * `pyKVFinder.export <export.html>`_
    * `Box configuration file template <../_cfg_files/box_file_template.html>`_

.. raw:: html

    <h4><u>Example</u></h4>

First, define a box configuration file (see `Box configuration file template <../_cfg_files/box_file_template.html>`_).

.. code-block:: python

    >>> import os
    >>> fn = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', 'custom-box.toml')
    >>> with open(fn, 'r') as f:
    ...     print(f.read())
    [box]
    p1 = [3.11, 7.34, 1.59]
    p2 = [11.51, 7.34, 1.59]
    p3 = [3.11, 10.74, 1.59]
    p4 = [3.11, 7.34, 6.19]

With the atomic information and coordinates read with ``pyKVFinder.read_pdb`` and a box configuration file, we can get the coordinates of grid vertices and select atoms inside custom 3D grid.

.. code-block:: python

    >>> from pyKVFinder import get_vertices_from_file
    >>> vertices, atomic = pyKVFinder.get_vertices_from_file(fn, atomic)

.. warning::

    Custom box coordinates adds Probe Out size in each direction to create the coordinates of grid vertices.
