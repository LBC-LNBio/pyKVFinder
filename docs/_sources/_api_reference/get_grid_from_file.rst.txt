pyKVFinder.get_grid_from_file
=============================

.. autofunction:: pyKVFinder.get_grid_from_file

.. seealso::

    * `pyKVFinder.read_pdb <read_pdb.html>`_
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

With the atomic information and coordinates read with ``pyKVFinder.read_pdb`` and a box configuration file, we can get the coordinates of grid vertices, select atoms inside custom 3D grid, define sine and cosine of rotation angles and define xyz grid units.

.. code-block:: python

    >>> from pyKVFinder import get_grid_from_file
    >>> vertices, atominfo, xyzr, sincos, nx, ny, nz = pyKVFinder.get_grid_from_file(fn, atominfo, xyzr)

.. warning::

    Custom box coordinates adds Probe Out size in each direction to create the coordinates of grid vertices.
