pyKVFinder.get_grid_from_file
=============================

Gets 3D grid vertices from box configuration file or parKVFinder parameters file, selects atoms inside custom 3D grid, define sine and cosine of 3D grid angles and define xyz grid units.

.. code-block:: python
    
    pyKVFinder.get_grid_from_file(fn, atominfo, xyzr, step = 0.6, probe_in = 1.4, probe_out = 4.0, nthreads = os.cpu_count() - 1)

:Args:

    ``fn`` : *str*
        A path to a box configuration file (TOML-formatted)
    ``atominfo`` : *numpy.ndarray*
        A numpy array with atomic information (residue number, chain, residue name, atom name)
    ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz atomic coordinates and radii values (x, y, z, radius)
    ``step`` : *float, default 0.6*
        Grid spacing (A)
    ``probe_in`` : *float, default 1.4*
        Probe In size (A)
    ``probe_out`` : *float, default 4.0*
        Probe Out size (A)
    ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads

:Returns:

    ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis) of the custom box
    ``atominfo`` : *numpy.ndarray*
        A numpy array with atomic information (residue number, chain, residue name, atom name) of atoms inside the custom box
    ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz atomic coordinates and radii values (x, y, z, radius) of atoms inside the custom box
    ``sincos`` : *numpy.ndarray*
        A numpy array with sine and cossine of the custom box rotation angles (sina, cosa, sinb, cosb)
    ``nx`` : *int*
        x grid units
    ``ny`` : *int*
        y grid units 
    ``nz`` : *int*
        z grid units 

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
