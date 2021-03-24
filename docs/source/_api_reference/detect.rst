pyKVFinder.detect
=================

Detects biomolecular cavities. 

Cavity points that belongs to the same cavity are assigned with an integer in the grid. Biomolecule points = 0, Unsigned cavity points = 1, Assigned cavity points >= 2 and Bulk points = -1.

.. code-block:: python
    
    pyKVFinder.detect(nx, ny, nz, xyzr, vetices, sincos, step = 0.6, probe_in = 1.4, probe_out = 4.0, removal_distance = 2.4, volume_cutoff = 5.0, lxyzr = None, ligand_cutoff = 5.0, box_adjustment = False, surface = 'SES', nthreads = os.cpu_count() - 1, verbose = False)

:Args:
    ``nx`` : *int*
        x 3D grid units
    ``nx`` : *int*
        y 3D grid units
    ``nx`` : *int*
        z 3D grid units
    ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz atomic coordinates and radii values (x, y, z, radius) 
    ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)
    ``sincos`` : *numpy.ndarray*
        A numpy array with sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb)
    ``step`` : *float, default 0.6*
        Grid spacing (A)
    ``probe_in`` : *float, default 1.4*
        Probe In size (A)
    ``probe_out`` : *float, default 4.0*
        Probe Out size (A)
    ``removal_distance`` : *float, default 2.4*
        Length to be removed from the cavity-bulk frontier (A)
    ``volume_cutoff`` : *float, default 5.0*
        Cavities volume filter (A3)
    ``lxyzr`` : *numpy.ndarray*
        A numpy array with xyz atomic coordinates and radii values (x, y, z, radius) of ligand atoms
    ``ligand_cutoff`` : *float, default 5.0*
        Radius value to limit a space around a ligand (A)
    ``box_adjustment`` :  *bool, default False*
        Whether a custom 3D grid is applied
    ``surface`` : *str, default 'SES'*
        Keyword options passed to surface representation: 
            * SES: Solvent Excluded Surface.
            * SAS: Solvent Accessible Surface
    ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads
    ``verbose`` : *bool, default False*
        Print extra information to standard output

:Returns:
    
    ``ncav``: *int*
        Number of detected cavities
    ``cavities``: *numpy.ndarray*
        A numpy array with integer labels in each position, that are:
            * -1: bulk point.
            * 0: biomolecule point.
            * 1: empty space point.
            * >=2: cavity point.

.. note::

    Points belonging to the same cavity receive the same integer label.

    The empty space points are regions that do not meet the volume cutoff defined in ``detect`` function.   

.. warning::

    If you are using box adjusment mode, do not forget to set ``box_adjustment`` flag to ``True`` and read box configuration file with `get_grid_from_file <get_grid_from_file.html>`_.

    If you are using ligand adjustment mode, do not forget to read ligand atom coordinates with `read_pdb <read_pdb.html>`_.

.. seealso::

    * `pyKVFinder.read_pdb <read_pdb.html>`_
    * `pyKVFinder.get_vertices <get_vertices.html>`_
    * `pyKVFinder.get_dimensions <get_dimensions.html>`_
    * `pyKVFinder.get_sincos <get_sincos.html>`_
    * `pyKVFinder.get_grid_from_file <get_grid_from_file.html>`_
    * `pyKVFinder.spatial <spatial.html>`_
    * `pyKVFinder.depth <depth.html>`_
    * `pyKVFinder.constitutional <constitutional.html>`_
    * `pyKVFinder.export <export.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

With the grid defined with ``pyKVFinder.get_vertices``, ``pyKVFinder.get_dimensions`` and ``pyKVFinder.get_sincos``, and atomic coordiantes loaded with ``pyKVFinder.read_pdb``, we can detect cavities on the whole target biomolecule:

.. code-block:: python

    >>> from pyKVFinder import detect
    >>> ncav, cavities = detect(nx, ny, nz, xyzr, vertices, sincos)
    >>> ncav
    18
    >>> cavities
    array([[[-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        ...,
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1]],

       ...,

       [[-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        ...,
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)

However, users may opt to perform cavity detection in a segmented space through ligand adjustment and/or box adjustment modes.

The cavity detection can be limited around the target ligand(s), which will be passed to pyKVFinder through a *.pdb* file. Thus, the detected cavities are limited within a radius (``ligand_cutoff``) of the target ligand(s).

.. code-block:: python

    >>> import os
    >>> ligand = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', 'ADN.pdb')
    >>> from pyKVFinder import read_pdb
    >>> _, lxyzr = read_pdb(ligand, vdw)
    >>> ncav, cavities = detect(nx, ny, nz, xyzr, vertices, sincos, lxyzr=lxyzr, ligand_cutoff=5.0)
    >>> ncav
    2
    >>> cavities
    array([[[-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        ...,
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1]],

       ...,

       [[-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        ...,
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)

Further, we can also perform cavity detection on a custom 3D grid, where we can explore closed regions with a custom box, which can be defined by a *.toml* file (see `Box configuration file template <box_file_template.html>`_).

.. code-block:: python

    >>> import os
    >>> fn = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', 'custom-box.toml')
    >>> with open(fn, 'r') as f:
    ...     print(f.read())

With this box adjustment mode, we must defined the 3D grid with ``pyKVFinder.get_grid_from_file``. 

.. code-block:: python

    >>> from pyKVFinder import get_grid_from_file
    >>> vertices, atominfo, xyzr, sincos, nx, ny, nz = pyKVFinder.get_grid_from_file(fn, atominfo, xyzr)

Then, we can perform cavity detection:

.. code-block:: python

    >>> ncav, cavities = pyKVFinder.detect(nx, ny, nz, xyzr, vertices, sincos, box_adjustment=True)
    >>> ncav
    1
    >>> cavities
    array([[[-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        ...,
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1]],

       ...,

       [[-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        ...,
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)

.. warning::

    If you are using box adjusment mode, do not forget to set ``box_adjustment`` flag to ``True``.
