pyKVFinder.detect
=================

.. autofunction:: pyKVFinder.detect

.. seealso::

    * `pyKVFinder.read_pdb <read_pdb.html>`_
    * `pyKVFinder.read_xyz <read_xyz.html>`_
    * `pyKVFinder.get_vertices <get_vertices.html>`_
    * `pyKVFinder.get_vertices_from_file <get_vertices_from_file.html>`_
    * `pyKVFinder.spatial <spatial.html>`_
    * `pyKVFinder.depth <depth.html>`_
    * `pyKVFinder.constitutional <constitutional.html>`_
    * `pyKVFinder.export <export.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

With the grid vertices defined with ``pyKVFinder.get_vertices`` and atomic data loaded with ``pyKVFinder.read_pdb`` or ``pyKVFinder.read_xyz``, we can detect cavities on the whole target biomolecule:

.. code-block:: python

    >>> from pyKVFinder import detect
    >>> ncav, cavities = detect(atomic, vertices)
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

The cavity detection can be limited around the target ligand(s), which will be passed to pyKVFinder through a *.pdb* or a *.xyz* files. Thus, the detected cavities are limited within a radius (``ligand_cutoff``) of the target ligand(s).

.. code-block:: python

    >>> import os
    >>> ligand = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', 'ADN.pdb')
    >>> from pyKVFinder import read_pdb
    >>> latomic = read_pdb(ligand)
    >>> ncav, cavities = detect(atomic, vertices, latomic=latomic, ligand_cutoff=5.0)
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
    [box]
    p1 = [3.11, 7.34, 1.59]
    p2 = [11.51, 7.34, 1.59]
    p3 = [3.11, 10.74, 1.59]
    p4 = [3.11, 7.34, 6.19]


With this box adjustment mode, we must defined the 3D grid with ``pyKVFinder.get_vertices_from_file``. 

.. code-block:: python

    >>> from pyKVFinder import get_vertices_from_file
    >>> vertices, atomic = get_vertices_from_file(fn, atomic)

Then, we can perform cavity detection:

.. code-block:: python

    >>> ncav, cavities = detect(atomic, vertices, box_adjustment=True)
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
