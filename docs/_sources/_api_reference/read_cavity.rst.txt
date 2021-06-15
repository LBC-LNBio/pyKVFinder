pyKVFinder.read_cavity
======================

.. autofunction:: pyKVFinder.read_cavity

.. seealso::

    * `pyKVFinder.read_pdb <read_pdb.html>`_
    * `pyKVFinder.get_vertices <get_vertices.html>`_
    * `pyKVFinder.get_grid_from_file <get_grid_from_file.html>`_
    * `pyKVFinder.spatial <spatial.html>`_
    * `pyKVFinder.depth <depth.html>`_
    * `pyKVFinder.constitutional <constitutional.html>`_
    * `pyKVFinder.hydropathy <hydropathy.html>`_
    * `pyKVFinder.export <export.html>`_
    * `van der Waals file template <../_cfg_files/vdw_file_template.html>`_

.. raw:: html

    <h4><u>Example</u></h4>

With a previously calculated cavity, that can be manually curated in a molecular visualization software, such as PyMOL, we can read it with its respective receptor back to pyKVFinder:

.. code-block:: python

    >>> import os
    >>> import pyKVFinder
    >>> from pyKVFinder import read_cavity
    >>> cavity = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.KVFinder.output.pdb')
    >>> receptor = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
    >>> grid = read_cavity(cavity, receptor)
    >>> grid
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

  The function takes the `built-in dictionary <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/vdw.dat>`_ when the ``vdw`` argument is not specified. If you wish to use a custom van der Waals radii file, you must set ``vdw`` argument to the custom file path, e.g. ``pyKVFinder.read_cavity(cavity, receptor, vdw='path/to/vdw.dat')``.
