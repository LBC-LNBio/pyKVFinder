pyKVFinder.read_pdb
===================

.. autofunction:: pyKVFinder.read_pdb

.. seealso::

    * `pyKVFinder.read_vdw <read_vdw.html>`_
    * `pyKVFinder.get_vertices <get_vertices.html>`_
    * `pyKVFinder.get_vertices_from_file <get_vertices_from_file.html>`_
    * `pyKVFinder.detect <detect.html>`_
    * `pyKVFinder.constitutional <constitutional.html>`_
    * `pyKVFinder.hydropathy <hydropathy.html>`_
    * `van der Waals file template <../_cfg_files/vdw_file_template.html>`_

.. raw:: html

    <h4><u>Example</u></h4>

With the vdW radii dictionary loaded with ``pyKVFinder.read_vdw``, we can read a target PDB file into Numpy arrays (atomic information and atomic coordinates):

.. code-block:: python

    >>> import os
    >>> import pyKVFinder
    >>> from pyKVFinder import read_pdb
    >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
    >>> atominfo, xyzr = read_pdb(pdb)
    >>> atominfo
    array([['13_E_GLU', 'N'],
       ['13_E_GLU', 'CA'],
       ['13_E_GLU', 'C'],
       ...,
       ['350_E_PHE', 'CE2'],
       ['350_E_PHE', 'CZ'],
       ['350_E_PHE', 'OXT']], dtype='<U9')
    >>> xyzr
    array([[ -6.693, -15.642, -14.858,   1.824],
       [ -6.73 , -14.62 , -15.897,   1.908],
       [ -7.49 , -13.357, -15.508,   1.908],
       ...,
       [  7.216,  18.878,  -9.885,   1.908],
       [  7.735,  17.624,  -9.558,   1.908],
       [  5.767,  19.234, -13.442,   1.69 ]])

.. warning::

  The function takes the `built-in dictionary <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/vdw.dat>`_ when the ``vdw`` argument is not specified. If you wish to use a custom van der Waals radii file, you must read it with ``pyKVFinder.read_vdw`` as shown earlier and pass it as ``pyKVFinder.read_pdb(pdb, vdw=vdw)``.
