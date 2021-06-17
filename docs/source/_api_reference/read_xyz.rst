pyKVFinder.read_xyz
===================

.. autofunction:: pyKVFinder.read_xyz

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

With the vdW radii dictionary loaded with ``pyKVFinder.read_vdw``, we can read a target XYZ file into Numpy arrays (atomic information and atomic coordinates):

.. code-block:: python

    >>> import os
    >>> import pyKVFinder
    >>> from pyKVFinder import read_xyz
    >>> xyz = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.xyz')
    >>> atomic = read_xyz(xyz)
    >>> atominfo
    array([['1', 'A', 'UNK', ..., '-15.642', '-14.858', '1.97'],
       ['2', 'A', 'UNK', ..., '-14.62', '-15.897', '1.66'],
       ['3', 'A', 'UNK', ..., '-13.357', '-15.508', '1.66'],
       ...,
       ['2790', 'A', 'UNK', ..., '18.878', '-9.885', '1.66'],
       ['2791', 'A', 'UNK', ..., '17.624001', '-9.558', '1.66'],
       ['2792', 'A', 'UNK', ..., '19.233999', '-13.442', '1.69']],
      dtype='<U32')

.. warning::

  The function takes the `built-in dictionary <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/vdw.dat>`_ when the ``vdw`` argument is not specified. If you wish to use a custom van der Waals radii file, you must read it with ``pyKVFinder.read_vdw`` as shown earlier and pass it as ``pyKVFinder.read_xyz(xyz, vdw=vdw)``.
