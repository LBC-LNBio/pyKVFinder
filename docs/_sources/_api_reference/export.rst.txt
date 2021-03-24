pyKVFinder.export
=================

Exports cavitiy (H) and surface (HA) points to PDB-formatted file with a variable (B; optional) in B-factor column, and hydropathy to PDB-formatted file in B-factor column at surface points (HA).

.. code-block:: python
    
    pyKVFinder.export(fn, cavities, surface, vertices, sincos, ncav, step = 0.6, B = None, output_hydropathy = 'hydropathy.pdb', scales = None, nthreads = os.cpu_count() - 1, append = False, model = 0)

:Args:

    ``fn`` : *str*
        A path to PDB file for writing cavities
    ``cavities`` : *numpy.ndarray*
        A numpy array with integer labels in each position, that are:
            * -1: bulk point.
            * 0: biomolecule point.
            * 1: empty space.
            * >=2: cavity point.
    ``surface`` : *numpy.ndarray*
        A numpy array with integer labels in each position, that are:
            * -1: bulk point.
            * 0: biomolecule or empty space point.
            * >=2: surface point.
    ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)
    ``sincos`` : *numpy.ndarray*
        A numpy array with sine and cossine of 3D grid angles (a, b)
    ``ncav`` : *int*
        Number of cavities in ``cavities`` and ``surface`` numpy arrays
    ``step`` : *float, default 0.6*
        Grid spacing (A)
    ``B`` : *numpy.ndarray*
        Values to be mapped on B-factor column in cavity points (B[nx][ny][nz])
    ``output_hydropathy`` :  *str, default 'hydropathy.pdb'*
        A path to hydropathy PDB file (surface points mapped with a hydrophobicity scale)
    ``scales``: *numpy.ndarray, default None*
        Hydrophobicity scale values to be mapped on B-factor column in surface points (scales[nx][ny][nz])
    ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads
    ``append`` : *bool, default False*
        Append cavities to PDB file
    ``model`` : *int, default 0*
        Model number

:Returns:

  File with PDB-formatted data corresponding to cavity points (H), surface points (HA) and a target variable (B) in B-factor column.
  
  (Optional) File with PDB-formatted data corresponding to hydropathy mapped in B-factor column at surface points (HA).

.. note::

    The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on. 

.. seealso::

    * `pyKVFinder.detect <detect.html>`_
    * `pyKVFinder.spatial <spatial.html>`_
    * `pyKVFinder.depths <depth.html>`_
    * `pyKVFinder.hydropathy <hydropathy.html>`_
    * `pyKVFinder.write_results <write_results.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

With the cavity and surface points identified and depth and hydrophobicity scale mapped in the 3D grid, we can:

* Export cavity points

.. code-block:: python

    >>> from pyKVFinder import export
    >>> export('cavity_wo_surface.pdb', cavities, None, vertices, sincos, ncav)

* Export cavity and surface points

.. code-block:: python

    >>> export('cavities.pdb', cavities, surface, vertices, sincos, ncav)

* Export cavity and surface points with depth mapped on them

.. code-block:: python

    >>> export('cavities_with_depth.pdb', cavities, surface, vertices, sincos, ncav, B=depths)

* Export surface points with hydrophobicity_scale mapped on them

.. code-block:: python

    >>> export(None, None, surface, vertices, sincos, ncav, output_hydropathy='hydropathy.pdb', scales=scales)

* Export all

.. code-block:: python

    >>> export('cavities.pdb', cavities, surface, vertices, sincos, ncav, B=depths, output_hydropathy='hydropathy.pdb', scales=scales)
