pyKVFinder.depth
================

Characterization of the depth of the detected cavities, including depth per cavity point and maximum and average depths of detected cavities.

.. code-block:: python
    
    pyKVFinder.depth(cavities, ncav, step = 0.6, nthreads = os.cpu_count() - 1, verbose = False)
  
:Args:

    ``cavities`` : *numpy.ndarray*
        A numpy array with integer labels in each position, that are:
            * -1: bulk point.
            * 0: biomolecule point.
            * 1: empty space.
            * >=2: cavity point.
    ``ncav`` : *int*
        Number of cavities in ``cavities`` numpy array
    ``step`` : *float, default 0.6*
        Grid spacing (A)
    ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads
    ``verbose`` : *bool, default False*
        Print extra information to standard output

:Returns:
    
    ``depths``: *numpy.ndarray*
        A numpy array with depth of cavity points (depth[nx][ny][nz])
    ``max_depth``: *dict*
        A dictionary with maximum depth of each detected cavity
    ``avg_depth``: *dict*
        A dictionary with average depth of each detected cavity

.. note::

  The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

.. seealso::
  
  * `pyKVFinder.detect <detect.html>`_
  * `pyKVFinder.export <export.html>`_
  * `pyKVFinder.write_results <write_results.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

With the cavity points identified with ``pyKVFinder.detect``, we can perform a depth characterization, that includes maximum depth, average depth and defining depth of cavity points:

.. code-block:: python

    >>> from pyKVFinder import depth
    >>> depths, max_depth, avg_depth = depth(cavities, ncav)
    >>> depths
    array([[[0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            ...,
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.]],

        ...,

        [[0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            ...,
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.]]])
    >>> max_depth
    {'KAA': 3.79, 'KAB': 2.68, 'KAC': 2.62, 'KAD': 0.85, 'KAE': 3.0, 'KAF': 0.85, 'KAG': 0.6, 'KAH': 10.73, 'KAI': 2.55, 'KAJ': 2.24, 'KAK': 0.0, 'KAL': 3.0, 'KAM': 1.2, 'KAN': 0.0, 'KAO': 1.04, 'KAP': 2.08, 'KAQ': 0.85, 'KAR': 0.6}
    >>> avg_depth
    {'KAA': 1.28, 'KAB': 0.86, 'KAC': 0.67, 'KAD': 0.29, 'KAE': 0.98, 'KAF': 0.24, 'KAG': 0.1, 'KAH': 3.75, 'KAI': 1.5, 'KAJ': 0.96, 'KAK': 0.0, 'KAL': 1.0, 'KAM': 0.24, 'KAN': 0.0, 'KAO': 0.29, 'KAP': 0.7, 'KAQ': 0.22, 'KAR': 0.12}
