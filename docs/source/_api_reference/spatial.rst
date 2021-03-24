pyKVFinder.spatial
==================

Spatial characterization (volume and area) of the detected cavities.

Surface points that belongs to the same cavity are assigned with an integer in the grid. Biomolecule or unsigned cavity points = 0, Assigned cavity points >= 2 and Bulk points = -1.

.. code-block:: python

    pyKVFinder.spatial(cavities, ncav, step = 0.6, nthreads = os.cpu_count() - 1, verbose = False)

:Args:

    ``cavities`` : *numpy.ndarray*
        A numpy array with integer labels in each position, that are:
            * -1: bulk point.
            * 0: biomolecule point.
            * 1: empty space point.
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
    
    ``surface`` : *numpy.ndarray*
        A numpy array with integer labels in each position, that are:
            * -1: bulk point.
            * 0: biomolecule or empty space point.
            * >=2: surface point.
    ``volume`` : *dict*
        A dictionary with volume of each detected cavity
    ``area`` : *dict*
        A dictionary with area of each detected cavity


.. note::

    The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on. 

.. seealso::

    * `pyKVFinder.detect <spatial.html>`_
    * `pyKVFinder.hydropathy <hydropathy.html>`_
    * `pyKVFinder.export <export.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

With the cavity points identified with ``pyKVFinder.detect``, we can perform a spatial characterization, that includes volume, area and defining surface points:

.. code-block:: python

    >>> from pyKVFinder import spatial
    >>> surface, volume, area = spatial(cavities, ncav)
    >>> surface
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
    >>> volume
    {'KAA': 137.16, 'KAB': 47.52, 'KAC': 66.96, 'KAD': 8.21, 'KAE': 43.63, 'KAF': 12.53, 'KAG': 6.26, 'KAH': 520.13, 'KAI': 12.31, 'KAJ': 26.57, 'KAK': 12.31, 'KAL': 33.91, 'KAM': 23.11, 'KAN': 102.82, 'KAO': 6.05, 'KAP': 15.55, 'KAQ': 7.99, 'KAR': 7.78}
    >>> area
    {'KAA': 120.52, 'KAB': 58.76, 'KAC': 72.06, 'KAD': 17.62, 'KAE': 56.44, 'KAF': 22.53, 'KAG': 15.38, 'KAH': 489.25, 'KAI': 29.87, 'KAJ': 44.85, 'KAK': 30.58, 'KAL': 43.59, 'KAM': 45.25, 'KAN': 129.77, 'KAO': 11.57, 'KAP': 24.8, 'KAQ': 12.59, 'KAR': 15.97}
