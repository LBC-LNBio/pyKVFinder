pyKVFinder.spatial
==================

.. autofunction:: pyKVFinder.spatial

.. seealso::

    * `pyKVFinder.detect <spatial.html>`_
    * `pyKVFinder.hydropathy <hydropathy.html>`_
    * `pyKVFinder.export <export.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

With the cavity points identified with ``pyKVFinder.detect``, we can perform a spatial characterization, that includes volume, area and defining surface points:

.. code-block:: python

    >>> from pyKVFinder import spatial
    >>> surface, volume, area = spatial(cavities)
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
    {'KAA': 126.41, 'KAB': 62.37, 'KAC': 74.57, 'KAD': 19.06, 'KAE': 57.08, 'KAF': 22.77, 'KAG': 15.38, 'KAH': 496.97, 'KAI': 30.58, 'KAJ': 45.64, 'KAK': 30.58, 'KAL': 45.58, 'KAM': 45.25, 'KAN': 129.77, 'KAO': 12.28, 'KAP': 25.04, 'KAQ': 13.46, 'KAR': 16.6}
