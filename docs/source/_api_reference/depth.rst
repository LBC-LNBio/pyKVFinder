pyKVFinder.depth
================

.. autofunction:: pyKVFinder.depth

.. seealso::
  
  * `pyKVFinder.detect <detect.html>`_
  * `pyKVFinder.export <export.html>`_
  * `pyKVFinder.write_results <write_results.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

With the cavity points identified with ``pyKVFinder.detect``, we can perform a depth characterization, that includes maximum depth, average depth and defining depth of cavity points:

.. code-block:: python

    >>> from pyKVFinder import depth
    >>> depths, max_depth, avg_depth = depth(cavities)
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
