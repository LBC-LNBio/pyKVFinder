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
    {'KAA': 3.79, 'KAB': 2.68, 'KAC': 2.62, 'KAD': 0.85, 'KAE': 3.0, 'KAF': 0.85, 'KAG': 0.6, 'KAH': 10.73, 'KAI': 0.0, 'KAJ': 2.24, 'KAK': 0.0, 'KAL': 3.0, 'KAM': 1.2, 'KAN': 0.0, 'KAO': 1.04, 'KAP': 2.08, 'KAQ': 0.85, 'KAR': 0.6}
    >>> avg_depth
    {'KAA': 1.35, 'KAB': 0.91, 'KAC': 0.68, 'KAD': 0.32, 'KAE': 0.99, 'KAF': 0.24, 'KAG': 0.1, 'KAH': 3.91, 'KAI': 0.0, 'KAJ': 0.96, 'KAK': 0.0, 'KAL': 1.07, 'KAM': 0.24, 'KAN': 0.0, 'KAO': 0.29, 'KAP': 0.7, 'KAQ': 0.22, 'KAR': 0.12}
