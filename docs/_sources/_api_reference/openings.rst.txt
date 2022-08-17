pyKVFinder.openings
===================

.. autofunction:: pyKVFinder.openings

.. seealso::

    * `pyKVFinder.detect <spatial.html>`_
    * `pyKVFinder.depth <depth.html>`_
    * `pyKVFinder.export <export.html>`_

.. raw:: html

  <h4><u>Example</u></h4>


With the cavity points identified with ``pyKVFinder.detect``, we can characterize their openings, that includes number and area of openings and defining opening points:

.. code-block:: python

    >>> from pyKVFinder import openings
    >>> nopenings, openings, aopenings = openings(cavities)
    >>> nopenings
    19
    >>> openings
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
    >>> aopenings
    {'KAA': 70.18, 'KAB': 36.5, 'KAC': 85.98, 'KAD': 18.93, 'KAE': 41.56, 'KAF': 28.61, 'KAG': 1.08, 'KAH': 17.34, 'KAI': 30.58, 'KAJ': 22.8, 'KAK': 30.58, 'KAL': 27.17, 'KAM': 94.99, 'KAN': 53.6, 'KAO': 129.77, 'KAP': 13.78, 'KAQ': 17.89, 'KAR': 22.56, 'KAS': 27.54}

With the cavity and opening points identified, we can:

* Export opening points with same nomenclature from ``aopenings``:

.. code-block:: python
    
    >>> from pyKVFinder import export
    >>> export("openings.pdb", openings, None, vertices)

* Export cavity points with opening points mapped on them:

    >>> export("cavities_with_openings.pdb", cavities, None, vertices, B=openings)
