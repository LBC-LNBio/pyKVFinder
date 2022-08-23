pyKVFinder.openings
===================

.. autofunction:: pyKVFinder.openings

.. seealso::

    * `pyKVFinder.detect <detect.html>`_
    * `pyKVFinder.depth <depth.html>`_
    * `pyKVFinder.export <export.html>`_
    * `pyKVFinder.export_openings <export_openings.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

With the cavity points identified with ``pyKVFinder.detect``, we can characterize their openings, that includes number and area of openings and defining opening points:

.. code-block:: python

    >>> from pyKVFinder import openings
    >>> nopenings, openings, aopenings = openings(cavities)
    >>> nopenings
    16
    >>> openings
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
            [-1, -1, -1, ..., -1, -1, -1]]])   
    >>> aopenings
    {'KAA': {'OAA': 47.41, 'OAG': 3.6}, 'KAB': {'OAB': 25.84}, 'KAC': {'OAC': 53.62}, 'KAD': {'OAD': 12.59}, 'KAE': {'OAE': 26.3}, 'KAF': {'OAF': 18.46}, 'KAG': {'OAH': 12.83}, 'KAH': {'OAK': 59.96}, 'KAJ': {'OAI': 16.11}, 'KAL': {'OAJ': 17.3}, 'KAM': {'OAL': 35.27}, 'KAO': {'OAM': 8.49}, 'KAP': {'OAN': 13.71}, 'KAQ': {'OAO': 13.16}, 'KAR': {'OAP': 15.36}}


With the cavity and opening points identified, we can:

* Export cavity points with opening points mapped on them:

.. code-block:: python

    >>> from pyKVFinder import export
    >>> export("cavities_with_openings.pdb", cavities, None, vertices, B=openings)

* Export opening points with same nomenclature from ``aopenings``:

.. code-block:: python
    
    >>> from pyKVFinder import export_openings
    >>> export_openings("openings.pdb", openings, vertices)
