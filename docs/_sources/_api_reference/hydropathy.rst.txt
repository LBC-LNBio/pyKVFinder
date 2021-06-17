pyKVFinder.hydropathy
=====================

.. autofunction:: pyKVFinder.hydropathy

.. seealso::

    * `pyKVFinder.read_pdb <read_pdb.html>`_
    * `pyKVFinder.read_xyz <read_xyz.html>`_
    * `pyKVFinder.spatial <spatial.html>`_
    * `pyKVFinder.export <export.html>`_
    * `pyKVFinder.write_results <write_results.html>`_
    * `Hydrophobicity scale file template <../_cfg_files/hydrophobicity_scale_file_template.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

With the surface points identified with ``pyKVFinder.spatial`` and atomic coordinates and information read with ``pyKVFinder.read_pdb`` or ``pyKVFinder.read_xyz``, we can perform a hydropathy characterization, that maps a target hydrophobicity scale on surface points and calculate the average hydropathy

.. code-block:: python

    >>> from pyKVFinder import hydropathy
    >>> scales, avg_hydropathy = hydropathy(surface, atomic, vertices)
    >>> scales
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
    >>> avg_hydropathy
    {'KAA': -0.73, 'KAB': -0.05, 'KAC': -0.07, 'KAD': -0.62, 'KAE': -0.81, 'KAF': -0.14, 'KAG': -0.33, 'KAH': -0.16, 'KAI': -0.4, 'KAJ': 0.62, 'KAK': -0.99, 'KAL': 0.36, 'KAM': -0.33, 'KAN': 0.18, 'KAO': 0.88, 'KAP': -0.96, 'KAQ': 0.48, 'KAR': 0.24, 'EisenbergWeiss': [-1.42, 2.6]}

However, users may opt to ignore backbones contacts (C, CA, N, O) with the cavity when mapping hydrophobicity scales on surface points. Then, users must set ``ignore_backbone`` flag to ``True``.

.. code-block:: python

    >>> from pyKVFinder import hydropathy
    >>> scales, avg_hydropathy = hydropathy(surface, atomic, vertices, ignore_backbone=True)
    >>> scales
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
    >>> avg_hydropathy
    {'KAA': -0.7, 'KAB': 0.12, 'KAC': -0.08, 'KAD': -0.56, 'KAE': -0.28, 'KAF': -0.25, 'KAG': -0.28, 'KAH': -0.09, 'KAI': -0.4, 'KAJ': 0.96, 'KAK': -0.87, 'KAL': 0.23, 'KAM': 0.06, 'KAN': -0.1, 'KAO': 0.99, 'KAP': -1.04, 'KAQ': 0.48, 'KAR': -0.84, 'EisenbergWeiss': [-1.42, 2.6]}
