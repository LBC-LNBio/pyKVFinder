pyKVFinder.hydropathy
=====================

Hydropathy characterization of the detected cavities. 

Map a target hydrophobicity scale per surface point and calculate average hydropathy of detected cavities.

.. code-block:: python

    pyKVFinder.hydropathy(surface, atominfo, xyzr, vertices, sincos, ncav, step = 0.6, probe_in = 1.4, hydrophobicity_scale = 'EisenbergWeiss', ignore_backbone = False, nthreads = os.cpu_count() - 1, verbose = False)

:Args:

    ``surface`` : *numpy.ndarray*
        A numpy array with integer labels in each position, that are:
            * -1: bulk point.
            * 0: biomolecule or empty space point.
            * >=2: surface point.
    ``atominfo`` : *numpy.ndarray*
        A numpy array with atomic information (residue number, chain, residue name, atom name)
    ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz atomic coordinates and radii values (x, y, z, radius) 
    ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)
    ``sincos`` : *numpy.ndarray*
        A numpy array with sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb)
    ``ncav`` : *int*
        Number of cavities
    ``step`` : *float, default 0.6*
        Grid spacing (A)
    ``probe_in`` : *float, default 1.4*
        Probe In size (A)
    ``hydrophobicity_scale`` : *str, default 'EisenbergWeiss'*
        Name of a built-in hydrophobicity scale (EisenbergWeiss, HessaHeijne, KyteDoolittle, MoonFleming, WimleyWhite, ZhaoLondon) or a path to a TOML-formatted file with a custom hydrophobicity scale file.
    ``ignore_backbone`` :  *bool, default False*
        Whether to ignore backbone atoms (C, CA, N, O) when defining interface residues
    ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads
    ``verbose`` : *bool, default False*
        Print extra information to standard output

:Returns:

    ``scales``: *numpy.ndarray*
        A numpy array with hydrophobicity scale value mapped at surface points (scales[nx][ny][nz])
    ``avg_hydropathy``: *dict*
        A dictionary with average hydropathy of each detected cavity and the range of the hydrophobicity scale (min, max)

.. note::

    The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on. 

    The hydrophobicity scale file defines the name of the scale and the hydrophobicity value for each residue and when not defined, it assigns zero to the missing residues (see `Hydrophobicity scale file template <hydrophobicity_scale_file_template.html>`_). The package contains six built-in hydrophobicity scales: Eisenberg & Weiss [1], Hessa & Heijne [2], Kyte & Doolittle [3], Moon & Fleming [4], Wimley & White [5] and Zhao & London [6].

.. seealso::

    * `pyKVFinder.spatial <spatial.html>`_
    * `pyKVFinder.export <export.html>`_
    * `pyKVFinder.write_results <write_results.html>`_
    * `Hydrophobicity scale file template <../_cfg_files/hydrophobicity_scale_file_template.html>`_

.. raw:: html

    <h4><u>References</u></h4>

1. Eisenberg D, Weiss RM, Terwilliger TC. The hydrophobic moment detects periodicity in protein hydrophobicity. Proceedings of the National Academy of Sciences. 1984;81. 

2. Hessa T, Kim H, Bihlmaier K, Lundin C, Boekel J, Andersson H, et al. Recognition of transmembrane helices by the endoplasmic reticulum translocon. Nature. 2005;433. 

3. Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. Journal of Molecular Biology. 1982;157. 

4. Moon CP, Fleming KG. Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers. Proceedings of the National Academy of Sciences. 2011;108. 

5. Wimley WC, White SH. Experimentally determined hydrophobicity scale for proteins at membrane interfaces. Nature Structural & Molecular Biology. 1996;3. 

6. Zhao G, London E. An amino acid “transmembrane tendency” scale that approaches the theoretical limit to accuracy for prediction of transmembrane helices: Relationship to biological hydrophobicity. Protein Science. 2006;15. 

.. raw:: html

  <h4><u>Example</u></h4>

With the surface points identified with ``pyKVFinder.spatial`` and atomic coordinates and information read with ``pyKVFinder.read_pdb``, we can perform a hydropathy characterization, that maps a target hydrophobicity scale on surface points and calculate the average hydropathy

.. code-block:: python

    >>> from pyKVFinder import hydropathy
    >>> scales, avg_hydropathy = hydropathy(surface, atominfo, xyzr, vertices, sincos, ncav)
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
    {'KAA': -0.71, 'KAB': -0.06, 'KAC': -0.07, 'KAD': -0.62, 'KAE': -0.81, 'KAF': -0.14, 'KAG': -0.33, 'KAH': -0.16, 'KAI': -0.4, 'KAJ': 0.62, 'KAK': -0.99, 'KAL': 0.35, 'KAM': -0.33, 'KAN': 0.18, 'KAO': 0.88, 'KAP': -0.96, 'KAQ': 0.48, 'KAR': 0.24, 'EisenbergWeiss': [-1.42, 2.6]

However, users may opt to ignore backbones contacts (C, CA, N, O) with the cavity when mapping hydrophobicity scales on surface points. Then, users must set ``ignore_backbone`` flag to ``True``.

.. code-block:: python

    >>> from pyKVFinder import hydropathy
    >>> scales, avg_hydropathy = hydropathy(surface, atominfo, xyzr, vertices, sincos, ncav, ignore_backbone=True)
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
    {'KAA': -0.69, 'KAB': 0.11, 'KAC': -0.08, 'KAD': -0.56, 'KAE': -0.28, 'KAF': -0.25, 'KAG': -0.28, 'KAH': -0.14, 'KAI': -0.4, 'KAJ': 0.97, 'KAK': -0.87, 'KAL': 0.22, 'KAM': 0.06, 'KAN': -0.1, 'KAO': 0.99, 'KAP': -1.04, 'KAQ': 0.48, 'KAR': -0.84, 'EisenbergWeiss': [-1.42, 2.6]}
