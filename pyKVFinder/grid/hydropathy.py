import os
import pathlib
from importlib import resources

import numpy

from .cavity import _get_cavity_label, _get_cavity_name, _select_cavities
from .geometry import _get_sincos


def _process_hydropathy(
    raw_avg_hydropathy: numpy.ndarray, ncav: int, selection: list[int] | None = None
) -> dict[str, float]:
    """Processes array of average hydropathy.

    Parameters
    ----------
    raw_avg_hydropathy : numpy.ndarray
        A numpy.ndarray of average hydropathy.
    ncav : int
        Number of cavities.
    selection : list[int] | None, optional
        A list of integer labels of each cavity to be selected.

    Returns
    -------
    avg_hydropathy : dict[str, float]
        A dictionary with average hydropathy per detected cavity.
    """
    avg_hydropathy = {}

    for index in range(ncav):
        key = _get_cavity_name(index)
        avg_hydropathy[key] = float(round(raw_avg_hydropathy[index], 2))

    if selection is not None:
        # Get keys from selection
        all_keys = list(avg_hydropathy.keys())
        keys = [all_keys[sele - 2] for sele in selection]

        # Get average hydropathy of selection
        avg_hydropathy = {key: avg_hydropathy[key] for key in keys}

    return avg_hydropathy


def hydropathy(
    surface: numpy.ndarray,
    atomic: numpy.ndarray | list[list[str | float | int]],
    vertices: numpy.ndarray | list[list[float]],
    step: float | int = 0.6,
    probe_in: float | int = 1.4,
    hydrophobicity_scale: str | pathlib.Path = "EisenbergWeiss",
    ignore_backbone: bool = False,
    selection: list[int] | list[str] | None = None,
    nthreads: int | None = None,
    verbose: bool = False,
) -> tuple[numpy.ndarray, dict[str, float]]:
    """Hydropathy characterization of the detected cavities.

    Map a target hydrophobicity scale per surface point and calculate average hydropathy of detected cavities.

    Parameters
    ----------
    surface : numpy.ndarray
        Surface points in the 3D grid (surface[nx][ny][nz]).
        Surface array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule or empty space points;

            * >=2: surface points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.
    atomic : Union[numpy.ndarray, List[List[Union[str, float, int]]]]
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    vertices : Union[numpy.ndarray, List[List[float]]]
        A numpy.ndarray or a list with xyz vertices coordinates (origin,
        X-axis, Y-axis, Z-axis).
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6.
    probe_in : Union[float, int], optional
        Probe In size (A), by default 1.4.
    hydrophobicity_scale : str, optional
        Name of a built-in hydrophobicity scale or a path to a
        TOML-formatted file with a custom hydrophobicity scale, by default
        `EisenbergWeiss`.
        The hydrophobicity scale file defines the name of the scale and the
        hydrophobicity value for each residue and when not defined, it assigns
        zero to the missing residues (see `Hydrophobicity scale file
        template`). The package contains seven built-in hydrophobicity scales:

            * EisenbergWeiss [1]_;
            * HessaHeijne [2]_;
            * KyteDoolittle [3]_;
            * MoonFleming [4]_;
            * RadzickaWolfenden [5]_;
            * WimleyWhite [6]_;
            * ZhaoLondon [7]_.

    ignore_backbone : bool, optional
        Whether to ignore backbone atoms (C, CA, N, O) when defining interface
        residues, by default False.
    selection : Union[List[int], List[str]], optional
        A list of integer labels or a list of cavity names to be selected, by default None.
    nthreads : int, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    scales : numpy.ndarray
        A numpy.ndarray with hydrophobicity scale value mapped at surface
        points (scales[nx][ny][nz]).
    avg_hydropathy : Dict[str, float]
        A dictionary with average hydropathy for each detected cavity and
        the range of the hydrophobicity scale (min, max).

    Raises
    ------
    TypeError
        `surface` must be a numpy.ndarray.
    ValueError
        `surface` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `atomic` must be a list or a numpy.ndarray.
    ValueError
        `atomic` has incorrect shape. It must be (n, 8).
    TypeError
        `vertices` must be a list or a numpy.ndarray.
    ValueError
        `vertices` has incorrect shape. It must be (4, 3).
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.
    TypeError
        `probe_in` must be a non-negative real number.
    ValueError
        `probe_in` must be a non-negative real number.
    TypeError
        `hydrophobicity_scale` must be a string or a pathlib.Path.
    TypeError
        `ignore_backbone` must be a boolean.
    TypeError
        `selection` must be a list of strings (cavity names) or integers (cavity labels).
    ValueError
        Invalid `selection`: {selection}.
    TypeError
        `nthreads` must be a positive integer.
    ValueError
        `nthreads` must be a positive integer.
    TypeError
        `verbose` must be a boolean.

    Note
    ----
    The cavity nomenclature is based on the integer label. The cavity marked
    with 2, the first integer corresponding to a cavity, is KAA, the cavity
    marked with 3 is KAB, the cavity marked with 4 is KAC and so on. After KZZ,
    the next cavity is Kaa, Kab, and so on. The naming convention supports up to
    1352 cavities.

    See Also
    --------
    read_pdb
    read_xyz
    spatial
    export
    write_results

    Example
    -------
    With the surface points identified with ``pyKVFinder.spatial`` and atomic coordinates and information read with ``pyKVFinder.read_pdb`` or ``pyKVFinder.read_xyz``, we can perform a hydropathy characterization, that maps a target hydrophobicity scale on surface points and calculate the average hydropathy

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

    References
    ----------
    .. [1] Eisenberg D, Weiss RM, Terwilliger TC. The hydrophobic moment
       detects periodicity in protein hydrophobicity. Proceedings of the
       National Academy of Sciences. 1984;81.

    .. [2] Hessa T, Kim H, Bihlmaier K, Lundin C, Boekel J, Andersson H, et al.
       Recognition of transmembrane helices by the endoplasmic reticulum
       translocon. Nature. 2005;433.

    .. [3] Kyte J, Doolittle RF. A simple method for displaying the hydropathic
       character of a protein. Journal of Molecular Biology. 1982;157.

    .. [4] Moon CP, Fleming KG. Side-chain hydrophobicity scale derived from
       transmembrane protein folding into lipid bilayers. Proceedings of the
       National Academy of Sciences. 2011;108.

    .. [5] Radzicka A, Wolfenden R. Comparing the polarities of the amino acids:
        Side-chain distribution coefficients between the vapor phase, cyclohexane,
        1-octanol, and neutral aqueous solution. Biochemistry. 1988;27.

    .. [6] Wimley WC, White SH. Experimentally determined hydrophobicity scale
       for proteins at membrane interfaces. Nature Structural & Molecular
       Biology. 1996;3.

    .. [7] Zhao G, London E. An amino acid “transmembrane tendency” scale that
       approaches the theoretical limit to accuracy for prediction of
       transmembrane helices: Relationship to biological hydrophobicity.
       Protein Science. 2006;15.
    """
    import tomlkit
    from pyKVFinder._pyKVFinder import _hydropathy

    # Check arguments
    if type(surface) not in [numpy.ndarray]:
        raise TypeError("`surface` must be a numpy.ndarray.")
    elif len(surface.shape) != 3:
        raise ValueError("`surface` has the incorrect shape. It must be (nx, ny, nz).")
    if type(atomic) not in [numpy.ndarray, list]:
        raise TypeError("`atomic` must be a list or a numpy.ndarray.")
    elif len(numpy.asarray(atomic).shape) != 2:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    elif numpy.asarray(atomic).shape[1] != 8:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    if type(vertices) not in [numpy.ndarray, list]:
        raise TypeError("`vertices` must be a list or a numpy.ndarray.")
    elif numpy.asarray(vertices).shape != (4, 3):
        raise ValueError("`vertices` has incorrect shape. It must be (4, 3).")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")
    if type(probe_in) not in [float, int]:
        raise TypeError("`probe_in` must be a non-negative real number.")
    elif probe_in < 0.0:
        raise ValueError("`probe_in` must be a non-negative real number.")
    if type(hydrophobicity_scale) not in [str, pathlib.Path]:
        raise TypeError("`hydrophobicity_scale` must be a string or a pathlib.Path.")
    if type(ignore_backbone) not in [bool]:
        raise TypeError("`ignore_backbone` must be a boolean.")
    if selection is not None:
        # Check selection types
        if all(isinstance(x, int) for x in selection):
            pass
        elif all(isinstance(x, str) for x in selection):
            selection = [_get_cavity_label(sele) for sele in selection]
        else:
            raise TypeError(
                "`selection` must be a list of strings (cavity names) or integers (cavity labels)."
            )
        # Check if selection includes valid cavity labels
        if any(x < 2 for x in selection):
            raise ValueError(f"Invalid `selection`: {selection}.")
    if nthreads is None:
        nthreads = os.cpu_count() - 1
    else:
        if type(nthreads) not in [int]:
            raise TypeError("`nthreads` must be a positive integer.")
        elif nthreads <= 0:
            raise ValueError("`nthreads` must be a positive integer.")
    if type(verbose) not in [bool]:
        raise TypeError("`verbose` must be a boolean.")

    # Convert types
    if isinstance(atomic, list):
        atomic = numpy.asarray(atomic)
    if isinstance(vertices, list):
        vertices = numpy.asarray(vertices)
    if isinstance(step, int):
        step = float(step)
    if isinstance(probe_in, int):
        probe_in = float(probe_in)

    # Convert numpy.ndarray data types
    surface = surface.astype("int32") if surface.dtype != "int32" else surface
    vertices = vertices.astype("float64") if vertices.dtype != "float64" else vertices

    # Extract atominfo from atomic
    atominfo = numpy.asarray(
        ([[f"{atom[0]}_{atom[1]}_{atom[2]}", atom[3]] for atom in atomic[:, :4]])
    )

    # Extract xyzr from atomic
    xyzr = atomic[:, 4:].astype(numpy.float64)

    # Get number of cavities
    ncav = int(surface.max() - 1)

    # Get sincos: sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb)
    sincos = _get_sincos(vertices)

    # Select cavities
    if selection is not None:
        surface = _select_cavities(surface, selection)

    # Get dimensions
    nx, ny, nz = surface.shape
    nvoxels = nx * ny * nz

    # Load hydrophobicity scales
    if hydrophobicity_scale in [
        "EisenbergWeiss",
        "HessaHeijne",
        "KyteDoolittle",
        "MoonFleming",
        "RadzickaWolfenden",
        "WimleyWhite",
        "ZhaoLondon",
    ]:
        hydrophobicity_scale = (
            resources.files("pyKVFinder.data") / f"{hydrophobicity_scale}.toml"
        )
    with open(hydrophobicity_scale, "r") as file:
        f = tomlkit.load(file)
    data, name = list(f.values())[0], list(f.keys())[0]
    resn, scale = list(data.keys()), numpy.asarray(list(data.values()))

    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Remove backbone from atominfo
    if ignore_backbone:
        mask = numpy.where(
            (atominfo[:, 1] != "C")
            & (atominfo[:, 1] != "CA")
            & (atominfo[:, 1] != "N")
            & (atominfo[:, 1] != "O")
        )
        atominfo = atominfo[mask[0],]
        xyzr = xyzr[mask[0],]

    # Get residue name from atominfo
    resname = list(map(lambda x: x.split("_")[2], atominfo[:, 0]))

    # Get hydrophobicity scales in 3D grid and average hydropathy
    scales, avg_hydropathy = _hydropathy(
        nvoxels,
        ncav,
        surface,
        xyzr,
        P1,
        sincos,
        resname,
        resn,
        scale,
        step,
        probe_in,
        nthreads,
        verbose,
    )
    avg_hydropathy = _process_hydropathy(avg_hydropathy, ncav, selection)
    avg_hydropathy[f"{name}"] = [float(scale.min()), float(scale.max())]

    return scales.reshape(nx, ny, nz), avg_hydropathy
