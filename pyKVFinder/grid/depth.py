import os

import numpy

from .cavity import _get_cavity_label, _get_cavity_name, _select_cavities


def _process_depth(
    raw_max_depth: numpy.ndarray,
    raw_avg_depth: numpy.ndarray,
    ncav: int,
    selection: list[int] | None = None,
) -> tuple[dict[str, float], dict[str, float]]:
    """Processes arrays of maximum and average depths.

    Parameters
    ----------
    raw_max_depth : numpy.ndarray
        A numpy.ndarray of maximum depth.
    raw_avg_depth : numpy.ndarray
        A numpy.ndarray of average depth.
    ncav : int
        Number of cavities.
    selection : list[int]| None, optional
        A list of integer labels of each cavity to be selected, by default None.

    Returns
    -------
    max_depth : dict[str, float]
        A dictionary with maximum depth of each detected cavity.
    avg_depth : dict[str, float]
        A dictionary with average depth of each detected cavity.
    """
    max_depth, avg_depth = {}, {}

    # Prepare maximum and average depth dictionary
    for index in range(ncav):
        key = _get_cavity_name(index)
        max_depth[key] = float(round(raw_max_depth[index], 2))
        avg_depth[key] = float(round(raw_avg_depth[index], 2))

    if selection is not None:
        # Get keys from selection
        all_keys = list(max_depth.keys())
        keys = [all_keys[sele - 2] for sele in selection]

        # Get volume and area of selection
        max_depth = {key: max_depth[key] for key in keys}
        avg_depth = {key: avg_depth[key] for key in keys}

    return max_depth, avg_depth


def depth(
    cavities: numpy.ndarray,
    step: float | int = 0.6,
    selection: list[int] | list[str] | None = None,
    nthreads: int | None = None,
    verbose: bool = False,
) -> tuple[numpy.ndarray, dict[str, float], dict[str, float]]:
    """Characterization of the depth of the detected cavities, including depth
    per cavity point and maximum and average depths of detected cavities.

    Parameters
    ----------
    cavities : numpy.ndarray
        Cavity points in the 3D grid (cavities[nx][ny][nz]).
        Cavities array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule points;

            * 1: empty space points;

            * >=2: cavity points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.
    step : float | int, optional
        Grid spacing (A).
    selection : list[int] | list[str] | None, optional
        A list of integer labels or a list of cavity names to be selected, by default None.
    nthreads : int | None, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    depths : numpy.ndarray
        A numpy.ndarray with depth of cavity points (depth[nx][ny][nz]).
    max_depth : dict[str, float]
        A dictionary with maximum depth of each detected cavity.
    avg_depth : dict[str, float]
        A dictionary with average depth of each detected cavity.

    Raises
    ------
    TypeError
        `cavities` must be a numpy.ndarray.
    ValueError
        `cavities` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.
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
    detect
    export
    write_results

    Example
    -------
    With the cavity points identified with ``detect``, we can perform a depth characterization, that includes maximum depth, average depth and defining depth of cavity points:

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
    """
    from pyKVFinder._pyKVFinder import _depth

    # Check arguments
    if type(cavities) not in [numpy.ndarray]:
        raise TypeError("`cavities` must be a numpy.ndarray.")
    elif len(cavities.shape) != 3:
        raise ValueError("`cavities` has the incorrect shape. It must be (nx, ny, nz).")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")
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
    if isinstance(step, int):
        step = float(step)

    # Check numpy.ndarray data types
    cavities = cavities.astype("int32") if cavities.dtype != "int32" else cavities

    # Get number of cavities
    ncav = int(cavities.max() - 1)

    # Select cavities
    if selection is not None:
        cavities = _select_cavities(cavities, selection)

    # Get cavities shape
    nx, ny, nz = cavities.shape

    # Get depth of cavity points, maximum depth and average depth
    depths, max_depth, avg_depth = _depth(
        cavities, nx * ny * nz, ncav, ncav, step, nthreads, verbose
    )
    max_depth, avg_depth = _process_depth(max_depth, avg_depth, ncav, selection)

    return depths.reshape(nx, ny, nz), max_depth, avg_depth
