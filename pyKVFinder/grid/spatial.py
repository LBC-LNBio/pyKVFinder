import os

import numpy

from .cavity import _get_cavity_label, _get_cavity_name, _select_cavities


def _process_spatial(
    raw_volume: numpy.ndarray,
    raw_area: numpy.ndarray,
    ncav: int,
    selection: list[int] | None = None,
) -> tuple[dict[str, float], dict[str, float]]:
    """Processes arrays of volumes and areas.

    Parameters
    ----------
    raw_volume : numpy.ndarray
        A numpy.ndarray of volumes.
    raw_area : numpy.ndarray
        A numpy.ndarray of areas.
    ncav : int
        Number of cavities.
    selection : list[int] | None, optional
        A list of integer labels of each cavity to be selected, by default None.

    Returns
    -------
    volume : dict[str, float]
        A dictionary with volume of each detected cavity.
    area : dict[str, float]
        A dictionary with area of each detected cavity.
    """
    volume, area = {}, {}

    # Prepare volume and area dictionary
    for index in range(ncav):
        key = _get_cavity_name(index)
        volume[key] = float(round(raw_volume[index], 2))
        area[key] = float(round(raw_area[index], 2))

    if selection is not None:
        # Get keys from selection
        all_keys = list(volume.keys())
        keys = [all_keys[sele - 2] for sele in selection]

        # Get volume and area of selection
        volume = {key: volume[key] for key in keys}
        area = {key: area[key] for key in keys}

    return volume, area


def spatial(
    cavities: numpy.ndarray,
    step: float | int = 0.6,
    selection: list[int] | list[str] | None = None,
    nthreads: int | None = None,
    verbose: bool = False,
) -> tuple[numpy.ndarray, dict[str, float], dict[str, float]]:
    """Spatial characterization (volume and area) of the detected cavities.

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
        Grid spacing (A), by default 0.6.
    selection : list[int] | list[str] | None, optional
        A list of integer labels or a list of cavity names to be selected, by default None.
    nthreads : int | None, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    surface : numpy.ndarray
        Surface points in the 3D grid (surface[nx][ny][nz]).
        Surface array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule or empty space points;

            * >=2: surface points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.
    volume : dict[str, float]
        A dictionary with volume of each detected cavity.
    area : dict[str, float]
        A dictionary with area of each detected cavity.

    Note
    ----
    The cavity nomenclature is based on the integer label. The cavity marked
    with 2, the first integer corresponding to a cavity, is KAA, the cavity
    marked with 3 is KAB, the cavity marked with 4 is KAC and so on. After KZZ,
    the next cavity is Kaa, Kab, and so on. The naming convention supports up to
    1352 cavities.

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

    See Also
    --------
    detect
    hydropathy
    export

    Example
    -------
    With the cavity points identified with ``detect``, we can perform a spatial characterization, that includes volume, area and defining surface points:

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
    """
    from pyKVFinder._pyKVFinder import _spatial

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

    # Convert numpy.ndarray data types
    cavities = cavities.astype("int32") if cavities.dtype != "int32" else cavities

    # Get number of cavities
    ncav = int(cavities.max() - 1)

    # Select cavities
    if selection is not None:
        cavities = _select_cavities(cavities, selection)

    # Get cavities shape
    nx, ny, nz = cavities.shape

    # Get surface points, volume and area
    surface, volume, area = _spatial(
        cavities, nx * ny * nz, ncav, ncav, step, nthreads, verbose
    )
    volume, area = _process_spatial(volume, area, ncav)

    return surface.reshape(nx, ny, nz), volume, area
