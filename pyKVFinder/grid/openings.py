import os

import numpy

from .cavity import _get_cavity_label, _get_cavity_name, _select_cavities
from .depth import depth


def _get_opening_name(index: int) -> str:
    """Get opening name, eg OAA, OAB, and so on, based on the index.

    Naming convention:
    - 0-675   -> OAA ... OZZ
    - 676-1351 -> Oaa ... Ozz

    Parameters
    ----------
    index : int
        Index in the dictionary.

    Returns
    -------
    opening_name : str
        Opening name
    """
    # Get block and offset
    block = index // (26 * 26)
    offset = index % (26 * 26)

    # Choose ASCII base: uppercase or lowercase
    base = 65 if block == 0 else 97  # 'A' or 'a'

    # Get first and second characters
    first = chr(base + (offset // 26))
    second = chr(base + (offset % 26))

    # Build opening name
    opening_name = f"O{first}{second}"

    return opening_name


def _get_opening_label(opening_name: str) -> int:
    """Get opening label, eg 2, 3, and so on, based on the opening name.

    Naming convention:
    - OAA-OZZ -> 2-677
    - Oaa-Ozz -> 678-1353

    Parameters
    ----------
    opening_name : str
        Opening name

    Returns
    -------
    copening_label : int
        Integer label of each opening.

    Raises
    ------
    ValueError
        Invalid opening name: `opening_name`.
    """
    # Check opening name
    if opening_name[0] != "O":
        raise ValueError(f"Invalid opening name: {opening_name}.")

    # Get opening label
    o1, o2 = opening_name[1], opening_name[2]

    # Uppercase block: KAA–KZZ
    if o1.isupper() and o2.isupper():
        base = 65  # 'A'
        block_offset = 0

    # Lowercase block: Kaa–Kzz
    elif o1.islower() and o2.islower():
        base = 97  # 'a'
        block_offset = 26 * 26

    # +2 preserves the original labeling convention
    opening_label = (block_offset + (ord(o1) - base) * 26 + (ord(o2) - base)) + 2

    return opening_label


def _process_openings(
    raw_openings: numpy.ndarray,
    opening2cavity: numpy.ndarray,
) -> dict[str, dict[str, float]]:
    """Processes arrays of openings' areas.

    Parameters
    ----------
    raw_openings : numpy.ndarray
        A numpy.ndarray of openings' areas.
    openings2cavity: numpy.ndarray
        A numpy.ndarray of openings as indexes and cavities as values.

    Returns
    -------
    area : dict[str, dict[str, float]]
        A dictionary with area of each detected opening.
    """
    area = {}

    # Get number of openings
    nopenings = raw_openings.shape[0]

    for index in range(nopenings):
        # Get opening name
        opening = _get_opening_name(index)

        # Get cavity name
        cavity = _get_cavity_name(opening2cavity[index])

        # Save opening area
        if cavity not in area.keys():
            area[cavity] = {}
        area[cavity][opening] = float(round(raw_openings[index], 2))

    # Sort keys
    area = dict(sorted(area.items()))

    return area


def openings(
    cavities: numpy.ndarray,
    depths: numpy.ndarray | None = None,
    step: float | int = 0.6,
    openings_cutoff: int = 1,
    selection: list[int] | list[str] | None = None,
    nthreads: int | None = None,
    verbose: bool = False,
) -> tuple[int, numpy.ndarray, dict[str, dict[str, float]]]:
    """Identify openings of the detected cavities and calculate their areas.

    Parameters
    ----------
    cavities : numpy.ndarray
        Cavity points in the 3D grid (cavities[nx][ny][nz]).
        Cavities array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule points;

            * 1: empty space points;

            * >=2: cavity points.
    depths : numpy.ndarray | None, optional
        A numpy.ndarray with depth of cavity points (depth[nx][ny][nz]), by default None. If None, depths is calculated from cavities.
    step : float | int, optional
        Grid spacing (A), by default 0.6.
    openings_cutoff : int, optional
        The minimum number of voxels an opening must have, by default 1.
    selection : list[int] | list[str], optional
        A list of integer labels or a list of cavity names to be selected, by default None.
    nthreads : int, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    nopenings : int
        Total number of openings.
    openings : numpy.ndarray
        Openings points in the 3D grid (openings[nx][ny][nz]).
        Openings array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: cavity or biomolecule points;

            * 1: empty space points;

            * >=2: Opening points.

        The empty space points are regions that do not meet the chosen
        openings cutoff to be considered an opening.
    aopenings : dict[str, dict[str, float]]
        A dictionary with area of each detected opening.

    Raises
    ------
    TypeError
        `cavities` must be a numpy.ndarray.
    ValueError
        `cavities` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `depths` must be a numpy.ndarray.
    ValueError
        `depths` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.
    TypeError
        `openings_cutoff` must be an integer.
    ValueError
        `openings_cutoff` must be a positive integer.
    TypeError
        `selection` must be a list of strings (cavity names) or integers (cavity labels).
    ValueError
        Invalid `selection`: `selection`.
    TypeError
        `nthreads` must be a positive integer.
    ValueError
        `nthreads` must be a positive integer.
    TypeError
        `verbose` must be a boolean

    Example
    -------
    With the cavity points identified with ``detect``, we can characterize their openings, that includes number and area of openings and defining opening points:

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

    >>> from pyKVFinder import export
    >>> export("cavities_with_openings.pdb", cavities, None, vertices, B=openings)

    * Export opening points with same nomenclature from ``aopenings``:

    >>> from pyKVFinder import export_openings
    >>> export_openings("openings.pdb", openings, vertices)

    See Also
    --------
    detect
    depth
    export
    export_openings
    """
    from pyKVFinder._pyKVFinder import _area, _openings, _openings2cavities

    # Check arguments
    if type(cavities) not in [numpy.ndarray]:
        raise TypeError("`cavities` must be a numpy.ndarray.")
    elif len(cavities.shape) != 3:
        raise ValueError("`cavities` has the incorrect shape. It must be (nx, ny, nz).")
    if depths is None:
        depths, _, _ = depth(cavities, step, selection, nthreads, verbose)
    elif type(depths) not in [numpy.ndarray]:
        raise TypeError("`depths` must be a numpy.ndarray.")
    elif len(depths.shape) != 3:
        raise ValueError("`depths` has the incorrect shape. It must be (nx, ny, nz).")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")
    if type(openings_cutoff) not in [int]:
        raise TypeError("`openings_cutoff` must be an integer.")
    elif openings_cutoff < 0:
        raise ValueError("`openings_cutoff` must be a positive integer.")
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

    # Select cavities
    if selection is not None:
        cavities = _select_cavities(cavities, selection)

    # Get number of cavities
    ncav = int(cavities.max() - 1)

    # Get cavities shape
    nx, ny, nz = cavities.shape

    # Find openings
    nopenings, openings = _openings(
        nx * ny * nz, cavities, depths, ncav, openings_cutoff, step, nthreads, verbose
    )
    if nopenings > 1352:
        print("Warning: The number of openings exceeds the maximum supported (1352). ")

    # Reshape openings
    openings = openings.reshape(nx, ny, nz)

    # Calculate openings areas
    if verbose:
        print("> Estimating openings area")
    aopenings = _area(openings, step, nopenings, nthreads)

    # Find which openings belongs to each cavity
    opening2cavity = _openings2cavities(nopenings, cavities, openings, nthreads)

    # Process openings
    aopenings = _process_openings(aopenings, opening2cavity)

    return nopenings, openings, aopenings
