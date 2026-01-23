import os
from typing import Optional

import numpy

__all__ = ["detect"]

from .geometry import _get_dimensions, _get_sincos


def detect(
    atomic: numpy.ndarray | list[list[str | float | int]],
    vertices: numpy.ndarray | list[list[float]],
    step: float | int = 0.6,
    probe_in: float | int = 1.4,
    probe_out: float | int = 4.0,
    removal_distance: float | int = 2.4,
    volume_cutoff: float | int = 5.0,
    latomic: Optional[numpy.ndarray | list[list[float]]] = None,
    ligand_cutoff: float | int = 5.0,
    box_adjustment: bool = False,
    surface: str = "SES",
    nthreads: int | None = None,
    verbose: bool = False,
) -> tuple[int, numpy.ndarray]:
    """Detects biomolecular cavities.

    Cavity points that belongs to the same cavity are assigned with an integer
    in the grid.

    Parameters
    ----------
    atomic : numpy.ndarray | list[list[Union[str, float, int]]]
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    vertices : numpy.ndarray | list[list[float]]
        A numpy.ndarray or a list with xyz vertices coordinates (origin,
        X-axis, Y-axis, Z-axis).
    step : float | int, optional
        Grid spacing (A), by default 0.6.
    probe_in : float | int, optional
        Probe In size (A), by default 1.4.
    probe_out : float | int, optional
        Probe Out size (A), by default 4.0.
    removal_distance : float | int, optional
        A length to be removed from the cavity-bulk frontier (A), by
        default 2.4.
    volume_cutoff : float | int, optional
        Volume filter for detected cavities (A3), by default 5.0.
    latomic : numpy.ndarray | list[list[Union[str, float, int]]], optional
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom of a target ligand, by default None.
    ligand_cutoff : float | int, optional
        A radius to limit a space around a ligand (A), by default 5.0.
    box_adjustment : bool, optional
        Whether a custom 3D grid is applied, by default False.
    surface : str, optional
        Surface representation. Keywords options are SES (Solvent Excluded Surface) or SAS (Solvent
        Accessible Surface), by default SES.
    nthreads : int | None, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    ncav : int
        Number of cavities.
    cavities : numpy.ndarray
        Cavity points in the 3D grid (cavities[nx][ny][nz]).
        Cavities array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule points;

            * 1: empty space points;

            * >=2: cavity points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.

    Raises
    ------
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
        `probe_out` must be a non-negative real number.
    ValueError
        `probe_out` must be a non-negative real number.
    ValueError
        `probe_out` must be greater than `probe_in`.
    TypeError
        `removal_distance` must be a non-negative real number.
    ValueError
        `removal_distance` must be a non-negative real number.
    TypeError
        `volume_cutoff` must be a non-negative real number.
    ValueError
        `volume_cutoff` must be a non-negative real number.
    TypeError
        `latomic` must be a list, a numpy.ndarray or None.
    ValueError
        `latomic` has incorrect shape. It must be (n, 8).
    TypeError
        `ligand_cutoff` must be a positive real number.
    ValueError
        `ligand_cutoff` must be a positive real number.
    TypeError
        `box_adjustment` must be a boolean.
    TypeError
        `surface` must be a string.
    TypeError
        `nthreads` must be a positive integer.
    ValueError
        `nthreads` must be a positive integer.
    ValueError
        `surface` must be SAS or SES, not `surface`.

    See also
    --------
    read_pdb
    read_xyz
    get_vertices
    get_vertices_from_file
    spatial
    depth
    constitutional
    export

    Warning
    -------
    If you are using box adjustment mode, do not forget to set box_adjustment
    flag to True and read the box configuration file with 'get_vertices_from_file'
    function.

    Warning
    -------
    If you are using ligand adjustment mode, do not forget to read ligand atom
    coordinates with 'read_pdb' or 'read_xyz' functions.

    Example
    -------
    With the grid vertices defined with ``get_vertices`` and atomic data loaded with ``read_pdb`` or ``read_xyz``, we can detect cavities on the whole target biomolecule:

    >>> from pyKVFinder import detect
    >>> ncav, cavities = detect(atomic, vertices)
    >>> ncav
    18
    >>> cavities
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

    However, users may opt to perform cavity detection in a segmented space through ligand adjustment and/or box adjustment modes.

    The cavity detection can be limited around the target ligand(s), which will be passed to pyKVFinder through a *.pdb* or a *.xyz* files. Thus, the detected cavities are limited within a radius (``ligand_cutoff``) of the target ligand(s).

    >>> import os
    >>> ligand = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', 'ADN.pdb')
    >>> from pyKVFinder import read_pdb
    >>> latomic = read_pdb(ligand)
    >>> ncav, cavities = detect(atomic, vertices, latomic=latomic, ligand_cutoff=5.0)
    >>> ncav
    2
    >>> cavities
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

    Further, we can also perform cavity detection on a custom 3D grid, where we can explore closed regions with a custom box, which can be defined by a *.toml* file (see `Box configuration file template`).

    >>> import os
    >>> fn = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', 'custom-box.toml')
    >>> with open(fn, 'r') as f:
    ...     print(f.read())
    [box]
    p1 = [3.11, 7.34, 1.59]
    p2 = [11.51, 7.34, 1.59]
    p3 = [3.11, 10.74, 1.59]
    p4 = [3.11, 7.34, 6.19]

    With this box adjustment mode, we must defined the 3D grid with ``get_vertices_from_file``.

    >>> from pyKVFinder import get_vertices_from_file
    >>> vertices, atomic = get_vertices_from_file(fn, atomic)

    Then, we can perform cavity detection:

    >>> ncav, cavities = detect(atomic, vertices, box_adjustment=True)
    >>> ncav
    1
    >>> cavities
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

    .. warning::
        If you are using box adjusment mode, do not forget to set ``box_adjustment`` flag to ``True``.
    """
    from pyKVFinder._pyKVFinder import _detect, _detect_ladj

    # Check arguments
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
    if type(probe_out) not in [float, int]:
        raise TypeError("`probe_out` must be a non-negative real number.")
    elif probe_out < 0.0:
        raise ValueError("`probe_out` must be a non-negative real number.")
    elif probe_out < probe_in:
        raise ValueError("`probe_out` must be greater than `probe_in`.")
    if type(removal_distance) not in [float, int]:
        raise TypeError("`removal_distance` must be a non-negative real number.")
    elif removal_distance < 0.0:
        raise ValueError("`removal_distance` must be a non-negative real number.")
    if type(volume_cutoff) not in [float, int]:
        raise TypeError("`volume_cutoff` must be a non-negative real number.")
    elif volume_cutoff < 0.0:
        raise ValueError("`volume_cutoff` must be a non-negative real number.")
    if latomic is not None:
        if type(latomic) not in [numpy.ndarray, list]:
            raise TypeError("`latomic` must be a list, a numpy.ndarray or None.")
        if len(numpy.asarray(latomic).shape) != 2:
            raise ValueError("`latomic` has incorrect shape. It must be (n, 8).")
        elif numpy.asarray(latomic).shape[1] != 8:
            raise ValueError("`latomic` has incorrect shape. It must be (n, 8).")
    if type(ligand_cutoff) not in [float, int]:
        raise TypeError("`ligand_cutoff` must be a positive real number.")
    elif ligand_cutoff <= 0.0:
        raise ValueError("`ligand_cutoff` must be a positive real number.")
    if type(box_adjustment) not in [bool]:
        raise TypeError("`box_adjustment` must be a boolean.")
    if type(surface) not in [str]:
        raise TypeError("`surface` must be a string.")
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
    if isinstance(probe_out, int):
        probe_out = float(probe_out)
    if isinstance(removal_distance, int):
        removal_distance = float(removal_distance)
    if isinstance(volume_cutoff, int):
        volume_cutoff = float(volume_cutoff)
    if isinstance(ligand_cutoff, int):
        ligand_cutoff = float(ligand_cutoff)
    if isinstance(latomic, list):
        latomic = numpy.asarray(latomic)

    # Convert numpy.ndarray data types
    vertices = vertices.astype("float64") if vertices.dtype != "float64" else vertices

    # Extract xyzr from atomic
    xyzr = atomic[:, 4:].astype(numpy.float64)

    # Extract lxyzr from latomic
    if latomic is not None:
        lxyzr = latomic[:, 4:].astype(numpy.float64)

    # Define ligand adjustment mode
    ligand_adjustment = True if latomic is not None else False

    # Unpack vertices
    P1, P2, P3, P4 = vertices

    if surface == "SES":
        if verbose:
            print("> Surface representation: Solvent Excluded Surface (SES)")
        surface = True
    elif surface == "SAS":
        if verbose:
            print("> Surface representation: Solvent Accessible Surface (SAS)")
        surface = False
    else:
        raise ValueError(f"`surface` must be SAS or SES, not {surface}.")

    # Get sincos: sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb)
    sincos = _get_sincos(vertices)

    # Get dimensions
    nx, ny, nz = _get_dimensions(vertices, step)

    # Calculate number of voxels
    nvoxels = nx * ny * nz

    # Detect cavities
    if ligand_adjustment:
        ncav, cavities = _detect_ladj(
            nvoxels,
            nx,
            ny,
            nz,
            xyzr,
            lxyzr,
            P1,
            sincos,
            step,
            probe_in,
            probe_out,
            removal_distance,
            volume_cutoff,
            ligand_adjustment,
            ligand_cutoff,
            box_adjustment,
            P2,
            surface,
            nthreads,
            verbose,
        )
    else:
        ncav, cavities = _detect(
            nvoxels,
            nx,
            ny,
            nz,
            xyzr,
            P1,
            sincos,
            step,
            probe_in,
            probe_out,
            removal_distance,
            volume_cutoff,
            box_adjustment,
            P2,
            surface,
            nthreads,
            verbose,
        )

    return ncav, cavities.reshape(nx, ny, nz)
