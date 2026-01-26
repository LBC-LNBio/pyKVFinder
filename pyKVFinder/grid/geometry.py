import os
import pathlib

import numpy


def get_vertices(
    atomic: numpy.ndarray | list[list[str | float | int]],
    probe_out: float | int = 4.0,
    step: float | int = 0.6,
) -> numpy.ndarray:
    """Gets 3D grid vertices.

    Parameters
    ----------
    atomic : numpy.ndarray | list[list[str | float | int]]
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    probe_out : float | int, optional
        Probe Out size (A), by default 4.0.
    step : float | int, optional
        Grid spacing (A), by default 0.6.

    Returns
    -------
    vertices : numpy.ndarray
        A numpy.ndarray with xyz vertices coordinates
        (origin, X-axis, Y-axis, Z-axis).

    Raises
    ------
    TypeError
        `atomic` must be a list or a numpy.ndarray.
    ValueError
        `atomic` has incorrect shape. It must be (n, 8).
    TypeError
        `probe_out` must be a non-negative real number.
    ValueError
        `probe_out` must be a non-negative real number.
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.

    See Also
    --------
    read_pdb
    read_xyz
    detect
    constitutional
    hydropathy
    export

    Example
    -------
    With the atomic data read with ``read_pdb`` or ``read_xyz``, we can get the coordinates of 3D grid vertices (origin, X-axis, Y-axis, Z-axis):

    >>> from pyKVFinder import get_vertices
    >>> vertices = get_vertices(atomic)
    >>> vertices
    array([[-19.911, -32.125, -30.806],
        [ 40.188, -32.125, -30.806],
        [-19.911,  43.446, -30.806],
        [-19.911, -32.125,  27.352]])
    """
    # Check arguments types
    if type(atomic) not in [numpy.ndarray, list]:
        raise TypeError("`atomic` must be a list or a numpy.ndarray.")
    elif len(numpy.asarray(atomic).shape) != 2:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    elif numpy.asarray(atomic).shape[1] != 8:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    if type(probe_out) not in [int, float, numpy.float64]:
        raise TypeError("`probe_out` must be a non-negative real number.")
    elif probe_out < 0.0:
        raise ValueError("`probe_out` must be a non-negative real number.")
    if type(step) not in [int, float, numpy.float64]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")

    # Convert atomic type
    if isinstance(atomic, list):
        atomic = numpy.asarray(atomic)

    # Extract xyzr from atomic
    xyzr = atomic[:, 4:].astype(numpy.float64)

    # Prepare vertices
    P1 = numpy.min(xyzr[:, 0:3], axis=0) - probe_out - step
    xmax, ymax, zmax = numpy.max(xyzr[:, 0:3], axis=0) + probe_out + step
    P2 = numpy.array([xmax, P1[1], P1[2]])
    P3 = numpy.array([P1[0], ymax, P1[2]])
    P4 = numpy.array([P1[0], P1[1], zmax])

    # Pack vertices
    vertices = numpy.array([P1, P2, P3, P4])

    return vertices


def get_vertices_from_file(
    fn: str | pathlib.Path,
    atomic: numpy.ndarray | list[list[str | float | int]],
    step: float | int = 0.6,
    probe_in: float | int = 1.4,
    probe_out: float | int = 4.0,
    nthreads: int | None = None,
) -> tuple[numpy.ndarray, numpy.ndarray | list[list[str | float | int]]]:
    """Gets 3D grid vertices from box configuration file or parKVFinder
    parameters file, selects atoms inside custom 3D grid, define sine
    and cosine of 3D grid angles and define xyz grid units.

    Parameters
    ----------
    fn : str | pathlib.Path
        A path to box configuration file (TOML-formatted).
    atomic : numpy.ndarray | list[list[str | float | int]]
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    step : float | int, optional
        Grid spacing (A), by default 0.6.
    probe_in : float | int, optional
        Probe In size (A), by default 1.4.
    probe_out : float | int, optional
        Probe Out size (A), by default 4.0.
    nthreads : int | None, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.

    Returns
    -------
    vertices : numpy.ndarray
        A numpy.ndarray with xyz vertices coordinates (origin, X-axis, Y-axis,
        Z-axis) of the custom box.
    atomic : numpy.ndarray | list[list[str | float | int]]
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom inside the custom box.

    Raises
    ------
    TypeError
        `fn` must be a string or a pathlib.Path.
    TypeError
        `atomic` must be a list or a numpy.ndarray.
    ValueError
        `atomic` has incorrect shape. It must be (n, 8).
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
        `nthreads` must be a positive integer.
    ValueError
        You must define (p1, p2, p3, p4) or (residues, padding) keys in `fn`.
    ValueError
        Box not properly defined in `fn`.

    Note
    ----
    The box configuration scale file defines the vertices of the 3D grid used by pyKVFinder to detect and characterize cavities. There are three methods for defining a custom 3D grid in pyKVFinder. The first directly defines four vertices of the 3D grid (origin, X-axis, Y-axis and Z-axis), the second defines a list of residues and a padding, and the the third uses parKVFinder parameters file created by its PyMOL plugin. For more details, see `Box configuration file template`.

    See Also
    --------
    read_pdb
    read_xyz
    detect
    constitutional
    hydropathy
    export

    Example
    -------
    First, define a box configuration file (see ``Box configuration file template``).

    >>> import os
    >>> fn = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', 'custom-box.toml')
    >>> with open(fn, 'r') as f:
    ...     print(f.read())
    [box]
    p1 = [3.11, 7.34, 1.59]
    p2 = [11.51, 7.34, 1.59]
    p3 = [3.11, 10.74, 1.59]
    p4 = [3.11, 7.34, 6.19]

    With the atomic information and coordinates read with ``pyKVFinder.read_pdb`` and a box configuration file, we can get the coordinates of grid vertices and select atoms inside custom 3D grid.

    >>> from pyKVFinder import get_vertices_from_file
    >>> vertices, atomic = pyKVFinder.get_vertices_from_file(fn, atomic)

    .. warning::
        Custom box coordinates adds Probe Out size in each direction to create the coordinates of grid vertices.
    """
    import tomlkit
    from pyKVFinder._pyKVFinder import _filter_pdb

    # Check arguments types
    if type(fn) not in [str, pathlib.Path]:
        raise TypeError("`fn` must be a string or a pathlib.Path.")
    if type(atomic) not in [numpy.ndarray, list]:
        raise TypeError("`atomic` must be a list or a numpy.ndarray.")
    elif len(numpy.asarray(atomic).shape) != 2:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    elif numpy.asarray(atomic).shape[1] != 8:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
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
    if nthreads is None:
        nthreads = os.cpu_count() - 1
    else:
        if type(nthreads) not in [int]:
            raise TypeError("`nthreads` must be a positive integer.")
        elif nthreads <= 0:
            raise ValueError("`nthreads` must be a positive integer.")

    # Convert type
    if isinstance(atomic, list):
        atomic = numpy.asarray(atomic)

    # Extract atominfo from atomic
    atominfo = numpy.asarray(
        ([[f"{atom[0]}_{atom[1]}_{atom[2]}", atom[3]] for atom in atomic[:, :4]])
    )

    # Extract xyzr from atomic
    xyzr = atomic[:, 4:].astype(numpy.float64)

    # Read box file
    with open(fn, "r") as file:
        tmp = tomlkit.load(file)
    if "SETTINGS" in tmp.keys():
        if "visiblebox" in tmp["SETTINGS"].keys():
            box = tmp["SETTINGS"]["visiblebox"]
            box = {key: list(box[key].values()) for key in box.keys()}
            vertices = _get_vertices_from_box(box, probe_out)
    else:
        box = tmp["box"] if "box" in tmp.keys() else tmp

        # Check conditions
        if all([key in box.keys() for key in ["p1", "p2", "p3", "p4"]]):
            if any([key in box.keys() for key in ["padding", "residues"]]):
                raise ValueError(
                    f"You must define (p1, p2, p3, p4) or (residues, padding) keys in {fn}."
                )
            vertices = _get_vertices_from_box(box, probe_out)
        elif "residues" in box.keys():
            if any([key in box.keys() for key in ["p1", "p2", "p3", "p4"]]):
                raise ValueError(
                    f"You must define (p1, p2, p3, p4) or (residues, padding) keys in {fn}."
                )
            if "padding" not in box.keys():
                box["padding"] = 3.5
            vertices = _get_vertices_from_residues(box, atominfo, xyzr, probe_out)
        else:
            raise ValueError(f"Box not properly defined in {fn}.")

    # Get sincos
    sincos = numpy.round(_get_sincos(vertices), 4)

    # Get dimensions
    nx, ny, nz = _get_dimensions(vertices, step)

    # Get atoms inside box only
    _filter_pdb(nx, ny, nz, xyzr, vertices[0], sincos, step, probe_in, nthreads)

    # Get indexes of the atoms inside the box
    indexes = xyzr[:, 3] != 0

    # Slice atominfo and xyzr
    atomic = atomic[indexes, :]

    return vertices, atomic


def _get_vertices_from_box(
    box: dict[str, list[float]], probe_out: float = 4.0
) -> numpy.ndarray:
    """Gets 3D grid vertices from box coordinates.

    Parameters
    ----------
    box : dict[str, list[float]]
        A dictionary with xyz coordinates of p1, p2, p3 and p4.
    probe_out : float, optional
        Probe Out size (A), by default 4.0.

    Returns
    -------
    vertices : numpy.ndarray
        A numpy.ndarray of vertices coordinates (origin, Xmax, Ymax, Zmax).

    Example
    -------
    >>> [box]
    >>> p1 = [x1, y1, z1]
    >>> p2 = [x2, y2, z2]
    >>> p3 = [x3, y3, z3]
    >>> p4 = [x4, y4, z4]
    """
    # Get vertices
    P1 = numpy.asarray(box["p1"])
    P2 = numpy.asarray(box["p2"])
    P3 = numpy.asarray(box["p3"])
    P4 = numpy.asarray(box["p4"])

    # Get sincos
    sincos = numpy.round(_get_sincos(numpy.asarray([P1, P2, P3, P4])), 4)

    # Get probe out additions
    # p1 = (x1, y1, z1)
    x1 = (
        -(probe_out * sincos[3])
        - (probe_out * sincos[0] * sincos[2])
        + (probe_out * sincos[1] * sincos[2])
    )
    y1 = -(probe_out * sincos[1]) - (probe_out * sincos[0])
    z1 = (
        -(probe_out * sincos[2])
        + (probe_out * sincos[0] * sincos[3])
        - (probe_out * sincos[1] * sincos[3])
    )
    # p2 = (x2, y2, z2)
    x2 = (
        (probe_out * sincos[3])
        - (probe_out * sincos[0] * sincos[2])
        + (probe_out * sincos[1] * sincos[2])
    )
    y2 = -(probe_out * sincos[1]) - (probe_out * sincos[0])
    z2 = (
        (probe_out * sincos[2])
        + (probe_out * sincos[0] * sincos[3])
        - (probe_out * sincos[1] * sincos[3])
    )
    # p3 = (x3, y3, z3)
    x3 = (
        -(probe_out * sincos[3])
        + (probe_out * sincos[0] * sincos[2])
        + (probe_out * sincos[1] * sincos[2])
    )
    y3 = (probe_out * sincos[1]) - (probe_out * sincos[0])
    z3 = (
        -(probe_out * sincos[2])
        - (probe_out * sincos[0] * sincos[3])
        - (probe_out * sincos[1] * sincos[3])
    )
    # p4 = (x4, y4, z4)
    x4 = (
        -(probe_out * sincos[3])
        - (probe_out * sincos[0] * sincos[2])
        - (probe_out * sincos[1] * sincos[2])
    )
    y4 = -(probe_out * sincos[1]) + (probe_out * sincos[0])
    z4 = (
        -(probe_out * sincos[2])
        + (probe_out * sincos[0] * sincos[3])
        + (probe_out * sincos[1] * sincos[3])
    )

    # Adjust vertices
    P1 += numpy.asarray([x1, y1, z1])
    P2 += numpy.asarray([x2, y2, z2])
    P3 += numpy.asarray([x3, y3, z3])
    P4 += numpy.asarray([x4, y4, z4])

    vertices = numpy.round([P1, P2, P3, P4], 3)

    return vertices


def _get_vertices_from_residues(
    box: dict[str, list[float]],
    atominfo: numpy.ndarray,
    xyzr: numpy.ndarray,
    probe_out: float = 4.0,
) -> numpy.ndarray:
    """Gets 3D grid vertices based on a list of residues (name and chain)
    and a padding value.

    Parameters
    ----------
    box : dict[str, list[float]]
        A dictionary with a list of residues (name and chain) and a
        padding value.
    atominfo : numpy.ndarray
        A numpy.ndarray with residue number, chain, residue name and
        atom name.
    xyzr : numpy.ndarray
        A numpy.ndarray with xyz coordinates and radius of input atoms.
    probe_out : float, optional
        Probe Out size (A), by default 4.0.

    Returns
    -------
    vertices : numpy.ndarray
        A numpy.ndarray of vertices coordinates (origin, Xmax, Ymax, Zmax).

    Box File
    --------
    >>> [box]
    >>> residues = [ ["resnum", "chain", "resname",], ["resnum", "chain",
    ... "resname"], ]
    >>> padding =  3.5
    """
    # Prepare residues list
    residues_list = numpy.array(["_".join(map(str, item)) for item in box["residues"]])

    # Get coordinates of residues
    indexes = numpy.isin(atominfo[:, 0], residues_list)
    xyzr = xyzr[indexes, 0:3]

    # Calculate vertices
    P1 = numpy.min(xyzr[:, 0:3], axis=0) - probe_out - box["padding"]
    xmax, ymax, zmax = numpy.max(xyzr[:, 0:3], axis=0) + probe_out + box["padding"]
    P2 = numpy.array([xmax, P1[1], P1[2]])
    P3 = numpy.array([P1[0], ymax, P1[2]])
    P4 = numpy.array([P1[0], P1[1], zmax])

    vertices = numpy.round([P1, P2, P3, P4], 3)

    return vertices


def _get_dimensions(
    vertices: numpy.ndarray | list[list[float]], step: float | int = 0.6
) -> tuple[int, int, int]:
    """Gets dimensions of 3D grid from vertices.

    Parameters
    ----------
    vertices : numpy.ndarray | list[list[float]]
        A numpy.ndarray or a list with xyz vertices coordinates
        (origin, X-axis, Y-axis, Z-axis).
    step : float | int, optional
        Grid spacing (A), by default 0.6.

    Returns
    -------
    nx : int
        x grid units.
    ny : int
        y grid units.
    nz : int
        z grid units.

    Raises
    ------
    TypeError
        `vertices` must be a list or a numpy.ndarray.
    ValueError
        `vertices` has incorrect shape.
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.
    """
    # Check arguments
    if type(vertices) not in [numpy.ndarray, list]:
        raise TypeError("`vertices` must be a list or a numpy.ndarray.")
    if numpy.asarray(vertices).shape != (4, 3):
        raise ValueError("`vertices` has incorrect shape. It must be (4, 3).")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")

    # Convert type
    if isinstance(vertices, list):
        vertices = numpy.asarray(vertices)
    if isinstance(step, int):
        step = float(step)

    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Calculate distance between points
    norm1 = numpy.linalg.norm(P2 - P1)
    norm2 = numpy.linalg.norm(P3 - P1)
    norm3 = numpy.linalg.norm(P4 - P1)

    # Calculate grid dimensions
    nx = int(norm1 / step) + 1 if norm1 % step != 0 else int(norm1 / step)
    ny = int(norm2 / step) + 1 if norm1 % step != 0 else int(norm1 / step)
    nz = int(norm3 / step) + 1 if norm1 % step != 0 else int(norm1 / step)

    return nx, ny, nz


def _get_sincos(vertices: numpy.ndarray | list[list[float]]) -> numpy.ndarray:
    """Gets sine and cossine of the grid rotation angles from a list of vertices
    coordinates.

    Parameters
    ----------
    vertices : numpy.ndarray | list[list[float]]
        A numpy.ndarray or a list of xyz vertices coordinates (origin, X-axis,
        Y-axis, Z-axis).

    Returns
    -------
    sincos : numpy.ndarray
        A numpy.ndarray with sine and cossine of the grid rotation
        angles (sina, cosa, sinb, cosb).

    Raises
    ------
    TypeError
        `vertices` must be a list or a numpy.ndarray.
    ValueError
        `vertices` has incorrect shape. It must be (4, 3).
    """
    # Check arguments
    if type(vertices) not in [numpy.ndarray, list]:
        raise TypeError("`vertices` must be a list or a numpy.ndarray.")
    if numpy.asarray(vertices).shape != (4, 3):
        raise ValueError("`vertices` has incorrect shape. It must be (4, 3).")

    # Convert type
    if isinstance(vertices, list):
        vertices = numpy.asarray(vertices)

    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Calculate distance between points
    norm1 = numpy.linalg.norm(P2 - P1)
    norm2 = numpy.linalg.norm(P3 - P1)
    norm3 = numpy.linalg.norm(P4 - P1)

    # Calculate sin and cos of angles a and b
    sincos = numpy.array(
        [
            (P4[1] - P1[1]) / norm3,  # sin a
            (P3[1] - P1[1]) / norm2,  # cos a
            (P2[2] - P1[2]) / norm1,  # sin b
            (P2[0] - P1[0]) / norm1,  # cos b
        ]
    )

    return sincos
