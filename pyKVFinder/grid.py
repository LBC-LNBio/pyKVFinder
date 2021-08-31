import os
import pathlib
import numpy
from typing import Union, Tuple, Optional, Dict, List

__all__ = [
    "get_vertices",
    "get_vertices_from_file",
    "detect",
    "spatial",
    "depth",
    "constitutional",
    "hydropathy",
    "export",
]


def get_vertices(
    atomic: Union[numpy.ndarray, List[List[Union[str, float, int]]]],
    probe_out: Union[float, int] = 4.0,
    step: Union[float, int] = 0.6,
) -> numpy.ndarray:
    """Gets 3D grid vertices.

    Parameters
    ----------
    atomic : Union[numpy.ndarray, List[List[Union[str, float, int]]]]
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    probe_out : Union[float, int], optional
        Probe Out size (A), by default 4.0.
    step : Union[float, int], optional
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
    """
    # Check arguments types
    if type(atomic) not in [numpy.ndarray, list]:
        raise TypeError("`atomic` must be a list or a numpy.ndarray.")
    elif len(numpy.asarray(atomic).shape) != 2:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    elif numpy.asarray(atomic).shape[1] != 8:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    if type(probe_out) not in [int, float]:
        raise TypeError("`probe_out` must be a non-negative real number.")
    elif probe_out < 0.0:
        raise ValueError("`probe_out` must be a non-negative real number.")
    if type(step) not in [int, float]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")

    # Convert atomic type
    if type(atomic) == list:
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
    fn: Union[str, pathlib.Path],
    atomic: Union[numpy.ndarray, List[List[Union[str, float, int]]]],
    step: Union[float, int] = 0.6,
    probe_in: Union[float, int] = 1.4,
    probe_out: Union[float, int] = 4.0,
    nthreads: Optional[int] = None,
) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Gets 3D grid vertices from box configuration file or parKVFinder
    parameters file, selects atoms inside custom 3D grid, define sine
    and cosine of 3D grid angles and define xyz grid units.

    Parameters
    ----------
    fn : Union[str, pathlib.Path]
        A path to box configuration file (TOML-formatted).
    atomic : Union[numpy.ndarray, List[List[Union[str, float, int]]]]
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6.
    probe_in : Union[float, int], optional
        Probe In size (A), by default 1.4.
    probe_out : Union[float, int], optional
        Probe Out size (A), by default 4.0.
    nthreads : int, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.

    Returns
    -------
    vertices : numpy.ndarray
        A numpy.ndarray with xyz vertices coordinates (origin, X-axis, Y-axis,
        Z-axis) of the custom box.
    atomic : Union[numpy.ndarray, List[List[Union[str, float, int]]]]
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
    The box configuration scale file defines the vertices of the 3D grid used
    by pyKVFinder to detect and characterize cavities. There are three methods
    for defining a custom 3D grid in pyKVFinder. The first directly defines
    four vertices of the 3D grid (origin, X-axis, Y-axis and Z-axis), the
    second defines a list of residues and a padding, and the the third uses
    parKVFinder parameters file created by its PyMOL plugin. For more details,
    see `Box configuration file template`.
    """
    from _pyKVFinder import _filter_pdb
    from toml import load

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
    if type(atomic) == list:
        atomic = numpy.asarray(atomic)

    # Extract atominfo from atomic
    atominfo = numpy.asarray(
        ([[f"{atom[0]}_{atom[1]}_{atom[2]}", atom[3]] for atom in atomic[:, :4]])
    )

    # Extract xyzr from atomic
    xyzr = atomic[:, 4:].astype(numpy.float64)

    # Read box file
    tmp = load(fn)
    if "SETTINGS" in tmp.keys():
        if "visiblebox" in tmp["SETTINGS"].keys():
            box = tmp["SETTINGS"]["visiblebox"]
            box = {key: list(box[key].values()) for key in box.keys()}
            vertices = _get_vertices_from_box(box, probe_out)
    else:
        box = tmp["box"] if "box" in tmp.keys() else tmp

        # Check conditions
        if all([key in box.keys() for key in ["p1", "p2", "p3", "p4"]]):
            if all([key in box.keys() for key in ["padding", "residues"]]):
                raise ValueError(
                    f"You must define (p1, p2, p3, p4) or (residues, padding) keys in {fn}."
                )
            vertices = _get_vertices_from_box(box, probe_out)
        elif "residues" in box.keys():
            if all([key in box.keys() for key in ["p1", "p2", "p3", "p4"]]):
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
    box: Dict[str, List[float]], probe_out: float = 4.0
) -> numpy.ndarray:
    """Gets 3D grid vertices from box coordinates.

    Parameters
    ----------
    box : Dict[str, List[float]]
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
    box: Dict[str, List[float]],
    atominfo: numpy.ndarray,
    xyzr: numpy.ndarray,
    probe_out: float = 4.0,
) -> numpy.ndarray:
    """Gets 3D grid vertices based on a list of residues (name and chain)
    and a padding value.

    Parameters
    ----------
    box : Dict[str, List[float]]
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
    box["residues"] = numpy.array(["_".join(item[0:3]) for item in box["residues"]])

    # Get coordinates of residues
    indexes = numpy.in1d(atominfo[:, 0], box["residues"])
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
    vertices: Union[numpy.ndarray, List[List[float]]], step: Union[float, int] = 0.6
) -> Tuple[int, int, int]:
    """Gets dimensions of 3D grid from vertices.

    Parameters
    ----------
    vertices : Union[numpy.ndarray, List[List[float]]]
        A numpy.ndarray or a list with xyz vertices coordinates
        (origin, X-axis, Y-axis, Z-axis).
    step : Union[float, int], optional
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
    if type(vertices) == list:
        vertices = numpy.asarray(vertices)
    if type(step) == int:
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


def _get_sincos(vertices: Union[numpy.ndarray, List[List[float]]]) -> numpy.ndarray:
    """Gets sine and cossine of the grid rotation angles from a list of vertices
    coordinates.

    Parameters
    ----------
    vertices : Union[numpy.ndarray, List[List[float]]]
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
    if type(vertices) == list:
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


def detect(
    atomic: Union[numpy.ndarray, List[List[Union[str, float, int]]]],
    vertices: Union[numpy.ndarray, List[List[float]]],
    step: Union[float, int] = 0.6,
    probe_in: Union[float, int] = 1.4,
    probe_out: Union[float, int] = 4.0,
    removal_distance: Union[float, int] = 2.4,
    volume_cutoff: Union[float, int] = 5.0,
    latomic: Optional[Union[numpy.ndarray, List[List[float]]]] = None,
    ligand_cutoff: Union[float, int] = 5.0,
    box_adjustment: bool = False,
    surface: str = "SES",
    nthreads: Optional[int] = None,
    verbose: bool = False,
) -> Tuple[int, numpy.ndarray]:
    """Detects biomolecular cavities.

    Cavity points that belongs to the same cavity are assigned with an integer
    in the grid.

    Parameters
    ----------
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
    probe_out : Union[float, int], optional
        Probe Out size (A), by default 4.0.
    removal_distance : Union[float, int], optional
        A length to be removed from the cavity-bulk frontier (A), by
        default 2.4.
    volume_cutoff : Union[float, int], optional
        Volume filter for detected cavities (A3), by default 5.0.
    latomic : Union[numpy.ndarray, List[List[Union[str, float, int]]]], optional
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom of a target ligand, by default None.
    ligand_cutoff : Union[float, int], optional
        A radius to limit a space around a ligand (A), by default 5.0.
    box_adjustment : bool, optional
        Whether a custom 3D grid is applied, by default False.
    surface : str, optional
        Surface representation. Keywords options are SES (Solvent Excluded Surface) or SAS (Solvent
        Accessible Surface), by default SES.
    nthreads : int, optional
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

    Warning
    -------
    If you are using box adjustment mode, do not forget to set box_adjustment
    flag to True and read the box configuration file with 'get_vertices_from_file'
    function.

    Warning
    -------
    If you are using ligand adjustment mode, do not forget to read ligand atom
    coordinates with 'read_pdb' or 'read_xyz' functions.
    """
    from _pyKVFinder import _detect, _detect_ladj

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
        if len(latomic.shape) != 2:
            raise ValueError("`latomic` has incorrect shape. It must be (n, 8).")
        elif latomic.shape[1] != 8:
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
    if type(atomic) == list:
        atomic = numpy.asarray(atomic)
    if type(vertices) == list:
        vertices = numpy.asarray(vertices)
    if type(step) == int:
        step = float(step)
    if type(probe_in) == int:
        probe_in = float(probe_in)
    if type(probe_out) == int:
        probe_out = float(probe_out)
    if type(removal_distance) == int:
        removal_distance = float(removal_distance)
    if type(volume_cutoff) == int:
        volume_cutoff = float(volume_cutoff)
    if type(ligand_cutoff) == int:
        ligand_cutoff = float(ligand_cutoff)
    if type(latomic) == list:
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
            print("> Surface representation: Solvent Excluded Surface (SES).")
        surface = True
    elif surface == "SAS":
        if verbose:
            print("> Surface representation: Solvent Accessible Surface (SAS).")
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


def _select_cavities(cavities: numpy.ndarray, selection: List[int]) -> numpy.ndarray:
    """Select cavities in the 3D grid by cavity labels.

    Parameters
    ----------
    cavities : numpy.ndarray
        Cavity points in the 3D grid (cavities[nx][ny][nz]).
        Cavities array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule points;

            * 1: empty space points;

            * >=2: cavity points.
    selection : List[int]
        A list of integer labels of each cavity to be selected.

    Returns
    -------
    selected : numpy.ndarray
        Selected cavity points in the 3D grid.
        Selected cavities array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule points;

            * 1: empty space points or unselected cavity points;

            * >=2: selected cavity points.
    """
    # Copy cavities
    selected = numpy.copy(cavities)

    # When outside selection, change cavities tags to 1
    for cav in range(2, cavities.max() + 1):
        if cav not in selection:
            selected[cavities == cav] = 1

    return selected


def _get_cavity_name(index: int) -> str:
    """Get cavity name, eg KAA, KAB, and so on, based on the index.

    Parameters
    ----------
    index : int
        Index in the dictionary.

    Returns
    -------
    cavity_name : str
        Cavity name
    """
    cavity_name = f"K{chr(65 + int(index / 26) % 26)}{chr(65 + (index % 26))}"
    return cavity_name


def _get_cavity_label(cavity_name: str) -> int:
    """Get cavity label, eg 2, 3, and so on, based on the cavity name.

    Parameters
    ----------
    cavity_name : str
        Cavity name

    Returns
    -------
    cavity_label : int
        Integer label of each cavity.

    Raises
    ------
    ValueError
        Invalid cavity name: `cavity_name`.
    """
    # Check cavity name
    if cavity_name[0] != "K":
        raise ValueError(f"Invalid cavity name: {cavity_name}.")

    # Get cavity label
    cavity_label = (ord(cavity_name[1]) - 65) * 26 + (ord(cavity_name[2]) - 65) + 2

    return cavity_label


def _process_spatial(
    raw_volume: numpy.ndarray,
    raw_area: numpy.ndarray,
    ncav: int,
    selection: Optional[List[int]] = None,
) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Processes arrays of volumes and areas.

    Parameters
    ----------
    raw_volume : numpy.ndarray
        A numpy.ndarray of volumes.
    raw_area : numpy.ndarray
        A numpy.ndarray of areas.
    ncav : int
        Number of cavities.
    selection : List[int], optional
        A list of integer labels of each cavity to be selected, by default None.

    Returns
    -------
    volume : Dict[str, float]
        A dictionary with volume of each detected cavity.
    area : Dict[str, float]
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
    step: Union[float, int] = 0.6,
    selection: Optional[Union[List[int], List[str]]] = None,
    nthreads: Optional[int] = None,
    verbose: bool = False,
) -> Tuple[numpy.ndarray, Dict[str, float], Dict[str, float]]:
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
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6.
    selection : Union[List[int], List[str]], optional
        A list of integer labels or a list of cavity names to be selected, by default None.
    nthreads : int, optional
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
    volume : Dict[str, float]
        A dictionary with volume of each detected cavity.
    area : Dict[str, float]
        A dictionary with area of each detected cavity.

    Note
    ----
    The cavity nomenclature is based on the integer label. The cavity marked
    with 2, the first integer corresponding to a cavity, is KAA, the cavity
    marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

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
    """
    from _pyKVFinder import _spatial

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
    if type(step) == int:
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


def _process_depth(
    raw_max_depth: numpy.ndarray,
    raw_avg_depth: numpy.ndarray,
    ncav: int,
    selection: Optional[List[int]] = None,
) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Processes arrays of maximum and average depths.

    Parameters
    ----------
    raw_max_depth : numpy.ndarray
        A numpy.ndarray of maximum depth.
    raw_avg_depth : numpy.ndarray
        A numpy.ndarray of average depth.
    ncav : int
        Number of cavities.
    selection : List[int], optional
        A list of integer labels of each cavity to be selected, by default None.

    Returns
    -------
    max_depth : Dict[str, float]
        A dictionary with maximum depth of each detected cavity.
    avg_depth : Dict[str, float]
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
    step: Union[float, int] = 0.6,
    selection: Optional[Union[List[int], List[str]]] = None,
    nthreads: Optional[int] = None,
    verbose: bool = False,
) -> Tuple[numpy.ndarray, Dict[str, float], Dict[str, float]]:
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
    step : Union[float, int], optional
        Grid spacing (A).
    selection : Union[List[int], List[str]], optional
        A list of integer labels or a list of cavity names to be selected, by default None.
    nthreads : int, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    depths : numpy.ndarray
        A numpy.ndarray with depth of cavity points (depth[nx][ny][nz]).
    max_depth : Dict[str, float]
        A dictionary with maximum depth of each detected cavity.
    avg_depth : Dict[str, float]
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
    marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

    Warning
    -------
    Cavities array has integer labels in each position, that are:

        * -1: bulk points;

        * 0: biomolecule points;

        * 1: empty space points;

        * >=2: cavity points.
    """
    from _pyKVFinder import _depth

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
    if type(step) == int:
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


def _process_residues(
    raw: List[str], ncav: int, selection: Optional[List[int]] = None
) -> Dict[str, List[List[str]]]:
    """Processes raw list of residues from _constitutional to a list of
    residue information per cavity name.

    Parameters
    ----------
    raw : List[str]
        A list of residues with cavities separated by '-1'.
    ncav : int
        Number of cavities.
    selection : List[int], optional
        A list of integer labels of each cavity to be selected, by default None.

    Returns
    -------
    residues : Dict[str, List[List[str]]]
        A dictionary with a list of interface residues for each detected
        cavity.
    """
    from itertools import groupby

    residues = {}

    if selection is None:
        selection = list(range(ncav + 1))
    else:
        selection = [sele - 2 for sele in selection]

    index = 0
    for flag, cavity_residues in groupby(raw, lambda res: res == "-1"):
        if not flag:
            key = _get_cavity_name(selection[index])
            residues[key] = [
                item.split("_") for item in list(dict.fromkeys(cavity_residues))
            ]
            index += 1

    return residues


def constitutional(
    cavities: numpy.ndarray,
    atomic: Union[numpy.ndarray, List[List[Union[str, float, int]]]],
    vertices: Union[numpy.ndarray, List[List[float]]],
    step: Union[float, int] = 0.6,
    probe_in: Union[float, int] = 1.4,
    ignore_backbone: bool = False,
    selection: Optional[Union[List[int], List[str]]] = None,
    nthreads: Optional[int] = None,
    verbose: bool = False,
) -> Dict[str, List[List[str]]]:
    """Constitutional characterization (interface residues) of the detected
    cavities.

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
    ignore_backbone : bool, optional
        Whether to ignore backbone atoms (C, CA, N, O) when defining interface
        residues, by default False.
    selection : Union[List[int], List[str]], optional
        A list of integer labels or a list of cavity names to be selectedA list of integer labels of each cavity to be selected, by default None.
    nthreads : int, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    residues: Dict[str, List[List[str]]]
        A dictionary with a list of interface residues for each detected
        cavity.

    Raises
    ------
    TypeError
        `cavities` must be a numpy.ndarray.
    ValueError
        `cavities` has the incorrect shape. It must be (nx, ny, nz).
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
    marked with 3 is KAB, the cavity marked with 4 is KAC and so on.


    Note
    ----
    The classes of residues are:

    * Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine.

    * Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine.

    * Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine.

    * Negatively charged (R4): Aspartate, Glutamate.

    * Positively charged (R5): Arginine, Histidine, Lysine.

    * Non-standard (RX): Non-standard residues.
    """
    from _pyKVFinder import _constitutional

    # Check arguments
    if type(cavities) not in [numpy.ndarray]:
        raise TypeError("`cavities` must be a numpy.ndarray.")
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
    if type(atomic) == list:
        atomic = numpy.asarray(atomic)
    if type(vertices) == list:
        vertices = numpy.asarray(vertices)
    if type(step) == int:
        step = float(step)
    if type(probe_in) == int:
        probe_in = float(probe_in)

    # Convert numpy.ndarray data types
    cavities = cavities.astype("int32") if cavities.dtype != "int32" else cavities
    vertices = vertices.astype("float64") if vertices.dtype != "float64" else vertices

    # Extract atominfo from atomic
    atominfo = numpy.asarray(
        ([[f"{atom[0]}_{atom[1]}_{atom[2]}", atom[3]] for atom in atomic[:, :4]])
    )

    # Extract xyzr from atomic
    xyzr = atomic[:, 4:].astype(numpy.float64)

    # Get number of cavities
    ncav = int(cavities.max() - 1)

    # Get sincos: sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb)
    sincos = _get_sincos(vertices)

    # Select cavities
    if selection is not None:
        cavities = _select_cavities(cavities, selection)

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
        atominfo = atominfo[
            mask[0],
        ]
        xyzr = xyzr[
            mask[0],
        ]

    # Prepare atominfo
    atominfo = atominfo[:, 0].tolist()

    # Get interface residues
    residues = _constitutional(
        cavities, atominfo, xyzr, P1, sincos, step, probe_in, ncav, nthreads, verbose
    )
    residues = _process_residues(residues, ncav, selection)

    return residues


def _process_hydropathy(
    raw_avg_hydropathy: numpy.ndarray, ncav: int, selection: Optional[List[int]] = None
) -> Dict[str, float]:
    """Processes array of average hydropathy.

    Parameters
    ----------
    raw_avg_hydropathy : numpy.ndarray
        A numpy.ndarray of average hydropathy.
    ncav : int
        Number of cavities.
    selection : List[int], optional
        A list of integer labels of each cavity to be selected.

    Returns
    -------
    avg_hydropathy : Dict[str, float]
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
    atomic: Union[numpy.ndarray, List[List[Union[str, float, int]]]],
    vertices: Union[numpy.ndarray, List[List[float]]],
    step: Union[float, int] = 0.6,
    probe_in: Union[float, int] = 1.4,
    hydrophobicity_scale: Union[str, pathlib.Path] = "EisenbergWeiss",
    ignore_backbone: bool = False,
    selection: Optional[Union[List[int], List[str]]] = None,
    nthreads: Optional[int] = None,
    verbose: bool = False,
) -> Tuple[numpy.ndarray, Dict[str, float]]:
    """Hydropathy characterization of the detected cavities.

    Map a target hydrophobicity scale per surface point and calculate average
    hydropathy of detected cavities.

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
        Name of a built-in hydrophobicity scale (EisenbergWeiss [1]_,
        HessaHeijne [2]_, KyteDoolitte [3]_, MoonFleming [4]_,
        WimleyWhite [5]_, ZhaoLondon [6]_) or a path to a
        TOML-formatted file with a custom hydrophobicity scale, by default
        `EisenbergWeiss`.
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
    marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

    Note
    ----
    The hydrophobicity scale file defines the name of the scale and the
    hydrophobicity value for each residue and when not defined, it assigns
    zero to the missing residues (see `Hydrophobicity scale file
    template`). The package contains six built-in hydrophobicity scales:
    EisenbergWeiss [1]_, HessaHeijne [2]_, KyteDoolittle [3]_,
    MoonFleming [4]_, WimleyWhite [5]_ and ZhaoLondon [6]_.

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

    .. [5] Wimley WC, White SH. Experimentally determined hydrophobicity scale
       for proteins at membrane interfaces. Nature Structural & Molecular
       Biology. 1996;3.

    .. [6] Zhao G, London E. An amino acid “transmembrane tendency” scale that
       approaches the theoretical limit to accuracy for prediction of
       transmembrane helices: Relationship to biological hydrophobicity.
       Protein Science. 2006;15.
    """
    import toml
    from _pyKVFinder import _hydropathy

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
    if type(atomic) == list:
        atomic = numpy.asarray(atomic)
    if type(vertices) == list:
        vertices = numpy.asarray(vertices)
    if type(step) == int:
        step = float(step)
    if type(probe_in) == int:
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
        "WimleyWhite",
        "ZhaoLondon",
    ]:
        hydrophobicity_scale = os.path.join(
            os.path.abspath(os.path.dirname(__file__)),
            f"data/{hydrophobicity_scale}.toml",
        )
    f = toml.load(hydrophobicity_scale)
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
        atominfo = atominfo[
            mask[0],
        ]
        xyzr = xyzr[
            mask[0],
        ]

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


def export(
    fn: Union[str, pathlib.Path],
    cavities: numpy.ndarray,
    surface: Optional[numpy.ndarray],
    vertices: Union[numpy.ndarray, List[List[float]]],
    step: Union[float, int] = 0.6,
    B: Union[numpy.ndarray, None] = None,
    output_hydropathy: Union[str, pathlib.Path] = "hydropathy.pdb",
    scales: Union[numpy.ndarray, None] = None,
    selection: Optional[Union[List[int], List[str]]] = None,
    nthreads: Optional[int] = None,
    append: bool = False,
    model: int = 0,
) -> None:
    """Exports cavitiy (H) and surface (HA) points to PDB-formatted file with
    a variable (B; optional) in B-factor column, and hydropathy to
    PDB-formatted file in B-factor column at surface points (HA).

    Parameters
    ----------
    fn : Union[str, pathlib.Path]
        A path to PDB file for writing cavities.
    cavities : numpy.ndarray
        Cavity points in the 3D grid (cavities[nx][ny][nz]).
        Cavities array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule points;

            * 1: empty space points;

            * >=2: cavity points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.
    surface : numpy.ndarray, optional
        Surface points in the 3D grid (surface[nx][ny][nz]). If None, surface
        is a numpy.zeros array with same shape of cavities.
        Surface array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule or empty space points;

            * >=2: surface points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.
    vertices : Union[numpy.ndarray, List[List[float]]]
        A numpy.ndarray or a list with xyz vertices coordinates (origin,
        X-axis, Y-axis, Z-axis).
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6.
    B : Union[numpy.ndarray, None], optional
        A numpy.ndarray with values to be mapped on B-factor column in cavity
        points (B[nx][ny][nz]), by default None.
    output_hydropathy : Union[str, pathlib.Path], optional
        A path to hydropathy PDB file (surface points mapped with a
        hydrophobicity scale), by default `hydropathy.pdb`.
    scales : Union[numpy.ndarray, None], optional
        A numpy.ndarray with hydrophobicity scale values to be mapped on
        B-factor column in surface points (scales[nx][ny][nz]), by default
        None.
    selection : Union[List[int], List[str]], optional
        A list of integer labels or a list of cavity names to be selectedA list of integer labels of each cavity to be selected, by default None.
    nthreads : int, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    append : bool, optional
        Whether to append cavities to the PDB file, by default False.
    model : int, optional
        Model number, by default 0.

    Raises
    ------
    TypeError
        `cavities` must be a numpy.ndarray.
    ValueError
        `cavities` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `surface` must be a numpy.ndarray.
    ValueError
        `surface` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `vertices` must be a list or a numpy.ndarray.
    ValueError
        `vertices` has incorrect shape. It must be (4, 3).
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.
    TypeError
        `B` must be a numpy.ndarray.
    ValueError
        `B` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `scales` must be a numpy.ndarray.
    ValueError
        `scales` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `selection` must be a list of strings (cavity names) or integers (cavity labels).
    ValueError
        Invalid `selection`: {selection}.
    TypeError
        `nthreads` must be a positive integer.
    ValueError
        `nthreads` must be a positive integer.
    TypeError
        `append` must be a boolean.
    TypeError
        `model` must be a integer.
    TypeError
        `fn` must be a string or pathlib.Path.
    TypeError
        `output_hydropathy` must be a string.
    Exception
        User must define `surface` when not defining `cavities`.

    Note
    ----
    The cavity nomenclature is based on the integer label. The cavity marked
    with 2, the first integer corresponding to a cavity, is KAA, the cavity
    marked with 3 is KAB, the cavity marked with 4 is KAC and so on.
    """
    from _pyKVFinder import _export, _export_b

    # Check arguments
    if cavities is not None:
        if type(cavities) not in [numpy.ndarray]:
            raise TypeError("`cavities` must be a numpy.ndarray.")
        elif len(cavities.shape) != 3:
            raise ValueError(
                "`cavities` has the incorrect shape. It must be (nx, ny, nz)."
            )
    if surface is not None:
        if type(surface) not in [numpy.ndarray]:
            raise TypeError("`surface` must be a numpy.ndarray.")
        elif len(surface.shape) != 3:
            raise ValueError(
                "`surface` has the incorrect shape. It must be (nx, ny, nz)."
            )
    if type(vertices) not in [numpy.ndarray, list]:
        raise TypeError("`vertices` must be a list or a numpy.ndarray.")
    elif numpy.asarray(vertices).shape != (4, 3):
        raise ValueError("`vertices` has incorrect shape. It must be (4, 3).")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")
    if B is not None:
        if type(B) not in [numpy.ndarray]:
            raise TypeError("`B` must be a numpy.ndarray.")
        elif len(B.shape) != 3:
            raise ValueError("`B` has the incorrect shape. It must be (nx, ny, nz).")
    if scales is not None:
        if type(scales) not in [numpy.ndarray]:
            raise TypeError("`scales` must be a numpy.ndarray.")
        elif len(scales.shape) != 3:
            raise ValueError(
                "`scales` has the incorrect shape. It must be (nx, ny, nz)."
            )
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
    if type(append) not in [bool]:
        raise TypeError("`append` must be a boolean.")
    if type(model) not in [int]:
        raise TypeError("`model` must be a integer.")

    # Convert types
    if type(vertices) == list:
        vertices = numpy.asarray(vertices)
    if type(step) == int:
        step = float(step)

    # Convert numpy.ndarray data types
    vertices = vertices.astype("float64") if vertices.dtype != "float64" else vertices
    if B is not None:
        B = B.astype("float64") if B.dtype != "float64" else B
    if scales is not None:
        scales = scales.astype("float64") if scales.dtype != "float64" else scales

    # Get sincos: sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb)
    sincos = _get_sincos(vertices)

    # Create base directories of results
    if fn is not None:
        if type(fn) not in [str, pathlib.Path]:
            raise TypeError("`fn` must be a string or a pathlib.Path.")
        os.makedirs(os.path.abspath(os.path.dirname(fn)), exist_ok=True)
    if output_hydropathy is not None:
        if type(output_hydropathy) not in [str, pathlib.Path]:
            raise TypeError("`output_hydropathy` must be a string or a pathlib.Path.")
        os.makedirs(os.path.abspath(os.path.dirname(output_hydropathy)), exist_ok=True)

    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # If surface is None, create an empty grid
    if surface is None:
        surface = numpy.zeros(cavities.shape, dtype="int32")
    else:
        surface = surface.astype("int32") if surface.dtype != "int32" else surface

    # Select cavities
    if selection is not None:
        surface = _select_cavities(surface, selection)
        if cavities is not None:
            cavities = _select_cavities(cavities, selection)

    if cavities is None:

        # Get number of cavities
        ncav = int(surface.max() - 1)

        if surface is None:
            raise Exception(f"User must define `surface` when not defining `cavities`.")
        else:
            _export_b(
                output_hydropathy,
                surface,
                surface,
                scales,
                P1,
                sincos,
                step,
                ncav,
                nthreads,
                append,
                model,
            )
    else:
        # Check and convert cavities dtype
        cavities = cavities.astype("int32") if cavities.dtype != "int32" else cavities

        # Get number of cavities
        ncav = int(cavities.max() - 1)

        # Export cavities
        if B is None:
            _export(
                fn, cavities, surface, P1, sincos, step, ncav, nthreads, append, model
            )
        else:
            _export_b(
                fn,
                cavities,
                surface,
                B,
                P1,
                sincos,
                step,
                ncav,
                nthreads,
                append,
                model,
            )

        # Export hydropathy surface points
        if scales is None:
            pass
        else:
            _export_b(
                output_hydropathy,
                surface,
                surface,
                scales,
                P1,
                sincos,
                step,
                ncav,
                nthreads,
                append,
                model,
            )
