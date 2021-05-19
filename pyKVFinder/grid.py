import os
import numpy
from typing import Union, Tuple, Dict, List

__all__ = [
    "get_vertices",
    "get_grid_from_file",
    "get_dimensions",
    "get_sincos",
    "detect",
    "spatial",
    "depth",
    "constitutional",
    "hydropathy",
    "export",
]


def get_vertices(
    xyzr: Union[numpy.ndarray, List],
    probe_out: Union[float, int] = 4.0,
    step: Union[float, int] = 0.6,
) -> numpy.ndarray:
    """
    Gets 3D grid vertices

    Parameters
    ----------
    xyzrc : Union[numpy.ndarray, List]
        A numpy array with xyz atomic coordinates and radii
        values for each atom
    probe_out : Union[float, int], optional
        Probe Out size (A)
    step : Union[float, int], optional
        Grid spacing (A)

    Returns
    -------
    vertices : numpy.ndarray
        A numpy array with xyz vertices coordinates
        (origin, X-axis, Y-axis, Z-axis)

    Raises
    ------
    TypeError
        `xyzr` must be a list or a numpy.ndarray
    TypeError
        `probe_out` must be a non-negative real number
    TypeError
        `step` must be a non-negative real number
    ValueError
        `xyzr` has incorrect shape. It must be (, 5)
    """
    if type(xyzr) not in [numpy.ndarray, list]:
        raise TypeError("`xyzr` must be a list or a numpy.ndarray.")
    if type(probe_out) not in [int, float]:
        raise TypeError("`probe_out` must be a non-negative real number.")
    if type(step) not in [int, float]:
        raise TypeError("`step` must be a non-negative real number.")
    if len(numpy.asarray(xyzr).shape) != 2:
        raise ValueError("`xyzr` has incorrect shape. It must be (n, 4)")
    elif numpy.asarray(xyzr).shape[1] != 4:
        raise ValueError("`xyzr` has incorrect shape. It must be (n, 4).")

    # Convert xyzr type
    if type(xyzr) == list:
        xyzr = numpy.asarray(xyzr)

    # Prepare vertices
    P1 = numpy.min(xyzr[:, 0:3], axis=0) - probe_out - step
    xmax, ymax, zmax = numpy.max(xyzr[:, 0:3], axis=0) + probe_out + step
    P2 = numpy.array([xmax, P1[1], P1[2]])
    P3 = numpy.array([P1[0], ymax, P1[2]])
    P4 = numpy.array([P1[0], P1[1], zmax])

    # Pack vertices
    vertices = numpy.array([P1, P2, P3, P4])

    return vertices


def get_grid_from_file(
    fn: str,
    atominfo: Union[numpy.ndarray, List],
    xyzr: Union[numpy.ndarray, List],
    step: Union[float, int] = 0.6,
    probe_in: Union[float, int] = 1.4,
    probe_out: Union[float, int] = 4.0,
    nthreads: int = os.cpu_count() - 1,
) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray, int, int, int]:
    f"""Gets 3D grid vertices from box configuration file or parKVFinder
    parameters file, selects atoms inside custom 3D grid, define sine
    and cosine of 3D grid angles and define xyz grid units.

    Parameters
    ----------
    fn : str
        A path to box configuration file (TOML-formatted)
    atominfo : Union[numpy.ndarray, List]
        A numpy array with atomic information (residue number, chain, residue
        name, atom name)
    xyzr : Union[numpy.ndarray, List]
        A numpy array with xyz atomic coordinates and radii values for each
        atom (x, y, z, radius)
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6
    probe_in : Union[float, int], optional
        Probe In size (A), by default 1.4
    probe_out : Union[float, int], optional
        Probe Out size (A), by default 4.0
    nthreads : int, optional
        Number of threads, by default {os.cpu_count()-1}

    Returns
    -------
    vertices : numpy.ndarray
        A numpy array with xyz vertices coordinates (origin, X-axis, Y-axis,
        Z-axis) of the custom box
    atominfo : numpy.ndarray
        A numpy array with atomic information (residue number, chain, residue
        name, atom name) of atoms inside the custom box
    xyzr : numpy.ndarray
        A numpy array with xyz atomic coordinates and radii values (x, y, z,
        radius) of atoms inside the custom boxx
    sincos : numpy.ndarray
        A numpy array with sine and cossine of the custom box rotation angles
        (sina, cosa, sinb, cosb)
    nx : int
        x grid units
    ny : int
        y grid units
    nz : int
        z grid units

    Raises
    ------
    TypeError
        `fn` must be a str
    TypeError
        `atominfo` must be a list or a numpy.ndarray
    TypeError
        `xyzr` must be a list or a numpy.ndarray
    TypeError
        `step` must be a non-negative real number
    TypeError
        `probe_in` must be a non-negative real number
    TypeError
        `probe_out` must be a non-negative real number
    TypeError
        `nthreads` must be a positive integer
    ValueError
        `xyzr` has incorrect shape. It must be (n, 4)
    ValueError
        `atominfo` has incorrect shape. It must be (n, 2)
    ValueError
        You must define (p1, p2, p3, p4) or (residues, padding) keys in `fn`
    ValueError
        Box not properly defined in `fn`

    Box Configuration File Template
    -------------------------------
    [box]
    p1 = [x1, y1, z1]
    p2 = [x2, y2, z2]
    p3 = [x3, y3, z3]
    p4 = [x4, y4, z4]

    or

    [box]
    residues = [ ["resnum", "chain", "resname", ], ["resnum", "chain",
    "resname"], ]
    padding =  3.5

    ParKVFinder Parameters File
    ---------------------------
    [SETTINGS.visiblebox.p1]
    x = x1
    y = y1
    z = z1

    [SETTINGS.visiblebox.p2]
    x = x2
    y = y2
    z = z2

    [SETTINGS.visiblebox.p3]
    x = x3
    y = y3
    z = z3

    [SETTINGS.visiblebox.p4]
    x = x4
    y = y4
    z = z4
    """
    from _grid import _filter_pdb
    from toml import load

    # Check arguments types
    if type(fn) not in [str]:
        raise TypeError("`fn` must be a str.")
    if type(atominfo) not in [numpy.ndarray, List]:
        raise TypeError("`atominfo` must be a list or a numpy.ndarray.")
    if type(xyzr) not in [numpy.ndarray, List]:
        raise TypeError("`xyzr` must be a list or a numpy.ndarray.")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a non-negative real number.")
    if type(probe_in) not in [float, int]:
        raise TypeError("`probe_in` must be a non-negative real number.")
    if type(probe_out) not in [float, int]:
        raise TypeError("`probe_out` must be a non-negative real number.")
    if type(nthreads) not in [int]:
        raise TypeError("`nthreads` must be a positive integer.")
    if len(numpy.asarray(xyzr).shape) != 2:
        raise ValueError("`xyzr` has incorrect shape. It must be (n, 4).")
    elif numpy.asarray(xyzr).shape[1] != 4:
        raise ValueError("`xyzr` has incorrect shape. It must be (n, 4).")
    if len(numpy.asarray(atominfo).shape) != 2:
        raise ValueError("`atominfo` has incorrect shape. It must be (n, 2).")
    elif numpy.asarray(atominfo).shape[1] != 4:
        raise ValueError("`atominfo` has incorrect shape. It must be (n, 2).")

    # Convert type
    if type(atominfo) == list:
        atominfo = numpy.asarray(atominfo)
    if type(xyzr) == list:
        xyzr = numpy.asarray(xyzr)

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

    # Get atoms inside box only
    sincos = numpy.round(get_sincos(vertices), 4)
    nx, ny, nz = get_dimensions(vertices, step)
    _filter_pdb(nx, ny, nz, xyzr, vertices[0], sincos, step, probe_in, nthreads)
    indexes = xyzr[:, 3] != 0

    # Slice atominfo and xyzr
    atominfo = atominfo[indexes, :]
    xyzr = xyzr[indexes, :]

    return vertices, atominfo, xyzr, sincos, nx, ny, nz


def _get_vertices_from_box(
    box: Dict[str, List[float]], probe_out: Union[float, int] = 4.0
) -> numpy.ndarray:
    """Gets 3D grid vertices from box coordinates

    Parameters
    ----------
    box : Dict[str, List[float]]
        A dictionary with xyz coordinates of p1, p2, p3 and p4
    probe_out : Union[float, int], optional
        Probe Out size (A), by default 4.0

    Returns
    -------
    vertices : numpy.ndarray
        A numpy.ndarray of vertices coordinates (origin, Xmax, Ymax, Zmax)

    Box File
    --------
    [box]
    p1 = [x1, y1, z1]
    p2 = [x2, y2, z2]
    p3 = [x3, y3, z3]
    p4 = [x4, y4, z4]
    """
    # Get vertices
    P1 = numpy.array(box["p1"])
    P2 = numpy.array(box["p2"])
    P3 = numpy.array(box["p3"])
    P4 = numpy.array(box["p4"])

    # Get sincos
    sincos = numpy.round(get_sincos(numpy.array([P1, P2, P3, P4])), 4)

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
    P1 += numpy.array([x1, y1, z1])
    P2 += numpy.array([x2, y2, z2])
    P3 += numpy.array([x3, y3, z3])
    P4 += numpy.array([x4, y4, z4])

    vertices = numpy.array([P1, P2, P3, P4])

    return vertices


def _get_vertices_from_residues(
    box: Dict[str, List[float]],
    atominfo: numpy.ndarray,
    xyzr: numpy.ndarray,
    probe_out: Union[float, int] = 4.0,
) -> numpy.ndarray:
    """Gets 3D grid vertices based on a list of residues (name and chain)
    and a padding value

    Parameters
    ----------
    box : Dict[str, List[float]]
        A dictionary with a list of residues (name and chain) and a
        padding value
    atominfo : numpy.ndarray
        A numpy.ndarray with residue number, chain, residue name and
        atom name
    xyzr : numpy.ndarray
        A numpy.ndarray with xyz coordinates and radius of input atoms
    probe_out : Union[float, int], optional
        Probe Out size (A), by default 4.0

    Returns
    -------
    vertices : numpy.ndarray
        A numpy.ndarray of vertices coordinates (origin, Xmax, Ymax, Zmax)

    Box File
    --------
    [box]
    residues = [ ["resnum", "chain", "resname",], ["resnum", "chain",
    "resname"], ]
    padding =  3.5
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

    vertices = numpy.array([P1, P2, P3, P4])

    return vertices


def get_dimensions(
    vertices: Union[numpy.ndarray, List], step: Union[float, int] = 0.6
) -> Tuple[int, int, int]:
    """Gets dimensions of 3D grid from vertices

    Parameters
    ----------
    vertices : Union[numpy.ndarray, List]
        A numpy.ndarray with xyz vertices coordinates
        (origin, X-axis, Y-axis, Z-axis)
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6

    Returns
    -------
    nx : int
        x grid units
    ny : int
        y grid units
    nz : int
        z grid units

    Raises
    ------
    TypeError
        `vertices` must be a list or a numpy.ndarray
    ValueError
        `vertices` has incorrect shape
    TypeError
        `step` must be a non-negative real number
    """
    # Check arguments
    if type(vertices) not in [numpy.ndarray, list]:
        raise TypeError("`vertices` must be a list or a numpy.ndarray.")
    if numpy.asarray(vertices).shape != (4, 3):
        raise ValueError("`vertices` has incorrect shape. It must be (4, 3).")
    if type(step) not in [float, int] and step > 0.0:
        raise TypeError("`step` must be a non-negative real number.")

    # Convert type
    if type(vertices) == list:
        vertices = numpy.asarray(vertices)

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


def get_sincos(vertices: Union[numpy.ndarray, List]) -> numpy.ndarray:
    """Gets sine and cossine of the grid rotation angles from a list of vertices
    coordinates

    Parameters
    ----------
    vertices : Union[numpy.ndarray, List]
        A list of xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)

    Returns
    -------
    numpy.ndarray
        A numpy.ndarray with sine and cossine of the grid rotation
        angles (sina, cosa, sinb, cosb)

    Raises
    ------
    TypeError
        `vertices` must be a list or a numpy.ndarray
    ValueError
        `vertices` has incorrect shape
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
    nx: int,
    ny: int,
    nz: int,
    xyzr: Union[numpy.ndarray, List],
    vertices: Union[numpy.ndarray, List],
    sincos: Union[numpy.ndarray, List],
    step: Union[float, int] = 0.6,
    probe_in: Union[float, int] = 1.4,
    probe_out: Union[float, int] = 4.0,
    removal_distance: Union[float, int] = 2.4,
    volume_cutoff: Union[float, int] = 5.0,
    lxyzr: Union[numpy.ndarray, List, None] = None,
    ligand_cutoff: Union[float, int] = 5.0,
    box_adjustment: bool = False,
    surface: str = "SES",
    nthreads: int = os.cpu_count() - 1,
    verbose: bool = False,
) -> Tuple[int, numpy.ndarray]:
    """Detects biomolecular cavities

    Parameters
    ----------
    nx : int
        x grid units
    ny : int
        y grid units
    nz : int
        z grid units
    xyzr : Union[numpy.ndarray, List]
        A numpy.ndarray with xyz atomic coordinates and radii values
        (x, y, z, radius) for each atom
    vertices : Union[numpy.ndarray, List]
        A numpy.ndarray with xyz vertices coordinates (origin, X-axis,
        Y-axis, Z-axis)
    sincos : Union[numpy.ndarray, List]
        A numpy.ndarray with sine and cossine of the grid rotation
        angles (sina, cosa, sinb, cosb)
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6
    probe_in : Union[float, int], optional
        Probe In size (A), by default 1.4
    probe_out : Union[float, int], optional
        Probe Out size (A), by default 4.0
    removal_distance : Union[float, int], optional
        A length to be removed from the cavity-bulk frontier (A), by
        default 2.4
    volume_cutoff : Union[float, int], optional
        Volume filter for detected cavities (A3), by default 5.0
    lxyzr : Union[numpy.ndarray, List, None]
        A numpy.ndarray with xyz atomic coordinates and radii values (x, y, z,
        radius) for each atom of a target ligand, by default None
    ligand_cutoff : Union[float, int], optional
        A radius to limit a space around a ligand (A), by default 5.0
    box_adjustment : bool, optional
        Whether a custom 3D grid is applied, by default False
    surface : str, optional
        Keywords options are SES (Solvent Excluded Surface) or SAS (Solvent
        Accessible Surface), by default SES
    nthreads : int
        Number of threads
    verbose : bool
        Print extra information to standard output

    Returns
    -------
    ncav : int
        Number of cavities
    cavities : numpy.ndarray
        Cavity points in the 3D grid (cavities[nx][ny][nz])

    Notes
    -----
    Cavities array has integer labels in each position, that are:
        * -1: bulk points
        * 0: biomolecule points
        * 1: empty space points
        * >=2: cavity points

    The empty space points are regions that do not meet the chosen volume cutoff

    Warnings
    --------
    If you are using box adjustment mode, do not forget to set box_adjustment
    flag to True and read the box configuration file with 'get_grid_from_file'
    function

    If you are using ligand adjustment mode, do not forget to read ligand atom
    coordinates with 'read_pdb' function

    Raises
    ------
    TypeError
        `nx` must be a positive integer
    TypeError
        `ny` must be a positive integer
    TypeError
        `nz` must be a positive integer
    TypeError
        `xyzr` must be a list or a numpy.ndarray
    TypeError
        `vertices` must be a list or a numpy.ndarray
    TypeError
        `sincos` must be a list or a numpy.ndarray
    TypeError
        `step` must be a non-negative real number
    TypeError
        `probe_in` must be a non-negative real number
    TypeError
        `probe_out` must be a non-negative real number
    TypeError
        `removal_distance` must be a non-negative real number
    TypeError
        `volume_cutoff` must be a non-negative real number
    TypeError
        `lxyzr` must be a list, a numpy.ndarray or None
    TypeError
        `ligand_cutoff` must be a non-negative real number
    TypeError
        `box_adjustment` must be a boolean
    TypeError
        `surface` must be a str
    TypeError
        `nthreads` must be a positive integer
    ValueError
        `xyzr` has incorrect shape. It must be (n, 4)
    ValueError
        `vertices` has incorrect shape. It must be (4, 3)
    ValueError
        `sincos` has incorrect shape. It must be (4,)
    ValueError
        `lxyzr` has incorrect shape. It must be (n, 4)
    ValueError
        `surface` must be SAS or SES, not `surface`
    """
    from _grid import _detect, _detect_ladj

    # Check types
    if type(nx) not in [int]:
        raise TypeError("`nx` must be a positive integer.")
    if type(ny) not in [int]:
        raise TypeError("`ny` must be a positive integer.")
    if type(nz) not in [int]:
        raise TypeError("`nz` must be a positive integer.")
    if type(xyzr) not in [numpy.ndarray, list]:
        raise TypeError("`xyzr` must be a list or a numpy.ndarray.")
    if type(vertices) not in [numpy.ndarray, list]:
        raise TypeError("`vertices` must be a list or a numpy.ndarray.")
    if type(sincos) not in [numpy.ndarray, list]:
        raise TypeError("`sincos` must be a list or a numpy.ndarray.")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a non-negative real number.")
    if type(probe_in) not in [float, int]:
        raise TypeError("`probe_in` must be a non-negative real number.")
    if type(probe_out) not in [float, int]:
        raise TypeError("`probe_out` must be a non-negative real number.")
    if type(removal_distance) not in [float, int]:
        raise TypeError("`removal_distance` must be a non-negative real number.")
    if type(volume_cutoff) not in [float, int]:
        raise TypeError("`volume_cutoff` must be a non-negative real number.")
    if type(lxyzr) not in [numpy.ndarray, list, None]:
        raise TypeError("`lxyzr` must be a list, a numpy.ndarray or None.")
    if type(ligand_cutoff) not in [float, int]:
        raise TypeError("`ligand_cutoff` must be a non-negative real number.")
    if type(box_adjustment) not in [bool]:
        raise TypeError("`box_adjustment` must be a boolean.")
    if type(surface) not in [str]:
        raise TypeError("`surface` must be a str.")
    if type(nthreads) not in [int]:
        raise TypeError("`nthreads` must be a positive integer.")
    if len(numpy.asarray(xyzr).shape) != 2:
        raise ValueError("`xyzr` has incorrect shape. It must be (n, 4).")
    elif numpy.asarray(xyzr).shape[1] != 4:
        raise ValueError("`xyzr` has incorrect shape. It must be (n, 4).")
    if numpy.asarray(vertices).shape != (4, 3):
        raise ValueError("`vertices` has incorrect shape. It must be (4, 3).")
    if len(numpy.asarray(sincos).shape) != 1:
        raise ValueError("`sincos` has incorrect shape. It must be (4,).")
    elif numpy.asarray(sincos).shape[0] != 4:
        raise ValueError("`sincos` has incorrect shape. It must be (4,).")

    # Convert types
    if type(xyzr) == list:
        xyzr = numpy.asarray(xyzr)
    if type(vertices) == list:
        vertices = numpy.asarray(vertices)
    if type(sincos) == list:
        sincos = numpy.asarray(sincos)
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
    if type(lxyzr) == list:
        lxyzr = numpy.asarray(lxyzr)

    # Convert numpy.ndarray data types
    xyzr = xyzr.astype("float64") if xyzr.dtype != "float64" else xyzr
    vertices = vertices.astype("float64") if vertices.dtype != "float64" else vertices
    sincos = sincos.astype("float64") if sincos.dtype != "float64" else sincos
    if lxyzr is not None:
        if len(lxyzr.shape) != 2:
            raise ValueError("`lxyzr` has incorrect shape. It must be (n, 4).")
        elif lxyzr.shape[1] != 4:
            raise ValueError("`lxyzr` has incorrect shape. It must be (n, 4).")
        lxyzr = lxyzr.astype("float64") if lxyzr.dtype != "float64" else lxyzr

    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Define ligand adjustment mode
    ligand_adjustment = True if lxyzr is not None else False

    # Calculate number of voxels
    nvoxels = nx * ny * nz

    if surface == "SES":
        if verbose:
            print("> Surface representation: Solvent Excluded Surface (SES)")
        surface = True
    elif surface == "SAS":
        if verbose:
            print("> Surface representation: Solvent Accessible Surface (SAS)")
        surface = False
    else:
        raise ValueError(f"`surface` must be SAS or SES, not {surface}")

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


def _process_spatial(
    raw_volume: numpy.ndarray, raw_area: numpy.ndarray, ncav: int
) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Processes arrays of volumes and areas

    Parameters
    ----------
    raw_volume : numpy.ndarray
        A numpy.ndarray of volumes
    raw_area : numpy.ndarray
        A numpy.ndarray of areas
    ncav : int
        Number of cavities

    Returns
    -------
    volume : Dict[str, float]
        A dictionary with cavity name/volume pairs
    area : Dict[str, float]
        A dictionary with cavity name/area pairs
    """
    volume, area = {}, {}

    for index in range(ncav):
        key = f"K{chr(65 + int(index / 26) % 26)}{chr(65 + (index % 26))}"
        volume[key] = float(round(raw_volume[index], 2))
        area[key] = float(round(raw_area[index], 2))

    return volume, area


def spatial(
    cavities: numpy.ndarray,
    ncav: int,
    step: float = 0.6,
    nthreads: int = os.cpu_count() - 1,
    verbose: bool = False,
) -> tuple:
    """
    Spatial characterization (volume and area) of the detected cavities

    Parameters
    ----------
        cavities (numpy.ndarray): cavity points in the 3D grid (cavities[nx][ny][nz])
        ncav (int): number of cavities
        step (float): grid spacing (A)
        nthreads (int): number of threads
        verbose: print extra information to standard output

    Returns
    -------
        surface (numpy.ndarray): surface points in the 3D grid (surface[nx][ny][nz])
        volume (dict): a dictionary with volume of each detected cavity
        area (dict): a dictionary with area of each detected cavity

    Notes
    -----
        The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

        Surface array has integer labels in each position, that are:
            * -1: bulk points
            * 0: biomolecule or empty space points
            * >=2: cavity points
    """
    from _grid import _spatial

    # Check and convert data types
    cavities = cavities.astype("int32") if cavities.dtype != "int32" else cavities

    # Get cavities shape
    nx, ny, nz = cavities.shape

    # Get surface points, volume and area
    surface, volume, area = _spatial(
        cavities, nx * ny * nz, ncav, ncav, step, nthreads, verbose
    )
    volume, area = _process_spatial(volume, area, ncav)

    return surface.reshape(nx, ny, nz), volume, area


def _process_depth(
    raw_max_depth: numpy.ndarray, raw_avg_depth: numpy.ndarray, ncav: int
) -> tuple:
    """
    Processes arrays of maximum and average depths

    Parameters
    ----------
        raw_max_depth (numpy.ndarray): an array of maximum depth
        raw_avg_depth (numpy.ndarray): an array of average depth
        ncav (int): number of cavities

    Returns
    -------
        max_depth (dict): dictionary with cavity name/maximum depth pairs
        avg_depth (dict): dictionary with cavity name/average depth pairs
    """
    max_depth, avg_depth = {}, {}
    for index in range(ncav):
        key = f"K{chr(65 + int(index / 26) % 26)}{chr(65 + (index % 26))}"
        max_depth[key] = float(round(raw_max_depth[index], 2))
        avg_depth[key] = float(round(raw_avg_depth[index], 2))
    return max_depth, avg_depth


def depth(
    cavities: numpy.ndarray,
    ncav: int,
    step: float = 0.6,
    nthreads: int = os.cpu_count() - 1,
    verbose: bool = False,
) -> tuple:
    """
    Characterization of the depth of the detected cavities, including depth per cavity point and maximum and average depths of detected cavities.

    Parameters
    ----------
        cavities (numpy.ndarray): cavity points in the 3D grid (cavities[nx][ny][nz])
        ncav (int): number of cavities
        step (float): grid spacing (A)
        nthreads (int): number of threads
        verbose: print extra information to standard output

    Returns
    -------
        depths (numpy.ndarray): a numpy array with depth of cavity points (depth[nx][ny][nz])
        max_depth (dict): a dictionary with maximum depth of each detected cavity
        avg_depth (dict): a dictionary with average depth of each detected cavity

    Notes
    -----
        The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

        Cavities array has integer labels in each position, that are:
            * -1: bulk points
            * 0: biomolecule points
            * 1: empty space points
            * >=2: cavity points
    """
    from _grid import _depth

    # Check and convert data types
    cavities = cavities.astype("int32") if cavities.dtype != "int32" else cavities

    # Get cavities shape
    nx, ny, nz = cavities.shape

    # Get depth of cavity points, maximum depth and average depth
    depths, max_depth, avg_depth = _depth(
        cavities, nx * ny * nz, ncav, ncav, step, nthreads, verbose
    )
    max_depth, avg_depth = _process_depth(max_depth, avg_depth, ncav)

    return depths.reshape(nx, ny, nz), max_depth, avg_depth


def _process_residues(raw: list) -> dict:
    """
    Processes raw list of residues from _constitutional to a list of residue information per cavity name

    Parameters
    ----------
        raw (list): a list of residues with cavities separated by '-1'

    Returns
    -------
        residues (dict): a dictionary with cavity name/list of interface residues pairs
    """
    from itertools import groupby

    residues = {}
    index = 0
    for flag, cavity_residues in groupby(raw, lambda res: res == "-1"):
        if not flag:
            key = f"K{chr(65 + int(index / 26) % 26)}{chr(65 + (index % 26))}"
            residues[key] = [
                item.split("_") for item in list(dict.fromkeys(cavity_residues))
            ]
            index += 1
    return residues


def constitutional(
    cavities: numpy.ndarray,
    atominfo: numpy.ndarray,
    xyzr: numpy.ndarray,
    vertices: numpy.ndarray,
    sincos: numpy.ndarray,
    ncav: int,
    step: float = 0.6,
    probe_in: float = 1.4,
    ignore_backbone: bool = False,
    nthreads: int = os.cpu_count() - 1,
    verbose: bool = False,
) -> dict:
    """
    Constitutional characterization (interface residues) of the detected cavities

    Parameters
    ----------
        cavities (numpy.ndarray): cavity points in the 3D grid (cavities[nx][ny][nz])
        atominfo (numpy.ndarray): a numpy array with atomic information (residue number, chain, residue name, atom name)
        xyzr (numpy.ndarray): a numpy array with xyz atomic coordinates and radii values (x, y, z, radius)
        vertices (numpy.ndarray): a numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)
        sincos (numpy.ndarray): a numpy array with sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb)
        ncav (int): number of cavities
        step (float): grid spacing (A)
        probe_in (float): Probe In size (A)
        ignore_backbone (bool): whether to ignore backbone atoms (C, CA, N, O) when defining interface residues
        nthreads (int): number of threads
        verbose: print extra information to standard output

    Returns
    -------
        residues (dict): a dictionary of list of interface residues pairs of each detected cavity

    Notes
    -----
        The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

        Cavities array has integer labels in each position, that are:
            * -1: bulk points
            * 0: biomolecule points
            * 1: empty spaces points
            * >=2: cavity points

    Classes
    -------
        Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine
        Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine
        Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine
        Negatively charged (R4): Aspartate, Glutamate
        Positively charged (R5): Arginine, Histidine, Lysine
        Non-standard (RX): Non-standard residues
    """
    from _grid import _constitutional

    # Check and convert data types
    cavities = cavities.astype("int32") if cavities.dtype != "int32" else cavities
    xyzr = xyzr.astype("float64") if xyzr.dtype != "float64" else xyzr
    vertices = vertices.astype("float64") if vertices.dtype != "float64" else vertices
    sincos = sincos.astype("float64") if sincos.dtype != "float64" else sincos

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
    residues = _process_residues(residues)

    return residues


def _process_hydropathy(raw_avg_hydropathy: numpy.ndarray, ncav: int) -> tuple:
    """
    Processes array of average hydropathy

    Parameters
    ----------
        raw_avg_hydropathy (numpy.ndarray): an array of average hydropathy
        ncav (int): number of cavities

    Returns
    -------
        avg_hydropathy (dict): dictionary with cavity name/average hydropathy pairs
    """
    avg_hydropathy = {}
    for index in range(ncav):
        key = f"K{chr(65 + int(index / 26) % 26)}{chr(65 + (index % 26))}"
        avg_hydropathy[key] = float(round(raw_avg_hydropathy[index], 2))
    return avg_hydropathy


def hydropathy(
    surface: numpy.ndarray,
    atominfo: numpy.ndarray,
    xyzr: numpy.ndarray,
    vertices: numpy.ndarray,
    sincos: numpy.ndarray,
    ncav: int,
    step: float = 0.6,
    probe_in: float = 1.4,
    hydrophobicity_scale: str = "EisenbergWeiss",
    ignore_backbone: bool = False,
    nthreads: int = os.cpu_count() - 1,
    verbose: bool = False,
) -> tuple:
    """
    Hydropathy characterization of the detected cavities.

    Map a target hydrophobicity scale per surface point and calculate average hydropathy of detected cavities.

    Parameters
    ----------
        surface (numpy.ndarray): surface points in the 3D grid (surface[nx][ny][nz])
        atominfo (numpy.ndarray): a numpy array with atomic information (residue number, chain, residue name, atom name)
        xyzr (numpy.ndarray): a numpy array with xyz atomic coordinates and radii values (x, y, z, radius)
        vertices (numpy.ndarray): a numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)
        sincos (numpy.ndarray): a numpy array with sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb)
        ncav (int): number of cavities
        step (float): grid spacing (A)
        probe_in (float): Probe In size (A)
        hydrophobicity_scale (str): name of a built-in hydrophobicity scale (EisenbergWeiss, HessaHeijne, KyteDoolittle, MoonFleming, WimleyWhite, ZhaoLondon) or a path to a TOML-formatted file with a custom hydrophobicity scale.
        ignore_backbone (bool): whether to ignore backbone atoms (C, CA, N, O) when defining interface residues
        nthreads (int): number of threads
        verbose: print extra information to standard output

    Returns
    -------
        scales (numpy.ndarray): a numpy array with hydrophobicity scale value mapped at surface points (scales[nx][ny][nz])
        avg_hydropathy (dict): a dictionary with average hydropathy of each detected cavity and the range of the hydrophobicity scale (min, max)

    Notes
    -----
        The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

        The hydrophobicity scale file defines the name of the scale and the hydrophobicity value for each residue and when not defined, it assigns zero to the missing residues (check Hydrophobicity Scale File Template below). The package contains six built-in hydrophobicity scales: EisenbergWeiss, HessaHeijne, KyteDoolittle, MoonFleming, WimleyWhite and ZhaoLondon.

    Hydrophobicity Scale File Template
    ----------------------------------
        [EisenbergWeiss]
        ALA = -0.64
        ARG = 2.6
        ASN = 0.8
        ASP = 0.92
        CYS = -0.3
        GLN = 0.87
        GLU = 0.76
        GLY = -0.49
        HIS = 0.41
        ILE = -1.42
        LEU = -1.09
        LYS = 1.54
        MET = -0.66
        PHE = -1.22
        PRO = -0.12
        SER = 0.18
        THR = 0.05
        TRP = -0.83
        TYR = -0.27
        VAL = -1.11
    """
    import toml
    from _grid import _hydropathy

    # Check and convert data types
    surface = surface.astype("int32") if surface.dtype != "int32" else surface
    xyzr = xyzr.astype("float64") if xyzr.dtype != "float64" else xyzr
    vertices = vertices.astype("float64") if vertices.dtype != "float64" else vertices
    sincos = sincos.astype("float64") if sincos.dtype != "float64" else sincos

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
    avg_hydropathy = _process_hydropathy(avg_hydropathy, ncav)
    avg_hydropathy[f"{name}"] = [float(scale.min()), float(scale.max())]

    return scales.reshape(nx, ny, nz), avg_hydropathy


def export(
    fn: str,
    cavities: numpy.ndarray,
    surface: numpy.ndarray,
    vertices: numpy.ndarray,
    sincos: numpy.ndarray,
    ncav: int,
    step: float = 0.6,
    B: numpy.ndarray = None,
    output_hydropathy: str = "hydropathy.pdb",
    scales: numpy.ndarray = None,
    nthreads: int = os.cpu_count() - 1,
    append: bool = False,
    model: int = 0,
) -> None:
    """
    Exports cavitiy (H) and surface (HA) points to PDB-formatted file with a variable (B; optional) in B-factor column, and hydropathy to PDB-formatted file in B-factor column at surface points (HA).

    Parameters
    ----------
        fn (str): a path to PDB file for writing cavities
        cavities (numpy.ndarray): cavity points in the 3D grid (cavities[nx][ny][nz])
        surface (numpy.ndarray): surface points in the 3D grid (surface[nx][ny][nz])
        vertices (numpy.ndarray): a numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)
        sincos (numpy.ndarray): a numpy array with sine and cossine of 3D grid angles (a, b)
        ncav (int): number of cavities
        step (float): grid spacing (A)
        probe_in (float): Probe In size (A)
        B (numpy.ndarray): values to be mapped on B-factor column in cavity points (B[nx][ny][nz])
        output_hydropathy (str): a path to hydropathy PDB file (surface points mapped with a hydrophobicity scale)
        scales (numpy.ndarray): hydrophobicity scale values to be mapped on B-factor column in surface points (scales[nx][ny][nz])
        nthreads (int): number of threads
        append (bool): append cavities to PDB file
        model (int): model number

    Returns
    -------
        None

    Note
    ----
        The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.
    """
    from _grid import _export, _export_b

    # Check and convert data types

    vertices = vertices.astype("float64") if vertices.dtype != "float64" else vertices
    sincos = sincos.astype("float64") if sincos.dtype != "float64" else sincos
    if B is not None:
        B = B.astype("float64") if B.dtype != "float64" else B
    if scales is not None:
        scales = scales.astype("float64") if scales.dtype != "float64" else scales

    # Create base directories of results
    if fn is not None:
        os.makedirs(os.path.abspath(os.path.dirname(fn)), exist_ok=True)
    if output_hydropathy is not None:
        os.makedirs(os.path.abspath(os.path.dirname(output_hydropathy)), exist_ok=True)

    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # If surface is None, create an empty grid
    if surface is None:
        surface = numpy.zeros(cavities.shape, dtype="int32")
    else:
        surface = surface.astype("int32") if surface.dtype != "int32" else surface

    if cavities is None:
        if surface is None:
            raise Exception(f"You must define surface when not defining cavities.")
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
