import os
import numpy

__all__ = ["get_vertices", "get_vertices_from_file", "get_dimensions", "get_sincos", "detect", "spatial", "depth", "constitutional", "export"]


def get_vertices(xyzr: numpy.ndarray, probe_out: float = 4.0, step: float = 0.6) -> numpy.ndarray:
    """
    Gets 3D grid vertices

    Parameters
    ----------
        xyzr (numpy.ndarray): an array with xyz coordinates and radius of input atoms
        probe_out (float): Probe Out size (A)
        step (float): grid spacing (A)

    Returns
    -------
        vertices (numpy.ndarray): an array of vertices coordinates (origin, Xmax, Ymax, Zmax)
    """
    P1 = numpy.min(xyzr[:, 0:3], axis=0) - probe_out - step
    xmax, ymax, zmax = numpy.max(xyzr[:, 0:3], axis=0) + probe_out + step
    P2 = numpy.array([xmax, P1[1], P1[2]])
    P3 = numpy.array([P1[0], ymax, P1[2]])
    P4 = numpy.array([P1[0], P1[1], zmax])

    vertices = numpy.array([P1, P2, P3, P4])

    return vertices


def get_vertices_from_file(fn: str, resinfo: numpy.ndarray, xyzr: numpy.ndarray, step: float = 0.6, probe_in: float = 1.4, probe_out: float = 4.0, nthreads: int = os.cpu_count() - 1) -> tuple:
    """
    Gets 3D grid vertices from box configuration file, selects atoms inside custom 3D grid, define sine and cosine of 3D grid angles and define xyz grid units

    Parameters
    ----------
        fn (str): path to box configuration file (TOML-formatted)
        resinfo (numpy.ndarray): an array with resnum, chain, resname and atom name
        xyzr (numpy.ndarray): an array with xyz coordinates and radius of input atoms
        step (float): grid spacing (A)
        probe_in (float): Probe In size (A)
        probe_out (float): Probe Out size (A)
        nthreads (int): number of threads for OpenMP

    Returns
    -------
        vertices (numpy.ndarray): an array of vertices coordinates (origin, Xmax, Ymax, Zmax)
        resinfo (numpy.ndarray): an array with resnum, chain, resname and atom name of atoms inside the box
        xyzr (numpy.ndarray): an array with xyz coordinates and radius of atoms inside the box
        sincos (numpy.ndarray): an array with sine and cossine of 3D grid angles (a, b)
        nx (int): x 3D grid units
        ny (int): y 3D grid units
        nz (int): z 3D grid units

    Box File Template
    -----------------
    [box]
    p1 = [x1, y1, z1]
    p2 = [x2, y2, z2]
    p3 = [x3, y3, z3]
    p4 = [x4, y4, z4]

    or

    [box]
    residues = [ ["resname", "chain",], ["resname", "chain",], ]
    padding =  3.5
    """
    from _grid import _filter_pdb
    from toml import load

    # Read box file
    tmp = load(fn)
    box = tmp['box'] if 'box' in tmp.keys() else tmp

    # Check conditions
    if all([key in box.keys() for key in ['p1', 'p2', 'p3', 'p4']]):
        if all([key in box.keys() for key in ['padding', 'residues']]):
            raise Exception(f"You must define (p1, p2, p3, p4) or (residues, padding) keys in {fn}.")
        vertices = _get_vertices_from_box(box, probe_out)
    elif 'residues' in box.keys():
        if all([key in box.keys() for key in ['p1', 'p2', 'p3', 'p4']]):
            raise Exception(f"You must define (p1, p2, p3, p4) or (residues, padding) keys in {fn}.")
        if 'padding' not in box.keys():
            box['padding'] = 3.5
        vertices = _get_vertices_from_residues(box, resinfo, xyzr, probe_out)
    else:
        raise Exception(f"Box not properly defined in {fn}")

    # Get atoms inside box only
    sincos = numpy.round(get_sincos(vertices), 4)
    nx, ny, nz = get_dimensions(vertices, step)
    _filter_pdb(nx, ny, nz, xyzr, vertices[0], sincos, step, probe_in, nthreads)
    indexes = (xyzr[:, 3] != 0)

    # Slice resinfo and xyzr
    resinfo = resinfo[indexes, :]
    xyzr = xyzr[indexes, :]

    return vertices, resinfo, xyzr, sincos, nx, ny, nz


def _get_vertices_from_box(box: dict, probe_out: float = 4.0) -> numpy.ndarray:
    """
    Gets 3D grid vertices from box coordinates

    Parameters
    ----------
        box (dict): an dict with xyzr coordinates of p1, p2, p3 and p4 keys
        probe_out (float): Probe Out size (A)

    Returns
    -------
        vertices (numpy.ndarray): an array of vertices coordinates (origin, Xmax, Ymax, Zmax)

    Box File
    --------
    [box]
    p1 = [x1, y1, z1]
    p2 = [x2, y2, z2]
    p3 = [x3, y3, z3]
    p4 = [x4, y4, z4]
    """
    # Get vertices
    P1 = numpy.array(box['p1'])
    P2 = numpy.array(box['p2'])
    P3 = numpy.array(box['p3'])
    P4 = numpy.array(box['p4'])

    # Get sincos
    sincos = numpy.round(get_sincos(numpy.array([P1, P2, P3, P4])), 4)

    # Get probe out additions
    # p1 = (x1, y1, z1)
    x1 = - (probe_out * sincos[3]) - (probe_out * sincos[0] * sincos[2]) + (probe_out * sincos[1] * sincos[2])
    y1 = - (probe_out * sincos[1]) - (probe_out * sincos[0])
    z1 = - (probe_out * sincos[2]) + (probe_out * sincos[0] * sincos[3]) - (probe_out * sincos[1] * sincos[3])
    # p2 = (x2, y2, z2)
    x2 = (probe_out * sincos[3]) - (probe_out * sincos[0] * sincos[2]) + (probe_out * sincos[1] * sincos[2])
    y2 = - (probe_out * sincos[1]) - (probe_out * sincos[0])
    z2 = (probe_out * sincos[2]) + (probe_out * sincos[0] * sincos[3]) - (probe_out * sincos[1] * sincos[3])
    # p3 = (x3, y3, z3)
    x3 = - (probe_out * sincos[3]) + (probe_out * sincos[0] * sincos[2]) + (probe_out * sincos[1] * sincos[2])
    y3 = (probe_out * sincos[1]) - (probe_out * sincos[0])
    z3 = - (probe_out * sincos[2]) - (probe_out * sincos[0] * sincos[3]) - (probe_out * sincos[1] * sincos[3])
    # p4 = (x4, y4, z4)
    x4 = - (probe_out * sincos[3]) - (probe_out * sincos[0] * sincos[2]) - (probe_out * sincos[1] * sincos[2])
    y4 = - (probe_out * sincos[1]) + (probe_out * sincos[0])
    z4 = - (probe_out * sincos[2]) + (probe_out * sincos[0] * sincos[3]) + (probe_out * sincos[1] * sincos[3])

    # Adjust vertices
    P1 += numpy.array([x1, y1, z1])
    P2 += numpy.array([x2, y2, z2])
    P3 += numpy.array([x3, y3, z3])
    P4 += numpy.array([x4, y4, z4])

    vertices = numpy.array([P1, P2, P3, P4])

    return vertices


def _get_vertices_from_residues(box: dict, resinfo: numpy.ndarray, xyzr: numpy.ndarray, probe_out: float = 4.0) -> numpy.ndarray:
    """
    Gets 3D grid vertices based on a list of residues (name and chain) and a padding value

    Parameters
    ----------
        box (dict): dictionary with a list of residues (name and chain) and a padding value
        resinfo (numpy.ndarray): an array with residue number, chain, residue name and atom name
        xyzr (numpy.ndarray): an array with xyz coordinates and radius of input atoms
        probe_out (float): Probe Out size (A)

    Returns
    -------
        vertices (numpy.ndarray): an array of vertices coordinates (origin, Xmax, Ymax, Zmax)

    Box File
    --------
    [box]
    residues = [ ["resname", "chain",], ["resname", "chain",], ]
    padding =  3.5
    """
    # Prepare residues list
    box['residues'] = numpy.array(['_'.join(item[0:2]) for item in box['residues']])

    # Get coordinates of residues
    indexes = numpy.in1d(resinfo[:, 0], box['residues'])
    xyzr = xyzr[indexes, 0:3]

    # Calculate vertices
    P1 = numpy.min(xyzr[:, 0:3], axis=0) - probe_out - box['padding']
    xmax, ymax, zmax = numpy.max(xyzr[:, 0:3], axis=0) + probe_out + box['padding']
    P2 = numpy.array([xmax, P1[1], P1[2]])
    P3 = numpy.array([P1[0], ymax, P1[2]])
    P4 = numpy.array([P1[0], P1[1], zmax])

    vertices = numpy.array([P1, P2, P3, P4])

    return vertices


def get_dimensions(vertices: numpy.ndarray, step: float = 0.6) -> tuple:
    """
    Gets dimensions of 3D grid from vertices

    Parameters
    ----------
        vertices (numpy.ndarray): an array of vertices coordinates (origin, Xmax, Ymax, Zmax)
        step (float): grid spacing (A)

    Returns
    -------
        nx (int): x 3D grid units
        ny (int): y 3D grid units
        nz (int): z 3D grid units
    """
    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Calculate distance between points
    norm1 = numpy.linalg.norm(P2 - P1)
    norm2 = numpy.linalg.norm(P3 - P1)
    norm3 = numpy.linalg.norm(P4 - P1)

    # Calculate grid dimensions
    nx = int(norm1 / step) + 1
    ny = int(norm2 / step) + 1
    nz = int(norm3 / step) + 1
    return nx, ny, nz


def get_sincos(vertices: numpy.ndarray):
    """
    Gets sine and cossine of 3D grid angles (a, b)

    Parameters
    ----------
        vertices (numpy.ndarray): an array of vertices coordinates (origin, Xmax, Ymax, Zmax)

    Returns
    -------
        sincos (numpy.ndarray): an array with sine and cossine of 3D grid angles (a, b)
    """
    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Calculate distance between points
    norm1 = numpy.linalg.norm(P2 - P1)
    norm2 = numpy.linalg.norm(P3 - P1)
    norm3 = numpy.linalg.norm(P4 - P1)

    # Calculate sin and cos of angles a and b
    sincos = numpy.array([
        (P4[1] - P1[1]) / norm3,  # sin a
        (P3[1] - P1[1]) / norm2,  # cos a
        (P2[2] - P1[2]) / norm1,  # sin b
        (P2[0] - P1[0]) / norm1   # cos b
    ])
    return sincos


def _process_spatial(raw_volume: numpy.ndarray, raw_area: numpy.ndarray, ncav: int) -> tuple:
    """
    Processes arrays of volumes and areas

    Parameters
    ----------
        raw_volume (numpy.ndarray): an array of volumes
        raw_area (numpy.ndarray): an array of areas
        ncav (int): number of cavities

    Returns
    -------
        volume (dict): dictionary with cavity name/volume pairs
        area (dict): dictionary with cavity name/area pairs
    """
    volume, area = {}, {}
    for index in range(ncav):
        key = f"K{chr(65 + int(index / 26) % 26)}{chr(65 + (index % 26))}"
        volume[key] = float(round(raw_volume[index], 2))
        area[key] = float(round(raw_area[index], 2))
    return volume, area


def detect(nx: int, ny: int, nz: int, xyzr: numpy.ndarray, vertices: numpy.ndarray, sincos: numpy.ndarray, step: float = 0.6, probe_in: float = 1.4, probe_out: float = 4.0, removal_distance: float = 2.4, volume_cutoff: float = 5.0, lxyzr: numpy.ndarray = None, ligand_cutoff: float = 5.0, box_adjustment: bool = False, surface: str = 'SES', nthreads: int = os.cpu_count() - 1, verbose: bool = False) -> tuple:
    """
    Detects biomolecular cavities

    Parameters
    ----------
        nx (int): x 3D grid units
        ny (int): y 3D grid units
        nz (int): z 3D grid units
        xyzr (numpy.ndarray): an array with xyz coordinates and radius of input pdb atoms
        vertices (numpy.ndarray): an array of vertices coordinates (origin, Xmax, Ymax, Zmax)
        sincos (numpy.ndarray): an array with sine and cossine of 3D grid angles (a, b)
        step (float): grid spacing (A)
        probe_in (float): Probe In size (A)
        probe_out (float): Probe Out size (A)
        removal_distance (float): length to be removed from the cavity-bulk frontier (A)
        volume_cutoff (float): cavities volume filter (A3)
        lxyzr (numpy.ndarray): an array with xyz coordinates and radius of ligand atoms
        ligand_cutoff (float): radius value to limit a space around a ligand (A)
        box_adjustment (bool): whether a custom 3D grid is applied
        surface (str): SES (Solvent Excluded Surface) or SAS (Solvent Accessible Surface)
        nthreads (int): number of threads
        verbose: print extra information to standard output

    Returns
    -------
        ncav (int): number of cavities
        cavities (numpy.ndarray): cavities 3D grid (cavities[nx][ny][nz])
    """
    from _grid import _detect, _detect_ladj
    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Define ligand adjustment mode
    ligand_adjustment = True if lxyzr is not None else False

    # Calculate number of voxels
    nvoxels = nx * ny * nz

    if surface == 'SES':
        if verbose:
            print("> Surface representation: Solvent Excluded Surface (SES)")
        surface = True
    elif surface == 'SAS':
        if verbose:
            print("> Surface representation: Solvent Accessible Surface (SAS)")
        surface = False
    else:
        raise Exception(f'surface variable should be SAS or SES, not {surface}')

    # Detect cavities
    if ligand_adjustment:
        ncav, cavities = _detect_ladj(nvoxels, nx, ny, nz, xyzr, lxyzr, P1, sincos, step, probe_in, probe_out, removal_distance, volume_cutoff, ligand_adjustment, ligand_cutoff, box_adjustment, P2, surface, nthreads, verbose)
    else:
        ncav, cavities = _detect(nvoxels, nx, ny, nz, xyzr, P1, sincos, step, probe_in, probe_out, removal_distance, volume_cutoff, box_adjustment, P2, surface, nthreads, verbose)

    return ncav, cavities.reshape(nx, ny, nz)


def spatial(cavities: numpy.ndarray, ncav: int, step: float = 0.6, nthreads: int = os.cpu_count() - 1, verbose: bool = False) -> tuple:
    """
    Spatial characterization (volume and area) of the detected cavities

    Parameters
    ----------
        cavities (numpy.ndarray): cavities 3D grid (cavities[nx][ny][nz])
        ncav (int): number of cavities
        step (float): grid spacing (A)
        nthreads (int): number of threads
        verbose: print extra information to standard output

    Returns
    -------
        surface (numpy.ndarray): surface points 3D grid (surface[nx][ny][nz])
        volume (dict): dictionary with cavity name/volume pairs
        area (dict): dictionary with cavity name/area pairs
    """
    from _grid import _spatial
    nx, ny, nz = cavities.shape
    # Get surface points, volume and area
    surface, volume, area = _spatial(cavities, nx * ny * nz, ncav, ncav, step, nthreads, verbose)
    volume, area = _process_spatial(volume, area, ncav)

    return surface.reshape(nx, ny, nz), volume, area


def _process_depth(raw_max_depth: numpy.ndarray, raw_avg_depth: numpy.ndarray, ncav: int) -> tuple:
    """
    Processes arrays of maximum and average depths

    Parameters
    ----------
        raw_max_depth (numpy.ndarray): an array of volumes
        raw_avg_depth (numpy.ndarray): an array of areas
        ncav (int): number of cavities

    Returns
    -------
        max_depth (dict): dictionary with cavity name/maximum depth pairs
        avg_depth (dict): dictionary with cavity name/average depth pairss
    """
    max_depth, avg_depth = {}, {}
    for index in range(ncav):
        key = f"K{chr(65 + int(index / 26) % 26)}{chr(65 + (index % 26))}"
        max_depth[key] = float(round(raw_max_depth[index], 2))
        avg_depth[key] = float(round(raw_avg_depth[index], 2))
    return max_depth, avg_depth


def depth(cavities: numpy.ndarray, ncav: int, step: float = 0.6, nthreads: int = os.cpu_count() - 1, verbose: bool = False) -> tuple:
    """
    Characterization of the depth of the detected cavities, including depth per cavity point and maximum and average depths of detected cavities.

    Parameters
    ----------
        cavities (numpy.ndarray): cavities 3D grid (cavities[nx][ny][nz])
        ncav (int): number of cavities
        step (float): grid spacing (A)
        nthreads (int): number of threads
        verbose: print extra information to standard output

    Returns
    -------
        depths (numpy.ndarray): depth of cavities 3D grid points (depth[nx][ny][nz])
        max_depth (dict): dictionary with cavity name/maximum depth pairs
        avg_depth (dict): dictionary with cavity name/average depth pairs
    """
    from _grid import _depth
    nx, ny, nz = cavities.shape
    # Get depth of cavity points, maximum depth and average depth
    depths, max_depth, avg_depth = _depth(cavities, nx * ny * nz, ncav, ncav, step, nthreads, verbose)
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
            residues[key] = [item.split('_') for item in list(dict.fromkeys(cavity_residues))]
            index += 1
    return residues


def constitutional(cavities: numpy.ndarray, resinfo: numpy.ndarray, xyzr: numpy.ndarray, vertices: numpy.ndarray, sincos: numpy.ndarray, ncav: int, step: float = 0.6, probe_in: float = 1.4, ignore_backbone: bool = False, nthreads: int = os.cpu_count() - 1, verbose: bool = False) -> dict:
    """
    Constitutional characterization (interface residues) of the detected cavities

    Parameters
    ----------
        cavities (numpy.ndarray): cavities 3D grid (cavities[nx][ny][nz])
        resinfo (numpy.ndarray): an array with residue number, chain, residue name and atom name
        xyzr (numpy.ndarray): an array with xyz coordinates and radius of input atoms
        vertices (numpy.ndarray): an array of vertices coordinates (origin, Xmax, Ymax, Zmax)
        sincos (numpy.ndarray): an array with sine and cossine of 3D grid angles (a, b)
        ncav (int): number of cavities
        step (float): grid spacing (A)
        probe_in (float): Probe In size (A)
        ignore_backbone (bool): whether to ignore backbone atoms (C, CA, N, O) when defining interface residues
        nthreads (int): number of threads
        verbose: print extra information to standard output

    Returns
    -------
        residues (dict): dictionary with cavity name/list of interface residues pairs
    """
    from _grid import _constitutional
    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Remove backbone from resinfo
    if ignore_backbone:
        mask = numpy.where(resinfo[:, 1] != 'C') and numpy.where(resinfo[:, 1] != 'CA') and numpy.where(resinfo[:, 1] != 'N') and numpy.where(resinfo[:, 1] != 'O')
        resinfo = resinfo[mask[0], ]
        xyzr = xyzr[mask[0]]

    # Prepare resinfo
    resinfo = resinfo[:, 0].tolist()

    # Get interface residues
    residues = _constitutional(cavities, resinfo, xyzr, P1, sincos, step, probe_in, ncav, nthreads, verbose)
    residues = _process_residues(residues)

    return residues


def hydropathy(surface: numpy.ndarray, resinfo: numpy.ndarray, xyzr: numpy.ndarray, vertices: numpy.ndarray, sincos: numpy.ndarray, hydrophobicity_scale: str, ncav: int, step: float = 0.6, probe_in = 1.4, ignore_backbone: bool = False, nthreads: int = os.cpu_count() - 1, verbose: bool = False) -> tuple:
    import toml
    from _grid import _hydropathy

    # Get dimensions
    nx, ny, nz = surface.shape
    nvoxels = nx * ny * nz

    # Load hydrophobicity scales
    data = list(toml.load(hydrophobicity_scale).values())[0]
    resn, scale = list(data.keys()), numpy.asarray(list(data.values()))

    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Remove backbone from resinfo
    if ignore_backbone:
        mask = numpy.where(resinfo[:, 1] != 'C') and numpy.where(resinfo[:, 1] != 'CA') and numpy.where(resinfo[:, 1] != 'N') and numpy.where(resinfo[:, 1] != 'O')
        resinfo = resinfo[mask[0], ]
        xyzr = xyzr[mask[0]]

    # Get residue name from resinfo
    resname = list(map(lambda x: x.split("_")[2], resinfo[:,0]))

    # Get hydrophobicity scales
    scales = _hydropathy(nvoxels, surface, xyzr, P1, sincos, resname, resn, scale, step, probe_in, ncav, nthreads, verbose).reshape(nx, ny, nz)

    print((scales != 0.0).sum())
    
    return scales


def export(fn: str, cavities: numpy.ndarray, surface: numpy.ndarray, vertices: numpy.ndarray, sincos: numpy.ndarray, ncav: int, step: float = 0.6, B: numpy.ndarray = None, nthreads: int = os.cpu_count() - 1, append: bool = False) -> None:
    """
    Exports cavities to PDB file, with variable as B-factor (optional)

    Parameters
    ----------
        fn (str): path to cavity pdb
        cavities (numpy.ndarray): cavities 3D grid (cavities[nx][ny][nz])
        surface (numpy.ndarray): surface points 3D grid (surface[nx][ny][nz])
        vertices (numpy.ndarray): an array of vertices coordinates (origin, Xmax, Ymax, Zmax)
        sincos (numpy.ndarray): an array with sine and cossine of 3D grid angles (a, b)
        ncav (int): number of cavities
        step (float): grid spacing (A)
        probe_in (float): Probe In size (A)
        B (numpy.ndarray): B-factor for cavity points (B[nx][ny][nz])
        nthreads (int): number of threads
        append (bool): append cavities PDB to `fn`

    Returns
    -------
        None
    """
    from _grid import _export, _export_b
    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Export cavities
    if B is None:
        _export(fn, cavities, surface, P1, sincos, step, ncav, nthreads, append)
    else:
        _export_b(fn, cavities, surface, B, P1, sincos, step, ncav, nthreads, append)
