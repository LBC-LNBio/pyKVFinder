import os as _os
import numpy as _np

__all__ = ["calculate_vertices", "get_vertices", "calculate_dimensions", "calculate_sincos", "detect", "spatial", "constitutional", "export"]


def calculate_vertices(xyzr: _np.ndarray, probe_out: float = 4.0, step: float = 0.6):
    P1 = _np.min(xyzr[:, 0:3], axis=0) - probe_out - step
    xmax, ymax, zmax = _np.max(xyzr[:, 0:3], axis=0) + probe_out + step
    P2 = _np.array([xmax, P1[1], P1[2]])
    P3 = _np.array([P1[0], ymax, P1[2]])
    P4 = _np.array([P1[0], P1[1], zmax])
    return _np.array([P1, P2, P3, P4])


def get_vertices(fn: str, pdb: _np.ndarray, xyzr: _np.ndarray, probe_in: float = 1.4, probe_out: float = 4.0, step: float = 0.6, nthreads: int = _os.cpu_count() - 1):
    from _grid import _filter_pdb
    from toml import load

    # Read box file
    tmp = load(fn)
    box = tmp['box'] if 'box' in tmp.keys() else tmp

    # Check conditions
    if all([key in box.keys() for key in ['p1', 'p2', 'p3', 'p4']]):
        if all([key in box.keys() for key in ['padding', 'residues']]):
            raise Exception(f"You must define (p1, p2, p3, p4) or (residues, padding) keys in {fn}.")
        vertices = _get_vertices_box(box, probe_out)
    elif 'residues' in box.keys():
        if all([key in box.keys() for key in ['p1', 'p2', 'p3', 'p4']]):
            raise Exception(f"You must define (p1, p2, p3, p4) or (residues, padding) keys in {fn}.")
        if 'padding' not in box.keys():
            box['padding'] = 3.5
        vertices = _get_residues_box(box, pdb, xyzr, probe_out)
    else:
        raise Exception(f"Box not properly defined in {fn}")

    # Get atoms inside box only
    sincos = _np.round(calculate_sincos(vertices), 4)
    nx, ny, nz = calculate_dimensions(vertices, step)
    _filter_pdb(nx, ny, nz, xyzr, vertices[0], sincos, step, probe_in, nthreads)
    indexes = (xyzr[:, 3] != 0)

    return vertices, pdb[indexes, :], xyzr[indexes, :], sincos, nx, ny, nz


def _get_vertices_box(box: dict, probe_out: float = 4.0):
    # Get vertices
    P1 = _np.array(box['p1'])
    P2 = _np.array(box['p2'])
    P3 = _np.array(box['p3'])
    P4 = _np.array(box['p4'])

    # Get sincos
    sincos = _np.round(calculate_sincos(_np.array([P1, P2, P3, P4])), 4)

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
    P1 += _np.array([x1, y1, z1])
    P2 += _np.array([x2, y2, z2])
    P3 += _np.array([x3, y3, z3])
    P4 += _np.array([x4, y4, z4])

    return _np.array([P1, P2, P3, P4])


def _get_residues_box(box: dict, pdb: _np.ndarray, xyzr: _np.ndarray, probe_out: float = 4.0):
    # Prepare residues list
    box['residues'] = _np.array(['_'.join(item[0:2]) for item in box['residues']])

    # Get coordinates of residues
    indexes = _np.in1d(pdb[:, 0], box['residues'])
    xyzr = xyzr[indexes, 0:3]

    # Calculate vertices
    P1 = _np.min(xyzr[:, 0:3], axis=0) - probe_out - box['padding']
    xmax, ymax, zmax = _np.max(xyzr[:, 0:3], axis=0) + probe_out + box['padding']
    P2 = _np.array([xmax, P1[1], P1[2]])
    P3 = _np.array([P1[0], ymax, P1[2]])
    P4 = _np.array([P1[0], P1[1], zmax])

    return _np.array([P1, P2, P3, P4])


def calculate_dimensions(vertices: _np.ndarray, step: float = 0.6):
    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Calculate distance between points
    norm1 = _np.linalg.norm(P2 - P1)
    norm2 = _np.linalg.norm(P3 - P1)
    norm3 = _np.linalg.norm(P4 - P1)

    # Calculate grid dimensions
    nx = int(norm1 / step) + 1
    ny = int(norm2 / step) + 1
    nz = int(norm3 / step) + 1
    return nx, ny, nz


def calculate_sincos(vertices: _np.ndarray):
    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Calculate distance between points
    norm1 = _np.linalg.norm(P2 - P1)
    norm2 = _np.linalg.norm(P3 - P1)
    norm3 = _np.linalg.norm(P4 - P1)

    # Calculate sin and cos of angles a and b
    sincos = _np.array([
        (P4[1] - P1[1]) / norm3,  # sin a
        (P3[1] - P1[1]) / norm2,  # cos a
        (P2[2] - P1[2]) / norm1,  # sin b
        (P2[0] - P1[0]) / norm1   # cos b
    ])
    return sincos


def _process_residues(raw: list) -> dict:
    from itertools import groupby
    residues = {}
    index = 0
    for flag, cavity_residues in groupby(raw, lambda res: res == "-1"):
        if not flag:
            key = f"K{chr(65 + int(index / 26) % 26)}{chr(65 + (index % 26))}"
            residues[key] = [item.split('_') for item in list(dict.fromkeys(cavity_residues))]
            index += 1
    return residues


def _process_spatial(raw_volume: _np.ndarray, raw_area: _np.ndarray, ncav: int) -> tuple:
    volume, area = {}, {}
    for index in range(ncav):
        key = f"K{chr(65 + int(index / 26) % 26)}{chr(65 + (index % 26))}"
        volume[key] = float(round(raw_volume[index], 2))
        area[key] = float(round(raw_area[index], 2))
    return volume, area


def detect(nx: int, ny: int, nz: int, xyzr: _np.ndarray, vertices: _np.ndarray, sincos: _np.ndarray, step: float = 0.6, probe_in: float = 1.4, probe_out: float = 4.0, removal_distance: float = 2.4, volume_cutoff: float = 5.0, lxyzr: _np.ndarray = None, ligand_cutoff: float = 5.0, box_adjustment: bool = False, surface: bool = True, nthreads: int = _os.cpu_count() - 1, verbose: bool = False):
    from _grid import _detect, _detect_ladj
    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Define ligand adjustment mode
    ligand_adjustment = True if lxyzr is not None else False

    # Calculate number of voxels
    nvoxels = nx * ny * nz

    # Detect cavities
    if ligand_adjustment:
        ncav, cavities = _detect_ladj(nvoxels, nx, ny, nz, xyzr, lxyzr, P1, sincos, step, probe_in, probe_out, removal_distance, volume_cutoff, ligand_adjustment, ligand_cutoff, box_adjustment, P2, surface, nthreads, verbose)
    else:
        ncav, cavities = _detect(nvoxels, nx, ny, nz, xyzr, P1, sincos, step, probe_in, probe_out, removal_distance, volume_cutoff, box_adjustment, P2, surface, nthreads, verbose)

    return ncav, cavities.reshape(nx, ny, nz)


def spatial(cavities: _np.ndarray, nx: int, ny: int, nz: int, ncav: int, step: float = 0.6, nthreads: int = _os.cpu_count() - 1, verbose: bool = False):
    from _grid import _spatial
    # Get surface points, volume and area
    surface, volume, area = _spatial(cavities, nx * ny * nz, ncav, ncav, step, nthreads, verbose)
    volume, area = _process_spatial(volume, area, ncav)

    return surface.reshape(nx, ny, nz), volume, area


def constitutional(cavities: _np.ndarray, pdb: _np.ndarray, xyzr: _np.ndarray, vertices: _np.ndarray, sincos: _np.ndarray, ncav: int, step: float = 0.6, probe_in: float = 1.4, ignore_backbone: bool = False, nthreads: int = _os.cpu_count() - 1, verbose: bool = False):
    from _grid import _constitutional
    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Remove backbone from pdb
    if ignore_backbone:
        mask = _np.where(pdb[:, 2] != 'C') and _np.where(pdb[:, 2] != 'CA') and _np.where(pdb[:, 2] != 'N') and _np.where(pdb[:, 2] != 'O')
        pdb = pdb[mask[0], ]
        xyzr = xyzr[mask[0]]

    # Prepare pdb
    pdb = pdb[:, 0].tolist()

    # Get interface residues
    residues = _constitutional(cavities, pdb, xyzr, P1, sincos, step, probe_in, ncav, nthreads, verbose)
    residues = _process_residues(residues)

    return residues


def export(fn: str, cavities: _np.ndarray, surface: _np.ndarray, vertices: _np.ndarray, sincos: _np.ndarray, ncav: int, step: float = 0.6, nthreads: int = _os.cpu_count() - 1):
    from _grid import _export
    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Export cavities
    _export(fn, cavities, surface, P1, sincos, step, ncav, nthreads)
