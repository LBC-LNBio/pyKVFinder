import os as _os
import numpy as _np
from typing import Optional
from _grid import _detect, _spatial, _constitutional, _export

__all__ = ["detect", "spatial", "constitutional", "export"]


def calculate_vertices(xyzr: _np.ndarray, probe_out: float = 4.0):
    P1 = _np.min(xyzr[:, 0:3], axis=0) - probe_out
    xmax, ymax, zmax = _np.max(xyzr[:, 0:3], axis=0) + probe_out
    P2 = _np.array([xmax, P1[1], P1[2]])
    P3 = _np.array([P1[0], ymax, P1[2]])
    P4 = _np.array([P1[0], P1[1], zmax])
    return P1, P2, P3, P4


def calculate_dimensions(P1: _np.ndarray, P2: _np.ndarray, P3: _np.ndarray, P4: _np.ndarray, step: float = 0.6):
    # Calculate distance between points
    norm1 = _np.linalg.norm(P2 - P1)
    norm2 = _np.linalg.norm(P3 - P1)
    norm3 = _np.linalg.norm(P4 - P1)

    # Calculate grid dimensions
    nx = int(norm1 / step) + 1
    ny = int(norm2 / step) + 1
    nz = int(norm3 / step) + 1
    return nx, ny, nz


def calculate_sincos(P1: _np.ndarray, P2: _np.ndarray, P3: _np.ndarray, P4: _np.ndarray):
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


def detect(nx:int, ny:int, nz:int, xyzr:_np.ndarray, lxyzr: Optional[_np.ndarray], reference:_np.ndarray, sincos: _np.ndarray, step:float=0.6, probe_in:float=1.4, probe_out:float=4.0, removal_distance:float=2.4, volume_cutoff:float=5.0, ligand_cutoff: float = 5.0, surface:bool=True, nthreads:int=_os.cpu_count()-1, verbose:bool=False):
    # Define ligand adjustment mode
    ligand_adjustment, lxyzr = (True, lxyzr) if lxyzr is not None else (False, _np.ndarray((0, 4)))
    # Calculate number of voxels
    nvoxels = nx * ny * nz
    # Detect cavities
    ncav, cavities = _detect(nvoxels, nx, ny, nz, xyzr, lxyzr, reference, sincos, step, probe_in, probe_out, removal_distance, volume_cutoff, ligand_adjustment, ligand_cutoff, surface, nthreads, verbose)
    return ncav, cavities.reshape(nx, ny, nz)


def spatial(cavities: _np.ndarray, nx: int, ny: int, nz: int, ncav: int, step: float = 0.6, nthreads: int = _os.cpu_count() - 1, verbose: bool = False):
    surface, volume, area = _spatial(cavities, nx * ny * nz, ncav, ncav, step, nthreads, verbose)
    volume, area = _process_spatial(volume, area, ncav)
    return surface.reshape(nx, ny, nz), volume, area


def constitutional(cavities: _np.ndarray, pdb: _np.ndarray, xyzr: _np.ndarray, reference: _np.ndarray, sincos: _np.ndarray, ncav: int, step: float = 0.6, probe_in: float = 1.4, ignore_backbone: bool = False, nthreads: int = _os.cpu_count() - 1, verbose: bool = False):
    if ignore_backbone:
        # Remove backbone from pdb
        mask = _np.where(pdb[:, 2] != 'C') and _np.where(pdb[:, 2] != 'CA') and _np.where(pdb[:, 2] != 'N') and _np.where(pdb[:, 2] != 'O')
        pdb = pdb[mask[0], ]
        xyzr = xyzr[mask[0]]
    pdb = pdb[:, 0].tolist()
    residues = _constitutional(cavities, pdb, xyzr, reference, sincos, step, probe_in, ncav, nthreads, verbose)
    residues = _process_residues(residues)
    return residues


def export(fn: str, cavities: _np.ndarray, surface: _np.ndarray, reference: _np.ndarray, sincos: _np.ndarray, ncav: int, step: float = 0.6, nthreads: int = _os.cpu_count() - 1):
    _export(fn, cavities, surface, reference, sincos, step, ncav, nthreads)
