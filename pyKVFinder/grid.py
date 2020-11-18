import os
import numpy as np
from itertools import groupby
from _grid import _detect, _spatial, _constitutional, _export

__all__ = ["detect", "spatial", "constitutional", "export"]

def process_residues(raw: list) -> dict:
    residues = {}
    index = 0
    for flag, cavity_residues in groupby(raw, lambda res: res == "-1"):
        if not flag:
            key = f"K{chr(65 + int(index / 26) % 26)}{chr(65 + (index % 26))}"
            residues[key] = [item.split('_') for item in list(dict.fromkeys(cavity_residues))]
            index += 1
    return residues


def process_spatial(raw_volume: np.ndarray, raw_area: np.ndarray, ncav: int) -> tuple:
    volume, area = {}, {}
    for index in range(ncav):
        key = f"K{chr(65 + int(index / 26) % 26)}{chr(65 + (index % 26))}"
        volume[key] = float(round(raw_volume[index], 2))
        area[key] = float(round(raw_area[index], 2))
    return volume, area

def detect(nx:int, ny:int, nz:int, xyzr:np.ndarray, reference:np.ndarray, sincos: np.ndarray, step:float=0.6, probe_in:float=1.4, probe_out:float=4.0, removal_distance:float=2.4, volume_cutoff:float=5.0, surface:bool=True, nthreads:int=os.cpu_count()-1, verbose:bool=False):
    nvoxels = nx * ny * nz
    ncav, cavities = _detect(nvoxels, nx, ny, nz, xyzr, reference, sincos, step, probe_in, probe_out, removal_distance, volume_cutoff, surface, nthreads, verbose)
    return ncav, cavities.reshape(nx, ny, nz)

def spatial(cavities:np.ndarray, nx:int, ny:int, nz:int, ncav:int, step:float=0.6, nthreads:int=os.cpu_count()-1, verbose:bool=False):
    surface, volume, area = _spatial(cavities, nx * ny * nz, ncav, ncav, step, nthreads, verbose)
    volume, area = process_spatial(volume, area, ncav)
    return surface.reshape(nx, ny, nz), volume, area

def constitutional(cavities:np.ndarray, pdb:np.ndarray, xyzr:np.ndarray, reference:np.ndarray, sincos:np.ndarray, ncav:int, step:float=0.6, probe_in:float=1.4, ignore_backbone:bool=False, nthreads:int=os.cpu_count()-1, verbose:bool=False):
    if ignore_backbone:
        # Remove backbone from pdb
        mask = np.where(pdb[:,2] != 'C') and np.where(pdb[:,2] != 'CA') and np.where(pdb[:,2] != 'N') and np.where(pdb[:,2] != 'O')
        pdb = pdb[mask[0],]
        xyzr = xyzr[mask[0]]
    pdb = pdb[:, 0].tolist()
    residues = _constitutional(cavities, pdb, xyzr, reference, sincos, step, probe_in, ncav, nthreads, verbose)
    residues = process_residues(residues)
    return residues

def export(fn:str, cavities:np.ndarray, surface:np.ndarray, reference:np.ndarray, sincos:np.ndarray, ncav:int, step:float=0.6, nthreads:int=os.cpu_count()-1):
    _export(fn, cavities, surface, reference, sincos, step, ncav, nthreads)
