import os
import sys
import numpy as np

here = os.path.abspath(os.path.dirname(__file__))
print(here)


def read_pdb(fn: str):
    pass


def read_vdw_dat(fn: str):
    """
    Read van der Waals radii from .dat format
    """
    vdw = {}
    
    with open(fn, 'r') as f:
        # Read line with data only (ignore empty lines)
        lines = [line.replace(' ', '') for line in f.read().splitlines() if line.replace('\t\t', '')]
        for line in lines:
            if line:
                if line.startswith('>'):
                    res = line.replace('>', '').replace('\t\t', '')
                    vdw[res] = {}
                else:
                    atom, radius = line.split('\t\t')
                    vdw[res][atom] = float(radius)
    
    return vdw


def read_vdw_json(fn: str):
    """
    Read van der Waals radii from .json format
    """
    import json
    with open(fn, "r") as f:
        return json.load(f)


def read_vdw_toml(fn: str):
    """
    Read van der Waals radii from .toml format
    """
    import toml
    with open(fn, "r") as f:
        return toml.load(f)
