import os
import sys
import numpy as np

here = os.path.abspath(os.path.dirname(__file__))
# print(here)


# def read_pdb(fn: str, vdw: dict):
#     """
#     Create a C function to read pdb and return a numpy array
#     """
    # # C reading function
    # numpy_char_vector = _read_pdb(fn: str)
    # Split strings by ,
    # np.char.split(numpy_char_vector, sep=' ')
    # return pdb, coords


def write_results(fn: str, volume: list, area: list):
    keys = [f"K{chr(65 + int(i / 26) % 26)}{chr( 65 + (i % 26))}" for i in range(len(volume))]
    results = {
        'RESULTS' : {
            'VOLUME': dict(zip(keys, volume)),
            'AREA': dict(zip(keys, area)), 
        }
    }
    import toml
    with open(fn, "w") as f:
        toml.dump(results, f)


def read_pdb(fn: str, vdw: dict) -> tuple:
    pdb = []
    coords = []
    with open(fn, "r") as f:
        for line in f.readlines():
            if line[:4] == 'ATOM' or line[:6] == 'HETATM':
                atom, xyzr = process_pdb_line(line, vdw)
                pdb.append(atom)
                coords.append(xyzr)
    return np.asarray(pdb), np.asarray(coords)


def process_pdb_line(line: str, vdw: dict) -> tuple:
    atom = line[12:16].strip()
    resname = line[17:20].strip()
    resnum = int(line[22:26])
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    atom_symbol = line[76:78].strip()
    if atom in vdw[resname].keys():
        radius = vdw[resname][atom]
    else:
        radius = vdw['GEN'][atom_symbol]
    return [resnum, resname, atom], [x, y, z, radius]


def read_vdw_dat(fn: str) -> dict:
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


def read_vdw_json(fn: str) -> dict:
    """
    Read van der Waals radii from .json format
    """
    import json
    with open(fn, "r") as f:
        return json.load(f)


def read_vdw_toml(fn: str) -> dict:
    """
    Read van der Waals radii from .toml format
    """
    import toml
    with open(fn, "r") as f:
        return toml.load(f)

if __name__ == "__main__":
    from time import time

    vdw = read_vdw_dat("data/vdw.dat")

    start = time()
    pdb = read_pdb("tests/1FMO.pdb", vdw)
    end = time()
    total = end-start

    start = time()
    # reference = np.array([np.min(coords[:, 0:2], axis=0), np.max(coords[:, 0:2], axis=0)])
    end = time()
    ref = end-start

    print(pdb)
    # print(coords)
    # print(reference)
    print(f"Time read_pdb:\t{total:.6f}")
    # print(f"Time reference:\t{ref:.6f}")