import os as _os
import toml as _toml
import numpy as _np
import logging as _logging

__all__ = ["read_pdb", "read_vdw", "write_results"]

here = _os.path.join(_os.path.abspath(_os.path.dirname(__file__)), "data/vdw.dat")

def read_pdb(fn: str, vdw: dict) -> tuple:
    pdb = []
    coords = []
    with open(fn, "r") as f:
        for line in f.readlines():
            if line[:4] == 'ATOM' or line[:6] == 'HETATM':
                atom, xyzr = process_pdb_line(line, vdw)
                pdb.append(atom)
                coords.append(xyzr)
    return _np.asarray(pdb), _np.asarray(coords)


def process_pdb_line(line: str, vdw: dict) -> tuple:
    atom = line[12:16].strip()
    resname = line[17:20].strip()
    resnum = int(line[22:26])
    chain = line[21]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    atom_symbol = line[76:78].strip()
    if atom in vdw[resname].keys():
        radius = vdw[resname][atom]
    else:
        radius = vdw['GEN'][atom_symbol]
        _logging.info(f"Warning: Atom {atom} of residue {resname} not found in dictionary")
        _logging.info(f"Warning: Using generic atom {atom_symbol} radius: {radius} \u00c5")
    return [f"{resnum}_{chain}", resname, atom], [x, y, z, radius]


def read_vdw(fn: str = here) -> dict:
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


def write_results(fn: str, pdb: str, ligand: str, output: str, volume: dict, area: dict, residues: dict, step: float):
    # Prepare paths
    pdb = _os.path.abspath(pdb)
    if ligand:
        ligand = _os.path.abspath(ligand)
    output = _os.path.abspath(output)

    # Create results dictionary
    results = {
        'FILES': {
            'INPUT': pdb,
            'LIGAND': ligand,
            'OUTPUT': output,
        },
        'PARAMETERS': {
            'STEP': step,
        },
        'RESULTS' : {
            'VOLUME': volume,
            'AREA': area,
            'RESIDUES': residues
        }
    }

    # Write results to toml file
    with open(fn, "w") as f:
        f.write("# pyKVFinder results\n\n")
        _toml.dump(results, f)
