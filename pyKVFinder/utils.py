import os
import logging
import argparse

__all__ = ["read_pdb", "read_vdw", "frequencies", "write_results"]

here = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/vdw.dat")


def read_vdw(fn: str = here) -> dict:
    """
    Reads van der Waals radii from .dat file

    Parameters
    ----------
        fn (str): van der Waals radii file

    Returns
    -------
        vdw (dict): A dictionary containing radii values with residue name and atom name as keys
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


def read_pdb(fn: str, vdw: dict) -> tuple:
    """
    Reads PDB file into numpy arrays

    Parameters
    ----------
        fn (str): Path to PDB file
        vdw (dict): Dictionary with radii values (vdw[resname][atom])

    Returns
    -------
        pdb (numpy.ndarray): an array with resnum, chain, resname and atom name
        xyzr (numpy.ndarray): an array with xyz coordinates and radii values
    """
    import numpy as _np
    pdb = []
    xyzr = []
    with open(fn, "r") as f:
        for line in f.readlines():
            if line[:4] == 'ATOM' or line[:6] == 'HETATM':
                atom, coords = _process_pdb_line(line, vdw)
                pdb.append(atom)
                xyzr.append(coords)
    return _np.asarray(pdb), _np.asarray(xyzr)


def _process_pdb_line(line: str, vdw: dict) -> tuple:
    """
    Extracts ATOM and HETATM information of PDB line

    Parameters
    ----------
        line (str): line of PDB file
        vdw (dict): Dictionary with radii values (vdw[resname][atom])

    Returns
    -------
        pdb (list): a list with resnum, chain, resname and atom name
        xyzr (list): a list with xyz coordinates and radius
    """
    atom = line[12:16].strip()
    resname = line[17:20].strip()
    resnum = int(line[22:26])
    chain = line[21]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    atom_symbol = line[76:78].strip()
    if resname in vdw.keys() and atom in vdw[resname].keys():
        radius = vdw[resname][atom]
    else:
        radius = vdw['GEN'][atom_symbol]
        logging.info(f"Warning: Atom {atom} of residue {resname} not found in dictionary")
        logging.info(f"Warning: Using generic atom {atom_symbol} radius: {radius} \u00c5")
    pdb = [f"{resnum}_{chain}_{resname}", atom]
    xyzr = [x, y, z, radius]
    return pdb, xyzr


def _process_box(args: argparse.Namespace):
    """
    Gets xyz coordinates of 3D grid vertices

    Parameters
    ----------
        args (argparse.Namespace): arguments passes by argparser CLI

    Returns
    -------
        box (list): a list of vertices coordinates (origin, Xmax, Ymax, Zmax)
    """
    import numpy as _np
    # Create box parameter
    box = {
        'p1': args.vertices[0],
        'p2': args.vertices[1],
        'p3': args.vertices[2],
        'p4': args.vertices[3],
    }

    # Adjust if box adjustment mode
    if args.box:
        # Get probe out additions
        # p1 = (x1, y1, z1)
        x1 = - (args.probe_out * args.sincos[3]) - (args.probe_out * args.sincos[0] * args.sincos[2]) + (args.probe_out * args.sincos[1] * args.sincos[2])
        y1 = - (args.probe_out * args.sincos[1]) - (args.probe_out * args.sincos[0])
        z1 = - (args.probe_out * args.sincos[2]) + (args.probe_out * args.sincos[0] * args.sincos[3]) - (args.probe_out * args.sincos[1] * args.sincos[3])
        # p2 = (x2, y2, z2)
        x2 = (args.probe_out * args.sincos[3]) - (args.probe_out * args.sincos[0] * args.sincos[2]) + (args.probe_out * args.sincos[1] * args.sincos[2])
        y2 = - (args.probe_out * args.sincos[1]) - (args.probe_out * args.sincos[0])
        z2 = (args.probe_out * args.sincos[2]) + (args.probe_out * args.sincos[0] * args.sincos[3]) - (args.probe_out * args.sincos[1] * args.sincos[3])
        # p3 = (x3, y3, z3)
        x3 = - (args.probe_out * args.sincos[3]) + (args.probe_out * args.sincos[0] * args.sincos[2]) + (args.probe_out * args.sincos[1] * args.sincos[2])
        y3 = (args.probe_out * args.sincos[1]) - (args.probe_out * args.sincos[0])
        z3 = - (args.probe_out * args.sincos[2]) - (args.probe_out * args.sincos[0] * args.sincos[3]) - (args.probe_out * args.sincos[1] * args.sincos[3])
        # p4 = (x4, y4, z4)
        x4 = - (args.probe_out * args.sincos[3]) - (args.probe_out * args.sincos[0] * args.sincos[2]) - (args.probe_out * args.sincos[1] * args.sincos[2])
        y4 = - (args.probe_out * args.sincos[1]) + (args.probe_out * args.sincos[0])
        z4 = - (args.probe_out * args.sincos[2]) + (args.probe_out * args.sincos[0] * args.sincos[3]) + (args.probe_out * args.sincos[1] * args.sincos[3])

        # Remove probe out addition
        box['p1'] -= _np.array([x1, y1, z1])
        box['p2'] -= _np.array([x2, y2, z2])
        box['p3'] -= _np.array([x3, y3, z3])
        box['p4'] -= _np.array([x4, y4, z4])

    # Prepare to dict to toml module
    box['p1'] = _np.around(box['p1'], 2).tolist()
    box['p2'] = _np.around(box['p2'], 2).tolist()
    box['p3'] = _np.around(box['p3'], 2).tolist()
    box['p4'] = _np.around(box['p4'], 2).tolist()

    return box


def _write_parameters(args: argparse.Namespace) -> None:
    """
    Writes parameters used in cavity detection and characterization of pyKVFinder
    to TOML-formatted file

    Parameters
    ----------
        args (argparse.Namespace): arguments passes by argparser CLI

    Returns
    -------
        None
    """
    import toml
    # Parameters filename
    fn = os.path.join(args.output_directory, f"{args.base_name}.parameters.toml")

    # Parameters dict
    parameters = {
        'FILES': {
            'INPUT': args.pdb,
            'LIGAND': args.ligand,
            'BASE_NAME': args.base_name,
            'OUTPUT_DIRECTORY': args.output_directory,
            'DICTIONARY': args.dictionary,
        },
        'SETTINGS': {
            'MODES': {
                'BOX_ADJUSTMENT': args.box,
                'LIGAND_ADJUSTMENT': True if args.ligand else False,
                'DEPTH': args.depth,
                'SURFACE': args.surface,
                'IGNORE_BACKBONE': args.ignore_backbone,
            },
            'STEP': args.step,
            'PROBES': {
                'PROBE_IN': args.probe_in,
                'PROBE_OUT': args.probe_out,
            },
            'CUTOFFS': {
                'VOLUME_CUTOFF': args.volume_cutoff,
                'LIGAND_CUTOFF': args.ligand_cutoff,
                'REMOVAL_DISTANCE': args.removal_distance
            },
            'BOX': _process_box(args),
        }
    }

    # Write to TOML file
    with open(fn, "w") as param:
        toml.dump(parameters, param)


def frequencies(residues: dict) -> dict:
    """
    Calculate frequencies of residues and class of residues (R1, R2, R3, R4 and R5) for detected cavities.

    Parameters
    ----------
        residues (dict): a dictionary with a list of interface residues of each detected cavity

    Returns
    -------
        frequency (dict): a dictionary with frequencies of residues and class of residues of each detected cavity
    
    Note
    ----
    R1: Aliphatic apolar
    R2: Aromatic
    R3: Polar Uncharged
    R4: Negatively charged
    R5: Positively charged
    RX: Not classified
    """
    # Create a dict for frequencies
    frequency = {}

    # Get cavity name and residues list for each detected cavity
    for name, reslist in residues.items():
        # Create a dict for cavity name
        frequency[name] = {
            'RESIDUES': {},
            'CLASS': {},
        }
        # Get unique residues names
        reslist = sorted(list(set([res[2] for res in reslist])))
        # Get residues frequency 
        for res in reslist:
            frequency[name]['RESIDUES'][res] = reslist.count(res)

        # Get class frequency
        frequency[name]['CLASS']['R1'] = frequency[name]['RESIDUES'].get('ALA', 0) + frequency[name]['RESIDUES'].get('GLY', 0) + frequency[name]['RESIDUES'].get('ILE', 0) + frequency[name]['RESIDUES'].get('LEU', 0) + frequency[name]['RESIDUES'].get('PRO', 0) + frequency[name]['RESIDUES'].get('VAL', 0)
        frequency[name]['CLASS']['R2'] = frequency[name]['RESIDUES'].get('PHE', 0) + frequency[name]['RESIDUES'].get('TRP', 0) + frequency[name]['RESIDUES'].get('TYR', 0)
        frequency[name]['CLASS']['R3'] = frequency[name]['RESIDUES'].get('ASN', 0) + frequency[name]['RESIDUES'].get('CYS', 0) + frequency[name]['RESIDUES'].get('GLN', 0) + frequency[name]['RESIDUES'].get('MET', 0) + frequency[name]['RESIDUES'].get('SER', 0) + frequency[name]['RESIDUES'].get('THR', 0)
        frequency[name]['CLASS']['R4'] = frequency[name]['RESIDUES'].get('ASP', 0) + frequency[name]['RESIDUES'].get('GLU', 0)
        frequency[name]['CLASS']['R5'] = frequency[name]['RESIDUES'].get('ARG', 0) + frequency[name]['RESIDUES'].get('HIS', 0) + frequency[name]['RESIDUES'].get('LYS', 0)
        frequency[name]['CLASS']['RX'] = len(reslist) - sum(frequency[name]['CLASS'].values())

    return frequency


def write_results(fn: str, pdb: str, ligand: str, output: str, volume: dict = None, area: dict = None, max_depth: dict = None, avg_depth: dict = None, residues: dict = None, step: float = 0.6) -> None:
    """
    Writes file paths and cavity characterization to TOML-formatted file

    Parameters
    ----------
        fn (str): path to KVFinder results TOML-formatted file
        pdb (str): path to input PDB file
        ligand (str): path to ligand PDB file
        output (str): path to cavity PDB file
        volume (dict): dictionary with cavity name/volume pairs
        area (dict): dictionary with cavity name/area pairs
        max_depth (dict): dictionary with cavity name/maximum depth pairs
        avg_depth (dict): dictionary with cavity name/average depth pairs
        residues (dict): dictionary with cavity name/list of interface residues pairs
        step (float): grid spacing (A)

    Returns
    -------
        None
    """
    import toml
    # Prepare paths
    pdb = os.path.abspath(pdb)
    if ligand:
        ligand = os.path.abspath(ligand)
    output = os.path.abspath(output)

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
        'RESULTS': {
            'VOLUME': volume,
            'AREA': area,
            'MAX_DEPTH': max_depth,
            'AVG_DEPTH': avg_depth,
            'RESIDUES': residues,
            'FREQUENCY': frequencies(residues),
        }
    }

    # Write results to toml file
    with open(fn, "w") as f:
        f.write("# pyKVFinder results\n\n")
        toml.dump(results, f)
