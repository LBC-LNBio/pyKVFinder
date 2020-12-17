import os
import logging
import argparse

__all__ = ["read_pdb", "read_vdw", "calculate_frequencies", "plot_frequencies", "write_results"]

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
        resinfo (numpy.ndarray): an array with resnum, chain, resname and atom name
        xyzr (numpy.ndarray): an array with xyz coordinates and radii values
    """
    import numpy as np
    resinfo = []
    xyzr = []
    with open(fn, "r") as f:
        for line in f.readlines():
            if line[:4] == 'ATOM' or line[:6] == 'HETATM':
                atom, coords = _process_pdb_line(line, vdw)
                resinfo.append(atom)
                xyzr.append(coords)
    return np.asarray(resinfo), np.asarray(xyzr)


def _process_pdb_line(line: str, vdw: dict) -> tuple:
    """
    Extracts ATOM and HETATM information of PDB line

    Parameters
    ----------
        line (str): line of PDB file
        vdw (dict): Dictionary with radii values (vdw[resname][atom])

    Returns
    -------
        resinfo (list): a list with resnum, chain, resname and atom name
        coords (list): a list with xyz coordinates and radius
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
    resinfo = [f"{resnum}_{chain}_{resname}", atom]
    coords = [x, y, z, radius]
    return resinfo, coords


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
    import numpy as np
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
        box['p1'] -= np.array([x1, y1, z1])
        box['p2'] -= np.array([x2, y2, z2])
        box['p3'] -= np.array([x3, y3, z3])
        box['p4'] -= np.array([x4, y4, z4])

    # Prepare to dict to toml module
    box['p1'] = np.around(box['p1'], 2).tolist()
    box['p2'] = np.around(box['p2'], 2).tolist()
    box['p3'] = np.around(box['p3'], 2).tolist()
    box['p4'] = np.around(box['p4'], 2).tolist()

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


def calculate_frequencies(residues: dict) -> dict:
    """
    Calculate frequencies of residues and class of residues (R1, R2, R3, R4 and R5) for detected cavities.

    Parameters
    ----------
        residues (dict): a dictionary with a list of interface residues of each detected cavity

    Returns
    -------
        frequencies (dict): a dictionary with frequencies of residues and class of residues of each detected cavity

    Note
    ----
    R1: Aliphatic apolar
    R2: Aromatic
    R3: Polar Uncharged
    R4: Negatively charged
    R5: Positively charged
    RX: Non-standard
    """
    # Create a dict for frequencies
    frequencies = {}

    # Get cavity name and residues list for each detected cavity
    for name, reslist in residues.items():
        # Create a dict for cavity name
        frequencies[name] = {
            'RESIDUES': {},
            'CLASS': {},
        }
        # Get unique residues names
        reslist = sorted(list(set([res[2] for res in reslist])))
        # Get residues frequencies
        for res in reslist:
            frequencies[name]['RESIDUES'][res] = reslist.count(res)

        # Get class frequencies
        frequencies[name]['CLASS']['R1'] = frequencies[name]['RESIDUES'].get('ALA', 0) + frequencies[name]['RESIDUES'].get('GLY', 0) + frequencies[name]['RESIDUES'].get('ILE', 0) + frequencies[name]['RESIDUES'].get('LEU', 0) + frequencies[name]['RESIDUES'].get('PRO', 0) + frequencies[name]['RESIDUES'].get('VAL', 0)
        frequencies[name]['CLASS']['R2'] = frequencies[name]['RESIDUES'].get('PHE', 0) + frequencies[name]['RESIDUES'].get('TRP', 0) + frequencies[name]['RESIDUES'].get('TYR', 0)
        frequencies[name]['CLASS']['R3'] = frequencies[name]['RESIDUES'].get('ASN', 0) + frequencies[name]['RESIDUES'].get('CYS', 0) + frequencies[name]['RESIDUES'].get('GLN', 0) + frequencies[name]['RESIDUES'].get('MET', 0) + frequencies[name]['RESIDUES'].get('SER', 0) + frequencies[name]['RESIDUES'].get('THR', 0)
        frequencies[name]['CLASS']['R4'] = frequencies[name]['RESIDUES'].get('ASP', 0) + frequencies[name]['RESIDUES'].get('GLU', 0)
        frequencies[name]['CLASS']['R5'] = frequencies[name]['RESIDUES'].get('ARG', 0) + frequencies[name]['RESIDUES'].get('HIS', 0) + frequencies[name]['RESIDUES'].get('LYS', 0)
        frequencies[name]['CLASS']['RX'] = len(reslist) - sum(frequencies[name]['CLASS'].values())

    return frequencies


def plot_frequencies(frequencies: dict, fn: str = 'histograms.pdf') -> None:
    """
    Plot histograms of calculated frequencies (residues and classes of residues) for each detected cavity in a target PDF file

    Parameters
    ----------
        frequencies (dict): A dictionary with frequencies of residues and class of residues of each detected cavity
        fn (str): A path to a PDF file

    Returns
    -------
        None

    Note
    ----
    R1: Aliphatic apolar
    R2: Aromatic
    R3: Polar Uncharged
    R4: Negatively charged
    R5: Positively charged
    RX: Non-standard
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    # Create base directories of output PDF file
    os.makedirs(os.path.abspath(os.path.dirname(pdf)), exist_ok=True)

    # Create a dictionary for standard amino acids
    tmp = {'ALA': 0, 'ARG': 0, 'ASN': 0, 'ASP': 0, 'CYS': 0, 'GLN': 0, 'GLU': 0, 'GLY': 0, 'HIS': 0, 'ILE': 0, 'LEU': 0, 'LYS': 0, 'MET': 0, 'PHE': 0, 'PRO': 0, 'SER': 0, 'THR': 0, 'TRP': 0, 'TYR': 0, 'VAL': 0}

    with PdfPages(fn) as pdf:
        # Standardize data
        ymax = 0
        for cavity_tag in frequencies.keys():
            # Include missing residues
            frequencies[cavity_tag]['RESIDUES'] = {**tmp, **frequencies[cavity_tag]['RESIDUES']}
            # Get y maximum
            if ymax < max(frequencies[cavity_tag]['CLASS'].values()):
                ymax = max(frequencies[cavity_tag]['CLASS'].values())
        ymax += 1

        # Pdf plots
        for cavity_tag in frequencies.keys():
            # Create page
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 9), dpi=300)
            fig.suptitle(r'Cavity ' + f'{cavity_tag}', fontsize=30)

            # Frequency residues
            x = list(frequencies[cavity_tag]['RESIDUES'].keys())
            y = frequencies[cavity_tag]['RESIDUES'].values()

            ax1.bar(x, y, align='center', edgecolor='black')
            ax1.set_xlabel(None)
            ax1.set_xlim(-1, len(x))
            ax1.tick_params(axis='x', labelsize=15, rotation=45)
            ax1.tick_params(axis='y', labelsize=20)
            ax1.set_ylabel(r'Frequency', fontsize=20)
            ax1.set_ylim(0, ymax)
            ax1.grid(which='major', axis='y', linestyle='--')

            # Frequency classes
            x = list(frequencies[cavity_tag]['CLASS'].keys())
            y = frequencies[cavity_tag]['CLASS'].values()
            colors = ['tab:cyan', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:gray']

            ax2.bar(x=x, height=y, align='center', edgecolor='black', color=colors)
            ax2.set_xlabel(None)
            ax2.set_xlim(-1, len(x))
            ax2.tick_params(axis='x', labelsize=20)
            ax2.tick_params(axis='y', labelsize=20)
            ax2.set_ylabel(None)
            ax2.set_ylim(0, ymax)
            ax2.grid(which='major', axis='y', linestyle='--')

            # Legend
            labels = [r'Aliphatic apolar', r'Aromatic', r'Polar uncharged', r'Negatively charged', r'Positively charged', r'Non-standard']
            handles = [plt.Rectangle((0, 0), 1, 1, facecolor=colors[label], edgecolor='black') for label in range(len(labels))]
            fig.legend(handles, labels, fontsize=15, fancybox=True, shadow=True, loc='lower center', ncol=6)

            # Adjust plots
            fig.tight_layout()
            fig.subplots_adjust(bottom=0.12)

            # Save page
            pdf.savefig()
            plt.close()


def write_results(fn: str, pdb: str, ligand: str, output: str, output_hydropathy: str = None, volume: dict = None, area: dict = None, max_depth: dict = None, avg_depth: dict = None, avg_hydropathy: dict = None, residues: dict = None, frequencies: dict = None, step: float = 0.6) -> None:
    """
    Writes file paths and cavity characterization to TOML-formatted file

    Parameters
    ----------
        fn (str): path to KVFinder results TOML-formatted file
        pdb (str): path to input PDB file
        ligand (str): path to ligand PDB file
        output (str): path to cavity PDB file
        output_hydropathy (str): path to hydropathy PDB file
        volume (dict): dictionary with cavity name/volume pairs
        area (dict): dictionary with cavity name/area pairs
        max_depth (dict): dictionary with cavity name/maximum depth pairs
        avg_depth (dict): dictionary with cavity name/average depth pairs
        avg_hydropapthy: dictionary with cavity name/average hydropathy pairs
        residues (dict): dictionary with cavity name/list of interface residues pairs
        frequencies (dict): dictionary with frequencies of residues and class of residues of each detected cavity
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
    if output:
        output = os.path.abspath(output)
    if output_hydropathy:
        output_hydropathy = os.path.abspath(output_hydropathy)

    # Create results dictionary
    results = {
        'FILES': {
            'INPUT': pdb,
            'LIGAND': ligand,
            'OUTPUT': output,
            'HYDROPATHY': output_hydropathy,
        },
        'PARAMETERS': {
            'STEP': step,
        },
        'RESULTS': {
            'VOLUME': volume,
            'AREA': area,
            'MAX_DEPTH': max_depth,
            'AVG_DEPTH': avg_depth,
            'AVG_HYDROPATHY': avg_hydropathy,
            'RESIDUES': residues,
            'FREQUENCY': frequencies,
        }
    }

    # Create base directories of results TOML file
    os.makedirs(os.path.abspath(os.path.dirname(fn)), exist_ok=True)

    # Write results to TOML file
    with open(fn, "w") as f:
        f.write("# pyKVFinder results\n\n")
        toml.dump(results, f)
