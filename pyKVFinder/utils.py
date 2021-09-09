import os
import logging
import argparse
import pathlib
import numpy
from typing import Dict, List, Union, Optional

__all__ = [
    "read_vdw",
    "read_pdb",
    "read_xyz",
    "read_cavity",
    "calculate_frequencies",
    "plot_frequencies",
    "write_results",
]

VDW = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/vdw.dat")


def read_vdw(
    fn: Optional[Union[str, pathlib.Path]] = None
) -> Dict[str, Dict[str, float]]:
    """Reads van der Waals radii from .dat file.

    Parameters
    ----------
    fn : Optional[Union[str, pathlib.Path]], optional
        A path to a van der Waals radii file, by default None. If None, apply the built-in van der
        Waals radii file: `vdw.dat`.

    Returns
    -------
    vdw : Dict[str, Dict[str, float]]
        A dictionary containing radii values.

    Raises
    ------
    TypeError
        `fn` must be a string or a pathlib.Path.
    ValueError
        A line in `vdw` has incorrect format. The values must be double
        tab-separated.
    ValueError
        A line in `vdw` has an incorrect radius type for an atom.

    Note
    ----
    The van der Waals radii file defines the radius values for each
    atom by residue and when not defined, it uses a generic value
    based on the atom type (see `van der Waals file template`).
    The package contains a built-in van der Waals radii file: `vdw.dat`.
    """
    # Check argument
    if fn is not None:
        if type(fn) not in [str, pathlib.Path]:
            raise TypeError("`fn` must be a string or a pathlib.Path.")
    else:
        # Define default vdw file
        fn = VDW

    # Create vdw dictionary
    vdw = {}

    # Open fn
    with open(fn, "r") as f:
        # Read line with data only (ignore empty lines)
        lines = [
            line.replace(" ", "")
            for line in f.read().splitlines()
            if line.replace("\t\t", "")
        ]
        for line in lines:
            if not line.startswith("#"):
                if line.startswith(">"):
                    res = line.replace(">", "").replace("\t\t", "").replace(" ", "")
                    vdw[res] = {}
                else:
                    try:
                        atom, radius = line.split("\t\t")
                    except ValueError:
                        if len(line.split("\t\t")) != 2:
                            raise ValueError(
                                "A line in `vdw` has incorrect format. \
The values must be double tab-separated."
                            )
                    try:
                        vdw[res][atom] = float(radius)
                    except ValueError:
                        raise ValueError(
                            "A line in `vdw` has an incorrect radius type for \
an atom."
                        )

    return vdw


def _process_pdb_line(
    line: str, vdw: Dict[str, Dict[str, float]]
) -> List[Union[str, float, int]]:
    """Extracts ATOM and HETATM information of PDB line.

    Parameters
    ----------
    line : str
        A line of a valid PDB file
    vdw : Dict[str, Dict[str, Dict[str, float]]]
        A dictionary containing radii values.

    Returns
    -------
    atomic : List[Union[str, float, int]]
        A list with resnum, chain, resname, atom name, xyz coordinates and radius.
    """
    # Get PDB infomation
    atom = line[12:16].strip()
    resname = line[17:20].strip()
    resnum = int(line[22:26])
    chain = line[21]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    atom_symbol = line[76:78].strip().upper()

    # Get atom and radius from vdw
    if resname in vdw.keys() and atom in vdw[resname].keys():
        radius = vdw[resname][atom]
    else:
        radius = vdw["GEN"][atom_symbol]
        logging.info(
            f"Warning: Atom {atom} of residue {resname} \
not found in dictionary."
        )
        logging.info(
            f"Warning: Using generic atom {atom_symbol} \
radius: {radius} \u00c5."
        )

    # Prepare output
    atomic = [resnum, chain, resname, atom, x, y, z, radius]

    return atomic


def read_pdb(
    fn: Union[str, pathlib.Path], vdw: Optional[Dict[str, Dict[str, float]]] = None
) -> numpy.ndarray:
    """Reads PDB file into numpy.ndarrays.

    Parameters
    ----------
    fn : Union[str, pathlib.Path]
        A path to PDB file.
    vdw : Dict[str, Dict[str, float]], optional
        A dictionary containing radii values, by default None. If None, use output of `pyKVFinder.read_vdw()`.

    Returns
    -------
    atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.

    Raises
    ------
    TypeError
        `fn` must be a string or a pathlib.Path.

    Note
    ----
    The van der Waals radii file defines the radius values for each atom
    by residue and when not defined, it uses a generic value based on the
    atom type. The function by default loads the built-in van der Waals radii
    file: `vdw.dat`.
    """
    # Check arguments
    if type(fn) not in [str, pathlib.Path]:
        raise TypeError("`fn` must be a string or a pathlib.Path.")

    # Define default vdw file
    if vdw is None:
        vdw = read_vdw(VDW)

    # Create lists
    atomic = []

    with open(fn, "r") as f:
        for line in f.readlines():
            if line[:4] == "ATOM" or line[:6] == "HETATM":
                atomic.append(_process_pdb_line(line, vdw))

    return numpy.asarray(atomic)


def read_xyz(
    fn: Union[str, pathlib.Path], vdw: Optional[Dict[str, Dict[str, float]]] = None
) -> numpy.ndarray:
    """Reads XYZ file into numpy.ndarrays.

    Parameters
    ----------
    fn : Union[str, pathlib.Path]
        A path to XYZ file.
    vdw : Dict[str, Dict[str, float]], optional
        A dictionary containing radii values, by default None. If None, use output of `pyKVFinder.read_vdw()`.

    Returns
    -------
    atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.

    Raises
    ------
    TypeError
        `fn` must be a string or a pathlib.Path.

    Note
    ----
    The van der Waals radii file defines the radius values for each atom
    by residue and when not defined, it uses a generic value based on the
    atom type. The function by default loads the built-in van der Waals radii
    file: `vdw.dat`.
    """
    # Check arguments
    if type(fn) not in [str, pathlib.Path]:
        raise TypeError("`fn` must be a string or a pathlib.Path.")

    # Define default vdw file
    if vdw is None:
        vdw = read_vdw(VDW)

    # Create lists
    atomic = []

    # Start resnum
    resnum = 0

    # Read XYZ file
    with open(fn, "r") as f:
        for line in f.readlines():
            line = line.split()
            if len(line) == 4:
                # Get PDB information
                atom_symbol = line[0]
                x = float(line[1])
                y = float(line[2])
                z = float(line[3])

                # Get radius (generic value)
                radius = vdw["GEN"][atom_symbol]

                # Get resnum
                resnum += 1

                # Append data
                atomic.append([resnum, "A", "UNK", atom_symbol, x, y, z, radius])

    return numpy.asarray(atomic)


def _read_cavity(cavity: Union[str, pathlib.Path]) -> numpy.ndarray:
    """Reads xyz coordinates and labels of a cavities file into numpy.ndarray.

    Parameters
    ----------
    cavity : Union[str, pathlib.Path]
        A path to a PDB-formatted file of cavities.

    Returns
    -------
    xyzl : numpy.ndarray
        A numpy.ndarray with xyz coordinates and cavity label for each cavity point.
    """
    from .grid import _get_cavity_label

    # Create xyzl (xyz coordinates and cavity label)
    xyzl = []

    # Read cavity file into list
    with open(cavity, "r") as f:
        for line in f.readlines():
            if line[:4] == "ATOM" or line[:6] == "HETATM":
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                label = _get_cavity_label(line[17:20].strip())
                xyzl.append([x, y, z, label])

    return numpy.asarray(xyzl)


def read_cavity(
    cavity: Union[str, pathlib.Path],
    receptor: Union[str, pathlib.Path],
    step: Union[float, int] = 0.6,
    probe_in: Union[float, int] = 1.4,
    probe_out: Union[float, int] = 4.0,
    surface: str = "SES",
    vdw: Optional[Dict[str, Dict[str, float]]] = None,
    nthreads: Optional[int] = None,
    verbose: bool = False,
) -> numpy.ndarray:
    """Read cavities and receptor inside a 3D grid.

    Parameters
    ----------
    cavity : Union[str, pathlib.Path]
        A path to a PDB file of cavities.
    receptor : Union[str, pathlib.Path]
        A path to a PDB or XYZ file of the receptor.
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6.
    probe_in : Union[float, int], optional
        Probe In size (A), by default 1.4.
    probe_out : Union[float, int], optional
        Probe Out size (A), by default 4.0.
    surface : str, optional
        Surface representation. Keywords options are SES (Solvent Excluded Surface) or SAS (Solvent
        Accessible Surface), by default "SES".
    vdw : Dict[str, Dict[str, float]], optional
        A dictionary containing radii values, by default None. If None, use output of `pyKVFinder.read_vdw()`.
    nthreads : Optional[int], optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    grid : numpy.ndarray
        Cavity and receptor points in the 3D grid (grid[nx][ny][nz]).
        Grid array has integer labels in each position, that are:

            * -1: bulk points or empty space points;

            * 0: biomolecule points;

            * >=2: cavity points.

    Raises
    ------
    TypeError
        `cavity` must be a string or a pathlib.Path.
    TypeError
        `receptor` must be a string or a pathlib.Path.
    TypeError
        `target` must have .pdb or .xyz extension.
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.
    TypeError
        `probe_in` must be a non-negative real number.
    ValueError
        `probe_in` must be a non-negative real number.
    TypeError
        `probe_out` must be a non-negative real number.
    ValueError
        `probe_out` must be a non-negative real number.
    ValueError
        `probe_out` must be greater than `probe_in`.
    TypeError
        `surface` must be a str.
    TypeError
        `nthreads` must be a positive integer.
    ValueError
        `nthreads` must be a positive integer.
    TypeError
        `verbose` must be a boolean.
    ValueError
        `surface` must be SAS or SES, not {surface}.
    """
    from .grid import get_vertices, _get_sincos, _get_dimensions
    from _pyKVFinder import _fill_receptor, _fill_cavity

    # Check arguments
    if type(cavity) not in [str, pathlib.Path]:
        raise TypeError("`cavity` must be a string or a pathlib.Path.")
    if type(receptor) not in [str, pathlib.Path]:
        raise TypeError("`receptor` must be a string or a pathlib.Path.")
    elif not receptor.endswith(".pdb") and not receptor.endswith(".xyz"):
        raise TypeError("`receptor` must have .pdb or .xyz extension.")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")
    if type(probe_in) not in [float, int]:
        raise TypeError("`probe_in` must be a non-negative real number.")
    elif probe_in < 0.0:
        raise ValueError("`probe_in` must be a non-negative real number.")
    if type(probe_out) not in [float, int]:
        raise TypeError("`probe_out` must be a non-negative real number.")
    elif probe_out < 0.0:
        raise ValueError("`probe_out` must be a non-negative real number.")
    elif probe_out < probe_in:
        raise ValueError("`probe_out` must be greater than `probe_in`.")
    if type(surface) not in [str]:
        raise TypeError("`surface` must be a str.")
    if nthreads is None:
        nthreads = os.cpu_count() - 1
    else:
        if type(nthreads) not in [int]:
            raise TypeError("`nthreads` must be a positive integer.")
        elif nthreads <= 0:
            raise ValueError("`nthreads` must be a positive integer.")
    if type(verbose) not in [bool]:
        raise TypeError("`verbose` must be a boolean.")

    # Convert types
    if type(step) == int:
        step = float(step)
    if type(probe_in) == int:
        probe_in = float(probe_in)
    if type(probe_out) == int:
        probe_out = float(probe_out)

    # Insert receptor inside 3D grid
    if verbose:
        print(f"> Inserting {receptor} into 3D grid")

    # Define default vdw file
    if vdw is None:
        vdw = read_vdw(VDW)

    # Load receptor coordinates and radii
    if receptor.endswith(".pdb"):
        atomic = read_pdb(receptor, vdw)
    elif receptor.endswith(".xyz"):
        atomic = read_xyz(receptor, vdw)

    # Extract xyzr from atomic
    xyzr = atomic[:, 4:].astype(numpy.float64)

    # Get vertices
    vertices = get_vertices(atomic, probe_out, step)

    # Get sincos
    sincos = _get_sincos(vertices)

    # Get dimensions
    nx, ny, nz = _get_dimensions(vertices, step)

    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Calculate number of voxels
    nvoxels = nx * ny * nz

    if surface == "SES":
        if verbose:
            print("> Surface representation: Solvent Excluded Surface (SES).")
        surface = True
    elif surface == "SAS":
        if verbose:
            print("> Surface representation: Solvent Accessible Surface (SAS).")
        surface = False
    else:
        raise ValueError(f"`surface` must be SAS or SES, not {surface}.")

    # Fill grid with receptor
    grid = _fill_receptor(
        nvoxels,
        nx,
        ny,
        nz,
        xyzr,
        P1,
        sincos,
        step,
        probe_in,
        surface,
        nthreads,
        verbose,
    ).reshape(nx, ny, nz)

    # Insert cavities inside 3D grid
    if verbose:
        print(f"> Inserting {cavity} into 3D grid")

    # Load cavities coordinates and labels
    xyzl = _read_cavity(cavity)

    # Fill grid with cavities
    _fill_cavity(grid, xyzl, P1, sincos, step, nthreads)

    return grid


def _process_box(args: argparse.Namespace) -> Dict[str, List[float]]:
    """Gets xyz coordinates of 3D grid vertices.

    Parameters
    ----------
    args (argparse.Namespace)
        Arguments passes by argparser CLI.

    Returns
    -------
    box : Dict[str, List[float]]
        A dictionary with a xyz coordinates (p1: origin,
        p2: X-axis, p3: Y-axis, p4: Z-axis) for each point.
    """
    # Create box parameter
    box = {
        "p1": args.vertices[0],
        "p2": args.vertices[1],
        "p3": args.vertices[2],
        "p4": args.vertices[3],
    }

    # Adjust if box adjustment mode
    if args.box:
        # Get probe out additions
        # p1 = (x1, y1, z1)
        x1 = (
            -(args.probe_out * args.sincos[3])
            - (args.probe_out * args.sincos[0] * args.sincos[2])
            + (args.probe_out * args.sincos[1] * args.sincos[2])
        )
        y1 = -(args.probe_out * args.sincos[1]) - (args.probe_out * args.sincos[0])
        z1 = (
            -(args.probe_out * args.sincos[2])
            + (args.probe_out * args.sincos[0] * args.sincos[3])
            - (args.probe_out * args.sincos[1] * args.sincos[3])
        )
        # p2 = (x2, y2, z2)
        x2 = (
            (args.probe_out * args.sincos[3])
            - (args.probe_out * args.sincos[0] * args.sincos[2])
            + (args.probe_out * args.sincos[1] * args.sincos[2])
        )
        y2 = -(args.probe_out * args.sincos[1]) - (args.probe_out * args.sincos[0])
        z2 = (
            (args.probe_out * args.sincos[2])
            + (args.probe_out * args.sincos[0] * args.sincos[3])
            - (args.probe_out * args.sincos[1] * args.sincos[3])
        )
        # p3 = (x3, y3, z3)
        x3 = (
            -(args.probe_out * args.sincos[3])
            + (args.probe_out * args.sincos[0] * args.sincos[2])
            + (args.probe_out * args.sincos[1] * args.sincos[2])
        )
        y3 = (args.probe_out * args.sincos[1]) - (args.probe_out * args.sincos[0])
        z3 = (
            -(args.probe_out * args.sincos[2])
            - (args.probe_out * args.sincos[0] * args.sincos[3])
            - (args.probe_out * args.sincos[1] * args.sincos[3])
        )
        # p4 = (x4, y4, z4)
        x4 = (
            -(args.probe_out * args.sincos[3])
            - (args.probe_out * args.sincos[0] * args.sincos[2])
            - (args.probe_out * args.sincos[1] * args.sincos[2])
        )
        y4 = -(args.probe_out * args.sincos[1]) + (args.probe_out * args.sincos[0])
        z4 = (
            -(args.probe_out * args.sincos[2])
            + (args.probe_out * args.sincos[0] * args.sincos[3])
            + (args.probe_out * args.sincos[1] * args.sincos[3])
        )

        # Remove probe out addition
        box["p1"] -= numpy.array([x1, y1, z1])
        box["p2"] -= numpy.array([x2, y2, z2])
        box["p3"] -= numpy.array([x3, y3, z3])
        box["p4"] -= numpy.array([x4, y4, z4])

    # Prepare to dict to toml module
    box["p1"] = numpy.around(box["p1"], 2).tolist()
    box["p2"] = numpy.around(box["p2"], 2).tolist()
    box["p3"] = numpy.around(box["p3"], 2).tolist()
    box["p4"] = numpy.around(box["p4"], 2).tolist()

    return box


def _write_parameters(args: argparse.Namespace) -> None:
    """Writes parameters used in cavity detection and characterization of
    pyKVFinder to TOML-formatted file.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments passes by argparser CLI.
    """
    import toml

    # Parameters filename
    fn = os.path.join(args.output_directory, f"{args.base_name}.parameters.toml")

    # Parameters dict
    parameters = {
        "FILES": {
            "INPUT": args.input,
            "LIGAND": args.ligand,
            "BASE_NAME": args.base_name,
            "OUTPUT_DIRECTORY": args.output_directory,
            "DICTIONARY": args.dictionary,
        },
        "SETTINGS": {
            "MODES": {
                "BOX_ADJUSTMENT": args.box,
                "LIGAND_ADJUSTMENT": True if args.ligand else False,
                "DEPTH": args.depth,
                "SURFACE": args.surface,
                "IGNORE_BACKBONE": args.ignore_backbone,
            },
            "STEP": args.step,
            "PROBES": {
                "PROBE_IN": args.probe_in,
                "PROBE_OUT": args.probe_out,
            },
            "CUTOFFS": {
                "VOLUME_CUTOFF": args.volume_cutoff,
                "LIGAND_CUTOFF": args.ligand_cutoff,
                "REMOVAL_DISTANCE": args.removal_distance,
            },
            "BOX": _process_box(args),
        },
    }

    # Write to TOML file
    with open(fn, "w") as param:
        toml.dump(parameters, param)


def calculate_frequencies(
    residues: Dict[str, List[List[str]]]
) -> Dict[str, Dict[str, Dict[str, int]]]:
    """Calculate frequencies of residues and class of residues
    (R1, R2, R3, R4 and R5) for detected cavities.

    Parameters
    ----------
    residues : Dict[str, List[List[str]]]
        A dictionary with a list of interface residues for each detected
        cavity.

    Returns
    -------
    frequencies : Dict[str, Dict[str, Dict[str, int]]]
        A dictionary with frequencies of residues and class for
        residues of each detected cavity.

    Note
    ----
    The cavity nomenclature is based on the integer label. The cavity
    marked with 2, the first integer corresponding to a cavity, is KAA, the
    cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

    Note
    ----
    The classes of residues are:

    * Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine.

    * Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine.

    * Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine.

    * Negatively charged (R4): Aspartate, Glutamate.

    * Positively charged (R5): Arginine, Histidine, Lysine.

    * Non-standard (RX): Non-standard residues.
    """
    # Create a dict for frequencies
    frequencies = {}

    # Get cavity name and residues list for each detected cavity
    for name, reslist in residues.items():
        # Create a dict for cavity name
        frequencies[name] = {
            "RESIDUES": {},
            "CLASS": {},
        }
        # Get unique residues names
        residues = [res[2] for res in reslist]
        reslist = sorted(list(set(residues)))

        # Get residues frequencies
        for res in reslist:
            frequencies[name]["RESIDUES"][res] = residues.count(res)

        # Get class frequencies
        frequencies[name]["CLASS"]["R1"] = (
            frequencies[name]["RESIDUES"].get("ALA", 0)
            + frequencies[name]["RESIDUES"].get("GLY", 0)
            + frequencies[name]["RESIDUES"].get("ILE", 0)
            + frequencies[name]["RESIDUES"].get("LEU", 0)
            + frequencies[name]["RESIDUES"].get("PRO", 0)
            + frequencies[name]["RESIDUES"].get("VAL", 0)
        )
        frequencies[name]["CLASS"]["R2"] = (
            frequencies[name]["RESIDUES"].get("PHE", 0)
            + frequencies[name]["RESIDUES"].get("TRP", 0)
            + frequencies[name]["RESIDUES"].get("TYR", 0)
        )
        frequencies[name]["CLASS"]["R3"] = (
            frequencies[name]["RESIDUES"].get("ASN", 0)
            + frequencies[name]["RESIDUES"].get("CYS", 0)
            + frequencies[name]["RESIDUES"].get("GLN", 0)
            + frequencies[name]["RESIDUES"].get("MET", 0)
            + frequencies[name]["RESIDUES"].get("SER", 0)
            + frequencies[name]["RESIDUES"].get("THR", 0)
        )
        frequencies[name]["CLASS"]["R4"] = frequencies[name]["RESIDUES"].get(
            "ASP", 0
        ) + frequencies[name]["RESIDUES"].get("GLU", 0)
        frequencies[name]["CLASS"]["R5"] = (
            frequencies[name]["RESIDUES"].get("ARG", 0)
            + frequencies[name]["RESIDUES"].get("HIS", 0)
            + frequencies[name]["RESIDUES"].get("LYS", 0)
        )
        frequencies[name]["CLASS"]["RX"] = len(residues) - sum(
            frequencies[name]["CLASS"].values()
        )

    return frequencies


def plot_frequencies(
    frequencies: Dict[str, Dict[str, Dict[str, int]]],
    fn: Union[str, pathlib.Path] = "barplots.pdf",
) -> None:
    """Plot bar charts of calculated frequencies (residues and classes of
    residues) for each detected cavity in a target PDF file.

    Parameters
    ----------
    frequencies : Dict[str, Dict[str, Dict[str, int]]]
        A dictionary with frequencies of residues and class for
        residues of each detected cavity.
    fn : Union[str, pathlib.Path], optional
        A path to PDF file for plotting bar charts of frequencies, by
        default `barplots.pdf`.

    Raises
    ------
    TypeError
        `fn` must be a string or a pathlib.Path.

    Note
    ----
    The cavity nomenclature is based on the integer label. The cavity
    marked with 2, the first integer corresponding to a cavity, is KAA, the
    cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

    Note
    ----
    The classes of residues are:

    * Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine.

    * Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine.

    * Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine.

    * Negatively charged (R4): Aspartate, Glutamate.

    * Positively charged (R5): Arginine, Histidine, Lysine.

    * Non-standard (RX): Non-standard residues.
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    # Check arguments
    if type(fn) not in [str, pathlib.Path]:
        raise TypeError("`fn` must be a string or a pathlib.Path.")

    # Create base directories of output PDF file
    os.makedirs(os.path.abspath(os.path.dirname(fn)), exist_ok=True)

    # Create a dictionary for standard amino acids
    tmp = {
        "ALA": 0,
        "ARG": 0,
        "ASN": 0,
        "ASP": 0,
        "CYS": 0,
        "GLN": 0,
        "GLU": 0,
        "GLY": 0,
        "HIS": 0,
        "ILE": 0,
        "LEU": 0,
        "LYS": 0,
        "MET": 0,
        "PHE": 0,
        "PRO": 0,
        "SER": 0,
        "THR": 0,
        "TRP": 0,
        "TYR": 0,
        "VAL": 0,
    }

    with PdfPages(fn) as pdf:
        # Standardize data
        ymax = 0
        for cavity_tag in frequencies.keys():
            # Include missing residues
            frequencies[cavity_tag]["RESIDUES"] = {
                **tmp,
                **frequencies[cavity_tag]["RESIDUES"],
            }
            # Get y maximum
            if ymax < max(frequencies[cavity_tag]["CLASS"].values()):
                ymax = max(frequencies[cavity_tag]["CLASS"].values())
        ymax += 1

        # Pdf plots
        for cavity_tag in frequencies.keys():
            # Create page
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 9), dpi=300)
            fig.suptitle(r"Cavity " + f"{cavity_tag}", fontsize=30)

            # Frequency residues
            x = list(frequencies[cavity_tag]["RESIDUES"].keys())
            y = frequencies[cavity_tag]["RESIDUES"].values()
            colors = [
                "tab:cyan",
                "tab:purple",
                "tab:green",
                "tab:red",
                "tab:green",
                "tab:green",
                "tab:red",
                "tab:cyan",
                "tab:purple",
                "tab:cyan",
                "tab:cyan",
                "tab:purple",
                "tab:green",
                "tab:orange",
                "tab:cyan",
                "tab:green",
                "tab:green",
                "tab:orange",
                "tab:orange",
                "tab:cyan",
            ]
            for _ in range(len(x) - len(colors)):
                colors.append("tab:gray")

            ax1.bar(x, y, align="center", edgecolor="black", color=colors)
            ax1.set_xlabel(None)
            ax1.set_xlim(-1, len(x))
            ax1.tick_params(axis="x", labelsize=15, rotation=45)
            ax1.tick_params(axis="y", labelsize=20)
            ax1.set_ylabel(r"Frequency", fontsize=20)
            ax1.set_ylim(0, ymax)
            ax1.grid(which="major", axis="y", linestyle="--")

            # Frequency classes
            x = list(frequencies[cavity_tag]["CLASS"].keys())
            y = frequencies[cavity_tag]["CLASS"].values()
            colors = [
                "tab:cyan",
                "tab:orange",
                "tab:green",
                "tab:red",
                "tab:purple",
                "tab:gray",
            ]

            ax2.bar(x=x, height=y, align="center", edgecolor="black", color=colors)
            ax2.set_xlabel(None)
            ax2.set_xlim(-1, len(x))
            ax2.tick_params(axis="x", labelsize=20)
            ax2.tick_params(axis="y", labelsize=20)
            ax2.set_ylabel(None)
            ax2.set_ylim(0, ymax)
            ax2.grid(which="major", axis="y", linestyle="--")

            # Legend
            labels = [
                r"Aliphatic apolar",
                r"Aromatic",
                r"Polar uncharged",
                r"Negatively charged",
                r"Positively charged",
                r"Non-standard",
            ]
            handles = [
                plt.Rectangle((0, 0), 1, 1, facecolor=colors[label], edgecolor="black")
                for label in range(len(labels))
            ]
            fig.legend(
                handles,
                labels,
                fontsize=15,
                fancybox=True,
                shadow=True,
                loc="lower center",
                ncol=6,
            )

            # Adjust plots
            fig.tight_layout()
            fig.subplots_adjust(bottom=0.12)

            # Save page
            pdf.savefig()
            plt.close()


def write_results(
    fn: Union[str, pathlib.Path],
    input: Optional[Union[str, pathlib.Path]],
    ligand: Optional[Union[str, pathlib.Path]],
    output: Optional[Union[str, pathlib.Path]],
    output_hydropathy: Optional[Union[str, pathlib.Path]] = None,
    volume: Optional[Dict[str, float]] = None,
    area: Optional[Dict[str, float]] = None,
    max_depth: Optional[Dict[str, float]] = None,
    avg_depth: Optional[Dict[str, float]] = None,
    avg_hydropathy: Optional[Dict[str, float]] = None,
    residues: Optional[Dict[str, List[List[str]]]] = None,
    frequencies: Optional[Dict[str, Dict[str, Dict[str, int]]]] = None,
    step: Union[float, int] = 0.6,
) -> None:
    """Writes file paths and cavity characterization to TOML-formatted file.

    Parameters
    ----------
    fn : Union[str, pathlib.Path]
        A path to TOML-formatted file for writing file paths and
        cavity characterization (volume, area, depth [optional] and interface
        residues) per cavity detected.
    input : Union[str, pathlib.Path], optional
        A path to input PDB or XYZ file.
    ligand : Union[str, pathlib.Path], optional
        A path to ligand PDB or XYZ file.
    output : Union[str, pathlib.Path], optional
        A path to cavity PDB file.
    output_hydropathy : Union[str, pathlib.Path], optional
        A path to hydropathy PDB file (surface points mapped with a
        hydrophobicity scale), by default None.
    volume : Dict[str, float], optional
        A dictionary with volume of each detected cavity, by default None.
    area : Dict[str, float], optional
        A dictionary with area of each detected cavity, by default None.
    max_depth : Dict[str, float], optional
        A dictionary with maximum depth of each detected cavity, by default
        None.
    avg_depth : Dict[str, float], optional
        A dictionary with average depth of each detected cavity, by default
        None.
    avg_hydropapthy : Dict[str, float], optional
        A dictionary with average hydropathy of each detected cavity and range
        of the hydrophobicity scale mapped, by default None.
    residues : Dict[str, List[List[str]]], optional
        A dictionary with interface residues of each detected cavity, by
        default None.
    frequencies : Dict[str, Dict[str, Dict[str, int]]], optional
        A dictionary with frequencies of interface residues and classes of
        residues of each detected cavity, by default None.
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6.

    Raises
    ------
    TypeError
        `fn` must be a string or a pathlib.Path.
    TypeError
        `input` must be a string or a pathlib.Path.
    TypeError
        `ligand` must be a string or a pathlib.Path.
    TypeError
        `output` must be a string or a pathlib.Path.
    TypeError
        `output_hydropathy` must be a string or a pathlib.Path.
    TypeError
        `volume` must be a dictionary.
    TypeError
        `area` must be a dictionary.
    TypeError
        `max_depth` must be a dictionary.
    TypeError
        `avg_depth` must be a dictionary.
    TypeError
        `avg_hydropathy` must be a dictionary.
    TypeError
        `residues` must be a dictionary.
    TypeError
        `frequencies` must be a dictionary.
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.

    Note
    ----
    The cavity nomenclature is based on the integer label. The cavity
    marked with 2, the first integer corresponding to a cavity, is KAA, the
    cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.
    """
    import toml

    # Check arguments
    if type(fn) not in [str, pathlib.Path]:
        raise TypeError("`fn` must be a string or a pathlib.Path.")
    if input is not None:
        if type(input) not in [str, pathlib.Path]:
            raise TypeError("`input` must be a string or a pathlib.Path.")
    if ligand is not None:
        if type(ligand) not in [str, pathlib.Path]:
            raise TypeError("`ligand` must be a string or a pathlib.Path.")
    if output is not None:
        if type(output) not in [str, pathlib.Path]:
            raise TypeError("`output` must be a string or a pathlib.Path.")
    if output_hydropathy is not None:
        if type(output_hydropathy) not in [str, pathlib.Path]:
            raise TypeError("`output_hydropathy` must be a string or a pathlib.Path.")
    if volume is not None:
        if type(volume) not in [dict]:
            raise TypeError("`volume` must be a dictionary.")
    if area is not None:
        if type(area) not in [dict]:
            raise TypeError("`area` must be a dictionary.")
    if max_depth is not None:
        if type(max_depth) not in [dict]:
            raise TypeError("`max_depth` must be a dictionary.")
    if avg_depth is not None:
        if type(avg_depth) not in [dict]:
            raise TypeError("`avg_depth` must be a dictionary.")
    if avg_hydropathy is not None:
        if type(avg_hydropathy) not in [dict]:
            raise TypeError("`avg_hydropathy` must be a dictionary.")
    if residues is not None:
        if type(residues) not in [dict]:
            raise TypeError("`residues` must be a dictionary.")
    if frequencies is not None:
        if type(frequencies) not in [dict]:
            raise TypeError("`frequencies` must be a dictionary.")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")

    # Convert types
    if type(step) == int:
        step = float(step)

    # Create base directories of results
    os.makedirs(os.path.abspath(os.path.dirname(fn)), exist_ok=True)

    # Prepare paths
    input = os.path.abspath(input)
    if ligand:
        ligand = os.path.abspath(ligand)
    if output:
        output = os.path.abspath(output)
    if output_hydropathy:
        output_hydropathy = os.path.abspath(output_hydropathy)

    # Create results dictionary
    results = {
        "FILES": {
            "INPUT": input,
            "LIGAND": ligand,
            "OUTPUT": output,
            "HYDROPATHY": output_hydropathy,
        },
        "PARAMETERS": {
            "STEP": step,
        },
        "RESULTS": {
            "VOLUME": volume,
            "AREA": area,
            "MAX_DEPTH": max_depth,
            "AVG_DEPTH": avg_depth,
            "AVG_HYDROPATHY": avg_hydropathy,
            "RESIDUES": residues,
            "FREQUENCY": frequencies,
        },
    }

    # Create base directories of results TOML file
    os.makedirs(os.path.abspath(os.path.dirname(fn)), exist_ok=True)

    # Write results to TOML file
    with open(fn, "w") as f:
        f.write("# pyKVFinder results\n\n")
        toml.dump(results, f)
