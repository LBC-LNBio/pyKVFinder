import os

import numpy

from .cavity import _get_cavity_label, _get_cavity_name, _select_cavities
from .geometry import _get_sincos


def _process_residues(
    raw: list[str], ncav: int, selection: list[int] | None = None
) -> dict[str, list[list[str]]]:
    """Processes raw list of residues from _constitutional to a list of
    residue information per cavity name.

    Parameters
    ----------
    raw : list[str]
        A list of residues with cavities separated by '-1'.
    ncav : int
        Number of cavities.
    selection : list[int], optional
        A list of integer labels of each cavity to be selected, by default None.

    Returns
    -------
    residues : dict[str, list[list[str]]]
        A dictionary with a list of interface residues for each detected
        cavity.
    """
    from itertools import groupby

    residues = {}

    if selection is None:
        selection = list(range(ncav + 1))
    else:
        selection = [sele - 2 for sele in selection]

    index = 0
    for flag, cavity_residues in groupby(raw, lambda res: res == "-1"):
        if not flag:
            key = _get_cavity_name(selection[index])
            residues[key] = [
                item.split("_") for item in list(dict.fromkeys(cavity_residues))
            ]
            index += 1

    return residues


def constitutional(
    cavities: numpy.ndarray,
    atomic: numpy.ndarray | list[list[str | float | int]],
    vertices: numpy.ndarray | list[list[float]],
    step: float | int = 0.6,
    probe_in: float | int = 1.4,
    ignore_backbone: bool = False,
    selection: list[int] | list[str] | None = None,
    nthreads: int | None = None,
    verbose: bool = False,
) -> dict[str, list[list[str]]]:
    """Constitutional characterization (interface residues) of the detected
    cavities.

    Parameters
    ----------
    cavities : numpy.ndarray
        Cavity points in the 3D grid (cavities[nx][ny][nz]).
        Cavities array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule points;

            * 1: empty space points;

            * >=2: cavity points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.
    atomic : numpy.ndarray | list[list[str | float | int]]
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    vertices : numpy.ndarray | list[list[float]]
        A numpy.ndarray or a list with xyz vertices coordinates (origin,
        X-axis, Y-axis, Z-axis).
    step : float | int, optional
        Grid spacing (A), by default 0.6.
    probe_in : float | int, optional
        Probe In size (A), by default 1.4.
    ignore_backbone : bool, optional
        Whether to ignore backbone atoms (C, CA, N, O) when defining interface
        residues, by default False.
    selection : list[int] | list[str] | None, optional
        A list of integer labels or a list of cavity names to be selectedA list of integer labels of each cavity to be selected, by default None.
    nthreads : int | None, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    residues: dict[str, list[list[str]]]
        A dictionary with a list of interface residues for each detected
        cavity.

    Raises
    ------
    TypeError
        `cavities` must be a numpy.ndarray.
    ValueError
        `cavities` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `atomic` must be a list or a numpy.ndarray.
    ValueError
        `atomic` has incorrect shape. It must be (n, 8).
    TypeError
        `vertices` must be a list or a numpy.ndarray.
    ValueError
        `vertices` has incorrect shape. It must be (4, 3).
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.
    TypeError
        `probe_in` must be a non-negative real number.
    ValueError
        `probe_in` must be a non-negative real number.
    TypeError
        `ignore_backbone` must be a boolean.
    TypeError
        `selection` must be a list of strings (cavity names) or integers (cavity labels).
    ValueError
        Invalid `selection`: {selection}.
    TypeError
        `nthreads` must be a positive integer.
    ValueError
        `nthreads` must be a positive integer.
    TypeError
        `verbose` must be a boolean.

    Note
    ----
    The cavity nomenclature is based on the integer label. The cavity marked
    with 2, the first integer corresponding to a cavity, is KAA, the cavity
    marked with 3 is KAB, the cavity marked with 4 is KAC and so on. After KZZ,
    the next cavity is Kaa, Kab, and so on. The naming convention supports up to
    1352 cavities.

    Note
    ----
    The classes of residues are:

    * Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine.

    * Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine.

    * Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine.

    * Negatively charged (R4): Aspartate, Glutamate.

    * Positively charged (R5): Arginine, Histidine, Lysine.

    * Non-standard (RX): Non-standard residues.

    See Also
    --------
    read_pdb
    read_xyz
    detect
    calculate_frequencies
    write_results

    Example
    -------
    With the cavity points identified with ``detect``, atomic data read with ``read_pdb`` or ``read_xyz``, we can perform a constitutional characterization, that identifies the interface residues surrounding the cavities:

    >>> from pyKVFinder import constitutional
    >>> residues = constitutional(cavities, atomic, vertices)
    >>> residues
    {'KAA': [['14', 'E', 'SER'], ['15', 'E', 'VAL'], ['18', 'E', 'PHE'], ['19', 'E', 'LEU'], ['100', 'E', 'PHE'], ['152', 'E', 'LEU'], ['155', 'E', 'GLU'], ['156', 'E', 'TYR'], ['292', 'E', 'LYS'], ['302', 'E', 'TRP'], ['303', 'E', 'ILE'], ['306', 'E', 'TYR']], 'KAB': [['18', 'E', 'PHE'], ['22', 'E', 'ALA'], ['25', 'E', 'ASP'], ['26', 'E', 'PHE'], ['29', 'E', 'LYS'], ['97', 'E', 'ALA'], ['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['156', 'E', 'TYR']], 'KAC': [['141', 'E', 'PRO'], ['142', 'E', 'HIS'], ['144', 'E', 'ARG'], ['145', 'E', 'PHE'], ['148', 'E', 'ALA'], ['299', 'E', 'THR'], ['300', 'E', 'THR'], ['305', 'E', 'ILE'], ['310', 'E', 'VAL'], ['311', 'E', 'GLU'], ['313', 'E', 'PRO']], 'KAD': [['122', 'E', 'TYR'], ['124', 'E', 'ALA'], ['176', 'E', 'GLN'], ['318', 'E', 'PHE'], ['320', 'E', 'GLY'], ['321', 'E', 'PRO'], ['322', 'E', 'GLY'], ['323', 'E', 'ASP']], 'KAE': [['95', 'E', 'LEU'], ['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['100', 'E', 'PHE'], ['103', 'E', 'LEU'], ['104', 'E', 'VAL'], ['105', 'E', 'LYS'], ['106', 'E', 'LEU']], 'KAF': [['123', 'E', 'VAL'], ['124', 'E', 'ALA'], ['175', 'E', 'ASP'], ['176', 'E', 'GLN'], ['181', 'E', 'GLN']], 'KAG': [['34', 'E', 'SER'], ['37', 'E', 'THR'], ['96', 'E', 'GLN'], ['106', 'E', 'LEU'], ['107', 'E', 'GLU'], ['108', 'E', 'PHE'], ['109', 'E', 'SER']], 'KAH': [['49', 'E', 'LEU'], ['50', 'E', 'GLY'], ['51', 'E', 'THR'], ['52', 'E', 'GLY'], ['53', 'E', 'SER'], ['54', 'E', 'PHE'], ['55', 'E', 'GLY'], ['56', 'E', 'ARG'], ['57', 'E', 'VAL'], ['70', 'E', 'ALA'], ['72', 'E', 'LYS'], ['74', 'E', 'LEU'], ['84', 'E', 'GLN'], ['87', 'E', 'HIS'], ['88', 'E', 'THR'], ['91', 'E', 'GLU'], ['104', 'E', 'VAL'], ['120', 'E', 'MET'], ['121', 'E', 'GLU'], ['122', 'E', 'TYR'], ['123', 'E', 'VAL'], ['127', 'E', 'GLU'], ['166', 'E', 'ASP'], ['168', 'E', 'LYS'], ['170', 'E', 'GLU'], ['171', 'E', 'ASN'], ['173', 'E', 'LEU'], ['183', 'E', 'THR'], ['184', 'E', 'ASP'], ['186', 'E', 'GLY'], ['187', 'E', 'PHE'], ['201', 'E', 'THR'], ['327', 'E', 'PHE']], 'KAI': [['131', 'E', 'HIS'], ['138', 'E', 'PHE'], ['142', 'E', 'HIS'], ['146', 'E', 'TYR'], ['174', 'E', 'ILE'], ['314', 'E', 'PHE']], 'KAJ': [['33', 'E', 'PRO'], ['89', 'E', 'LEU'], ['92', 'E', 'LYS'], ['93', 'E', 'ARG'], ['96', 'E', 'GLN'], ['349', 'E', 'GLU'], ['350', 'E', 'PHE']], 'KAK': [['157', 'E', 'LEU'], ['162', 'E', 'LEU'], ['163', 'E', 'ILE'], ['164', 'E', 'TYR'], ['185', 'E', 'PHE'], ['188', 'E', 'ALA']], 'KAL': [['49', 'E', 'LEU'], ['127', 'E', 'GLU'], ['129', 'E', 'PHE'], ['130', 'E', 'SER'], ['326', 'E', 'ASN'], ['327', 'E', 'PHE'], ['328', 'E', 'ASP'], ['330', 'E', 'TYR']], 'KAM': [['51', 'E', 'THR'], ['55', 'E', 'GLY'], ['56', 'E', 'ARG'], ['73', 'E', 'ILE'], ['74', 'E', 'LEU'], ['75', 'E', 'ASP'], ['115', 'E', 'ASN'], ['335', 'E', 'ILE'], ['336', 'E', 'ARG']], 'KAN': [['165', 'E', 'ARG'], ['166', 'E', 'ASP'], ['167', 'E', 'LEU'], ['199', 'E', 'CYS'], ['200', 'E', 'GLY'], ['201', 'E', 'THR'], ['204', 'E', 'TYR'], ['205', 'E', 'LEU'], ['206', 'E', 'ALA'], ['209', 'E', 'ILE'], ['219', 'E', 'VAL'], ['220', 'E', 'ASP'], ['223', 'E', 'ALA']], 'KAO': [['48', 'E', 'THR'], ['51', 'E', 'THR'], ['56', 'E', 'ARG'], ['330', 'E', 'TYR'], ['331', 'E', 'GLU']], 'KAP': [['222', 'E', 'TRP'], ['238', 'E', 'PHE'], ['253', 'E', 'GLY'], ['254', 'E', 'LYS'], ['255', 'E', 'VAL'], ['273', 'E', 'LEU']], 'KAQ': [['207', 'E', 'PRO'], ['208', 'E', 'GLU'], ['211', 'E', 'LEU'], ['213', 'E', 'LYS'], ['275', 'E', 'VAL'], ['277', 'E', 'LEU']], 'KAR': [['237', 'E', 'PRO'], ['238', 'E', 'PHE'], ['249', 'E', 'LYS'], ['254', 'E', 'LYS'], ['255', 'E', 'VAL'], ['256', 'E', 'ARG']]}

    However, users may opt to ignore backbones contacts (C, CA, N, O) with the cavity when defining interface residues. Then, users must set ``ignore_backbone`` flag to ``True``.

    >>> residues = constitutional(cavities, atomic, vertices, ignore_backbone=True)
    >>> residues
    {'KAA': [['15', 'E', 'VAL'], ['18', 'E', 'PHE'], ['19', 'E', 'LEU'], ['100', 'E', 'PHE'], ['152', 'E', 'LEU'], ['155', 'E', 'GLU'], ['156', 'E', 'TYR'], ['292', 'E', 'LYS'], ['302', 'E', 'TRP'], ['303', 'E', 'ILE'], ['306', 'E', 'TYR']], 'KAB': [['18', 'E', 'PHE'], ['22', 'E', 'ALA'], ['25', 'E', 'ASP'], ['26', 'E', 'PHE'], ['29', 'E', 'LYS'], ['99', 'E', 'ASN'], ['156', 'E', 'TYR']], 'KAC': [['144', 'E', 'ARG'], ['145', 'E', 'PHE'], ['148', 'E', 'ALA'], ['299', 'E', 'THR'], ['300', 'E', 'THR'], ['305', 'E', 'ILE'], ['310', 'E', 'VAL'], ['311', 'E', 'GLU'], ['313', 'E', 'PRO']], 'KAD': [['122', 'E', 'TYR'], ['124', 'E', 'ALA'], ['176', 'E', 'GLN'], ['318', 'E', 'PHE']], 'KAE': [['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['103', 'E', 'LEU'], ['105', 'E', 'LYS'], ['106', 'E', 'LEU']], 'KAF': [['123', 'E', 'VAL'], ['175', 'E', 'ASP'], ['181', 'E', 'GLN']], 'KAG': [['34', 'E', 'SER'], ['37', 'E', 'THR'], ['96', 'E', 'GLN'], ['106', 'E', 'LEU'], ['109', 'E', 'SER']], 'KAH': [['49', 'E', 'LEU'], ['53', 'E', 'SER'], ['54', 'E', 'PHE'], ['57', 'E', 'VAL'], ['70', 'E', 'ALA'], ['72', 'E', 'LYS'], ['74', 'E', 'LEU'], ['84', 'E', 'GLN'], ['87', 'E', 'HIS'], ['88', 'E', 'THR'], ['91', 'E', 'GLU'], ['104', 'E', 'VAL'], ['120', 'E', 'MET'], ['122', 'E', 'TYR'], ['123', 'E', 'VAL'], ['127', 'E', 'GLU'], ['166', 'E', 'ASP'], ['168', 'E', 'LYS'], ['170', 'E', 'GLU'], ['171', 'E', 'ASN'], ['173', 'E', 'LEU'], ['183', 'E', 'THR'], ['184', 'E', 'ASP'], ['187', 'E', 'PHE'], ['201', 'E', 'THR'], ['327', 'E', 'PHE']], 'KAI': [['131', 'E', 'HIS'], ['138', 'E', 'PHE'], ['142', 'E', 'HIS'], ['146', 'E', 'TYR'], ['174', 'E', 'ILE'], ['314', 'E', 'PHE']], 'KAJ': [['33', 'E', 'PRO'], ['89', 'E', 'LEU'], ['92', 'E', 'LYS'], ['93', 'E', 'ARG'], ['96', 'E', 'GLN'], ['349', 'E', 'GLU'], ['350', 'E', 'PHE']], 'KAK': [['157', 'E', 'LEU'], ['162', 'E', 'LEU'], ['164', 'E', 'TYR'], ['185', 'E', 'PHE'], ['188', 'E', 'ALA']], 'KAL': [['127', 'E', 'GLU'], ['129', 'E', 'PHE'], ['130', 'E', 'SER'], ['327', 'E', 'PHE'], ['328', 'E', 'ASP'], ['330', 'E', 'TYR']], 'KAM': [['51', 'E', 'THR'], ['56', 'E', 'ARG'], ['73', 'E', 'ILE'], ['115', 'E', 'ASN'], ['335', 'E', 'ILE']], 'KAN': [['165', 'E', 'ARG'], ['166', 'E', 'ASP'], ['167', 'E', 'LEU'], ['201', 'E', 'THR'], ['204', 'E', 'TYR'], ['205', 'E', 'LEU'], ['206', 'E', 'ALA'], ['209', 'E', 'ILE'], ['219', 'E', 'VAL'], ['220', 'E', 'ASP'], ['223', 'E', 'ALA']], 'KAO': [['48', 'E', 'THR'], ['51', 'E', 'THR'], ['56', 'E', 'ARG'], ['330', 'E', 'TYR']], 'KAP': [['222', 'E', 'TRP'], ['238', 'E', 'PHE'], ['255', 'E', 'VAL'], ['273', 'E', 'LEU']], 'KAQ': [['207', 'E', 'PRO'], ['208', 'E', 'GLU'], ['211', 'E', 'LEU'], ['213', 'E', 'LYS'], ['277', 'E', 'LEU']], 'KAR': [['238', 'E', 'PHE'], ['249', 'E', 'LYS'], ['255', 'E', 'VAL'], ['256', 'E', 'ARG']]}
    """
    from pyKVFinder._pyKVFinder import _constitutional

    # Check arguments
    if type(cavities) not in [numpy.ndarray]:
        raise TypeError("`cavities` must be a numpy.ndarray.")
    elif len(cavities.shape) != 3:
        raise ValueError("`cavities` has the incorrect shape. It must be (nx, ny, nz).")
    if type(atomic) not in [numpy.ndarray, list]:
        raise TypeError("`atomic` must be a list or a numpy.ndarray.")
    elif len(numpy.asarray(atomic).shape) != 2:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    elif numpy.asarray(atomic).shape[1] != 8:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    if type(vertices) not in [numpy.ndarray, list]:
        raise TypeError("`vertices` must be a list or a numpy.ndarray.")
    elif numpy.asarray(vertices).shape != (4, 3):
        raise ValueError("`vertices` has incorrect shape. It must be (4, 3).")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")
    if type(probe_in) not in [float, int]:
        raise TypeError("`probe_in` must be a non-negative real number.")
    elif probe_in < 0.0:
        raise ValueError("`probe_in` must be a non-negative real number.")
    if type(ignore_backbone) not in [bool]:
        raise TypeError("`ignore_backbone` must be a boolean.")
    if selection is not None:
        # Check selection types
        if all(isinstance(x, int) for x in selection):
            pass
        elif all(isinstance(x, str) for x in selection):
            selection = [_get_cavity_label(sele) for sele in selection]
        else:
            raise TypeError(
                "`selection` must be a list of strings (cavity names) or integers (cavity labels)."
            )
        # Check if selection includes valid cavity labels
        if any(x < 2 for x in selection):
            raise ValueError(f"Invalid `selection`: {selection}.")
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
    if isinstance(atomic, list):
        atomic = numpy.asarray(atomic)
    if isinstance(vertices, list):
        vertices = numpy.asarray(vertices)
    if isinstance(step, int):
        step = float(step)
    if isinstance(probe_in, int):
        probe_in = float(probe_in)

    # Convert numpy.ndarray data types
    cavities = cavities.astype("int32") if cavities.dtype != "int32" else cavities
    vertices = vertices.astype("float64") if vertices.dtype != "float64" else vertices

    # Extract atominfo from atomic
    atominfo = numpy.asarray(
        ([[f"{atom[0]}_{atom[1]}_{atom[2]}", atom[3]] for atom in atomic[:, :4]])
    )

    # Extract xyzr from atomic
    xyzr = atomic[:, 4:].astype(numpy.float64)

    # Get number of cavities
    ncav = int(cavities.max() - 1)

    # Get sincos: sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb)
    sincos = _get_sincos(vertices)

    # Select cavities
    if selection is not None:
        cavities = _select_cavities(cavities, selection)

    # Unpack vertices
    P1, P2, P3, P4 = vertices

    # Remove backbone from atominfo
    if ignore_backbone:
        mask = numpy.where(
            (atominfo[:, 1] != "C")
            & (atominfo[:, 1] != "CA")
            & (atominfo[:, 1] != "N")
            & (atominfo[:, 1] != "O")
        )
        atominfo = atominfo[mask[0],]
        xyzr = xyzr[mask[0],]

    # Prepare atominfo
    atominfo = atominfo[:, 0].tolist()

    # Get interface residues
    # BUG: Segmentation Fault when running it after depth characterization
    residues = _constitutional(
        cavities, atominfo, xyzr, P1, sincos, step, probe_in, ncav, nthreads, verbose
    )
    residues = _process_residues(residues, ncav, selection)

    return residues
