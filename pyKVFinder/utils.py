import argparse
import logging
import os
import pathlib
from typing import Dict, List, Optional, Union

import numpy

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
    based on the atom type (see ``van der Waals file template``).
    The package contains a built-in van der Waals radii file: ``vdw.dat``.

    See Also
    --------
    read_pdb
    read_xyz
    Molecule

    Example
    -------
    The ``read_vdw`` function takes the `built-in dictionary <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/vdw.dat>`_ when a *.dat* file is not specified.

    >>> from pyKVFinder import read_vdw
    >>> vdw = read_vdw()
    >>> vdw
    {'ALA': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB1': 1.487, '1HB': 1.487, 'HB2': 1.487, '2HB': 1.487, 'HB3': 1.487, '3HB': 1.487, 'C': 1.908, 'O': 1.6612}, 'ARG': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'HD2': 1.387, '1HD': 1.387, '2HD': 1.387, 'HD3': 1.387, 'HD1': 1.387, 'NE': 1.75, 'HE': 0.6, 'CZ': 1.908, 'NH1': 1.75, 'HH11': 0.6, '1HH1': 0.6, 'HH12': 0.6, '2HH1': 0.6, 'NH2': 1.75, 'HH21': 0.6, '2HH2': 0.6, 'HH22': 0.6, '1HH2': 0.6, 'C': 1.908, 'O': 1.6612}, 'ASH': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'OD1': 1.6612, 'OD2': 1.721, 'HD2': 0.0001, 'C': 1.908, 'O': 1.6612}, 'ASN': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'OD1': 1.6612, 'ND2': 1.824, 'HD21': 0.6, '1HD2': 0.6, 'HD22': 0.6, '2HD2': 0.6, 'C': 1.908, 'O': 1.6612}, 'ASP': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'OD1': 1.6612, 'OD2': 1.6612, 'C': 1.908, 'O': 1.6612}, 'CYM': {'N': 1.824, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB3': 1.387, 'HB2': 1.387, 'SG': 2.0, 'C': 1.908, 'O': 1.6612}, 'CYS': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, '2HB': 1.387, '1HB': 1.387, 'HB3': 1.387, 'HB1': 1.387, 'SG': 2.0, 'HG': 0.6, 'C': 1.908, 'O': 1.6612}, 'CYX': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, 'HB3': 1.387, 'SG': 2.0, 'C': 1.908, 'O': 1.6612}, 'GLH': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'HG2': 1.487, 'HG3': 1.487, 'CD': 1.908, 'OE1': 1.6612, 'OE2': 1.721, 'HE2': 0.0001, 'C': 1.908, 'O': 1.6612}, 'GLN': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'OE1': 1.6612, 'NE2': 1.824, 'HE21': 0.6, '1HE2': 0.6, 'HE22': 0.6, '2HE2': 0.6, 'C': 1.908, 'O': 1.6612}, 'GLU': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'OE1': 1.6612, 'OE2': 1.6612, 'C': 1.908, 'O': 1.6612}, 'GLY': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA2': 1.387, 'HA1': 1.387, '1HA': 1.387, '2HA': 1.387, 'HA3': 1.387, 'C': 1.908, 'O': 1.6612}, 'HID': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'HD1': 0.6, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'HIE': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'HE2': 0.6, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'HIP': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'HD1': 0.6, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'HE2': 0.6, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'ILE': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.487, 'CG2': 1.908, 'HG21': 1.487, '1HG2': 1.487, 'HG22': 1.487, '2HG2': 1.487, 'HG23': 1.487, '3HG2': 1.487, 'CG1': 1.908, 'HG12': 1.487, '2HG1': 1.487, 'HG13': 1.487, 'HG11': 1.487, '1HG1': 1.487, 'CD1': 1.908, 'HD11': 1.487, '1HD1': 1.487, 'HD12': 1.487, '2HD1': 1.487, 'HD13': 1.487, '3HD1': 1.487, 'C': 1.908, 'O': 1.6612}, 'LEU': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG': 1.487, 'CD1': 1.908, 'HD11': 1.487, '1HD1': 1.487, 'HD12': 1.487, '2HD1': 1.487, 'HD13': 1.487, '3HD1': 1.487, 'CD2': 1.908, 'HD21': 1.487, '1HD2': 1.487, 'HD22': 1.487, '2HD2': 1.487, 'HD23': 1.487, '3HD2': 1.487, 'C': 1.908, 'O': 1.6612}, 'LYN': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'HG2': 1.487, 'HG3': 1.487, 'CD': 1.908, 'HD2': 1.487, 'HD3': 1.487, 'CE': 1.908, 'HE2': 1.1, 'HE3': 1.1, 'NZ': 1.824, 'HZ2': 0.6, 'HZ3': 0.6, 'C': 1.908, 'O': 1.6612}, 'LYS': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'HD2': 1.487, '1HD': 1.487, '2HD': 1.487, 'HD3': 1.487, 'HD1': 1.487, 'CE': 1.908, 'HE2': 1.1, '2HE': 1.1, 'HE3': 1.1, '1HE': 1.1, 'HE1': 1.1, 'NZ': 1.824, 'HZ1': 0.6, '1HZ': 0.6, 'HZ2': 0.6, '2HZ': 0.6, 'HZ3': 0.6, '3HZ': 0.6, 'C': 1.908, 'O': 1.6612}, 'MET': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.387, '2HG': 1.387, 'HG3': 1.387, 'HG1': 1.387, '1HG': 1.387, 'SD': 2.0, 'CE': 1.908, 'HE1': 1.387, '1HE': 1.387, 'HE2': 1.387, '2HE': 1.387, 'HE3': 1.387, '3HE': 1.387, 'C': 1.908, 'O': 1.6612}, 'PHE': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'HZ': 1.459, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'C': 1.908, 'O': 1.6612}, 'PRO': {'N': 1.824, 'CD': 1.908, 'HD2': 1.387, '1HD': 1.387, '2HD': 1.387, 'HD3': 1.387, 'HD1': 1.387, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CA': 1.908, 'HA': 1.387, 'C': 1.908, 'O': 1.6612}, 'SER': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, '2HB': 1.387, '1HB': 1.387, 'HB3': 1.387, 'HB1': 1.387, 'OG': 1.721, 'HG': 0.0001, 'C': 1.908, 'O': 1.6612}, 'THR': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.387, 'CG2': 1.908, 'HG21': 1.487, '1HG2': 1.487, 'HG22': 1.487, '2HG2': 1.487, 'HG23': 1.487, '3HG2': 1.487, 'OG1': 1.721, 'HG1': 0.0001, 'C': 1.908, 'O': 1.6612}, 'TRP': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.85, 'CD1': 2.0, 'HD1': 1.409, 'NE1': 1.75, 'HE1': 0.6, 'CE2': 1.85, 'CZ2': 1.908, 'HZ2': 1.459, 'CH2': 1.908, 'HH2': 1.459, 'CZ3': 1.908, 'HZ3': 1.459, 'CE3': 1.908, 'HE3': 1.459, 'CD2': 1.85, 'C': 1.908, 'O': 1.6612}, 'TYR': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'OH': 1.721, 'HH': 0.0001, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'C': 1.908, 'O': 1.6612}, 'VAL': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.487, 'CG1': 1.908, 'CG2': 1.908, 'HG11': 1.487, '1HG2': 1.487, '1HG1': 1.487, 'HG21': 1.487, 'HG12': 1.487, '2HG1': 1.487, 'HG22': 1.487, '2HG2': 1.487, 'HG13': 1.487, '3HG2': 1.487, '3HG1': 1.487, 'HG23': 1.487, 'C': 1.908, 'O': 1.6612}, 'HIS': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.85, 'ND1': 1.75, 'HD1': 0.6, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'PTR': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'OH': 1.6837, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'SEP': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, 'HB3': 1.387, '1HB': 1.387, '2HB': 1.387, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'TPO': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.387, 'CG2': 1.908, 'HG21': 1.487, 'HG22': 1.487, 'HG23': 1.487, '1HG2': 1.487, '2HG2': 1.487, '3HG2': 1.487, 'OG1': 1.6837, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'H2D': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'HE2': 0.6, 'CD2': 2.0, 'HD2': 1.409, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'Y1P': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.721, 'O2P': 1.6612, 'O3P': 1.6612, 'H1P': 0.0001, 'C': 1.908, 'O': 1.6612}, 'T1P': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.387, 'CG2': 1.908, 'HG21': 1.487, 'HG22': 1.487, 'HG23': 1.487, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.721, 'O2P': 1.6612, 'O3P': 1.6612, 'H1P': 0.0001, 'C': 1.908, 'O': 1.6612}, 'S1P': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, 'HB3': 1.387, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.721, 'O2P': 1.6612, 'O3P': 1.6612, 'H1P': 0.0001, 'C': 1.908, 'O': 1.6612}, 'GEN': {'AC': 2.0, 'AG': 1.72, 'AL': 2.0, 'AM': 2.0, 'AR': 1.88, 'AS': 1.85, 'AT': 2.0, 'AU': 1.66, 'B': 2.0, 'BA': 2.0, 'BE': 2.0, 'BH': 2.0, 'BI': 2.0, 'BK': 2.0, 'BR': 1.85, 'C': 1.66, 'CA': 2.0, 'CD': 1.58, 'CE': 2.0, 'CF': 2.0, 'CL': 1.75, 'CM': 2.0, 'CO': 2.0, 'CR': 2.0, 'CS': 2.0, 'CU': 1.4, 'DB': 2.0, 'DS': 2.0, 'DY': 2.0, 'ER': 2.0, 'ES': 2.0, 'EU': 2.0, 'F': 1.47, 'FE': 2.0, 'FM': 2.0, 'FR': 2.0, 'GA': 1.87, 'GD': 2.0, 'GE': 2.0, 'H': 0.91, 'HE': 1.4, 'HF': 2.0, 'HG': 1.55, 'HO': 2.0, 'HS': 2.0, 'I': 1.98, 'IN': 1.93, 'IR': 2.0, 'K': 2.75, 'KR': 2.02, 'LA': 2.0, 'LI': 1.82, 'LR': 2.0, 'LU': 2.0, 'MD': 2.0, 'MG': 1.73, 'MN': 2.0, 'MO': 2.0, 'MT': 2.0, 'N': 1.97, 'NA': 2.27, 'NB': 2.0, 'ND': 2.0, 'NE': 1.54, 'NI': 1.63, 'NO': 2.0, 'NP': 2.0, 'O': 1.69, 'OS': 2.0, 'P': 2.1, 'PA': 2.0, 'PB': 2.02, 'PD': 1.63, 'PM': 2.0, 'PO': 2.0, 'PR': 2.0, 'PT': 1.72, 'PU': 2.0, 'RA': 2.0, 'RB': 2.0, 'RE': 2.0, 'RF': 2.0, 'RH': 2.0, 'RN': 2.0, 'RU': 2.0, 'S': 2.09, 'SB': 2.0, 'SC': 2.0, 'SE': 1.9, 'SG': 2.0, 'SI': 2.1, 'SM': 2.0, 'SN': 2.17, 'SR': 2.0, 'TA': 2.0, 'TB': 2.0, 'TC': 2.0, 'TE': 2.06, 'TH': 2.0, 'TI': 2.0, 'TL': 1.96, 'TM': 2.0, 'U': 1.86, 'V': 2.0, 'W': 2.0, 'XE': 2.16, 'Y': 2.0, 'YB': 2.0, 'ZN': 1.39, 'ZR': 2.0}}

    The van der Waals radii can be define by:

        * creating a Python dictionary:

        >>> vdw = {'GEN': {'C': 1.66, 'CA': 2.0, 'N': 1.97, 'O': 1.69, 'H': 0.91}}

        * specifying a *.dat* file following template of `van der Waals radii file`.

        >>> with open('vdw.dat', 'w') as f:
        ...     f.write('>GEN\\nC\\t\\t1.66\\nCA\\t\\t2.00\\nN\\t\\t1.97\\nO\\t\\t1.69\\nH\\t\\t0.91\\n')
        >>> vdw = read_vdw('vdw.dat')
        >>> vdw
        {'GEN': {'C': 1.66, 'CA': 2.0, 'N': 1.97, 'O': 1.69, 'H': 0.91}}
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
    fn: Union[str, pathlib.Path],
    vdw: Optional[Dict[str, Dict[str, float]]] = None,
    model: Optional[int] = None,
) -> numpy.ndarray:
    """Reads PDB file into numpy.ndarrays.

    Parameters
    ----------
    fn : Union[str, pathlib.Path]
        A path to PDB file.
    vdw : Dict[str, Dict[str, float]], optional
        A dictionary containing radii values, by default None. If None, use output of ``read_vdw()``.
    model : int, optional
        The model number of a multi-model PDB file, by default None. If None, keep atoms from all models.

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
    The van der Waals radii file defines the radius values for each atom by residue and when not defined, it uses a generic value based on the atom type. The function by default loads the built-in van der Waals radii file: `vdw.dat`.

    See Also
    --------
    read_vdw
    read_xyz
    get_vertices
    get_vertices_from_file
    detect
    constitutional
    hydropathy

    Example
    -------
    With the vdW radii dictionary loaded with ``read_vdw``, we can read a target PDB file into Numpy array (atomic data):

    >>> import os
    >>> import pyKVFinder
    >>> from pyKVFinder import read_pdb
    >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
    >>> atomic = read_pdb(pdb)
    >>> atomic
    array([['13', 'E', 'GLU', ..., '-15.642', '-14.858', '1.824'],
       ['13', 'E', 'GLU', ..., '-14.62', '-15.897', '1.908'],
       ['13', 'E', 'GLU', ..., '-13.357', '-15.508', '1.908'],
       ...,
       ['350', 'E', 'PHE', ..., '18.878', '-9.885', '1.908'],
       ['350', 'E', 'PHE', ..., '17.624', '-9.558', '1.908'],
       ['350', 'E', 'PHE', ..., '19.234', '-13.442', '1.69']],
      dtype='<U32')

    .. warning::
        The function takes the `built-in dictionary <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/vdw.dat>`_ when the ``vdw`` argument is not specified. If you wish to use a custom van der Waals radii file, you must read it with ``read_vdw`` as shown earlier and pass it as ``read_pdb(pdb, vdw=vdw)``.
    """
    # Check arguments
    if type(fn) not in [str, pathlib.Path]:
        raise TypeError("`fn` must be a string or a pathlib.Path.")
    if model is not None:
        if type(model) not in [int]:
            raise TypeError("`model` must be an integer.")

    # Define default vdw file
    if vdw is None:
        vdw = read_vdw(VDW)

    # Create lists
    atomic = []

    # Keep all models
    keep = True if model is None else False

    # Read file and process atoms
    with open(fn, "r") as f:
        for line in f.readlines():
            if model is not None:
                if line[:5] == "MODEL":
                    nmodel = int(line[5:].replace(" ", "").rstrip("\n"))
                    keep = True if model == nmodel else False
            if keep:
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
        A dictionary containing radii values, by default None. If None, use output of ``read_vdw()``.

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

    See Also
    --------
    read_vdw
    read_pdb
    get_vertices
    get_vertices_from_file
    detect
    constitutional
    hydropathy

    Example
    -------
    With the vdW radii dictionary loaded with ``pyKVFinder.read_vdw``, we can read a target XYZ file into Numpy arrays (atomic information and atomic coordinates):

    >>> import os
    >>> import pyKVFinder
    >>> from pyKVFinder import read_xyz
    >>> xyz = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.xyz')
    >>> atomic = read_xyz(xyz)
    >>> atominfo
    array([['1', 'A', 'UNK', ..., '-15.642', '-14.858', '1.97'],
       ['2', 'A', 'UNK', ..., '-14.62', '-15.897', '1.66'],
       ['3', 'A', 'UNK', ..., '-13.357', '-15.508', '1.66'],
       ...,
       ['2790', 'A', 'UNK', ..., '18.878', '-9.885', '1.66'],
       ['2791', 'A', 'UNK', ..., '17.624001', '-9.558', '1.66'],
       ['2792', 'A', 'UNK', ..., '19.233999', '-13.442', '1.69']],
      dtype='<U32')

    .. warning::
        The function takes the `built-in dictionary <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/vdw.dat>`_ when the ``vdw`` argument is not specified. If you wish to use a custom van der Waals radii file, you must read it with ``read_vdw`` as shown earlier and pass it as ``read_xyz(xyz, vdw=vdw)``.

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
        A dictionary containing radii values, by default None. If None, use output of ``read_vdw()``.
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

    Note
    ----
    The function takes the `built-in dictionary <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/vdw.dat>`_ when the ``vdw`` argument is not specified. If you wish to use a custom van der Waals radii file, you must read it with ``read_vdw`` as shown earlier and pass it as ``read_cavity(cavity, receptor, vdw=vdw)``.

    See Also
    --------
    read_pdb
    read_xyz
    get_vertices
    get_vertices_from_file
    spatial
    depth
    constitutional
    hydropathy
    export

    Example
    -------
    With a previously calculated cavity, that can be manually curated in a molecular visualization software, such as PyMOL, we can read it with its respective receptor back to pyKVFinder:

    >>> import os
    >>> import pyKVFinder
    >>> from pyKVFinder import read_cavity
    >>> cavity = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.KVFinder.output.pdb')
    >>> receptor = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
    >>> grid = read_cavity(cavity, receptor)
    >>> grid
    array([[[-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        ...,
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1]],
       ...,
       [[-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        ...,
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1],
        [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)
    """
    from _pyKVFinder import _fill_cavity, _fill_receptor

    from .grid import _get_dimensions, _get_sincos, get_vertices

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
            print("> Surface representation: Solvent Excluded Surface (SES)")
        surface = True
    elif surface == "SAS":
        if verbose:
            print("> Surface representation: Solvent Accessible Surface (SAS)")
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

    See Also
    --------
    constitutional
    plot_frequencies
    write_results

    Example
    -------
    With the interface residues identified with ``constitutional``, we can calculate residues and classes of residues frequencies:

    >>> from pyKVFinder import calculate_frequencies
    >>> residues
    {'KAA': [['49', 'E', 'LEU'], ['50', 'E', 'GLY'], ['51', 'E', 'THR'], ['52', 'E', 'GLY'], ['53', 'E', 'SER'], ['54', 'E', 'PHE'], ['55', 'E', 'GLY'], ['56', 'E', 'ARG'], ['57', 'E', 'VAL'], ['70', 'E', 'ALA'], ['72', 'E', 'LYS'], ['74', 'E', 'LEU'], ['84', 'E', 'GLN'], ['87', 'E', 'HIS'], ['88', 'E', 'THR'], ['91', 'E', 'GLU'], ['104', 'E', 'VAL'], ['120', 'E', 'MET'], ['121', 'E', 'GLU'], ['122', 'E', 'TYR'], ['123', 'E', 'VAL'], ['127', 'E', 'GLU'], ['166', 'E', 'ASP'], ['168', 'E', 'LYS'], ['170', 'E', 'GLU'], ['171', 'E', 'ASN'], ['173', 'E', 'LEU'], ['183', 'E', 'THR'], ['184', 'E', 'ASP'], ['186', 'E', 'GLY'], ['187', 'E', 'PHE'], ['201', 'E', 'THR'], ['327', 'E', 'PHE']]}
    >>> frequencies = calculate_frequencies(residues)
    >>> frequencies
    {'KAA': {'RESIDUES': {'ALA': 1, 'ARG': 1, 'ASN': 1, 'ASP': 2, 'GLN': 1, 'GLU': 4, 'GLY': 4, 'HIS': 1, 'LEU': 3, 'LYS': 2, 'MET': 1, 'PHE': 3, 'SER': 1, 'THR': 4, 'TYR': 1, 'VAL': 3}, 'CLASS': {'R1': 11, 'R2': 4, 'R3': 8, 'R4': 6, 'R5': 4, 'RX': 0}}}
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

    See Also
    --------
    calculate_frequencies

    Example
    -------
    With the residues and classes of residues frequencies calculated with ``calculate_frequencies``, we can plot the bar charts of these frequencies in a PDF file.

    >>> from pyKVFinder import plot_frequencies
    >>> frequencies
    {'KAA': {'RESIDUES': {'ALA': 1, 'ARG': 1, 'ASN': 1, 'ASP': 2, 'GLN': 1, 'GLU': 4, 'GLY': 4, 'HIS': 1, 'LEU': 3, 'LYS': 2, 'MET': 1, 'PHE': 3, 'SER': 1, 'THR': 4, 'TYR': 1, 'VAL': 3}, 'CLASS': {'R1': 11, 'R2': 4, 'R3': 8, 'R4': 6, 'R5': 4, 'RX': 0}}}
    >>> plot_frequencies(frequencies, fn='barplots.pdf')

    .. image:: ../_images/barplots.png
        :width: 600
        :align: center
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

    See Also
    --------
    detect
    spatial
    depth
    constitutional
    hydropathy
    export

    Example
    -------
    With the cavity and surface points identified and depth and hydrophobicity scale mapped in the 3D grid, we can:

    * Write cavity detection:

    >>> from pyKVFinder import write_results
    >>> import os
    >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
    >>> write_results('results.toml', input=pdb, ligand=None, output='cavities.pdb')

    * Write spatial characterization

    >>> volume
    {'KAA': 137.16, 'KAB': 47.52, 'KAC': 66.96, 'KAD': 8.21, 'KAE': 43.63, 'KAF': 12.53, 'KAG': 6.26, 'KAH': 520.13, 'KAI': 12.31, 'KAJ': 26.57, 'KAK': 12.31, 'KAL': 33.91, 'KAM': 23.11, 'KAN': 102.82, 'KAO': 6.05, 'KAP': 15.55, 'KAQ': 7.99, 'KAR': 7.78}
    >>> area
    {'KAA': 120.52, 'KAB': 58.76, 'KAC': 72.06, 'KAD': 17.62, 'KAE': 56.44, 'KAF': 22.53, 'KAG': 15.38, 'KAH': 489.25, 'KAI': 29.87, 'KAJ': 44.85, 'KAK': 30.58, 'KAL': 43.59, 'KAM': 45.25, 'KAN': 129.77, 'KAO': 11.57, 'KAP': 24.8, 'KAQ': 12.59, 'KAR': 15.97}
    >>> write_results('results.toml', input=pdb, ligand=None, output=None, volume=volume, area=area)

    * Write constitutional characterization

    >>> residues
    {'KAA': [['14', 'E', 'SER'], ['15', 'E', 'VAL'], ['18', 'E', 'PHE'], ['19', 'E', 'LEU'], ['100', 'E', 'PHE'], ['152', 'E', 'LEU'], ['155', 'E', 'GLU'], ['156', 'E', 'TYR'], ['292', 'E', 'LYS'], ['302', 'E', 'TRP'], ['303', 'E', 'ILE'], ['306', 'E', 'TYR']], 'KAB': [['18', 'E', 'PHE'], ['22', 'E', 'ALA'], ['25', 'E', 'ASP'], ['26', 'E', 'PHE'], ['29', 'E', 'LYS'], ['97', 'E', 'ALA'], ['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['156', 'E', 'TYR']], 'KAC': [['141', 'E', 'PRO'], ['142', 'E', 'HIS'], ['144', 'E', 'ARG'], ['145', 'E', 'PHE'], ['148', 'E', 'ALA'], ['299', 'E', 'THR'], ['300', 'E', 'THR'], ['305', 'E', 'ILE'], ['310', 'E', 'VAL'], ['311', 'E', 'GLU'], ['313', 'E', 'PRO']], 'KAD': [['122', 'E', 'TYR'], ['124', 'E', 'ALA'], ['176', 'E', 'GLN'], ['318', 'E', 'PHE'], ['320', 'E', 'GLY'], ['321', 'E', 'PRO'], ['322', 'E', 'GLY'], ['323', 'E', 'ASP']], 'KAE': [['95', 'E', 'LEU'], ['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['100', 'E', 'PHE'], ['103', 'E', 'LEU'], ['104', 'E', 'VAL'], ['105', 'E', 'LYS'], ['106', 'E', 'LEU']], 'KAF': [['123', 'E', 'VAL'], ['124', 'E', 'ALA'], ['175', 'E', 'ASP'], ['176', 'E', 'GLN'], ['181', 'E', 'GLN']], 'KAG': [['34', 'E', 'SER'], ['37', 'E', 'THR'], ['96', 'E', 'GLN'], ['106', 'E', 'LEU'], ['107', 'E', 'GLU'], ['108', 'E', 'PHE'], ['109', 'E', 'SER']], 'KAH': [['49', 'E', 'LEU'], ['50', 'E', 'GLY'], ['51', 'E', 'THR'], ['52', 'E', 'GLY'], ['53', 'E', 'SER'], ['54', 'E', 'PHE'], ['55', 'E', 'GLY'], ['56', 'E', 'ARG'], ['57', 'E', 'VAL'], ['70', 'E', 'ALA'], ['72', 'E', 'LYS'], ['74', 'E', 'LEU'], ['84', 'E', 'GLN'], ['87', 'E', 'HIS'], ['88', 'E', 'THR'], ['91', 'E', 'GLU'], ['104', 'E', 'VAL'], ['120', 'E', 'MET'], ['121', 'E', 'GLU'], ['122', 'E', 'TYR'], ['123', 'E', 'VAL'], ['127', 'E', 'GLU'], ['166', 'E', 'ASP'], ['168', 'E', 'LYS'], ['170', 'E', 'GLU'], ['171', 'E', 'ASN'], ['173', 'E', 'LEU'], ['183', 'E', 'THR'], ['184', 'E', 'ASP'], ['186', 'E', 'GLY'], ['187', 'E', 'PHE'], ['201', 'E', 'THR'], ['327', 'E', 'PHE']], 'KAI': [['131', 'E', 'HIS'], ['138', 'E', 'PHE'], ['142', 'E', 'HIS'], ['146', 'E', 'TYR'], ['174', 'E', 'ILE'], ['314', 'E', 'PHE']], 'KAJ': [['33', 'E', 'PRO'], ['89', 'E', 'LEU'], ['92', 'E', 'LYS'], ['93', 'E', 'ARG'], ['96', 'E', 'GLN'], ['349', 'E', 'GLU'], ['350', 'E', 'PHE']], 'KAK': [['157', 'E', 'LEU'], ['162', 'E', 'LEU'], ['163', 'E', 'ILE'], ['164', 'E', 'TYR'], ['185', 'E', 'PHE'], ['188', 'E', 'ALA']], 'KAL': [['49', 'E', 'LEU'], ['127', 'E', 'GLU'], ['129', 'E', 'PHE'], ['130', 'E', 'SER'], ['326', 'E', 'ASN'], ['327', 'E', 'PHE'], ['328', 'E', 'ASP'], ['330', 'E', 'TYR']], 'KAM': [['51', 'E', 'THR'], ['55', 'E', 'GLY'], ['56', 'E', 'ARG'], ['73', 'E', 'ILE'], ['74', 'E', 'LEU'], ['75', 'E', 'ASP'], ['115', 'E', 'ASN'], ['335', 'E', 'ILE'], ['336', 'E', 'ARG']], 'KAN': [['165', 'E', 'ARG'], ['166', 'E', 'ASP'], ['167', 'E', 'LEU'], ['199', 'E', 'CYS'], ['200', 'E', 'GLY'], ['201', 'E', 'THR'], ['204', 'E', 'TYR'], ['205', 'E', 'LEU'], ['206', 'E', 'ALA'], ['209', 'E', 'ILE'], ['219', 'E', 'VAL'], ['220', 'E', 'ASP'], ['223', 'E', 'ALA']], 'KAO': [['48', 'E', 'THR'], ['51', 'E', 'THR'], ['56', 'E', 'ARG'], ['330', 'E', 'TYR'], ['331', 'E', 'GLU']], 'KAP': [['222', 'E', 'TRP'], ['238', 'E', 'PHE'], ['253', 'E', 'GLY'], ['254', 'E', 'LYS'], ['255', 'E', 'VAL'], ['273', 'E', 'LEU']], 'KAQ': [['207', 'E', 'PRO'], ['208', 'E', 'GLU'], ['211', 'E', 'LEU'], ['213', 'E', 'LYS'], ['275', 'E', 'VAL'], ['277', 'E', 'LEU']], 'KAR': [['237', 'E', 'PRO'], ['238', 'E', 'PHE'], ['249', 'E', 'LYS'], ['254', 'E', 'LYS'], ['255', 'E', 'VAL'], ['256', 'E', 'ARG']]}
    >>> frequencies
    {'KAA': {'RESIDUES': {'GLU': 1, 'ILE': 1, 'LEU': 2, 'LYS': 1, 'PHE': 2, 'SER': 1, 'TRP': 1, 'TYR': 2, 'VAL': 1}, 'CLASS': {'R1': 4, 'R2': 5, 'R3': 1, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAB': {'RESIDUES': {'ALA': 2, 'ASN': 1, 'ASP': 1, 'LYS': 1, 'PHE': 2, 'TYR': 1, 'VAL': 1}, 'CLASS': {'R1': 3, 'R2': 3, 'R3': 1, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAC': {'RESIDUES': {'ALA': 1, 'ARG': 1, 'GLU': 1, 'HIS': 1, 'ILE': 1, 'PHE': 1, 'PRO': 2, 'THR': 2, 'VAL': 1}, 'CLASS': {'R1': 5, 'R2': 1, 'R3': 2, 'R4': 1, 'R5': 2, 'RX': 0}}, 'KAD': {'RESIDUES': {'ALA': 1, 'ASP': 1, 'GLN': 1, 'GLY': 2, 'PHE': 1, 'PRO': 1, 'TYR': 1}, 'CLASS': {'R1': 4, 'R2': 2, 'R3': 1, 'R4': 1, 'R5': 0, 'RX': 0}}, 'KAE': {'RESIDUES': {'ASN': 1, 'LEU': 3, 'LYS': 1, 'PHE': 1, 'VAL': 2}, 'CLASS': {'R1': 5, 'R2': 1, 'R3': 1, 'R4': 0, 'R5': 1, 'RX': 0}}, 'KAF': {'RESIDUES': {'ALA': 1, 'ASP': 1, 'GLN': 2, 'VAL': 1}, 'CLASS': {'R1': 2, 'R2': 0, 'R3': 2, 'R4': 1, 'R5': 0, 'RX': 0}}, 'KAG': {'RESIDUES': {'GLN': 1, 'GLU': 1, 'LEU': 1, 'PHE': 1, 'SER': 2, 'THR': 1}, 'CLASS': {'R1': 1, 'R2': 1, 'R3': 4, 'R4': 1, 'R5': 0, 'RX': 0}}, 'KAH': {'RESIDUES': {'ALA': 1, 'ARG': 1, 'ASN': 1, 'ASP': 2, 'GLN': 1, 'GLU': 4, 'GLY': 4, 'HIS': 1, 'LEU': 3, 'LYS': 2, 'MET': 1, 'PHE': 3, 'SER': 1, 'THR': 4, 'TYR': 1, 'VAL': 3}, 'CLASS': {'R1': 11, 'R2': 4, 'R3': 8, 'R4': 6, 'R5': 4, 'RX': 0}}, 'KAI': {'RESIDUES': {'HIS': 2, 'ILE': 1, 'PHE': 2, 'TYR': 1}, 'CLASS': {'R1': 1, 'R2': 3, 'R3': 0, 'R4': 0, 'R5': 2, 'RX': 0}}, 'KAJ': {'RESIDUES': {'ARG': 1, 'GLN': 1, 'GLU': 1, 'LEU': 1, 'LYS': 1, 'PHE': 1, 'PRO': 1}, 'CLASS': {'R1': 2, 'R2': 1, 'R3': 1, 'R4': 1, 'R5': 2, 'RX': 0}}, 'KAK': {'RESIDUES': {'ALA': 1, 'ILE': 1, 'LEU': 2, 'PHE': 1, 'TYR': 1}, 'CLASS': {'R1': 4, 'R2': 2, 'R3': 0, 'R4': 0, 'R5': 0, 'RX': 0}}, 'KAL': {'RESIDUES': {'ASN': 1, 'ASP': 1, 'GLU': 1, 'LEU': 1, 'PHE': 2, 'SER': 1, 'TYR': 1}, 'CLASS': {'R1': 1, 'R2': 3, 'R3': 2, 'R4': 2, 'R5': 0, 'RX': 0}}, 'KAM': {'RESIDUES': {'ARG': 2, 'ASN': 1, 'ASP': 1, 'GLY': 1, 'ILE': 2, 'LEU': 1, 'THR': 1}, 'CLASS': {'R1': 4, 'R2': 0, 'R3': 2, 'R4': 1, 'R5': 2, 'RX': 0}}, 'KAN': {'RESIDUES': {'ALA': 2, 'ARG': 1, 'ASP': 2, 'CYS': 1, 'GLY': 1, 'ILE': 1, 'LEU': 2, 'THR': 1, 'TYR': 1, 'VAL': 1}, 'CLASS': {'R1': 7, 'R2': 1, 'R3': 2, 'R4': 2, 'R5': 1, 'RX': 0}}, 'KAO': {'RESIDUES': {'ARG': 1, 'GLU': 1, 'THR': 2, 'TYR': 1}, 'CLASS': {'R1': 0, 'R2': 1, 'R3': 2, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAP': {'RESIDUES': {'GLY': 1, 'LEU': 1, 'LYS': 1, 'PHE': 1, 'TRP': 1, 'VAL': 1}, 'CLASS': {'R1': 3, 'R2': 2, 'R3': 0, 'R4': 0, 'R5': 1, 'RX': 0}}, 'KAQ': {'RESIDUES': {'GLU': 1, 'LEU': 2, 'LYS': 1, 'PRO': 1, 'VAL': 1}, 'CLASS': {'R1': 4, 'R2': 0, 'R3': 0, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAR': {'RESIDUES': {'ARG': 1, 'LYS': 2, 'PHE': 1, 'PRO': 1, 'VAL': 1}, 'CLASS': {'R1': 2, 'R2': 1, 'R3': 0, 'R4': 0, 'R5': 3, 'RX': 0}}}
    >>> write_results('results.toml', input=pdb, ligand=None, output=output_cavity, volume=volume, area=area, residues=residues, frequencies=frequencies)

    * Write depth characterization

    >>> max_depth
    {'KAA': 3.79, 'KAB': 2.68, 'KAC': 2.62, 'KAD': 0.85, 'KAE': 3.0, 'KAF': 0.85, 'KAG': 0.6, 'KAH': 10.73, 'KAI': 2.55, 'KAJ': 2.24, 'KAK': 0.0, 'KAL': 3.0, 'KAM': 1.2, 'KAN': 0.0, 'KAO': 1.04, 'KAP': 2.08, 'KAQ': 0.85, 'KAR': 0.6}
    >>> avg_depth
    {'KAA': 1.28, 'KAB': 0.86, 'KAC': 0.67, 'KAD': 0.29, 'KAE': 0.98, 'KAF': 0.24, 'KAG': 0.1, 'KAH': 3.75, 'KAI': 1.5, 'KAJ': 0.96, 'KAK': 0.0, 'KAL': 1.0, 'KAM': 0.24, 'KAN': 0.0, 'KAO': 0.29, 'KAP': 0.7, 'KAQ': 0.22, 'KAR': 0.12}
    >>> write_results('results.toml', input=pdb, ligand=None, output=None, max_depth=max_depth, avg_depth=avg_depth)

    * Write hydropathy characterization

    >>> avg_hydropathy
    {'KAA': -0.73, 'KAB': -0.06, 'KAC': -0.07, 'KAD': -0.62, 'KAE': -0.81, 'KAF': -0.14, 'KAG': -0.33, 'KAH': -0.16, 'KAI': -0.4, 'KAJ': 0.62, 'KAK': -0.99, 'KAL': 0.35, 'KAM': -0.33, 'KAN': 0.18, 'KAO': 0.88, 'KAP': -0.96, 'KAQ': 0.48, 'KAR': 0.24, 'EisenbergWeiss': [-1.42, 2.6]}
    >>> write_results('results.toml', input=pdb, ligand=None, output=None, avg_hydropathy=avg_hydropathy)

    * Write all

    >>> write_results('results.toml', input=pdb, ligand=None, output='cavities.pdb', volume=volume, area=area, max_depth=max_depth, avg_depth=avg_depth, avg_hydropathy=avg_hydropathy, residues=residues, frequencies=frequencies)

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

    # Create output dictionary for results file
    files = {
        "INPUT": input,
        "LIGAND": ligand,
        "OUTPUT": output,
    }
    parameters = {
        "STEP": step,
    }
    results = (
        {
            "VOLUME": volume,
            "AREA": area,
            "MAX_DEPTH": max_depth,
            "AVG_DEPTH": avg_depth,
            "AVG_HYDROPATHY": avg_hydropathy,
            "RESIDUES": residues,
            "FREQUENCY": frequencies,
        }
        if (volume is not None)
        or (area is not None)
        or (max_depth is not None)
        or (avg_depth is not None)
        or (avg_hydropathy is not None)
        or (residues is not None)
        or (frequencies is not None)
        else None
    )

    # Create base directories of results TOML file
    os.makedirs(os.path.abspath(os.path.dirname(fn)), exist_ok=True)

    # Write results to TOML file
    with open(fn, "w") as f:
        f.write("# pyKVFinder results\n\n")
        toml.dump(
            {
                "FILES": files,
                "PARAMETERS": parameters,
                "RESULTS": results,
            },
            f,
        )
