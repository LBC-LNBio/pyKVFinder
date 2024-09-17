# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This is the source code for the visualization embedded in PyMOL 
pyKVFinder Tools. Changes in this file are not advised, as it controls
interactions with pyKVFinder and PyQt5.
"""

import os

from pymol import cmd
from PyQt5.QtWidgets import QListWidget

def _show_residues(results: dict, widget: QListWidget, input_pdb: str, cavity_pdb: str) -> None:
    """
    Show the residues for the selected cavity.

    Parameters
    ----------
    results : dict
        The results from the pyKVFinder calculation.
    widget : QListWidget
        The widget to display the results.
    """
    # Check if input pdb is loaded
    control = 0
    for item in cmd.get_names("all"):
        if item == input_pdb:
            control = 1
            break
    if control == 0:
        return

    # Get the selected cavity
    index = [item.text() for item in widget.selectedItems()]

    # Clean objects
    cmd.set("auto_zoom", 0)
    cmd.delete("res")
    cmd.delete("residues")

    # Return if no cavity is selected
    if len(index) < 1:
        return

    # Get residues from selected cavities
    residues = []
    for index in indexes:
        for residue in results["RESULTS"]["RESIDUES"][index]:
            if residue not in residues:
                residues.append(residue)

    # Select residues
    command = ""
    while len(residues) > 0:
        res, chain, _ = residues.pop(0)
        command = f"{command} (resid {res} and chain {chain}) or"
    command = f"obj {input_pdb} and ({command[:-3]})"
    cmd.select("res", command)

    # Create residues object and show
    cmd.create("residues", "res")
    cmd.delete("res")
    cmd.hide("everything", "residues")
    cmd.show("sticks", "residues")
    cmd.disable(cavity_pdb)
    cmd.enable(cavity_pdb)
    cmd.set("auto_zoom", 1)