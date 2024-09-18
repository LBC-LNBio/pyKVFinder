# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This is the source code for the visualization embedded in PyMOL 
pyKVFinder Tools. Changes in this file are not advised, as it controls
interactions with pyKVFinder and PyQt5.
"""

from pymol import cmd
from PyQt5.QtWidgets import QListWidget


def _show_residues(results: dict, widget: QListWidget, input_pdb: str) -> None:
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
    indexes = [item.text() for item in widget.selectedItems()]

    # Clean objects
    cmd.set("auto_zoom", 0)
    cmd.delete("res")
    cmd.delete("residues")

    # Return if no cavity is selected
    if len(indexes) < 1:
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
    cmd.set("auto_zoom", 1)


def _show_cavities(
    results: dict, widget1: QListWidget, widget2: QListWidget, cavity_pdb: str
):
    """
    Show the selected cavities in Volume and Area lists.

    Parameters
    ----------
    results : dict
        The results from the pyKVFinder calculation.
    widget1 : QListWidget
        The widget 1.
    widget2 : QListWidget
        The widget 2.
    cavity_pdb : str
        The name of the cavity pdb.
    """
    # Check if cavity file is loaded
    control = 0
    for item in cmd.get_names("all"):
        if item == cavity_pdb:
            control = 1
            break
    if control == 0:
        return

    # Get items from widget1
    indexes = [item.text()[0:3] for item in widget1.selectedItems()]

    # Select items from widget2
    for index in range(widget1.count()):
        if widget2.item(index).text()[0:3] in indexes:
            widget2.item(index).setSelected(True)
        else:
            widget2.item(index).setSelected(False)

    # Clean objects
    cmd.set("auto_zoom", 0)
    cmd.delete("cavs")
    cmd.delete("cavities")

    # Return if no cavity is selected
    if len(indexes) < 1:
        return

    # Color filling cavity points as blue nonbonded
    command = f"obj {cavity_pdb} and (resname "
    while len(indexes) > 0:
        command = f"{command}{indexes.pop(0)},"
    command = f"{command[:-1]})"
    cmd.select("cavs", command)

    # Create cavities object with blue nonbonded
    cmd.create("cavities", "cavs")
    cmd.delete("cavs")
    cmd.color("blue", "cavities")
    cmd.show("nonbonded", "cavities")

    # Color surface cavity points as red nb_spheres
    cmd.select("cavs", "cavities and name HA")
    cmd.color("red", "cavs")
    cmd.show("nb_spheres", "cavs")
    cmd.delete("cavs")

    # Reset cavities output object
    cmd.disable(cavity_pdb)
    cmd.enable(cavity_pdb)
    for item in cmd.get_names("all"):
        if item == "hydropathy":
            cmd.disable("hydropathy")
            cmd.enable("hydropathy")
        if item == "depths":
            cmd.disable("depths")
            cmd.enable("depths")
    cmd.set("auto_zoom", 1)
