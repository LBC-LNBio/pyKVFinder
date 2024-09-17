# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This is the source code for the actions performed by PyMOL pyKVFinder Tools.
Changes in this file are not advised, as it controls interactions with
pyKVFinder and PyQt5.
"""

import os

from pymol import cmd
from PyQt5.QtCore import QDir
from PyQt5.QtWidgets import QComboBox, QFileDialog, QLineEdit, QListWidget


def _load_area(results: dict, widget: QListWidget) -> None:
    """
    Load a area file into PyMOL.

    Parameters
    ----------
    results : dict
        The results from the pyKVFinder calculation.
    widget : QListWidget
        The widget to display the results.
    """
    # Get cavity indexes
    indexes = sorted(results["RESULTS"]["AREA"].keys())
    # Include Area
    for index in indexes:
        item = f"{index}: {results['RESULTS']['AREA'][index]}"
        widget.addItem(item)


def _load_avg_depth(results: dict, widget: QListWidget) -> None:
    """
    Load a avg_depth file into PyMOL.

    Parameters
    ----------
    results : dict
        The results from the pyKVFinder calculation.
    widget : QListWidget
        The widget to display the results.
    """
    # Get cavity indexes
    indexes = sorted(results["RESULTS"]["AVG_DEPTH"].keys())
    # Include Avg Depth
    for index in indexes:
        item = f"{index}: {results['RESULTS']['AVG_DEPTH'][index]}"
        widget.addItem(item)


def _load_avg_hydropathy(results: dict, widget: QListWidget) -> None:
    """
    Load a avg_hydropathy file into PyMOL.

    Parameters
    ----------
    results : dict
        The results from the pyKVFinder calculation.
    widget : QListWidget
        The widget to display the results.
    """
    # Get cavity indexes
    indexes = sorted(results["RESULTS"]["AVG_HYDROPATHY"].keys())
    # Include Avg Hydropathy
    for index in indexes:
        if not (
            (index == "EisenbergWeiss")
            or (index == "HessaHeijne")
            or (index == "KyteDoolittle")
            or (index == "RadzickaWolfenden")
            or (index == "MoonFleming")
            or (index == "WimleyWhite")
            or (index == "ZhaoLondon")
        ):
            item = f"{index}: {results['RESULTS']['AVG_HYDROPATHY'][index]}"
            widget.addItem(item)


def _load_cavity_file(filename: str, pymolname: str) -> None:
    """
    Load a cavity file into PyMOL.

    Parameters
    ----------
    filename : str
        The path to the cavity file.
    pymolname : str
        The name to assign to the object in PyMOL.
    """
    # Remove previous results in objects with same cavity name
    for obj in cmd.get_names("all"):
        if pymolname == obj:
            cmd.delete(obj)

    # Load the cavity file
    cmd.load(filename, pymolname, zoom=0)

    # Hide and show to update the display
    cmd.hide("everything", pymolname)
    cmd.show("nonbonded", pymolname)


def _load_max_depth(results: dict, widget: QListWidget) -> None:
    """
    Load a max_depth file into PyMOL.

    Parameters
    ----------
    results : dict
        The results from the pyKVFinder calculation.
    widget : QListWidget
        The widget to display the results.
    """
    # Get cavity indexes
    indexes = sorted(results["RESULTS"]["MAX_DEPTH"].keys())
    # Include Max Depth
    for index in indexes:
        item = f"{index}: {results['RESULTS']['MAX_DEPTH'][index]}"
        widget.addItem(item)


def _load_molecule_file(filename: str, pymolname: str) -> None:
    """
    Load a molecule file into PyMOL.

    Parameters
    ----------
    filename : str
        The path to the molecule file.
    pymolname : str
        The name to assign to the object in PyMOL.
    """
    # Remove previous results in objects with same cavity name
    for obj in cmd.get_names("all"):
        if pymolname == obj:
            cmd.delete(obj)

    # Load the molecule file
    cmd.load(filename, pymolname, zoom=0)


def _load_residues(results: dict, widget: QListWidget) -> None:
    """
    Load a residues file into PyMOL.

    Parameters
    ----------
    results : dict
        The results from the pyKVFinder calculation.
    widget : QListWidget
        The widget to display the results.
    """
    # Get cavity indexes
    indexes = sorted(results["RESULTS"]["RESIDUES"].keys())
    # Include Residues
    for index in indexes:
        widget.addItem(index)


def _load_volume(results: dict, widget: QListWidget) -> None:
    """
    Load a volume file into PyMOL.

    Parameters
    ----------
    results : dict
        The results from the pyKVFinder calculation.
    widget : QListWidget
        The widget to display the results.
    """
    # Get cavity indexes
    indexes = sorted(results["RESULTS"]["VOLUME"].keys())
    # Include Volume
    for index in indexes:
        item = f"{index}: {results['RESULTS']['VOLUME'][index]}"
        widget.addItem(item)


def _refresh_list(widget: QComboBox) -> None:
    """
    Refresh the list of objects in the given widget.

    Parameters
    ----------
    widget : QComboBox
        The widget to refresh.
    """
    # Clear the current list
    widget.clear()

    # Get the list of objects in the PyMOL session
    objects = cmd.get_names("all")
    for obj in objects:
        if (
            (cmd.get_type(obj) == "object:molecule")  # molecules
            and (cmd.count_atoms(obj) > 0)  # empty objects
            and (obj != "box")  # custom grid
            and (obj != "cavities")  # cavities highlight
            and (obj != "residues")  # residues highlight
            and (obj[-16:] != ".KVFinder.output")  # cavities output
            and (obj != "target_exclusive")
        ):
            widget.addItem(obj)


def _select_directory(widget: QLineEdit) -> None:
    """
    Open a dialog to select the base directory for the output files.
    Callback |for the "Browse ..." button.
    """
    # Open the dialog
    dirpath = QFileDialog.getExistingDirectory(
        caption="Select Directory", directory=os.getcwd()
    )

    # Set the selected directory
    if dirpath:
        dirpath = QDir.toNativeSeparators(dirpath)
        if os.path.isdir(dirpath):
            widget.setText(dirpath)


def _select_file(caption: str, filters: str, widget: QLineEdit) -> None:
    """
    Open a dialog to select a file. Callback for the "Browse ..." buttons.

    Parameters
    ----------
    caption : str
        The caption for the dialog.
    entry : str
        The entry to set the selected file.
    filters : str
        The filters for the dialog.
    """
    # Open the dialog
    filepath, _ = QFileDialog.getOpenFileName(
        caption=caption, filter=filters, directory=os.getcwd()
    )

    # Set the selected file
    if filepath:
        filepath = QDir.toNativeSeparators(filepath)
        if os.path.isfile(filepath):
            widget.setText(filepath)
