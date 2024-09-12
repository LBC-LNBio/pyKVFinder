#! /usr/bin/env python3
# -*- coding: utf-8 -*-

#####################################################################################
#    This is the PyMOL2 pyKVFinder Tools for PyMOL v2.x. It was developed using    #
#    Qt interface and Python. Changes in this file are not advised, as it controls  #
#    all interactions with pyKVFinder.                                             #
#                                                                                   #
#    PyMOL KVFinder Web Tools is free software: you can redistribute it and/or      #
#    modify it under the terms of the GNU General Public License as published       #
#    by the Free Software Foundation, either version 3 of the License, or           #
#    (at your option) any later version.                                            #
#                                                                                   #
#    PyMOL2 pyKVFinder Tools is distributed in the hope that it will be useful,    #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                 #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  #
#    GNU General Public License for more details.                                   #
#                                                                                   #
#    You should have received a copy of the GNU General Public License  along with  #
#    PyMOL2 pyKVFinder Tools.  If not, see <http://www.gnu.org/licenses/>.         #
#                                                                                   #
#####################################################################################


from __future__ import absolute_import, print_function, annotations

import os, sys
from typing import Optional, Any, Dict
from PyQt5.QtWidgets import QMainWindow, QDialog
from PyQt5.QtCore import QThread, pyqtSlot, pyqtSignal


# global reference to avoid garbage collection of our dialog
dialog = None


########## Relevant information ##########
# pyKVFinder executable (Ubuntu/macOS)  #
executable = "pyKVFinder"  #
# pyKVFinder executable (Windows)       #
# executable = 'pyKVFinder-win64.exe'   #
##########################################


class _Default(object):
    def __init__(self):
        super(_Default, self).__init__()
        #######################
        ### Main Parameters ###
        #######################
        self.step = 0.0
        self.resolution = "Low"
        self.probe_in = 1.4
        self.probe_out = 4.0
        self.removal_distance = 2.4
        self.volume_cutoff = 5.0
        self.surface = "Molecular Surface (VdW)"
        self.cavity_representation = "Filtered"
        self.base_name = "output"
        self.output_dir_path = os.getcwd()
        #######################
        ### File Locations  ###
        #######################
        self.pyKVFinder = None
        self.dictionary = None
        #######################
        ###  Search Space   ###
        #######################
        # Box Adjustment
        self.box_adjustment = False
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.min_x = 0.0
        self.max_x = 0.0
        self.min_y = 0.0
        self.max_y = 0.0
        self.min_z = 0.0
        self.max_z = 0.0
        self.angle1 = 0
        self.angle2 = 0
        self.padding = 3.5
        # Ligand Adjustment
        self.ligand_adjustment = False
        self.ligand_cutoff = 5.0


def __init_plugin__(app=None):
    """
    Add an entry to the PyMOL v2.x "Plugin" menu
    """
    from pymol.plugins import addmenuitemqt

    addmenuitemqt("PyMOL pyKVFinder Tools", run_plugin_gui)


def run_plugin_gui():
    """
    Open our custom dialog
    """
    global dialog

    if dialog is None:
        dialog = PyMOL2KVFinderwebTools()

    dialog.show()


class PyMOL2KVFinderwebTools(QMainWindow):
    """
    PyMOL KVFinder Web Tools

    - Create pyKVFinder client Graphical User Interface (GUI)
    with PyQt5 in PyMOL v2.x viewer
    - Define functions and callbacks for GUI
    """

    def __init__(self, server="http://localhost", port="8081"):
        super(PyMOL2KVFinderwebTools, self).__init__()
        from PyQt5.QtNetwork import QNetworkAccessManager

        # Get KVFinder_PATH
        KVFinder_PATH = self.get_KVFinder_PATH()

        # Define Default Parameters
        self._default = _Default()
        self._default.pyKVFinder = os.path.join(KVFinder_PATH, executable)
        self._default.dictionary = os.path.join(KVFinder_PATH, "dictionary")

        # Initialize PyMOL2KVFinderwebTools GUI
        self.initialize_gui()

        # Restore Default Parameters
        self.restore(is_startup=True)

        # Set box centers
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0

        # Results
        self.results = None
        self.input_pdb = None
        self.ligand_pdb = None
        self.cavity_pdb = None

    def initialize_gui(self) -> None:
        """
        Qt elements are located in self
        """
        # pymol.Qt provides the PyQt5 interface
        from PyQt5 import QtWidgets
        from PyQt5.uic import loadUi

        # from pymol.Qt.utils import loadUi

        # populate the QMainWindow from our *.ui file
        uifile = os.path.join(os.path.dirname(__file__), "PyMOL-pyKVFinder-Tools.ui")
        loadUi(uifile, self)

        # ScrollBars binded to QListWidgets in Descriptors
        scroll_bar_volume = QtWidgets.QScrollBar(self)
        self.volume_list.setVerticalScrollBar(scroll_bar_volume)
        scroll_bar_area = QtWidgets.QScrollBar(self)
        self.area_list.setVerticalScrollBar(scroll_bar_area)
        scroll_bar_residues = QtWidgets.QScrollBar(self)
        self.residues_list.setVerticalScrollBar(scroll_bar_residues)

        ########################
        ### Buttons Callback ###
        ########################

        # hook up QMainWindow buttons callbacks
        self.button_run.clicked.connect(self.run)
        self.button_exit.clicked.connect(self.close)
        self.button_restore.clicked.connect(self.restore)
        self.button_grid.clicked.connect(self.show_grid)
        self.button_save_parameters.clicked.connect(self.save_parameters)

        # hook up Browse buttons callbacks
        self.button_browse.clicked.connect(self.select_directory)
        self.button_browse2.clicked.connect(
            lambda: self.select_file(
                "Choose pyKVFinder executable", self.pyKVFinder, "*"
            )
        )
        self.button_browse3.clicked.connect(
            lambda: self.select_file(
                "Choose van der Waals radii dictionary", self.dictionary, "*"
            )
        )
        self.button_browse4.clicked.connect(
            lambda: self.select_file(
                "Choose KVFinder Results File",
                self.results_file_entry,
                "KVFinder Results File (*.toml);;All files (*)",
            )
        )

        # hook up Refresh buttons callback
        self.refresh_input.clicked.connect(lambda: self.refresh(self.input))

        # hook up resolution-step CheckBox callbacks
        self.resolution_label.clicked.connect(self.check_resolution)
        self.step_size_label.clicked.connect(self.check_step_size)

        # hook up Search Space button callbacks
        # Box Adjustment
        self.button_draw_box.clicked.connect(self.set_box)
        self.button_delete_box.clicked.connect(self.delete_box)
        self.button_redraw_box.clicked.connect(self.redraw_box)
        self.button_box_adjustment_help.clicked.connect(self.box_adjustment_help)
        # Ligand Adjustment
        self.refresh_ligand.clicked.connect(lambda: self.refresh(self.ligand))

        # hook up methods to results tab
        self.button_load_results.clicked.connect(self.load_results)
        self.volume_list.itemSelectionChanged.connect(
            lambda list1=self.volume_list, list2=self.area_list: self.show_cavities(
                list1, list2
            )
        )
        self.area_list.itemSelectionChanged.connect(
            lambda list1=self.area_list, list2=self.volume_list: self.show_cavities(
                list1, list2
            )
        )
        self.avg_depth_list.itemSelectionChanged.connect(
            lambda list1=self.avg_depth_list, list2=self.max_depth_list: self.show_depth(
                list1, list2
            )
        )
        self.max_depth_list.itemSelectionChanged.connect(
            lambda list1=self.max_depth_list, list2=self.avg_depth_list: self.show_depth(
                list1, list2
            )
        )
        self.avg_hydropathy_list.itemSelectionChanged.connect(
            lambda list1=self.avg_hydropathy_list: self.show_hydropathy(list1)
        )
        self.residues_list.itemSelectionChanged.connect(self.show_residues)
        self.default_view.toggled.connect(self.show_default_view)
        self.depth_view.toggled.connect(self.show_depth_view)
        self.hydropathy_view.toggled.connect(self.show_hydropathy_view)

    def check_resolution(self):
        if self.resolution_label.isChecked():
            self.resolution.setEnabled(True)
            self.resolution.setCurrentText(self._default.resolution)
            self.step_size_label.setChecked(False)
            self.step_size.setEnabled(False)
            self.step_size.setValue(self._default.step)
        else:
            self.resolution.setEnabled(False)
            self.resolution.setCurrentText("Off")
            self.step_size_label.setChecked(True)
            self.step_size.setEnabled(True)
            self.step_size.setValue(0.6)

    def check_step_size(self):
        if self.step_size_label.isChecked():
            self.resolution_label.setChecked(False)
            self.resolution.setEnabled(False)
            self.resolution.setCurrentText("Off")
            self.step_size.setEnabled(True)
            self.step_size.setValue(0.6)
        else:
            self.resolution_label.setChecked(True)
            self.resolution.setEnabled(True)
            self.resolution.setCurrentText(self._default.resolution)
            self.step_size.setEnabled(False)
            self.step_size.setValue(self._default.step)

    def get_KVFinder_PATH(self) -> str:
        """
        Get KVFinder_PATH environment variable
        """
        from PyQt5.QtWidgets import QMessageBox

        # Get KVFinder_PATH
        KVFinder_PATH = os.getenv("KVFinder_PATH")

        # Check if KVFinder_PATH was found by os.getenv()
        if KVFinder_PATH is None:
            # Check configuration files
            for fn in [".bash_profile", ".bashrc", ".zshrc"]:
                fn = os.path.join(os.getenv("HOME"), fn)
                if os.path.exists(fn):
                    with open(fn, "r") as envs:
                        for line in envs:
                            if line.find("export KVFinder_PATH") == 0:
                                KVFinder_PATH = line.split("=")[1].rstrip("\n")
                                QMessageBox.warning(
                                    self,
                                    "Warning",
                                    f"Check File Locations!\nKVFinder_PATH was foun in {fn}.",
                                )
                                return KVFinder_PATH
            # KVFinder_PATH was not found
            QMessageBox.warning(
                self,
                "Warning",
                f"KVFinder_PATH was not found!\nSet paths on File Locations!\nOtherwise, pyKVFinder cannot be executed in PyMOL2 pyKVFinder Tools.",
            )
            return ""
        else:
            return KVFinder_PATH

    def run(self) -> None:
        import subprocess, time

        # Create parameters.toml
        if self.save_parameters():
            # Running pyKVFinder
            print(
                f"\n[==> Running pyKVFinder for: {os.path.join(self.output_dir_path.text(), 'KV_Files', self.base_name.text(), f'{self.input.currentText()}.pdb')}"
            )
            start = time.time()
            subprocess.call(
                self.pyKVFinder.text().replace(" ", "\\ "), stdout=subprocess.PIPE
            )
            ncavs = self.get_number_of_cavities()
            elapsed_time = time.time() - start
            print(f"> Cavities detected: {ncavs}")
            print(f"> Elapsed time: {elapsed_time:.2f} seconds")

            # Copy parameters file
            if os.path.exists(
                f"{os.path.join(self.output_dir_path.text(), 'KV_Files', f'parameters_{self.base_name.text()}.toml')}"
            ):
                os.remove(
                    f"{os.path.join(self.output_dir_path.text(), 'KV_Files', f'parameters_{self.base_name.text()}.toml')}"
                )
            os.rename(
                "parameters.toml",
                f"{os.path.join(self.output_dir_path.text(), 'KV_Files', f'parameters_{self.base_name.text()}.toml')}",
            )

            # Load successfull run
            self.results_file_entry.setText(
                f"{os.path.join(self.output_dir_path.text(), 'KV_Files', self.base_name.text(), f'{self.base_name.text()}.KVFinder.results.toml')}"
            )
            if ncavs > 0:
                self.tabs.setCurrentIndex(2)
                self.load_results()
            elif ncavs == 0:
                from PyQt5.QtWidgets import QMessageBox

                QMessageBox.warning(self, "Warning!", "No cavities found!")
            else:
                from PyQt5.QtWidgets import QMessageBox

                QMessageBox.critical(
                    self, "Error!", "An error occurred during cavity detection!"
                )
                return False
        else:
            from PyQt5.QtWidgets import QMessageBox

            QMessageBox.critical(
                self,
                "Error",
                "An error occurred while creating the parameters file! Check the pyKVFinder parameters!",
            )
            return False

        return True

    def get_number_of_cavities(self):
        # Read results (Ubuntu/macOS)
        results = toml.load(
            f"{os.path.join(self.output_dir_path.text(), 'KV_Files', self.base_name.text(), f'{self.base_name.text()}.KVFinder.results.toml')}"
        )

        return len(results["RESULTS"]["VOLUME"].keys())

    def show_grid(self) -> None:
        """
        Callback for the "Show Grid" button
        - Get minimum and maximum coordinates of the KVFinder-web 3D-grid, dependent on selected parameters.
        :return: Call draw_grid function with minimum and maximum coordinates or return Error.
        """
        from pymol import cmd

        global x, y, z

        if self.input.count() > 0:
            # Get minimum and maximum dimensions of target PDB
            pdb = self.input.currentText()
            ([min_x, min_y, min_z], [max_x, max_y, max_z]) = cmd.get_extent(pdb)

            # Get Probe Out value
            probe_out = self.probe_out.value()
            probe_out = round(probe_out - round(probe_out, 4) % round(0.6, 4), 1)

            # Prepare dimensions
            min_x = round(min_x - (min_x % 0.6), 1) - probe_out
            min_y = round(min_y - (min_y % 0.6), 1) - probe_out
            min_z = round(min_z - (min_z % 0.6), 1) - probe_out
            max_x = round(max_x - (max_x % 0.6) + 0.6, 1) + probe_out
            max_y = round(max_y - (max_y % 0.6) + 0.6, 1) + probe_out
            max_z = round(max_z - (max_z % 0.6) + 0.6, 1) + probe_out

            # Get center of each dimension (x, y, z)
            x = (min_x + max_x) / 2
            y = (min_y + max_y) / 2
            z = (min_z + max_z) / 2

            # Draw Grid
            self.draw_grid(min_x, max_x, min_y, max_y, min_z, max_z)
        else:
            from PyQt5.QtWidgets import QMessageBox

            QMessageBox.critical(self, "Error", "Select an input PDB!")
            return

    def draw_grid(self, min_x, max_x, min_y, max_y, min_z, max_z) -> None:
        """
        Draw Grid in PyMOL.
        :param min_x: minimum X coordinate.
        :param max_x: maximum X coordinate.
        :param min_y: minimum Y coordinate.
        :param max_y: maximum Y coordinate.
        :param min_z: minimum Z coordinate.
        :param max_z: maximum Z coordinate.
        :return: grid object in PyMOL.
        """
        from pymol import cmd
        from math import sin, cos

        # Prepare dimensions
        angle1 = 0.0
        angle2 = 0.0
        min_x = x - min_x
        max_x = max_x - x
        min_y = y - min_y
        max_y = max_y - y
        min_z = z - min_z
        max_z = max_z - z

        # Get positions of grid vertices
        # P1
        x1 = (
            -min_x * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + x
        )

        y1 = -min_y * cos(angle1) + (-min_z) * sin(angle1) + y

        z1 = (
            min_x * sin(angle2)
            + min_y * sin(angle1) * cos(angle2)
            - min_z * cos(angle1) * cos(angle2)
            + z
        )

        # P2
        x2 = (
            max_x * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + x
        )

        y2 = (-min_y) * cos(angle1) + (-min_z) * sin(angle1) + y

        z2 = (
            (-max_x) * sin(angle2)
            - (-min_y) * sin(angle1) * cos(angle2)
            + (-min_z) * cos(angle1) * cos(angle2)
            + z
        )

        # P3
        x3 = (
            (-min_x) * cos(angle2)
            - max_y * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + x
        )

        y3 = max_y * cos(angle1) + (-min_z) * sin(angle1) + y

        z3 = (
            -(-min_x) * sin(angle2)
            - max_y * sin(angle1) * cos(angle2)
            + (-min_z) * cos(angle1) * cos(angle2)
            + z
        )

        # P4
        x4 = (
            (-min_x) * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + max_z * cos(angle1) * sin(angle2)
            + x
        )

        y4 = (-min_y) * cos(angle1) + max_z * sin(angle1) + y

        z4 = (
            -(-min_x) * sin(angle2)
            - (-min_y) * sin(angle1) * cos(angle2)
            + max_z * cos(angle1) * cos(angle2)
            + z
        )

        # P5
        x5 = (
            max_x * cos(angle2)
            - max_y * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + x
        )

        y5 = max_y * cos(angle1) + (-min_z) * sin(angle1) + y

        z5 = (
            (-max_x) * sin(angle2)
            - max_y * sin(angle1) * cos(angle2)
            + (-min_z) * cos(angle1) * cos(angle2)
            + z
        )

        # P6
        x6 = (
            max_x * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + max_z * cos(angle1) * sin(angle2)
            + x
        )

        y6 = (-min_y) * cos(angle1) + max_z * sin(angle1) + y

        z6 = (
            (-max_x) * sin(angle2)
            - (-min_y) * sin(angle1) * cos(angle2)
            + max_z * cos(angle1) * cos(angle2)
            + z
        )

        # P7
        x7 = (
            (-min_x) * cos(angle2)
            - max_y * sin(angle1) * sin(angle2)
            + max_z * cos(angle1) * sin(angle2)
            + x
        )

        y7 = max_y * cos(angle1) + max_z * sin(angle1) + y

        z7 = (
            -(-min_x) * sin(angle2)
            - max_y * sin(angle1) * cos(angle2)
            + max_z * cos(angle1) * cos(angle2)
            + z
        )

        # P8
        x8 = (
            max_x * cos(angle2)
            - max_y * sin(angle1) * sin(angle2)
            + max_z * cos(angle1) * sin(angle2)
            + x
        )

        y8 = max_y * cos(angle1) + max_z * sin(angle1) + y

        z8 = (
            (-max_x) * sin(angle2)
            - max_y * sin(angle1) * cos(angle2)
            + max_z * cos(angle1) * cos(angle2)
            + z
        )

        # Create box object
        if "grid" in cmd.get_names("objects"):
            cmd.delete("grid")

        # Create vertices
        cmd.pseudoatom("grid", name="v2", pos=[x2, y2, z2], color="white")
        cmd.pseudoatom("grid", name="v3", pos=[x3, y3, z3], color="white")
        cmd.pseudoatom("grid", name="v4", pos=[x4, y4, z4], color="white")
        cmd.pseudoatom("grid", name="v5", pos=[x5, y5, z5], color="white")
        cmd.pseudoatom("grid", name="v6", pos=[x6, y6, z6], color="white")
        cmd.pseudoatom("grid", name="v7", pos=[x7, y7, z7], color="white")
        cmd.pseudoatom("grid", name="v8", pos=[x8, y8, z8], color="white")

        # Connect vertices
        cmd.select("vertices", "(name v3,v7)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v2,v6)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v5,v8)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v2,v5)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v4,v6)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v4,v7)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v3,v5)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v6,v8)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v7,v8)")
        cmd.bond("vertices", "vertices")
        cmd.pseudoatom("grid", name="v1x", pos=[x1, y1, z1], color="white")
        cmd.pseudoatom("grid", name="v2x", pos=[x2, y2, z2], color="white")
        cmd.select("vertices", "(name v1x,v2x)")
        cmd.bond("vertices", "vertices")
        cmd.pseudoatom("grid", name="v1y", pos=[x1, y1, z1], color="white")
        cmd.pseudoatom("grid", name="v3y", pos=[x3, y3, z3], color="white")
        cmd.select("vertices", "(name v1y,v3y)")
        cmd.bond("vertices", "vertices")
        cmd.pseudoatom("grid", name="v4z", pos=[x4, y4, z4], color="white")
        cmd.pseudoatom("grid", name="v1z", pos=[x1, y1, z1], color="white")
        cmd.select("vertices", "(name v1z,v4z)")
        cmd.bond("vertices", "vertices")
        cmd.delete("vertices")

    def restore(self, is_startup=False) -> None:
        """
        Callback for the "Restore Default Values" button
        """
        from pymol import cmd
        from PyQt5.QtWidgets import QMessageBox, QCheckBox

        # Restore Results Tab
        if not is_startup:
            reply = QMessageBox(self)
            reply.setText("Also restore Results Visualization tab?")
            reply.setWindowTitle("Restore Values")
            reply.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
            reply.setIcon(QMessageBox.Information)
            reply.checkbox = QCheckBox("Also remove input and ligand PDBs?")
            reply.layout = reply.layout()
            reply.layout.addWidget(reply.checkbox, 1, 2)
            if reply.exec_() == QMessageBox.Yes:
                # Remove cavities, residues and pdbs (input, ligand, cavity)
                cmd.delete("cavities")
                cmd.delete("residues")
                if self.input_pdb and reply.checkbox.isChecked():
                    cmd.delete(self.input_pdb)
                if self.ligand_pdb and reply.checkbox.isChecked():
                    cmd.delete(self.ligand_pdb)
                if self.cavity_pdb:
                    cmd.delete(self.cavity_pdb)
                results = self.input_pdb = self.ligand_pdb = self.cavity_pdb = None
                cmd.frame(1)

                # Clean results
                self.clean_results()
                self.results_file_entry.clear()

        # Restore PDB and ligand input
        self.refresh(self.input)
        self.refresh(self.ligand)

        # Delete grid
        cmd.delete("grid")

        ### Main tab ###
        self.step_size_label.setChecked(False)
        self.step_size.setValue(self._default.step)
        self.step_size.setEnabled(False)
        self.resolution_label.setChecked(True)
        self.resolution.setCurrentText(self._default.resolution)
        self.resolution.setEnabled(True)
        self.base_name.setText(self._default.base_name)
        self.probe_in.setValue(self._default.probe_in)
        self.probe_out.setValue(self._default.probe_out)
        self.volume_cutoff.setValue(self._default.volume_cutoff)
        self.removal_distance.setValue(self._default.removal_distance)
        self.surface.setCurrentText(self._default.surface)
        self.cavity_representation.setCurrentText(self._default.cavity_representation)
        self.output_dir_path.setText(self._default.output_dir_path)
        self.pyKVFinder.setText(self._default.pyKVFinder)
        self.dictionary.setText(self._default.dictionary)

        ### Search Space Tab ###
        # Box Adjustment
        self.box_adjustment.setChecked(self._default.box_adjustment)
        self.padding.setValue(self._default.padding)
        self.delete_box()
        # Ligand Adjustment
        self.ligand_adjustment.setChecked(self._default.ligand_adjustment)
        self.ligand.clear()
        self.ligand_cutoff.setValue(self._default.ligand_cutoff)

    def refresh(self, combo_box) -> None:
        """
        Callback for the "Refresh" button
        """
        from pymol import cmd

        combo_box.clear()
        for item in cmd.get_names("all"):
            if (
                cmd.get_type(item) == "object:molecule"
                and item != "box"
                and item != "grid"
                and item != "cavities"
                and item != "residues"
                and item[-16:] != ".KVFinder.output"
                and item != "target_exclusive"
            ):
                combo_box.addItem(item)

        return

    def select_directory(self) -> None:
        """
        Callback for the "Browse ..." button
        Open a QFileDialog to select a directory.
        """
        from PyQt5.QtWidgets import QFileDialog
        from PyQt5.QtCore import QDir

        fname = QFileDialog.getExistingDirectory(
            caption="Choose Output Directory", directory=os.getcwd()
        )

        if fname:
            fname = QDir.toNativeSeparators(fname)
            if os.path.isdir(fname):
                self.output_dir_path.setText(fname)

        return

    def select_file(self, caption, entry, filters) -> None:
        """
        Callback for the "Browse ..." button
        Open a QFileDialog to select a file.
        """
        from PyQt5.QtWidgets import QFileDialog
        from PyQt5.QtCore import QDir

        # Get results file
        fname, _ = QFileDialog.getOpenFileName(
            self, caption=caption, directory=os.getcwd(), filter=filters
        )

        if fname:
            fname = QDir.toNativeSeparators(fname)
            if os.path.exists(fname):
                entry.setText(fname)

        return

    def set_box(self) -> None:
        """
        Create box coordinates, enable 'Delete Box' and 'Redraw Box' buttons and call draw_box function.
        :param padding: box padding value.
        """
        from pymol import cmd

        # Delete Box object in PyMOL
        if "box" in cmd.get_names("all"):
            cmd.delete("box")
        # Get dimensions of selected residues
        selection = "sele"
        if selection in cmd.get_names("selections"):
            ([min_x, min_y, min_z], [max_x, max_y, max_z]) = cmd.get_extent(selection)
        else:
            ([min_x, min_y, min_z], [max_x, max_y, max_z]) = cmd.get_extent("")

        # Get center of each dimension (x, y, z)
        self.x = (min_x + max_x) / 2
        self.y = (min_y + max_y) / 2
        self.z = (min_z + max_z) / 2

        # Set Box variables in interface
        self.min_x.setValue(round(self.x - (min_x - self.padding.value()), 1))
        self.max_x.setValue(round((max_x + self.padding.value()) - self.x, 1))
        self.min_y.setValue(round(self.y - (min_y - self.padding.value()), 1))
        self.max_y.setValue(round((max_y + self.padding.value()) - self.y, 1))
        self.min_z.setValue(round(self.z - (min_z - self.padding.value()), 1))
        self.max_z.setValue(round((max_z + self.padding.value()) - self.z, 1))
        self.angle1.setValue(0)
        self.angle2.setValue(0)

        # Setting background box values
        self.min_x_set = self.min_x.value()
        self.max_x_set = self.max_x.value()
        self.min_y_set = self.min_y.value()
        self.max_y_set = self.max_y.value()
        self.min_z_set = self.min_z.value()
        self.max_z_set = self.max_z.value()
        self.angle1_set = self.angle1.value()
        self.angle2_set = self.angle2.value()
        self.padding_set = self.padding.value()

        # Draw box
        self.draw_box()

        # Enable/Disable buttons
        self.button_draw_box.setEnabled(False)
        self.button_redraw_box.setEnabled(True)
        self.min_x.setEnabled(True)
        self.min_y.setEnabled(True)
        self.min_z.setEnabled(True)
        self.max_x.setEnabled(True)
        self.max_y.setEnabled(True)
        self.max_z.setEnabled(True)
        self.angle1.setEnabled(True)
        self.angle2.setEnabled(True)

    def draw_box(self) -> None:
        """
        Draw box in PyMOL interface.
        :return: box object.
        """
        from math import pi, sin, cos
        import pymol
        from pymol import cmd

        # Convert angle
        angle1 = (self.angle1.value() / 180.0) * pi
        angle2 = (self.angle2.value() / 180.0) * pi

        # Get positions of box vertices
        # P1
        x1 = (
            -self.min_x.value() * cos(angle2)
            - (-self.min_y.value()) * sin(angle1) * sin(angle2)
            + (-self.min_z.value()) * cos(angle1) * sin(angle2)
            + self.x
        )

        y1 = (
            -self.min_y.value() * cos(angle1)
            + (-self.min_z.value()) * sin(angle1)
            + self.y
        )

        z1 = (
            self.min_x.value() * sin(angle2)
            + self.min_y.value() * sin(angle1) * cos(angle2)
            - self.min_z.value() * cos(angle1) * cos(angle2)
            + self.z
        )

        # P2
        x2 = (
            self.max_x.value() * cos(angle2)
            - (-self.min_y.value()) * sin(angle1) * sin(angle2)
            + (-self.min_z.value()) * cos(angle1) * sin(angle2)
            + self.x
        )

        y2 = (
            (-self.min_y.value()) * cos(angle1)
            + (-self.min_z.value()) * sin(angle1)
            + self.y
        )

        z2 = (
            (-self.max_x.value()) * sin(angle2)
            - (-self.min_y.value()) * sin(angle1) * cos(angle2)
            + (-self.min_z.value()) * cos(angle1) * cos(angle2)
            + self.z
        )

        # P3
        x3 = (
            (-self.min_x.value()) * cos(angle2)
            - self.max_y.value() * sin(angle1) * sin(angle2)
            + (-self.min_z.value()) * cos(angle1) * sin(angle2)
            + self.x
        )

        y3 = (
            self.max_y.value() * cos(angle1)
            + (-self.min_z.value()) * sin(angle1)
            + self.y
        )

        z3 = (
            -(-self.min_x.value()) * sin(angle2)
            - self.max_y.value() * sin(angle1) * cos(angle2)
            + (-self.min_z.value()) * cos(angle1) * cos(angle2)
            + self.z
        )

        # P4
        x4 = (
            (-self.min_x.value()) * cos(angle2)
            - (-self.min_y.value()) * sin(angle1) * sin(angle2)
            + self.max_z.value() * cos(angle1) * sin(angle2)
            + self.x
        )

        y4 = (
            (-self.min_y.value()) * cos(angle1)
            + self.max_z.value() * sin(angle1)
            + self.y
        )

        z4 = (
            -(-self.min_x.value()) * sin(angle2)
            - (-self.min_y.value()) * sin(angle1) * cos(angle2)
            + self.max_z.value() * cos(angle1) * cos(angle2)
            + self.z
        )

        # P5
        x5 = (
            self.max_x.value() * cos(angle2)
            - self.max_y.value() * sin(angle1) * sin(angle2)
            + (-self.min_z.value()) * cos(angle1) * sin(angle2)
            + self.x
        )

        y5 = (
            self.max_y.value() * cos(angle1)
            + (-self.min_z.value()) * sin(angle1)
            + self.y
        )

        z5 = (
            (-self.max_x.value()) * sin(angle2)
            - self.max_y.value() * sin(angle1) * cos(angle2)
            + (-self.min_z.value()) * cos(angle1) * cos(angle2)
            + self.z
        )

        # P6
        x6 = (
            self.max_x.value() * cos(angle2)
            - (-self.min_y.value()) * sin(angle1) * sin(angle2)
            + self.max_z.value() * cos(angle1) * sin(angle2)
            + self.x
        )

        y6 = (
            (-self.min_y.value()) * cos(angle1)
            + self.max_z.value() * sin(angle1)
            + self.y
        )

        z6 = (
            (-self.max_x.value()) * sin(angle2)
            - (-self.min_y.value()) * sin(angle1) * cos(angle2)
            + self.max_z.value() * cos(angle1) * cos(angle2)
            + self.z
        )

        # P7
        x7 = (
            (-self.min_x.value()) * cos(angle2)
            - self.max_y.value() * sin(angle1) * sin(angle2)
            + self.max_z.value() * cos(angle1) * sin(angle2)
            + self.x
        )

        y7 = (
            self.max_y.value() * cos(angle1) + self.max_z.value() * sin(angle1) + self.y
        )

        z7 = (
            -(-self.min_x.value()) * sin(angle2)
            - self.max_y.value() * sin(angle1) * cos(angle2)
            + self.max_z.value() * cos(angle1) * cos(angle2)
            + self.z
        )

        # P8
        x8 = (
            self.max_x.value() * cos(angle2)
            - self.max_y.value() * sin(angle1) * sin(angle2)
            + self.max_z.value() * cos(angle1) * sin(angle2)
            + self.x
        )

        y8 = (
            self.max_y.value() * cos(angle1) + self.max_z.value() * sin(angle1) + self.y
        )

        z8 = (
            (-self.max_x.value()) * sin(angle2)
            - self.max_y.value() * sin(angle1) * cos(angle2)
            + self.max_z.value() * cos(angle1) * cos(angle2)
            + self.z
        )

        # Create box object
        pymol.stored.list = []
        if "box" in cmd.get_names("all"):
            cmd.iterate("box", "stored.list.append((name, color))", quiet=1)
        list_color = pymol.stored.list
        cmd.delete("box")
        if len(list_color) > 0:
            for item in list_color:
                at_name = item[0]
                at_c = item[1]
                cmd.set_color(at_name + "color", cmd.get_color_tuple(at_c))
        else:
            for at_name in [
                "v2",
                "v3",
                "v4",
                "v5",
                "v6",
                "v7",
                "v8",
                "v1x",
                "v1y",
                "v1z",
                "v2x",
                "v3y",
                "v4z",
            ]:
                cmd.set_color(at_name + "color", [0.86, 0.86, 0.86])

        # Create vertices
        cmd.pseudoatom("box", name="v2", pos=[x2, y2, z2], color="v2color")
        cmd.pseudoatom("box", name="v3", pos=[x3, y3, z3], color="v3color")
        cmd.pseudoatom("box", name="v4", pos=[x4, y4, z4], color="v4color")
        cmd.pseudoatom("box", name="v5", pos=[x5, y5, z5], color="v5color")
        cmd.pseudoatom("box", name="v6", pos=[x6, y6, z6], color="v6color")
        cmd.pseudoatom("box", name="v7", pos=[x7, y7, z7], color="v7color")
        cmd.pseudoatom("box", name="v8", pos=[x8, y8, z8], color="v8color")

        # Connect vertices
        cmd.select("vertices", "(name v3,v7)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v2,v6)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v5,v8)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v2,v5)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v4,v6)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v4,v7)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v3,v5)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v6,v8)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v7,v8)")
        cmd.bond("vertices", "vertices")
        cmd.pseudoatom("box", name="v1x", pos=[x1, y1, z1], color="red")
        cmd.pseudoatom("box", name="v2x", pos=[x2, y2, z2], color="red")
        cmd.select("vertices", "(name v1x,v2x)")
        cmd.bond("vertices", "vertices")
        cmd.pseudoatom("box", name="v1y", pos=[x1, y1, z1], color="forest")
        cmd.pseudoatom("box", name="v3y", pos=[x3, y3, z3], color="forest")
        cmd.select("vertices", "(name v1y,v3y)")
        cmd.bond("vertices", "vertices")
        cmd.pseudoatom("box", name="v4z", pos=[x4, y4, z4], color="blue")
        cmd.pseudoatom("box", name="v1z", pos=[x1, y1, z1], color="blue")
        cmd.select("vertices", "(name v1z,v4z)")
        cmd.bond("vertices", "vertices")
        cmd.delete("vertices")

    def delete_box(self) -> None:
        """
        Delete box object, disable 'Delete Box' and 'Redraw Box' buttons and enable 'Draw Box' button.
        """
        from pymol import cmd

        # Reset all box variables
        self.x = 0
        self.y = 0
        self.z = 0
        # self.min_x_set = 0.0
        # self.max_x_set = 0.0
        # self.min_y_set = 0.0
        # self.max_y_set = 0.0
        # self.min_z_set = 0.0
        # self.max_z_set = 0.0
        # self.angle1_set = 0.0
        # self.angle2_set = 0.0
        # self.padding_set = 3.5

        # Delete Box and Vertices objects in PyMOL
        cmd.delete("vertices")
        cmd.delete("box")

        # Set Box variables in the interface
        self.min_x.setValue(self._default.min_x)
        self.max_x.setValue(self._default.max_x)
        self.min_y.setValue(self._default.min_y)
        self.max_y.setValue(self._default.max_y)
        self.min_z.setValue(self._default.min_z)
        self.max_z.setValue(self._default.max_z)
        self.angle1.setValue(self._default.angle1)
        self.angle2.setValue(self._default.angle2)

        # Change state of buttons in the interface
        self.button_draw_box.setEnabled(True)
        self.button_redraw_box.setEnabled(False)
        self.min_x.setEnabled(False)
        self.min_y.setEnabled(False)
        self.min_z.setEnabled(False)
        self.max_x.setEnabled(False)
        self.max_y.setEnabled(False)
        self.max_z.setEnabled(False)
        self.angle1.setEnabled(False)
        self.angle2.setEnabled(False)

    def redraw_box(self) -> None:
        """
        Redraw box in PyMOL interface.
        :param padding: box padding.
        :return: box object.
        """
        from pymol import cmd

        # Provided a selection
        if "sele" in cmd.get_names("selections"):
            # Get dimensions of selected residues
            ([min_x, min_y, min_z], [max_x, max_y, max_z]) = cmd.get_extent("sele")

            if (
                self.min_x.value() != self.min_x_set
                or self.max_x.value() != self.max_x_set
                or self.min_y.value() != self.min_y_set
                or self.max_y.value() != self.max_y_set
                or self.min_z.value() != self.min_z_set
                or self.max_z.value() != self.max_z_set
                or self.angle1.value() != self.angle1_set
                or self.angle2.value() != self.angle2_set
            ):
                self.min_x_set = self.min_x.value()
                self.max_x_set = self.max_x.value()
                self.min_y_set = self.min_y.value()
                self.max_y_set = self.max_y.value()
                self.min_z_set = self.min_z.value()
                self.max_z_set = self.max_z.value()
                self.angle1_set = self.angle1.value()
                self.angle2_set = self.angle2.value()
            # Padding or selection altered
            else:
                # Get center of each dimension (x, y, z)
                self.x = (min_x + max_x) / 2
                self.y = (min_y + max_y) / 2
                self.z = (min_z + max_z) / 2

                # Set background box values
                self.min_x_set = (
                    round(self.x - (min_x - self.padding.value()), 1)
                    + self.min_x.value()
                    - self.min_x_set
                )
                self.max_x_set = (
                    round((max_x + self.padding.value()) - self.x, 1)
                    + self.max_x.value()
                    - self.max_x_set
                )
                self.min_y_set = (
                    round(self.y - (min_y - self.padding.value()), 1)
                    + self.min_y.value()
                    - self.min_y_set
                )
                self.max_y_set = (
                    round((max_y + self.padding.value()) - self.y, 1)
                    + self.max_y.value()
                    - self.max_y_set
                )
                self.min_z_set = (
                    round(self.z - (min_z - self.padding.value()), 1)
                    + self.min_z.value()
                    - self.min_z_set
                )
                self.max_z_set = (
                    round((max_z + self.padding.value()) - self.z, 1)
                    + self.max_z.value()
                    - self.max_z_set
                )
                self.angle1_set = 0 + self.angle1.value()
                self.angle2_set = 0 + self.angle2.value()
                self.padding_set = self.padding.value()
        # Not provided a selection
        else:
            if (
                self.min_x.value() != self.min_x_set
                or self.max_x.value() != self.max_x_set
                or self.min_y.value() != self.min_y_set
                or self.max_y.value() != self.max_y_set
                or self.min_z.value() != self.min_z_set
                or self.max_z.value() != self.max_z_set
                or self.angle1.value() != self.angle1_set
                or self.angle2.value() != self.angle2_set
            ):
                self.min_x_set = self.min_x.value()
                self.max_x_set = self.max_x.value()
                self.min_y_set = self.min_y.value()
                self.max_y_set = self.max_y.value()
                self.min_z_set = self.min_z.value()
                self.max_z_set = self.max_z.value()
                self.angle1_set = self.angle1.value()
                self.angle2_set = self.angle2.value()

            if self.padding_set != self.padding.value():
                # Prepare dimensions without old padding
                min_x = self.padding_set - self.min_x_set
                max_x = self.max_x_set - self.padding_set
                min_y = self.padding_set - self.min_y_set
                max_y = self.max_y_set - self.padding_set
                min_z = self.padding_set - self.min_z_set
                max_z = self.max_z_set - self.padding_set

                # Get center of each dimension (x, y, z)
                self.x = (min_x + max_x) / 2
                self.y = (min_y + max_y) / 2
                self.z = (min_z + max_z) / 2

                # Set background box values
                self.min_x_set = round(self.x - (min_x - self.padding.value()), 1)
                self.max_x_set = round((max_x + self.padding.value()) - self.x, 1)
                self.min_y_set = round(self.y - (min_y - self.padding.value()), 1)
                self.max_y_set = round((max_y + self.padding.value()) - self.y, 1)
                self.min_z_set = round(self.z - (min_z - self.padding.value()), 1)
                self.max_z_set = round((max_z + self.padding.value()) - self.z, 1)
                self.angle1_set = self.angle1.value()
                self.angle2_set = self.angle2.value()
                self.padding_set = self.padding.value()

        # Set Box variables in the interface
        self.min_x.setValue(self.min_x_set)
        self.max_x.setValue(self.max_x_set)
        self.min_y.setValue(self.min_y_set)
        self.max_y.setValue(self.max_y_set)
        self.min_z.setValue(self.min_z_set)
        self.max_z.setValue(self.max_z_set)
        self.angle1.setValue(self.angle1_set)
        self.angle2.setValue(self.angle2_set)

        # Redraw box
        self.draw_box()

    def box_adjustment_help(self) -> None:
        from PyQt5 import QtWidgets, QtCore

        text = QtCore.QCoreApplication.translate(
            "pyKVFinder",
            '<html><head/><body><p align="justify"><span style=" font-weight:600; text-decoration: underline;">Box Adjustment mode:</span></p><p align="justify">- Create a selection (optional);</p><p align="justify">- Define a <span style=" font-weight:600;">Padding</span> (optional);</p><p align="justify">- Click on <span style=" font-weight:600;">Draw Box</span> button.</p><p align="justify"><br/><span style="text-decoration: underline;">Customize your <span style=" font-weight:600;">box</span></span>:</p><p align="justify">- Change one item at a time (e.g. <span style=" font-style:italic;">Padding</span>, <span style=" font-style:italic;">Minimum X</span>, <span style=" font-style:italic;">Maximum X</span>, ...);</p><p align="justify">- Click on <span style=" font-weight:600;">Redraw Box</span> button.<br/></p><p><span style=" font-weight:400; text-decoration: underline;">Delete </span><span style=" text-decoration: underline;">box</span><span style=" font-weight:400; text-decoration: underline;">:</span></p><p align="justify">- Click on <span style=" font-weight:600;">Delete Box</span> button.<br/></p><p align="justify"><span style="text-decoration: underline;">Colors of the <span style=" font-weight:600;">box</span> object:</span></p><p align="justify">- <span style=" font-weight:600;">Red</span> corresponds to <span style=" font-weight:600;">X</span> axis;</p><p align="justify">- <span style=" font-weight:600;">Green</span> corresponds to <span style=" font-weight:600;">Y</span> axis;</p><p align="justify">- <span style=" font-weight:600;">Blue</span> corresponds to <span style=" font-weight:600;">Z</span> axis.</p></body></html>',
            None,
        )
        help_information = QtWidgets.QMessageBox(self)
        help_information.setText(text)
        help_information.setWindowTitle("Help")
        help_information.setStyleSheet("QLabel{min-width:500 px;}")
        help_information.exec_()

    def save_parameters(self) -> None:
        from pymol import cmd

        # Create base directory
        basedir = os.path.join(self.output_dir_path.text(), "KV_Files")
        if not os.path.isdir(basedir):
            os.mkdir(basedir)

        # Create base_name directory
        basedir = os.path.join(basedir, self.base_name.text())
        if not os.path.isdir(basedir):
            os.mkdir(basedir)

        # Save input pdb
        if self.input.currentText() != "":
            for x in cmd.get_names("all"):
                if x == self.input.currentText():
                    pdb = os.path.join(
                        os.path.join(basedir, f"{self.input.currentText()}.pdb")
                    )
                    cmd.save(pdb, self.input.currentText(), 0, "pdb")
        else:
            from PyQt5.QtWidgets import QMessageBox

            QMessageBox.critical(self, "Error", "Select an input PDB!")
            return False

        # Save ligand pdb
        if self.ligand_adjustment.isChecked():
            if self.ligand.currentText() != "":
                for x in cmd.get_names("all"):
                    if x == self.ligand.currentText():
                        ligand = os.path.join(
                            os.path.join(
                                basedir, f"{self.ligand.currentText()}.ligand.pdb"
                            )
                        )
                        cmd.save(ligand, self.ligand.currentText(), 0, "pdb")
            else:
                from PyQt5.QtWidgets import QMessageBox

                QMessageBox.critical(self, "Error", "Select an ligand PDB!")
                return False
        else:
            ligand = "-"

        if self.box_adjustment.isChecked():
            if "box" not in cmd.get_names("all"):
                from PyQt5.QtWidgets import QMessageBox

                QMessageBox.critical(self, "Error", "Draw a box in PyMOL!")
                return False

        with open("parameters.toml", "w") as f:
            f.write("# TOML configuration file for pyKVFinder software.\n")
            f.write('\ntitle = "pyKVFinder parameters file"\n')

            f.write("\n[FILES_PATH]\n")
            f.write("# The path of van der Waals radii dictionary for pyKVFinder.\n")
            f.write(f'dictionary = "{self.dictionary.text()}"\n')
            f.write("# The path of the input PDB file.\n")
            f.write(f'pdb = "{pdb}"\n')
            f.write("# The path of the output directory.\n")
            f.write(f'output = "{self.output_dir_path.text()}"\n')
            f.write("# Base name for output files.\n")
            f.write(f'base_name = "{self.base_name.text()}"\n')
            f.write("# Path for the ligand's PDB file.\n")
            f.write(f'ligand = "{ligand}"\n')

            f.write("\n[SETTINGS]\n")
            f.write("# Settings for pyKVFinder software.\n")
            f.write("\n[SETTINGS.modes]\n")
            f.write(
                "# Whole protein mode defines the search space as the whole protein.\n"
            )
            f.write(
                f"whole_protein_mode = {'true' if not self.box_adjustment.isChecked() else 'false'}\n"
            )
            f.write(
                "# Box adjustment mode defines the search space as the box drawn in PyMOL.\n"
            )
            f.write(
                f"box_mode = {'true' if self.box_adjustment.isChecked() else 'false'}\n"
            )
            f.write(
                "# Resolution mode implicitly sets the step size (grid spacing) of the 3D grid.\n"
            )
            f.write(
                "# If set to High, sets a voxel volume of 0.2. If set to Medium, sets a voxel volume of 0.1. If set to Low, it sets a voxel volume of 0.01. If set to Off, the step size must be set explicitly.\n"
            )
            f.write(
                f"resolution_mode = \"{self.resolution.currentText() if self.resolution_label.isChecked() else 'Off'}\"\n"
            )
            f.write(
                "# Surface mode defines the type of surface representation to be applied, van der Waals molecular surface (true) or solvent accessible surface (false).\n"
            )
            f.write(
                f"surface_mode = {'true' if self.surface.currentText() == 'Molecular Surface (VdW)' else 'false'}\n"
            )
            f.write(
                "# Cavity representation defines whether cavities are exported to the output PDB file as filled cavities (true) or filtered cavities (false).\n"
            )
            f.write(
                f"kvp_mode = {'true' if self.cavity_representation.currentText() == 'Full' else 'false'}\n"
            )
            f.write(
                "# Ligand adjustment mode defines the search space around the ligand.\n"
            )
            f.write(
                f"ligand_mode = {'true' if self.ligand_adjustment.isChecked() else 'false'}\n"
            )

            f.write("\n[SETTINGS.step_size]\n")
            f.write(
                "# Sets the 3D grid spacing. It directly affects accuracy and runtime.\n"
            )
            step = self.step_size.value() if self.step_size_label.isChecked() else 0.0
            f.write(f"step_size = {step:.2f}\n")

            f.write("\n[SETTINGS.probes]\n")
            f.write(
                "# pyKVFinder works with a dual probe system. A smaller probe, called Probe In, and a bigger one, called Probe Out, rolls around the protein.\n"
            )
            f.write(
                "# Points reached by the Probe In, but not the Probe Out are considered cavity points.\n"
            )
            f.write("# Sets the Probe In diameter. Default: 1.4 angstroms.\n")
            f.write(f"probe_in = {self.probe_in.value():.2f}\n")
            f.write("# Sets the Probe Out diameter. Default: 4.0 angstroms.\n")
            f.write(f"probe_out = {self.probe_out.value():.2f}\n")

            f.write("\n[SETTINGS.cutoffs]\n")
            f.write(
                "# Sets a volume cutoff for the detected cavities. Default: 5.0 angstroms.\n"
            )
            f.write(f"volume_cutoff = {self.volume_cutoff.value():.2f}\n")
            f.write(
                "# Sets a distance cutoff for a search space around the ligand in ligand adjustment mode. Default: 5.0 angstroms.\n"
            )
            f.write(f"ligand_cutoff = {self.ligand_cutoff.value():.2f}\n")
            f.write(
                "# Sets a removal distance for the cavity frontier, which is defined by comparing Probe In and Probe Out surfaces. Default: 2.4 angstroms.\n"
            )
            f.write(f"removal_distance = {self.removal_distance.value():.2f}\n")

            f.write("\n[SETTINGS.visiblebox]\n")
            f.write(
                "# Coordinates of the vertices that define the visible 3D grid. Only four points are required to define the search space.\n\n"
            )
            box = self.create_box_parameters()
            d = {"SETTINGS": {"visiblebox": box}}
            toml.dump(o=d, f=f)

            f.write("\n[SETTINGS.internalbox]\n")
            f.write("# Coordinates of the internal 3D grid. Used for calculations.\n\n")
            box = self.create_box_parameters(is_internal_box=True)
            d = {"SETTINGS": {"internalbox": box}}
            toml.dump(o=d, f=f)

        return True

    def create_box_parameters(
        self, is_internal_box=False
    ) -> Dict[str, Dict[str, float]]:
        from math import pi, cos, sin

        # Get box parameters
        if self.box_adjustment.isChecked():
            min_x = self.min_x_set
            max_x = self.max_x_set
            min_y = self.min_y_set
            max_y = self.max_y_set
            min_z = self.min_z_set
            max_z = self.max_z_set
            angle1 = self.angle1_set
            angle2 = self.angle2_set
        else:
            min_x = 0.0
            max_x = 0.0
            min_y = 0.0
            max_y = 0.0
            min_z = 0.0
            max_z = 0.0
            angle1 = 0.0
            angle2 = 0.0

        # Add probe_out to internal box
        if is_internal_box:
            min_x += self.probe_out.value()
            max_x += self.probe_out.value()
            min_y += self.probe_out.value()
            max_y += self.probe_out.value()
            min_z += self.probe_out.value()
            max_z += self.probe_out.value()

        # Convert angle
        angle1 = (angle1 / 180.0) * pi
        angle2 = (angle2 / 180.0) * pi

        # Get positions of box vertices
        # P1
        x1 = (
            -min_x * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + self.x
        )

        y1 = -min_y * cos(angle1) + (-min_z) * sin(angle1) + self.y

        z1 = (
            min_x * sin(angle2)
            + min_y * sin(angle1) * cos(angle2)
            - min_z * cos(angle1) * cos(angle2)
            + self.z
        )

        # P2
        x2 = (
            max_x * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + self.x
        )

        y2 = (-min_y) * cos(angle1) + (-min_z) * sin(angle1) + self.y

        z2 = (
            (-max_x) * sin(angle2)
            - (-min_y) * sin(angle1) * cos(angle2)
            + (-min_z) * cos(angle1) * cos(angle2)
            + self.z
        )

        # P3
        x3 = (
            (-min_x) * cos(angle2)
            - max_y * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + self.x
        )

        y3 = max_y * cos(angle1) + (-min_z) * sin(angle1) + self.y

        z3 = (
            -(-min_x) * sin(angle2)
            - max_y * sin(angle1) * cos(angle2)
            + (-min_z) * cos(angle1) * cos(angle2)
            + self.z
        )

        # P4
        x4 = (
            (-min_x) * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + max_z * cos(angle1) * sin(angle2)
            + self.x
        )

        y4 = (-min_y) * cos(angle1) + max_z * sin(angle1) + self.y

        z4 = (
            -(-min_x) * sin(angle2)
            - (-min_y) * sin(angle1) * cos(angle2)
            + max_z * cos(angle1) * cos(angle2)
            + self.z
        )

        # Create points
        p1 = {"x": x1, "y": y1, "z": z1}
        p2 = {"x": x2, "y": y2, "z": z2}
        p3 = {"x": x3, "y": y3, "z": z3}
        p4 = {"x": x4, "y": y4, "z": z4}
        box = {"p1": p1, "p2": p2, "p3": p3, "p4": p4}

        return box

    def closeEvent(self, event) -> None:
        """
        Add one step to closeEvent from QMainWindow
        """
        global dialog
        dialog = None

    def load_results(self) -> None:
        from pymol import cmd

        # Get results file
        results_file = self.results_file_entry.text()

        # Check if results file exist
        if os.path.exists(results_file) and results_file.endswith(".toml"):
            print(f"> Loading results from: {self.results_file_entry.text()}")
        else:
            from PyQt5.QtWidgets import QMessageBox

            error_msg = QMessageBox.critical(
                self, "Error", "Results file cannot be opened! Check results file path."
            )
            return False

        # Create global variable for results
        global results

        # Read results (Ubuntu/macOS)
        results = toml.load(results_file)

        if "FILES" in results.keys():
            results["FILES_PATH"] = results.pop("FILES")
        elif "FILES_PATH" in results.keys():
            pass
        else:
            from PyQt5.QtWidgets import QMessageBox

            error_msg = QMessageBox.critical(
                self,
                "Error",
                "Results file has incorrect format! Please check your file.",
            )
            return False

        # Clean results
        self.clean_results()

        # Refresh information
        self.refresh_information()

        # Refresh volume
        self.refresh_volume()

        # Refresh area
        self.refresh_area()

        # Refresh depth
        self.refresh_avg_depth()
        self.refresh_max_depth()
        self.refresh_avg_hydropathy()

        # Refresh residues
        self.refresh_residues()

        # Set default view in results
        self.default_view.setChecked(True)

        # Load files as PyMOL objects
        cmd.delete("cavities")
        cmd.delete("residues")
        cmd.frame(1)

        # Load input
        if "INPUT" in results["FILES_PATH"].keys():
            input_fn = results["FILES_PATH"]["INPUT"]
            self.input_pdb = os.path.basename(input_fn.replace(".pdb", ""))
            self.load_file(input_fn, self.input_pdb)
        else:
            self.input_pdb = None

        # Load ligand
        if "LIGAND" in results["FILES_PATH"].keys():
            ligand_fn = results["FILES_PATH"]["LIGAND"]
            self.ligand_pdb = os.path.basename(ligand_fn.replace(".pdb", ""))
            self.load_file(ligand_fn, self.ligand_pdb)
        else:
            self.ligand_pdb = None

        # Load cavity
        cavity_fn = results["FILES_PATH"]["OUTPUT"]
        self.cavity_pdb = os.path.basename(cavity_fn.replace(".pdb", ""))
        self.load_cavity(cavity_fn, self.cavity_pdb)

        return

    @staticmethod
    def load_cavity(fname, name) -> None:
        from pymol import cmd

        # Remove previous results in objects with same cavity name
        for obj in cmd.get_names("all"):
            if name == obj:
                cmd.delete(obj)

        # Load cavity filename
        if os.path.exists(fname):
            cmd.load(fname, name, zoom=0)
            cmd.hide("everything", name)
            cmd.show("nonbonded", name)

    @staticmethod
    def load_file(fname, name) -> None:
        from pymol import cmd

        # Remove previous results in objects with same cavity name
        for obj in cmd.get_names("all"):
            if name == obj:
                cmd.delete(obj)

        # Load cavity filename
        if os.path.exists(fname):
            cmd.load(fname, name, zoom=0)

    def refresh_information(self) -> None:
        # Input File
        if "INPUT" in results["FILES_PATH"].keys():
            self.input_file_entry.setText(f"{results['FILES_PATH']['INPUT']}")
        else:
            self.input_file_entry.setText(f"")

        # Ligand File
        if "LIGAND" in results["FILES_PATH"].keys():
            self.ligand_file_entry.setText(f"{results['FILES_PATH']['LIGAND']}")
        else:
            self.ligand_file_entry.setText(f"")

        # Cavities File
        self.cavities_file_entry.setText(f"{results['FILES_PATH']['OUTPUT']}")

        # Step Size
        if "PARAMETERS" in results.keys():
            if "STEP" in results["PARAMETERS"].keys():
                self.step_size_entry.setText(f"{results['PARAMETERS']['STEP']:.2f}")

        return

    def refresh_volume(self) -> None:
        # Get cavity indexes
        indexes = sorted(results["RESULTS"]["VOLUME"].keys())
        # Include Volume
        for index in indexes:
            item = f"{index}: {results['RESULTS']['VOLUME'][index]}"
            self.volume_list.addItem(item)
        return

    def refresh_area(self) -> None:
        # Get cavity indexes
        indexes = sorted(results["RESULTS"]["AREA"].keys())
        # Include Area
        for index in indexes:
            item = f"{index}: {results['RESULTS']['AREA'][index]}"
            self.area_list.addItem(item)
        return

    def refresh_avg_depth(self) -> None:
        # Get cavity indexes
        indexes = sorted(results["RESULTS"]["AVG_DEPTH"].keys())
        # Include Average Depth
        for index in indexes:
            item = f"{index}: {results['RESULTS']['AVG_DEPTH'][index]}"
            self.avg_depth_list.addItem(item)
        return

    def refresh_max_depth(self) -> None:
        # Get cavity indexes
        indexes = sorted(results["RESULTS"]["MAX_DEPTH"].keys())
        # Include Maximum Depth
        for index in indexes:
            item = f"{index}: {results['RESULTS']['MAX_DEPTH'][index]}"
            self.max_depth_list.addItem(item)
        return

    def refresh_avg_hydropathy(self) -> None:
        # Get cavity indexes
        indexes = sorted(results["RESULTS"]["AVG_HYDROPATHY"].keys())
        # Include Average Hydropathy
        for index in indexes:
            if index != "EisenbergWeiss":
                item = f"{index}: {results['RESULTS']['AVG_HYDROPATHY'][index]}"
                self.avg_hydropathy_list.addItem(item)
        return

    def refresh_residues(self) -> None:
        # Get cavity indexes
        indexes = sorted(results["RESULTS"]["RESIDUES"].keys())
        # Include Interface Residues
        for index in indexes:
            self.residues_list.addItem(index)
        return

    def show_residues(self) -> None:
        from pymol import cmd

        # Get selected cavities from residues list
        cavs = [item.text() for item in self.residues_list.selectedItems()]

        # Clean objects
        cmd.set("auto_zoom", 0)
        cmd.delete("res")
        cmd.delete("residues")

        # Return if no cavity is selected
        if len(cavs) < 1:
            return

        # Get residues from cavities selected
        residues = []
        for cav in cavs:
            for residue in results["RESULTS"]["RESIDUES"][cav]:
                if residue not in residues:
                    residues.append(residue)

        # Check if input pdb is loaded
        control = 0
        for item in cmd.get_names("all"):
            if item == self.input_pdb:
                control = 1
        if control == 0:
            return

        # Select residues
        command = ""
        while len(residues) > 0:
            res, chain, _ = residues.pop(0)
            command = f"{command} (resid {res} and chain {chain}) or"
        command = f"obj {self.input_pdb} and ({command[:-3]})"
        cmd.select("res", command)

        # Create residues object
        cmd.create("residues", "res")
        cmd.delete("res")
        cmd.hide("everything", "residues")
        cmd.show("sticks", "residues")
        cmd.disable(self.cavity_pdb)
        cmd.enable(self.cavity_pdb)
        cmd.set("auto_zoom", 1)

    def show_cavities(self, list1, list2) -> None:
        from pymol import cmd

        # Get items from list1
        cavs = [item.text()[0:3] for item in list1.selectedItems()]

        # Select items of list2
        number_of_items = list1.count()
        for index in range(number_of_items):
            if list2.item(index).text()[0:3] in cavs:
                list2.item(index).setSelected(True)
            else:
                list2.item(index).setSelected(False)

        # Clean objects
        cmd.set("auto_zoom", 0)
        cmd.delete("cavs")
        cmd.delete("cavities")

        # Return if no cavity is selected
        if len(cavs) < 1:
            return

        # Check if cavity file is loaded
        control = 0
        for item in cmd.get_names("all"):
            if item == self.cavity_pdb:
                control = 1
        if control == 0:
            return

        # Color filling cavity points as blue nonbonded
        command = f"obj {self.cavity_pdb} and (resname "
        while len(cavs) > 0:
            command = f"{command}{cavs.pop(0)},"
        command = f"{command[:-1]})"
        cmd.select("cavs", command)

        # Create cavities object with blue nonbonded
        cmd.create("cavities", "cavs")
        cmd.delete("cavs")
        cmd.color("blue", "cavities")
        cmd.show("nonbonded", "cavities")

        # Color surface cavity points as red nb_spheres
        cmd.select("cavs", "cavities and name HS+HA")
        cmd.color("red", "cavs")
        cmd.show("nb_spheres", "cavs")
        cmd.delete("cavs")

        # Reset cavities output object
        cmd.disable(self.cavity_pdb)
        cmd.enable(self.cavity_pdb)
        for item in cmd.get_names("all"):
            if item == "hydropathy":
                cmd.disable("hydropathy")
                cmd.enable("hydropathy")
            if item == "depths":
                cmd.disable("depths")
                cmd.enable("depths")
        cmd.set("auto_zoom", 1)

    def show_depth(self, list1, list2) -> None:
        from pymol import cmd

        # Get items from list1
        cavs = [item.text()[0:3] for item in list1.selectedItems()]

        # Select items of list2
        number_of_items = list1.count()
        for index in range(number_of_items):
            if list2.item(index).text()[0:3] in cavs:
                list2.item(index).setSelected(True)
            else:
                list2.item(index).setSelected(False)

        # Clean objects
        cmd.set("auto_zoom", 0)
        cmd.delete("deps")
        cmd.delete("depths")

        # Return if no cavity is selected
        if len(cavs) < 1:
            return

        # Check if cavity file is loaded
        control = 0
        for item in cmd.get_names("all"):
            if item == self.cavity_pdb:
                control = 1
        if control == 0:
            return

        # Color filling cavity points as blue nonbonded
        command = f"obj {self.cavity_pdb} and (resname "
        while len(cavs) > 0:
            command = f"{command}{cavs.pop(0)},"
        command = f"{command[:-1]})"
        cmd.select("deps", command)

        # Create cavities object with blue nonbonded
        cmd.create("depths", "deps")
        cmd.delete("deps")
        cmd.spectrum("b", "rainbow", "depths")
        cmd.show("nb_spheres", "depths")

        # Reset cavities output object
        cmd.disable(self.cavity_pdb)
        for item in cmd.get_names("all"):
            if item == "cavities":
                cmd.disable("cavities")
                cmd.enable("cavities")
            if item == "depths":
                cmd.disable("hydropathy")
                cmd.enable("hydropathy")
        cmd.enable(self.cavity_pdb)
        cmd.set("auto_zoom", 1)

    def show_hydropathy(self, list1) -> None:
        from pymol import cmd

        # Get items from list1
        cavs = [item.text()[0:3] for item in list1.selectedItems()]

        # Clean objects
        cmd.set("auto_zoom", 0)
        cmd.delete("hyd")
        cmd.delete("hydropathy")

        # Return if no cavity is selected
        if len(cavs) < 1:
            return

        # Check if cavity file is loaded
        control = 0
        for item in cmd.get_names("all"):
            if item == self.cavity_pdb:
                control = 1
        if control == 0:
            return

        # Color filling cavity points as blue nonbonded
        command = f"obj {self.cavity_pdb} and (resname "
        while len(cavs) > 0:
            command = f"{command}{cavs.pop(0)},"
        command = f"{command[:-1]}) and (name HA+HS)"
        cmd.select("hyd", command)

        # Create cavities object with blue nonbonded
        cmd.create("hydropathy", "hyd")
        cmd.delete("hyd")
        cmd.spectrum("q", "yellow_white_blue", "hydropathy")
        cmd.show("nb_spheres", "hydropathy")

        # Reset cavities output object
        cmd.disable(self.cavity_pdb)
        for item in cmd.get_names("all"):
            if item == "cavities":
                cmd.disable("cavities")
                cmd.enable("cavities")
            if item == "depths":
                cmd.disable("depths")
                cmd.enable("depths")
        cmd.enable(self.cavity_pdb)
        cmd.set("auto_zoom", 1)

    def show_default_view(self) -> None:
        from pymol import cmd

        # Clean objects
        cmd.set("auto_zoom", 0)
        cmd.delete("view")

        # Check if cavity file is loaded
        control = 0
        for item in cmd.get_names("all"):
            if item == self.cavity_pdb:
                control = 1
        if control == 0:
            return

        # Color filling cavity points as blue nonbonded
        command = f"obj {self.cavity_pdb} and (name H+HA+HS)"
        command = f"{command[:-1]})"
        cmd.select("view", command)

        # Create cavities object with blue nonbonded
        cmd.hide("everything", self.cavity_pdb)
        cmd.show("nonbonded", "view")
        cmd.color("white", "view")
        cmd.delete("view")

    def show_depth_view(self) -> None:
        from pymol import cmd

        # Clean objects
        cmd.set("auto_zoom", 0)
        cmd.delete("view")

        # Check if cavity file is loaded
        control = 0
        for item in cmd.get_names("all"):
            if item == self.cavity_pdb:
                control = 1
        if control == 0:
            return

        # Color filling cavity points as blue nonbonded
        command = f"obj {self.cavity_pdb} and (name H+HA+HS)"
        command = f"{command[:-1]})"
        cmd.select("view", command)

        # Create cavities object with blue nonbonded
        cmd.hide("everything", self.cavity_pdb)
        cmd.show("nonbonded", "view")
        cmd.spectrum("b", "rainbow", "view")
        cmd.delete("view")

    def show_hydropathy_view(self) -> None:
        from pymol import cmd

        # Clean objects
        cmd.set("auto_zoom", 0)
        cmd.delete("view")

        # Check if cavity file is loaded
        control = 0
        for item in cmd.get_names("all"):
            if item == self.cavity_pdb:
                control = 1
        if control == 0:
            return

        # Color filling cavity points as blue nonbonded
        command = f"obj {self.cavity_pdb} and (name HA+HS)"
        command = f"{command[:-1]})"
        cmd.select("view", command)

        # Create cavities object with blue nonbonded
        cmd.hide("everything", self.cavity_pdb)
        cmd.show("nonbonded", "view")
        cmd.spectrum("q", "yellow_white_blue", "view")
        cmd.delete("view")

    def clean_results(self) -> None:
        # Input File
        self.input_file_entry.setText(f"")

        # Ligand File
        self.ligand_file_entry.setText(f"")

        # Cavities File
        self.cavities_file_entry.setText(f"")

        # Step Size
        self.step_size_entry.setText(f"")

        # Volume
        self.volume_list.clear()

        # Area
        self.area_list.clear()

        # Depth
        self.avg_depth_list.clear()
        self.max_depth_list.clear()

        # Hydropathy
        self.avg_hydropathy_list.clear()

        # Residues
        self.residues_list.clear()
