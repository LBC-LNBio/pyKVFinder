# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This is the source code for the graphical user interface of PyMOL pyKVFinder
Tools. The code is written in Python and uses the PyQt5 library for the
graphical user interface. Changes in this file are not advised, as it
controls all interactions with pyKVFinder.
"""

import os
import time

import tomlkit
from pymol import cmd
from PyQt5.QtWidgets import QCheckBox, QMainWindow, QMessageBox, QScrollBar
from PyQt5.uic import loadUi

import pyKVFinder

from .actions import (
    _load_area,
    _load_avg_depth,
    _load_avg_hydropathy,
    _load_cavity_file,
    _load_max_depth,
    _load_molecule_file,
    _load_residues,
    _load_volume,
    _refresh_list,
    _select_directory,
    _select_file,
)
from .box import Box
from .visualization import (
    _show_cavities,
    _show_depth,
    _show_descriptors,
    _show_hydropathy,
    _show_residues,
)


class PyMOLpyKVFinderTools(QMainWindow):
    """
    PyMOL pyKVFinder Tools. This class is the main window for the pyKVFinder
    plugin. It contains all the widgets and methods to interact with the
    plugin.
    """

    def __init__(self) -> None:
        super().__init__()
        # Initialize GUI
        self._initialize_gui()

        # Set case sensitive
        cmd.set("ignore_case", 0)

        # Restore default values
        self._startup()

        # Startup custom grid
        self.box = Box()

        # Results
        global results
        results = None
        self.pyKVFinderResults = None
        self.input_pdb = None
        self.ligand_pdb = None
        self.cavity_pdb = None

    def _initialize_gui(self) -> None:
        """
        Initialize the graphical user interface of the plugin.
        """
        # populate the QMainWindow from our *.ui file
        uifile = os.path.join(os.path.dirname(__file__), "PyMOL-pyKVFinder-Tools.ui")
        loadUi(uifile, self)

        # Dialog buttons binds
        self.button_run.clicked.connect(self.run)
        self.button_exit.clicked.connect(self.close)
        self.button_surface.clicked.connect(self.model_surface)
        self.button_restore.clicked.connect(self.restore)

        # ScrollBars binded to QListWidgets in Descriptors
        scroll_bar_volume = QScrollBar(self)
        self.volume_list.setVerticalScrollBar(scroll_bar_volume)
        scroll_bar_area = QScrollBar(self)
        self.area_list.setVerticalScrollBar(scroll_bar_area)
        scroll_bar_residues = QScrollBar(self)
        self.residues_list.setVerticalScrollBar(scroll_bar_residues)

        # Browse buttons
        self.button_browse_basedir.clicked.connect(
            lambda: _select_directory(self.basedir)
        )
        self.button_browse_vdw.clicked.connect(
            lambda: _select_file("Choose van der Waals radii dictionary", "*", self.vdw)
        )
        self.button_browse_results.clicked.connect(
            lambda: _select_file(
                "Choose KVFinder Results File",
                "KVFinder Results File (*.toml);;All files (*)",
                self.results_file_entry,
            )
        )

        # Surface representation bindings in Main Tab
        self.surface.currentIndexChanged.connect(self._surface_selection)

        # Box Adjustment
        self.button_draw_box.clicked.connect(
            lambda: self.box._set_box(
                self.button_draw_box,
                self.button_redraw_box,
                self.min_x,
                self.max_x,
                self.min_y,
                self.max_y,
                self.min_z,
                self.max_z,
                self.angle1,
                self.angle2,
                self.padding,
            )
        )
        self.button_delete_box.clicked.connect(
            lambda: self.box._delete_box(
                self.button_draw_box,
                self.button_redraw_box,
                self.min_x,
                self.max_x,
                self.min_y,
                self.max_y,
                self.min_z,
                self.max_z,
                self.angle1,
                self.angle2,
            )
        )
        self.button_redraw_box.clicked.connect(
            lambda: self.box._redraw_box(
                self.min_x,
                self.max_x,
                self.min_y,
                self.max_y,
                self.min_z,
                self.max_z,
                self.angle1,
                self.angle2,
                self.padding,
            )
        )
        self.button_box_adjustment_help.clicked.connect(lambda: self.box._help(self))

        # Refresh the list of objects in the input and ligand widgets
        self.refresh_input.clicked.connect(lambda: _refresh_list(self.input))
        self.refresh_ligand.clicked.connect(lambda: _refresh_list(self.ligand))

        # Visualization in results tab
        self.button_load_results.clicked.connect(self._load_run)
        self.volume_list.itemSelectionChanged.connect(
            lambda: _show_cavities(self.volume_list, self.area_list, self.cavity_pdb)
        )
        self.area_list.itemSelectionChanged.connect(
            lambda: _show_cavities(self.area_list, self.volume_list, self.cavity_pdb)
        )
        self.residues_list.itemSelectionChanged.connect(
            lambda: _show_residues(results, self.residues_list, self.input_pdb)
        )
        self.avg_depth_list.itemSelectionChanged.connect(
            lambda: _show_depth(
                self.avg_depth_list, self.max_depth_list, self.cavity_pdb
            )
        )
        self.max_depth_list.itemSelectionChanged.connect(
            lambda: _show_depth(
                self.max_depth_list, self.avg_depth_list, self.cavity_pdb
            )
        )
        self.avg_hydropathy_list.itemSelectionChanged.connect(
            lambda: _show_hydropathy(results, self.avg_hydropathy_list, self.cavity_pdb)
        )
        self.default_view.toggled.connect(
            lambda: _show_descriptors(self.cavity_pdb, "default")
        )
        self.depth_view.toggled.connect(
            lambda: _show_descriptors(self.cavity_pdb, "depth")
        )
        self.hydropathy_view.toggled.connect(
            lambda: _show_descriptors(self.cavity_pdb, "hydropathy", results)
        )

    def _surface_selection(self) -> None:
        """
        Get the surface selection from the GUI. This method is called when the
        user selects a surface type in the GUI.
        """
        # Get the surface type
        surface = self.surface.currentText()[:-1][-3:]
        if surface == "vdW":
            self.probe_in.setValue(0.0)
            self.probe_in.setEnabled(False)
        elif surface == "SES":
            self.probe_in.setValue(1.4)
            self.probe_in.setEnabled(True)
        elif surface == "SAS":
            self.probe_in.setValue(1.4)
            self.probe_in.setEnabled(True)

    def _restore_results(self, remove_inputs: bool) -> None:
        """
        Restore the results tab to the default values. This method is called
        when the user clicks the 'Restore Defaults' button in the GUI.
        """
        # Remove cavities, residues and pdbs (input, ligand, cavity)
        cmd.delete("cavities")
        cmd.delete("residues")
        if self.input_pdb and remove_inputs:
            cmd.delete(self.input_pdb)
        if self.ligand_pdb and remove_inputs:
            cmd.delete(self.ligand_pdb)
        if self.cavity_pdb:
            cmd.delete(self.cavity_pdb)
        results = self.input_pdb = self.ligand_pdb = self.cavity_pdb = None
        cmd.frame(1)

        # Clean results
        self._clean_results_tab()
        self.results_file_entry.clear()

    def _startup(self) -> None:
        """
        Set the default values for the widgets in the GUI. This method is
        called when the plugin is started.
        """
        # Startup Main tab
        _refresh_list(self.input)
        self.step_size.setValue(0.6)
        self.probe_in.setValue(1.4)
        self.probe_out.setValue(4.0)
        self.removal_distance.setValue(2.4)
        self.volume_cutoff.setValue(5.0)
        self.surface.setCurrentIndex(1)  # 0: vdW, 1: SES, 2: SAS
        self.basename.setText("output")
        self.basedir.setText(os.getcwd())
        self.vdw.setText(
            os.path.join(os.path.dirname(pyKVFinder.__file__), "data", "vdw.dat")
        )

        # Startup Box Adjustment tab
        self.box_adjustment.setChecked(False)
        self.padding.setValue(3.5)
        # self._delete_box() # TODO: Implement this method

        # Startup Ligand Adjustment tab
        self.ligand_adjustment.setChecked(False)
        self.ligand.clear()
        self.ligand_cutoff.setValue(5.0)

    def closeEvent(self, event) -> None:
        """
        Exit the plugin. This method is called when the user clicks the 'Exit'
        button in the GUI.
        """
        global dialog
        dialog = None

    def restore(self) -> None:
        """
        Restore the default values for the widgets in the GUI. This method is
        called when the user clicks the 'Restore Defaults' button in the GUI.
        """
        # Restore to startup values
        self._startup()

        # Restore results tab
        reply = QMessageBox(self)
        reply.setText("Also restore Results Visualization tab?")
        reply.setWindowTitle("Restore Values")
        reply.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        reply.setIcon(QMessageBox.Information)
        reply.checkbox = QCheckBox("Also remove input and ligand PDBs?")
        reply.layout = reply.layout()
        reply.layout.addWidget(reply.checkbox, 1, 2)
        if reply.exec_() == QMessageBox.Yes:
            self._restore_results(reply.checkbox.isChecked())

    def _get_parameters(self) -> dict:
        """
        Get the parameters from the GUI. This method is called when the user
        clicks the 'Run pyKVFinder' button in the GUI. It will return a
        dictionary with the parameters to run pyKVFinder.

        Returns
        -------
        dict
            The parameters to run pyKVFinder.
        """
        # Get the parameters from the GUI
        parameters = {
            "input": os.path.join(
                self.basedir.text(),
                self.basename.text(),
                f"{self.input.currentText()}.pdb",
            ),
            "ligand": (
                os.path.join(
                    self.basedir.text(),
                    self.basename.text(),
                    f"{self.ligand.currentText()}.pdb",
                )
                if self.ligand.currentText() != ""
                else None
            ),
            "step": self.step_size.value(),
            "probe_in": self.probe_in.value(),
            "probe_out": self.probe_out.value(),
            "removal_distance": self.removal_distance.value(),
            "volume_cutoff": self.volume_cutoff.value(),
            "surface": (
                "SES"
                if self.surface.currentText()[:-1][-3:] == "vdW"
                else self.surface.currentText()[:-1][-3:]
            ),
            "basename": self.basename.text(),
            "basedir": os.path.join(self.basedir.text(), self.basename.text()),
            "vdw": self.vdw.text(),
            "box": (
                self.box._create_box_file(
                    os.path.join(
                        self.basedir.text(),
                        self.basename.text(),
                        f"box.{self.basename.text()}.toml",
                    )
                )
                if self.box_adjustment.isChecked()
                else None
            ),
            "padding": self.padding.value(),
            "ligand_adjustment": self.ligand_adjustment.isChecked(),
            "ligand_cutoff": self.ligand_cutoff.value(),
            "hydrophobicity_scale": self.hscale.currentText().replace(" & ", ""),
            "ignore_backbone": self.ignore_backbone.isChecked(),
        }

        return parameters

    def _prepare_run(self, parameters: dict) -> bool:
        """
        Prepare the plugin to run pyKVFinder. This method is called when the
        user clicks the 'Run pyKVFinder' button in the GUI. It will check if
        the input files are valid and if the output directory is writable.

        Parameters
        ----------
        parameters : dict
            The parameters to run pyKVFinder.
        """
        # Create basedir if it does not exist
        os.makedirs(parameters["basedir"], exist_ok=True)

        # Save input pdb file
        if self.input.currentText() != "":
            for obj in cmd.get_names("all"):
                if obj == self.input.currentText():
                    cmd.save(parameters["input"], self.input.currentText(), 0, "pdb")
        else:
            QMessageBox.critical(
                self, "Error", "Please select an input.", QMessageBox.Ok
            )
            return False

        # Save ligand pdb file
        if self.ligand_adjustment.isChecked():
            if self.ligand.currentText() != "":
                for obj in cmd.get_names("all"):
                    if obj == self.ligand.currentText():
                        cmd.save(
                            parameters["ligand"], self.ligand.currentText(), 0, "pdb"
                        )
            else:
                QMessageBox.critical(
                    self, "Error", "Please select a ligand.", QMessageBox.Ok
                )
                return False

        # Check box adjustment mode
        if self.box_adjustment.isChecked():
            if "box" not in cmd.get_names("all"):
                QMessageBox.critical(self, "Error", "Draw a box in PyMOL!")
                return False

        # Remove None values from parameters
        parameters = {k: v for k, v in parameters.items() if v is not None}

        # Save parameters file
        with open(
            os.path.join(
                parameters["basedir"], f'{parameters["basename"]}.parameters.toml'
            ),
            "w",
        ) as f:
            tomlkit.dump(parameters, f)

        return True

    def _load_run(self) -> None:
        """
        Load the results of pyKVFinder into PyMOL. This method is called after
        running pyKVFinder. It will load the cavities and the results into
        PyMOL.
        """
        # Get results file
        results_file = self.results_file_entry.text()

        # Check if results file exist
        if os.path.exists(results_file) and results_file.endswith(".toml"):
            print(f"> Loading results from: {results_file}")
        else:
            QMessageBox.critical(
                self, "Error", "Results file cannot be opened! Check results file path."
            )

        # Create global variable for results
        global results

        # Read results (Ubuntu/macOS)
        with open(results_file, "r") as f:
            results = tomlkit.load(f)

        # Clean results tab
        self._clean_results_tab()

        # Load file information
        self._load_file_information()

        # Load characterization (volume, area, depth, residues and hydropathy)
        _load_volume(results, self.volume_list)
        _load_area(results, self.area_list)
        _load_avg_depth(results, self.avg_depth_list)
        _load_max_depth(results, self.max_depth_list)
        _load_avg_hydropathy(results, self.avg_hydropathy_list)
        _load_residues(results, self.residues_list)

        # Set default view in results
        self.default_view.setChecked(True)

        # Remove previous objects from PyMOL
        objects = cmd.get_names("all")
        for obj in ["cavities", "residues", "depth", "hydropathy"]:
            if obj in objects:
                cmd.delete(obj)
        cmd.frame(1)

        # Load input file
        if "INPUT" in results["FILES"].keys():
            self.input_pdb = os.path.basename(results["FILES"]["INPUT"]).replace(
                ".pdb", ""
            )
            _load_molecule_file(results["FILES"]["INPUT"], self.input_pdb)
        else:
            self.input_pdb = None

        # Load ligand file
        if "LIGAND" in results["FILES"].keys():
            self.ligand_pdb = os.path.basename(results["FILES"]["LIGAND"]).replace(
                ".pdb", ""
            )
            _load_molecule_file(results["FILES"]["LIGAND"], self.ligand_pdb)
        else:
            self.ligand_pdb = None

        # Load cavities file
        self.cavity_pdb = os.path.basename(results["FILES"]["OUTPUT"]).replace(
            ".pdb", ""
        )
        _load_cavity_file(results["FILES"]["OUTPUT"], self.cavity_pdb)

    def _load_file_information(self) -> None:
        """
        Load the file information into the results tab. This method is called
        after running pyKVFinder. It will load the file information into the
        results tab.
        """
        if "INPUT" in results["FILES"].keys():
            self.input_file_entry.setText(f"{results['FILES']['INPUT']}")
        else:
            self.input_file_entry.setText(f"")

        # Ligand File
        if "LIGAND" in results["FILES"].keys():
            self.ligand_file_entry.setText(f"{results['FILES']['LIGAND']}")
        else:
            self.ligand_file_entry.setText(f"")

        # Cavities File
        self.cavities_file_entry.setText(f"{results['FILES']['OUTPUT']}")

        # Step Size
        if "PARAMETERS" in results.keys():
            if "STEP" in results["PARAMETERS"].keys():
                self.step_size_entry.setText(f"{results['PARAMETERS']['STEP']:.2f}")

    def _clean_results_tab(self) -> None:
        """
        Clean the results tab. This method is called before loading the results
        of pyKVFinder into PyMOL. It will clear the lists in the results tab.
        """
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

    def model_surface(self) -> None:
        """
        Model the surface of the input. This method is called when the user
        clicks the 'Model Surface' button in the GUI. It will model the surface
        of the input structure in PyMOL.
        """
        # Get parameters from the GUI
        parameters = self._get_parameters()

        # Save input pdb file
        if self.input.currentText() != "":
            for obj in cmd.get_names("all"):
                if obj == self.input.currentText():
                    cmd.save(parameters["input"], self.input.currentText(), 0, "pdb")
        else:
            QMessageBox.critical(
                self, "Error", "Please select an input.", QMessageBox.Ok
            )
            return False

        # Load molecule
        surface = pyKVFinder.Molecule(
            molecule=parameters["input"], radii=parameters["vdw"]
        )

        # Model the surface
        if parameters["surface"] == "vdW":
            surface.vdw(step=parameters["step"])
        elif parameters["surface"] == "SES":
            surface.surface(
                step=parameters["step"],
                probe=parameters["probe_in"],
                surface=parameters["surface"],
            )
        elif parameters["surface"] == "SAS":
            surface.surface(
                step=parameters["step"],
                probe=parameters["probe_in"],
                surface=parameters["surface"],
            )

        # Save surface file
        pymolname = f"{parameters['basename']}.KVFinder.surface"
        filename = os.path.join(parameters["basedir"], f"{pymolname}.pdb")
        surface.export(filename)

        # Remove previous results in objects with same cavity name
        for obj in cmd.get_names("all"):
            if pymolname == obj:
                cmd.delete(obj)

        # Load the cavity file
        cmd.load(filename, pymolname, zoom=0)

        # Hide and show to update the display
        cmd.hide("everything", pymolname)
        cmd.show("sphere", pymolname)
        cmd.alter(pymolname, f"vdw={parameters['step'] / 2}")
        cmd.rebuild()

    def run(self) -> bool:
        """
        Run pyKVFinder with the current parameters. This method is called when
        the user clicks the 'Run pyKVFinder' button in the GUI. After
        completion, it will display the results in the GUI.
        """
        # Get parameters from the GUI
        parameters = self._get_parameters()

        # Prepare parameters to run pyKVFinder
        if not self._prepare_run(parameters):
            QMessageBox.critical(
                self,
                "Error!",
                "An error occurred while preparing the run. Please check the \
parameters.",
                QMessageBox.Ok,
            )
            return False

        # Print the current status for users
        print(f"[==> Running pyKVFinder for {self.input.currentText()} ...")
        start_time = time.time()

        # Run pyKVFinder
        self.pyKVFinderResults = pyKVFinder.run_workflow(
            input=parameters["input"],
            ligand=parameters["ligand"],
            vdw=parameters["vdw"],
            box=parameters["box"],
            step=parameters["step"],
            probe_in=parameters["probe_in"],
            probe_out=parameters["probe_out"],
            removal_distance=parameters["removal_distance"],
            volume_cutoff=parameters["volume_cutoff"],
            ligand_cutoff=parameters["ligand_cutoff"],
            include_depth=True,
            include_hydropathy=True,
            hydrophobicity_scale=parameters["hydrophobicity_scale"],
            surface=parameters["surface"],
            ignore_backbone=parameters["ignore_backbone"],
        )

        # Export results and cavities to file
        results_file = os.path.join(
            parameters["basedir"],
            f"{parameters['basename']}.KVFinder.results.toml",
        )
        cavity_file = os.path.join(
            parameters["basedir"],
            f"{parameters['basename']}.KVFinder.output.pdb",
        )

        # If pyKVFinder found cavities, export results
        if self.pyKVFinderResults is not None:
            self.pyKVFinderResults.export_all(
                fn=results_file,
                output=cavity_file,
                include_frequencies_pdf=False,
            )

        # Elapsed time
        elapsed_time = time.time() - start_time
        print(
            f"> Cavities detected: {self.pyKVFinderResults.ncav if self.pyKVFinderResults else 0}"
        )
        print(f"> Elapsed time: {elapsed_time:.2f} seconds")

        # Sucessful run
        if self.pyKVFinderResults is not None:
            self.results_file_entry.setText(results_file)
            self.tabs.setCurrentIndex(2)
            self._load_run()
        # Unsucessful run
        else:
            self._restore_results(remove_inputs=False)
            QMessageBox.warning(
                self,
                "Warning",
                "No cavities were detected. Please check the parameters.",
                QMessageBox.Ok,
            )
            return False

        return True
