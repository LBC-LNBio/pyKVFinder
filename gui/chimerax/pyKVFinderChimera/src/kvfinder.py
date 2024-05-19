 #vim: set expandtab shiftwidth=4 softtabstop=4:

# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies,
# including partial copies, of the software or any revisions
# or derivations thereof.
# === UCSF ChimeraX Copyright ===

from chimerax.core.objects import all_objects
from chimerax.core.tools import ToolInstance
from chimerax.atomic import StructureSeq, Structure, selected_atoms, all_atoms, structure_atoms, all_atomic_structures
from chimerax.core.commands import run
from chimerax.std_commands import style
from chimerax.atomic.structure import AtomicStructure, Model
from os.path import expanduser
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QDialog
from PyQt5.QtCore import QThread, pyqtSlot, pyqtSignal
import chimerax as cx
import os
import pyKVFinder
import numpy as np
import sys
import toml

dialog = None

class _Default(object):
    def __init__(self):
        super(_Default, self).__init__()
        #######################
        ### Main Parameters ###
        #######################
        self.step = 0.6
        self.resolution = "Low"
        self.probe_in = 1.4
        self.probe_out = 4.0
        self.removal_distance = 2.4
        self.volume_cutoff = 5.0
        self.surface = "Solvent Excluded Surface (SES)"
        # self.cavity_representation = "Filtered"
        self.base_name = "output"
        if os.name == "posix":
            self.output_dir_path = expanduser('~/')
        else:
            self.output_dir_path = expanduser('~\\')
        self.region_option = "Default"
        #######################
        ### File Locations  ###
        #######################
        self.parKVFinder = None
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


class KVFinder(ToolInstance):

    SESSION_ENDURING = False    # Does this instance persist when session closes
    SESSION_SAVE = True         # We do save/restore in sessions
    #help = "help:user/tools/tutorial.html"
                                # Let ChimeraX know about our help page
    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)

        _translate = QtCore.QCoreApplication.translate
        self.display_name = "KVFinder"
        self._default = _Default()
        self.region_option = self._default.region_option

        self.app = QtWidgets.QApplication(sys.argv)
        self.tool_window = QtWidgets.QMainWindow()
        self.tool_window.fill_context_menu = self.fill_context_menu
        self.tool_window.setWindowTitle("pyKVFinder ChimeraX Tool (Qt)")
        self.ui = Ui_pyKVFinder()
        self.ui.setupUi(self.tool_window)  

        self.ui.debug = SampleGUI(self)
        self.ui.debug.setObjectName("debug")
        self.ui.tabs.addTab(self.ui.debug, "")

        self.ui.tabs.setTabText(self.ui.tabs.indexOf(self.ui.debug), _translate("pyKVFinder", "Debug"))
        # self.tool_window.adjustSize()
        self.tool_window.show()

        # Set box centers
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0

        # Results
        self.results = None
        self.input_pdb = None
        self.ligand_pdb = None
        self.cavity_pdb = None

        # Lists Selected
        self.vs_selected = []
        self.am_selected = []
        self.hyd_selected = []
        self.res_selected = []



        self._connect_ui()

        # Restore Default Parameters
        self.restore(is_startup=True)

    def _connect_ui(self):
        
        # ScrollBars binded to QListWidgets in Descriptors
        scroll_bar_volume = QtWidgets.QScrollBar(self.tool_window)
        self.ui.volume_list.setVerticalScrollBar(scroll_bar_volume)
        scroll_bar_area = QtWidgets.QScrollBar(self.tool_window)
        self.ui.area_list.setVerticalScrollBar(scroll_bar_area)
        scroll_bar_residues = QtWidgets.QScrollBar(self.tool_window)
        self.ui.residues_list.setVerticalScrollBar(scroll_bar_residues)

        ########################
        ### Buttons Callback ###
        ########################

        # hook up QMainWindow buttons callbacks
        self.ui.button_run.clicked.connect(self.run)
        # ui.button_exit.clicked.connect(tw.close)
        self.ui.button_restore.clicked.connect(self.restore)
        # ui.button_grid.clicked.connect(self.show_grid)
        self.ui.button_save_parameters.clicked.connect(self.save_parameters)

        # hook up Refresh buttons callback
        self.ui.refresh_input.clicked.connect(lambda: self.refresh(self.ui.input))

        # hook up resolution-step CheckBox callbacks
        # self.ui.resolution_label.clicked.connect(self.check_resolution)


        self.ui.groupButton.buttonClicked.connect(self._optionCheck)

        # hook up Browse buttons callbacks
        self.ui.button_browse.clicked.connect(self.select_directory)

        self.ui.button_browse3.clicked.connect(
            lambda: self.select_file(
                "Choose van der Waals radii dictionary", self.ui.dictionary, "*"
            )
        )

        self.ui.button_browse4.clicked.connect(
            lambda: self.select_file(
                "Choose KVFinder Results File",
                self.ui.results_file_entry,
                "KVFinder Results File (*.toml);;All files (*)",
            )
        )

        self.ui.button_box_adjustment_help.clicked.connect(self.box_adjustment_help)
        self.ui.button_exit.clicked.connect(self.tool_window.close)
        self.ui.button_load_results.clicked.connect(self.load_results)
    

        self.ui.volume_list.itemSelectionChanged.connect(
            lambda list1=self.ui.volume_list, list2=self.ui.area_list: self.show_cavities(
                list1, list2
            )
        )
        self.ui.area_list.itemSelectionChanged.connect(
            lambda list1=self.ui.area_list, list2=self.ui.volume_list: self.show_cavities(
                list1, list2
            )
        )

        self.ui.avg_depth_list.itemSelectionChanged.connect(
            lambda list1=self.ui.avg_depth_list, list2=self.ui.max_depth_list: self.show_depth(
                list1, list2
            )
        )
        self.ui.max_depth_list.itemSelectionChanged.connect(
            lambda list1=self.ui.max_depth_list, list2=self.ui.avg_depth_list: self.show_depth(
                list1, list2
            )
        )

        self.ui.avg_hydropathy_list.itemSelectionChanged.connect(
            lambda list1=self.ui.avg_hydropathy_list: self.show_hydropathy(list1)
        )

        self.ui.residues_list.itemSelectionChanged.connect(self.show_residues)
        self.ui.default_view.toggled.connect(self.show_default_view)
        self.ui.depth_view.toggled.connect(self.show_depth_view)
        self.ui.hydropathy_view.toggled.connect(self.show_hydropathy_view)


        # hook up Search Space button callbacks
        # Box Adjustment
        self.ui.button_draw_box.clicked.connect(self.set_box)
        self.ui.button_delete_box.clicked.connect(self.delete_box)
        self.ui.button_redraw_box.clicked.connect(self.redraw_box)
        
        # Grid
        self.ui.button_grid.clicked.connect(self.show_grid)
        
        # Ligand Adjustment
        self.ui.refresh_ligand.clicked.connect(lambda: self.refresh(self.ui.ligand))

    def _optionCheck(self, btn):

        if btn.isChecked():
            # self.session.logger.info(f'You selected {btn.text()}')
            self.region_option = btn.text()

    def restore(self, is_startup=False) -> None:
        """
        Callback for the "Restore Default Values" button
        """
        #from pymol import cmd
        from PyQt5.QtWidgets import QMessageBox, QCheckBox

        # Restore Results Tab
        if not is_startup:
            reply = QMessageBox(self.tool_window)
            reply.setText("Also restore Results Visualization tab?")
            reply.setWindowTitle("Restore Values")
            reply.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
            reply.setIcon(QMessageBox.Information)
            reply.checkbox = QCheckBox("Also remove input and ligand PDBs?")
            reply.layout = reply.layout()
            reply.layout.addWidget(reply.checkbox, 1, 2)
            if reply.exec_() == QMessageBox.Yes:
                # Remove cavities, residues and pdbs (input, ligand, cavity)

                if self.input_pdb and reply.checkbox.isChecked():
                    input_pdb = self._get_model(self.input_pdb)
                    input_pdb.delete()
                if self.cavity_pdb:
                    cavity_pdb = self._get_model(self.cavity_pdb)
                    cavity_pdb.delete()
                results = self.input_pdb = self.ligand_pdb = self.cavity_pdb = None

                # Clean results
                self.clean_results()
                self.ui.results_file_entry.clear()

        # Restore PDB and ligand input
        self.refresh(self.ui.input)
        self.refresh(self.ui.ligand)

        # Delete grid
        #cmd.delete("grid")

        ### Main tab ###
        self.ui.step_size.setValue(self._default.step)
        self.ui.step_size.setEnabled(True)

        # self.ui.resolution_label.setChecked(True)
        # self.ui.resolution.setCurrentText(self._default.resolution)
        # self.ui.resolution.setEnabled(True)

        self.ui.base_name.setText(self._default.base_name)
        self.ui.probe_in.setValue(self._default.probe_in)
        self.ui.probe_out.setValue(self._default.probe_out)
        self.ui.volume_cutoff.setValue(self._default.volume_cutoff)
        self.ui.removal_distance.setValue(self._default.removal_distance)
        self.ui.surface.setCurrentText(self._default.surface)
        # self.ui.cavity_representation.setCurrentText(self._default.cavity_representation)
        self.ui.output_dir_path.setText(self._default.output_dir_path)
        # self.ui.parKVFinder.setText(self._default.parKVFinder)
        self.ui.dictionary.setText(self._default.dictionary)

        ### Search Space Tab ###
        # Box Adjustment
        self.ui.box_adjustment.setChecked(self._default.box_adjustment)
        self.ui.padding.setValue(self._default.padding)
        self.delete_box()
        # Ligand Adjustment
        self.ui.ligand_adjustment.setChecked(self._default.ligand_adjustment)
        self.ui.ligand.clear()
        self.ui.ligand_cutoff.setValue(self._default.ligand_cutoff)

    def refresh(self, combo_box ) -> None:
        """
        Callback for the "Refresh" button
        """
        combo_box.clear()

        if combo_box == self.ui.input:
            models = all_objects(self.session).models
            
            for model in models:
                name = model.name + " " + model.id_string
                combo_box.addItem(name)
        
        if combo_box == self.ui.ligand:

            models = all_objects(self.session).models
            
            for model in models:
                name = model.name + " " + model.id_string
                combo_box.addItem(name)
        
        return    
    
    def _run_pyKVFinder(self, atomic, box_adjustment = False):

        import time
        print(
            f"\n[==> Running pyKVFinder for: {os.path.join(self.ui.output_dir_path.text(), 'KV_Files', self.ui.base_name.text(), f'{self.ui.input.currentText()}')}"
        )
        start = time.time()
        probe_out = self.ui.probe_out.value()
        probe_in = self.ui.probe_in.value()
        removal_distance = self.ui.removal_distance.value()
        volume_cutoff = self.ui.volume_cutoff.value()
        step = self.ui.step_size.value()
        ignore_backbone = True if self.ui.ignore_backbone_checkbox.isChecked() else False
        surface =  'SES' if self.ui.surface.currentText() == 'Solvent Excluded Surface (SES)' else 'SAS'
        if self.ui.ligand_adjustment.isChecked() and self.ui.ligand.currentData() != self.ui.input.currentData() and len(self.ui.ligand.currentData()) > 0:
            ligand_model = self._get_model(self.ui.ligand.currentData())
            
            ligand = self.extract_pdb_session(list_models=ligand_model, selected=False)
            ligand_cutoff = self.ui.ligand_cutoff.value()
            print("Gerei o ligante")
        else:
            ligand = None
            ligand_cutoff = 5
        if not box_adjustment:
            vertices = pyKVFinder.get_vertices(atomic, probe_out=probe_out, step=step)
            ncavs, cavities = pyKVFinder.detect(atomic, vertices, step=step, latomic=ligand, ligand_cutoff = ligand_cutoff, probe_in=probe_in, probe_out=probe_out, removal_distance=removal_distance, volume_cutoff=volume_cutoff, surface=surface)
        else:
            fn = os.path.join(self.ui.output_dir_path.text(), 'KV_Files', self.ui.base_name.text(), "parameters.toml")
            
            print(f"fn: {fn}")
            print(f"atomic: {atomic}")
            vertices, atomic = pyKVFinder.get_vertices_from_file(fn, atomic, step=step, probe_in=probe_in, probe_out=probe_out)
            ncavs, cavities = pyKVFinder.detect(atomic, vertices, step=step, latomic=ligand, ligand_cutoff = ligand_cutoff, probe_in=probe_in, probe_out=probe_out, removal_distance=removal_distance, volume_cutoff=volume_cutoff, box_adjustment=True, surface=surface)
        elapsed_time = time.time() - start
        print(f"> Cavities detected: {ncavs}")
        print(f"> Elapsed time: {elapsed_time:.2f} seconds")

                    # Load successfull run
        self.ui.results_file_entry.setText(
            f"{os.path.join(self.ui.output_dir_path.text(), 'KV_Files', self.ui.base_name.text(), f'{self.ui.base_name.text()}.KVFinder.results.toml')}"
        )

        if ncavs > 0:
            self.ui.tabs.setCurrentIndex(2)
            surface, volume, area, residues, scales, avg_hydropathy, depths, max_depth, avg_depth, frequencies = self.characterization(cavities=cavities, step=step, atomic=atomic, vertices=vertices, probe_in=probe_in, ignore_backbone=ignore_backbone)
            
            if os.path.exists(
                f"{os.path.join(self.ui.output_dir_path.text(), 'KV_Files', self.ui.base_name.text(), f'{self.ui.base_name.text()}.KVFinder.results.toml')}"
            ):
                os.remove(
                    f"{os.path.join(self.ui.output_dir_path.text(), 'KV_Files', self.ui.base_name.text(), f'{self.ui.base_name.text()}.KVFinder.results.toml')}"
                )

            output_cavity = os.path.join(self.ui.output_dir_path.text(), 'KV_Files', self.ui.base_name.text(), f'{self.ui.base_name.text()}.KVFInder.cavity')
            pyKVFinder.export(output_cavity, cavities, surface, vertices, step=step, B=depths, Q=scales)

            pdb = os.path.join(self.ui.output_dir_path.text(), 'KV_Files', self.ui.base_name.text(), f'{self.ui.base_name.text()}.KVFinder.output')

            output_results = os.path.join(self.ui.output_dir_path.text(), 'KV_Files', self.ui.base_name.text(),f'results.toml')
            pyKVFinder.write_results(output_results, ligand=None, input=pdb, output=output_cavity, volume=volume, area=area, max_depth=max_depth, avg_depth=avg_depth, avg_hydropathy=avg_hydropathy, residues=residues, frequencies=frequencies, step=step)
            
            os.rename(
                output_results,
                f"{os.path.join(self.ui.output_dir_path.text(), 'KV_Files', self.ui.base_name.text(), f'{self.ui.base_name.text()}.KVFinder.results.toml')}"
            )


            self.load_results()
        elif ncavs == 0:
            QtWidgets.QMessageBox.warning(self.tool_window, "Warning!", "No cavities found!")
    
    # Modificando
    def run(self):
        """This function calls extract_pdb_session to obtain the desired atomic structure, 
        and then invokes _run_pyKVFinder to identify the cavities.

        """

        if self.save_parameters():
            
            models = self._get_model(self.ui.input.currentData())
            spec = ",".join([model.atomspec for model in models])
            if self.ui.box_adjustment.isChecked():
                atomic = self.extract_pdb_session(selected=False, list_models=models)
                self._run_pyKVFinder(atomic, box_adjustment = True)
            elif self.region_option == "Default":
                atomic = self.extract_pdb_session(selected=False, list_models=models)
                self._run_pyKVFinder(atomic)
            elif self.region_option == "Selected":
                atomic = self.extract_pdb_session(selected=True, list_models=models)
                self._run_pyKVFinder(atomic)              
            elif self.region_option == "Protein":
                run(self.session, f"sel {spec} & protein")
                atomic = self.extract_pdb_session(selected=True, list_models=models)
                self._run_pyKVFinder(atomic)    
            elif self.region_option == "All ligands without solvent":
                run(self.session, f"sel {spec} & ~solvent")
                atomic = self.extract_pdb_session(selected=True, list_models=models)
                self._run_pyKVFinder(atomic)    
            else:
                QtWidgets.QMessageBox.critical(
                    self.tool_window, "Error!", "An error occurred during cavity detection!"
                )             

        else:
            from PyQt5.QtWidgets import QMessageBox

            QMessageBox.critical(
                self.tool_window,
                "Error",
                "An error occurred while creating the parameters file! Check the parKVFinder parameters!",
            )

    
    def run_old(self):
        """This function calls extract_pdb_session to obtain the desired atomic structure, 
        and then invokes _run_pyKVFinder to identify the cavities.

        """
        if self.save_parameters():
            model = self._get_model(self.ui.input.currentData())
            spec = model.atomspec
            if self.ui.box_adjustment.isChecked():
                atomic = self.extract_pdb_session(selected=False, name=self.ui.input.currentText())
                self._run_pyKVFinder(atomic, box_adjustment = True)
            elif self.region_option == "Default":
                atomic = self.extract_pdb_session(selected=False, name=self.ui.input.currentText())
                self._run_pyKVFinder(atomic)
            elif self.region_option == "Selected":
                atomic = self.extract_pdb_session(selected=True, name=self.ui.input.currentText())
                self._run_pyKVFinder(atomic)              
            elif self.region_option == "Protein":
                run(self.session, f"sel {spec} & protein")
                atomic = self.extract_pdb_session(selected=True, name=self.ui.input.currentText())
                self._run_pyKVFinder(atomic)    
            elif self.region_option == "All ligands without solvent":
                run(self.session, f"sel {spec} & ~solvent")
                atomic = self.extract_pdb_session(selected=True, name=self.ui.input.currentText())
                self._run_pyKVFinder(atomic)    
            else:

                QtWidgets.QMessageBox.critical(
                    self.tool_window, "Error!", "An error occurred during cavity detection!"
                )             

        else:
            from PyQt5.QtWidgets import QMessageBox

            QMessageBox.critical(
                self.tool_window,
                "Error",
                "An error occurred while creating the parameters file! Check the parKVFinder parameters!",
            )

    def load_results(self) -> None:

        # Get results file
        results_file = self.ui.results_file_entry.text()

        # Check if results file exist
        if os.path.exists(results_file) and results_file.endswith(".toml"):
            print(f"> Loading results from: {self.ui.results_file_entry.text()}")
        else:
            from PyQt5.QtWidgets import QMessageBox

            error_msg = QMessageBox.critical(
                self.tool_window, "Error", "Results file cannot be opened! Check results file path."
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
                self.tool_window,
                "Error",
                "Results file has incorrect format! Please check your file.",
            )
            return False

        # # Clean results
        self.clean_results()

        # # Refresh information
        self.refresh_information()

        # # Refresh volume
        self.refresh_volume()

        # # Refresh area
        self.refresh_area()

        # # Refresh depth
        self.refresh_avg_depth()
        self.refresh_max_depth()
        self.refresh_avg_hydropathy()

        # # Refresh residues
        self.refresh_residues()

        # # Set default view in results
        self.ui.default_view.setChecked(True)


        # # Load files as PyMOL objects
        # cmd.delete("cavities")
        # cmd.delete("residues")
        # cmd.frame(1)

        # Load input
        models_loaded = [model.name for model in self.session.models]
        if "INPUT" in results["FILES_PATH"].keys():
            input_fn = results["FILES_PATH"]["INPUT"]
            self.input_pdb = os.path.basename(input_fn)
            if os.path.split(input_fn)[-1] not in models_loaded:
                self.load_file(input_fn, self.input_pdb)
        else:
            self.input_pdb = None

        # Load ligand
        if "LIGAND" in results["FILES_PATH"].keys():
            ligand_fn = results["FILES_PATH"]["LIGAND"]
            self.ligand_pdb = os.path.basename(ligand_fn)
            if os.path.split(ligand_fn)[-1] not in models_loaded:
                self.load_file(ligand_fn, self.ligand_pdb)
        else:
            self.ligand_pdb = None

        # Load cavity
        cavity_fn = results["FILES_PATH"]["OUTPUT"]
        self.cavity_pdb = os.path.basename(cavity_fn)
        if os.path.split(cavity_fn)[-1] not in models_loaded:
            self.load_file(cavity_fn, self.cavity_pdb)
        else:
            models = all_objects(self.session).models
            for model in models:
                if model.name == self.cavity_pdb:
                    model.delete()
                    self.load_file(cavity_fn, self.cavity_pdb)
        # style.style(self.session, model, atom_style='ball', dashes=5)
        # run(self.session, f"surface {model.atomspec}; transparency 50")

        return

    def characterization(self, cavities, step, atomic, vertices, probe_in, ignore_backbone):

        surface, volume, area = pyKVFinder.spatial(cavities, step=step)
        residues = pyKVFinder.constitutional(cavities, atomic, vertices, step=step, probe_in=probe_in, ignore_backbone=ignore_backbone)
        frequencies = pyKVFinder.calculate_frequencies(residues)
        scales, avg_hydropathy = pyKVFinder.hydropathy(surface, atomic, vertices, step=step, probe_in=probe_in, ignore_backbone=ignore_backbone)
        depths, max_depth, avg_depth = pyKVFinder.depth(cavities, step=step)

        return surface, volume, area, residues, scales, avg_hydropathy, depths, max_depth, avg_depth, frequencies
 
    def save_parameters(self) -> None:

        # Create base directory
        basedir = os.path.join(self.ui.output_dir_path.text(), "KV_Files")
        if not os.path.isdir(basedir):
            os.mkdir(basedir)

        # Create base_name directory
        basedir = os.path.join(basedir, self.ui.base_name.text())
        if not os.path.isdir(basedir):
            os.mkdir(basedir)    

        # Save input pdb
        if len(self.ui.input.currentData()) > 0:
            model = self._get_model(self.ui.input.currentData())
            print(f"model: {model}")
            pdb = os.path.join(
                os.path.join(basedir, f"{self.ui.base_name.text()}.KVFinder.output.pdb")
            )
            from chimerax.pdb import save_pdb
            save_pdb(self.session, pdb, models=model)
                    #cmd.save(pdb, self.input.currentText(), 0, "pdb")
        else:
            from PyQt5.QtWidgets import QMessageBox

            QMessageBox.critical(self.tool_window, "Error", "Select an input PDB!")
            return False   
        
        # Save ligand pdb
        if self.ui.ligand_adjustment.isChecked():
            if self.ui.ligand.currentText() != "":
                ligandModel = self._get_model(self.ui.ligand.currentData())
                
                if ligandModel:
                    ligand = os.path.join(
                        os.path.join(
                            basedir, f"{self.ui.base_name.text()}.KVFinder.ligand.pdb"
                        )
                    )
                    save_pdb(self.session, ligand, models=ligandModel)
            else:
                from PyQt5.QtWidgets import QMessageBox

                QMessageBox.critical(self.tool_window, "Error", "Select an ligand PDB!")
                return False
        else:
            ligand = ""

        if self.ui.box_adjustment.isChecked():
            box = self._get_model("box")
            if not box:
                from PyQt5.QtWidgets import QMessageBox

                QMessageBox.critical(self.tool_window, "Error", "Draw a box in ChimeraX!")
                return False
        pdb = pdb.replace("\\", "/")
        output = self.ui.output_dir_path.text().replace("\\", "/")
        ligand = ligand.replace("\\", "/")
        with open(os.path.join(self.ui.output_dir_path.text(), 'KV_Files', self.ui.base_name.text(), "parameters.toml"), "w") as f:
            f.write("# TOML configuration file for parKVFinder software.\n")
            f.write('\ntitle = "parKVFinder parameters file"\n')

            f.write("\n[FILES_PATH]\n")
            f.write("# The path of van der Waals radii dictionary for parKVFinder.\n")
            f.write(f'dictionary = "{self.ui.dictionary.text()}"\n')
            f.write("# The path of the input PDB file.\n")
            f.write(f'pdb = "{pdb}"\n')
            f.write("# The path of the output directory.\n")
            f.write(f"output = '{output}'\n")
            f.write("# Base name for output files.\n")
            f.write(f'base_name = "{self.ui.base_name.text()}"\n')
            f.write("# Path for the ligand's PDB file.\n")
            f.write(f'ligand = "{ligand}"\n')

            f.write("\n[SETTINGS]\n")
            f.write("# Settings for parKVFinder software.\n")
            f.write("\n[SETTINGS.modes]\n")
            f.write(
                "# Whole protein mode defines the search space as the whole protein.\n"
            )
            f.write(
                f"whole_protein_mode = {'true' if not self.ui.box_adjustment.isChecked() else 'false'}\n"
            )
            f.write(
                "# Box adjustment mode defines the search space as the box drawn in PyMOL.\n"
            )
            f.write(
                f"box_mode = {'true' if self.ui.box_adjustment.isChecked() else 'false'}\n"
            )
            f.write(
                "# Surface mode defines the type of surface representation to be applied, solvent excluded surface (true) or solvent accessible surface (false).\n"
            )
            f.write(
                f"surface_mode = {'true' if self.ui.surface.currentText() == 'Solvent Excluded Surface (SES)' else 'false'}\n"
            )
            # f.write(
            #     "# Cavity representation defines whether cavities are exported to the output PDB file as filled cavities (true) or filtered cavities (false).\n"
            # )
            # f.write(
            #     f"kvp_mode = {'true' if self.ui.cavity_representation.currentText() == 'Full' else 'false'}\n"
            # )
            f.write(
                "# Ligand adjustment mode defines the search space around the ligand.\n"
            )
            f.write(
                f"ligand_mode = {'true' if self.ui.ligand_adjustment.isChecked() else 'false'}\n"
            )

            f.write("\n[SETTINGS.step_size]\n")
            f.write(
                "# Sets the 3D grid spacing. It directly affects accuracy and runtime.\n"
            )
            step = self.ui.step_size.value()
            f.write(f"step_size = {step:.2f}\n")

            f.write("\n[SETTINGS.probes]\n")
            f.write(
                "# parKVFinder works with a dual probe system. A smaller probe, called Probe In, and a bigger one, called Probe Out, rolls around the protein.\n"
            )
            f.write(
                "# Points reached by the Probe In, but not the Probe Out are considered cavity points.\n"
            )
            f.write("# Sets the Probe In diameter. Default: 1.4 angstroms.\n")
            f.write(f"probe_in = {self.ui.probe_in.value():.2f}\n")
            f.write("# Sets the Probe Out diameter. Default: 4.0 angstroms.\n")
            f.write(f"probe_out = {self.ui.probe_out.value():.2f}\n")

            f.write("\n[SETTINGS.cutoffs]\n")
            f.write(
                "# Sets a volume cutoff for the detected cavities. Default: 5.0 angstroms.\n"
            )
            f.write(f"volume_cutoff = {self.ui.volume_cutoff.value():.2f}\n")
            f.write(
                "# Sets a distance cutoff for a search space around the ligand in ligand adjustment mode. Default: 5.0 angstroms.\n"
            )
            f.write(f"ligand_cutoff = {self.ui.ligand_cutoff.value():.2f}\n")
            f.write(
                "# Sets a removal distance for the cavity frontier, which is defined by comparing Probe In and Probe Out surfaces. Default: 2.4 angstroms.\n"
            )
            f.write(f"removal_distance = {self.ui.removal_distance.value():.2f}\n")

            f.write("\n[SETTINGS.visiblebox]\n")
            f.write(
                "# Coordinates of the vertices that define the visible 3D grid. Only four points are required to define the search space.\n\n"
            )
            box = self.create_box_parameters()
            d = {"SETTINGS": {"visiblebox": box}}
            toml.dump(o=d, f=f, encoder=toml.TomlNumpyEncoder())

            f.write("\n[SETTINGS.internalbox]\n")
            f.write("# Coordinates of the internal 3D grid. Used for calculations.\n\n")
            box = self.create_box_parameters(is_internal_box=True)
            d = {"SETTINGS": {"internalbox": box}}
            toml.dump(o=d, f=f, encoder=toml.TomlNumpyEncoder())

        return True

    def load_file(self, fname, name) -> None:

        model = self._get_model(name)

        if model:
            model.delete()

        # Load cavity filename
        if os.path.exists(fname):
            from chimerax.pdb import open_pdb

            models, status_message = open_pdb(self.session, fname)
            self.session.models.add(models)

    def clean_results(self) -> None:
        # Input File
        self.ui.input_file_entry.setText(f"")

        # Ligand File
        self.ui.ligand_file_entry.setText(f"")

        # Cavities File
        self.ui.cavities_file_entry.setText(f"")

        # Step Size
        self.ui.step_size_entry.setText(f"")

        # Volume
        self.ui.volume_list.clear()

        # Area
        self.ui.area_list.clear()

        # Depth
        self.ui.avg_depth_list.clear()
        self.ui.max_depth_list.clear()

        # Hydropathy
        self.ui.avg_hydropathy_list.clear()

        # Residues
        self.ui.residues_list.clear()

    def refresh_information(self) -> None:
        # Input File
        if "INPUT" in results["FILES_PATH"].keys():
            self.ui.input_file_entry.setText(f"{results['FILES_PATH']['INPUT']}")
        else:
            self.ui.input_file_entry.setText(f"")

        # Ligand File
        if "LIGAND" in results["FILES_PATH"].keys():
            self.ui.ligand_file_entry.setText(f"{results['FILES_PATH']['LIGAND']}")
        else:
            self.ui.ligand_file_entry.setText(f"")

        # Cavities File
        self.ui.cavities_file_entry.setText(f"{results['FILES_PATH']['OUTPUT']}")

        # Step Size
        if "PARAMETERS" in results.keys():
            if "STEP" in results["PARAMETERS"].keys():
                self.ui.step_size_entry.setText(f"{results['PARAMETERS']['STEP']:.2f}")

        return

    def box_adjustment_help(self) -> None:
        from PyQt5 import QtWidgets, QtCore

        text = QtCore.QCoreApplication.translate(
            "parKVFinder",
            '<html><head/><body><p align="justify"><span style=" font-weight:600; text-decoration: underline;">Box Adjustment mode:</span></p><p align="justify">- Create a selection (optional);</p><p align="justify">- Define a <span style=" font-weight:600;">Padding</span> (optional);</p><p align="justify">- Click on <span style=" font-weight:600;">Draw Box</span> button.</p><p align="justify"><br/><span style="text-decoration: underline;">Customize your <span style=" font-weight:600;">box</span></span>:</p><p align="justify">- Change one item at a time (e.g. <span style=" font-style:italic;">Padding</span>, <span style=" font-style:italic;">Minimum X</span>, <span style=" font-style:italic;">Maximum X</span>, ...);</p><p align="justify">- Click on <span style=" font-weight:600;">Redraw Box</span> button.<br/></p><p><span style=" font-weight:400; text-decoration: underline;">Delete </span><span style=" text-decoration: underline;">box</span><span style=" font-weight:400; text-decoration: underline;">:</span></p><p align="justify">- Click on <span style=" font-weight:600;">Delete Box</span> button.<br/></p><p align="justify"><span style="text-decoration: underline;">Colors of the <span style=" font-weight:600;">box</span> object:</span></p><p align="justify">- <span style=" font-weight:600;">Red</span> corresponds to <span style=" font-weight:600;">X</span> axis;</p><p align="justify">- <span style=" font-weight:600;">Green</span> corresponds to <span style=" font-weight:600;">Y</span> axis;</p><p align="justify">- <span style=" font-weight:600;">Blue</span> corresponds to <span style=" font-weight:600;">Z</span> axis.</p></body></html>',
            None,
        )
        help_information = QtWidgets.QMessageBox(self.tool_window)
        help_information.setText(text)
        help_information.setWindowTitle("Help")
        help_information.setStyleSheet("QLabel{min-width:500 px;}")
        help_information.exec_()

    def set_box(self) -> None:
        """
        Create box coordinates, enable 'Delete Box' and 'Redraw Box' buttons and call draw_box function.
        :param padding: box padding value.
        """

        models = all_objects(self.session).models

        for model in models:
            if model.name == "box":
                model.delete()

        sel_atoms = selected_atoms(self.session)           
        
        if len(sel_atoms) > 0 :
            minCoords = sel_atoms.coords.min(axis=0)
            maxCoords = sel_atoms.coords.max(axis=0)
            
            min_x, min_y, min_z = minCoords[0], minCoords[1], minCoords[2]
            max_x, max_y, max_z = maxCoords[0], maxCoords[1], maxCoords[2]
        else:
            model = self._get_model(self.input_pdb)

            if model:
                run(self.session, f"sel protein & {model.atomspec}")
                sel_atoms = selected_atoms(self.session)
            else:
                run(self.session, f"sel protein")
                sel_atoms = selected_atoms(self.session)

        print(f"Min coords (Box): {min_x, min_y, min_z}\nMax coords (Box): {max_x, max_y, max_z}")

        # Get center of each dimension (x, y, z)
        self.x = (min_x + max_x) / 2
        self.y = (min_y + max_y) / 2
        self.z = (min_z + max_z) / 2

        # Set Box variables in interface
        self.ui.min_x.setValue(round(self.x - (min_x - self.ui.padding.value()), 1))
        self.ui.max_x.setValue(round((max_x + self.ui.padding.value()) - self.x, 1))
        self.ui.min_y.setValue(round(self.y - (min_y - self.ui.padding.value()), 1))
        self.ui.max_y.setValue(round((max_y + self.ui.padding.value()) - self.y, 1))
        self.ui.min_z.setValue(round(self.z - (min_z - self.ui.padding.value()), 1))
        self.ui.max_z.setValue(round((max_z + self.ui.padding.value()) - self.z, 1))
        self.ui.angle1.setValue(0)
        self.ui.angle2.setValue(0)

        # Setting background box values
        self.min_x_set = self.ui.min_x.value()
        self.max_x_set = self.ui.max_x.value()
        self.min_y_set = self.ui.min_y.value()
        self.max_y_set = self.ui.max_y.value()
        self.min_z_set = self.ui.min_z.value()
        self.max_z_set = self.ui.max_z.value()

        self.angle1_set = self.ui.angle1.value()
        self.angle2_set = self.ui.angle2.value()
        
        self.padding_set = self.ui.padding.value()

        # Draw box
        self.draw_box()

        # Enable/Disable buttons
        self.ui.button_draw_box.setEnabled(False)
        self.ui.button_redraw_box.setEnabled(True)

        self.ui.min_x.setEnabled(True)
        self.ui.min_y.setEnabled(True)
        self.ui.min_z.setEnabled(True)
        self.ui.max_x.setEnabled(True)
        self.ui.max_y.setEnabled(True)
        self.ui.max_z.setEnabled(True)

        self.ui.angle1.setEnabled(True)
        self.ui.angle2.setEnabled(True)

    def delete_box(self) -> None:
        """
        Delete box object, disable 'Delete Box' and 'Redraw Box' buttons and enable 'Draw Box' button.
        """

        # Reset all box variables
        self.x = 0
        self.y = 0
        self.z = 0
        self.min_x_set = 0.0
        self.max_x_set = 0.0
        self.min_y_set = 0.0
        self.max_y_set = 0.0
        self.min_z_set = 0.0
        self.max_z_set = 0.0
        self.angle1_set = 0.0
        self.angle2_set = 0.0
        self.padding_set = 3.5

        # Delete Box and Vertices objects in PyMOL
        models = all_objects(self.session).models

        for model in models:
            if model.name == "box":
                model.delete()       

        # Set Box variables in the interface
        self.ui.min_x.setValue(self._default.min_x)
        self.ui.max_x.setValue(self._default.max_x)
        self.ui.min_y.setValue(self._default.min_y)
        self.ui.max_y.setValue(self._default.max_y)
        self.ui.min_z.setValue(self._default.min_z)
        self.ui.max_z.setValue(self._default.max_z)
        self.ui.angle1.setValue(self._default.angle1)
        self.ui.angle2.setValue(self._default.angle2)

        # Change state of buttons in the interface
        self.ui.button_draw_box.setEnabled(True)
        self.ui.button_redraw_box.setEnabled(False)
        self.ui.min_x.setEnabled(False)
        self.ui.min_y.setEnabled(False)
        self.ui.min_z.setEnabled(False)
        self.ui.max_x.setEnabled(False)
        self.ui.max_y.setEnabled(False)
        self.ui.max_z.setEnabled(False)
        self.ui.angle1.setEnabled(False)
        self.ui.angle2.setEnabled(False)

    def redraw_box(self) -> None:
        """
        Redraw box in ChimeraX interface.
        :param padding: box padding.
        :return: box object.
        """

        sel_atoms = selected_atoms(self.session)

        # Provided a selection
        if len(sel_atoms) > 0:
            # Get dimensions of selected residues
            minCoords = sel_atoms.coords.min(axis=0)
            maxCoords = sel_atoms.coords.max(axis=0)
            
            min_x, min_y, min_z = minCoords[0], minCoords[1], minCoords[2]
            max_x, max_y, max_z = maxCoords[0], maxCoords[1], maxCoords[2]

            if (
                self.ui.min_x.value() != self.min_x_set
                or self.ui.max_x.value() != self.max_x_set
                or self.ui.min_y.value() != self.min_y_set
                or self.ui.max_y.value() != self.max_y_set
                or self.ui.min_z.value() != self.min_z_set
                or self.ui.max_z.value() != self.max_z_set
                or self.ui.angle1.value() != self.angle1_set
                or self.ui.angle2.value() != self.angle2_set
            ):
                self.min_x_set = self.ui.min_x.value()
                self.max_x_set = self.ui.max_x.value()
                self.min_y_set = self.ui.min_y.value()
                self.max_y_set = self.ui.max_y.value()
                self.min_z_set = self.ui.min_z.value()
                self.max_z_set = self.ui.max_z.value()
                self.angle1_set = self.ui.angle1.value()
                self.angle2_set = self.ui.angle2.value()
            # Padding or selection altered
            else:
                # Get center of each dimension (x, y, z)
                self.x = (min_x + max_x) / 2
                self.y = (min_y + max_y) / 2
                self.z = (min_z + max_z) / 2

                # Set background box values
                self.min_x_set = (
                    round(self.x - (min_x - self.ui.padding.value()), 1)
                    + self.ui.min_x.value()
                    - self.min_x_set
                )
                self.max_x_set = (
                    round((max_x + self.ui.padding.value()) - self.x, 1)
                    + self.ui.max_x.value()
                    - self.max_x_set
                )
                self.min_y_set = (
                    round(self.y - (min_y - self.ui.padding.value()), 1)
                    + self.ui.min_y.value()
                    - self.min_y_set
                )
                self.max_y_set = (
                    round((max_y + self.ui.padding.value()) - self.y, 1)
                    + self.ui.max_y.value()
                    - self.max_y_set
                )
                self.min_z_set = (
                    round(self.z - (min_z - self.ui.padding.value()), 1)
                    + self.ui.min_z.value()
                    - self.min_z_set
                )
                self.max_z_set = (
                    round((max_z + self.ui.padding.value()) - self.z, 1)
                    + self.ui.max_z.value()
                    - self.max_z_set
                )
                self.angle1_set = 0 + self.ui.angle1.value()
                self.angle2_set = 0 + self.ui.angle2.value()
                self.padding_set = self.ui.padding.value()
        # Not provided a selection
        else:
            if (
                self.ui.min_x.value() != self.min_x_set
                or self.ui.max_x.value() != self.max_x_set
                or self.ui.min_y.value() != self.min_y_set
                or self.ui.max_y.value() != self.max_y_set
                or self.ui.min_z.value() != self.min_z_set
                or self.ui.max_z.value() != self.max_z_set
                or self.ui.angle1.value() != self.angle1_set
                or self.ui.angle2.value() != self.angle2_set
            ):
                self.min_x_set = self.ui.min_x.value()
                self.max_x_set = self.ui.max_x.value()
                self.min_y_set = self.ui.min_y.value()
                self.max_y_set = self.ui.max_y.value()
                self.min_z_set = self.ui.min_z.value()
                self.max_z_set = self.ui.max_z.value()
                self.angle1_set = self.ui.angle1.value()
                self.angle2_set = self.ui.angle2.value()

            if self.padding_set != self.ui.padding.value():
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
                self.min_x_set = round(self.x - (min_x - self.ui.padding.value()), 1)
                self.max_x_set = round((max_x + self.ui.padding.value()) - self.x, 1)
                self.min_y_set = round(self.y - (min_y - self.ui.padding.value()), 1)
                self.max_y_set = round((max_y + self.ui.padding.value()) - self.y, 1)
                self.min_z_set = round(self.z - (min_z - self.ui.padding.value()), 1)
                self.max_z_set = round((max_z + self.ui.padding.value()) - self.z, 1)
                self.angle1_set = self.ui.angle1.value()
                self.angle2_set = self.ui.angle2.value()
                self.padding_set = self.ui.padding.value()

        # Set Box variables in the interface
        self.ui.min_x.setValue(self.min_x_set)
        self.ui.max_x.setValue(self.max_x_set)
        self.ui.min_y.setValue(self.min_y_set)
        self.ui.max_y.setValue(self.max_y_set)
        self.ui.min_z.setValue(self.min_z_set)
        self.ui.max_z.setValue(self.max_z_set)
        self.ui.angle1.setValue(self.angle1_set)
        self.ui.angle2.setValue(self.angle2_set)

        models = all_objects(self.session).models

        for model in models:
            if model.name == "box":
                model.delete()    

        # Redraw box
        self.draw_box()
        
    def box_geometry(self, p1, p2, p3, p4, p5, p6, p7, p8):

        #       v2 ---- v3
        #        |\      |\
        #        | v6 ---- v7 = urf
        #        |  |    | |
        #        |  |    | |
        # llb = v0 -|---v1 |
        #         \ |     \|
        #          v4 ---- v5
        
        vertices = np.array([

            # -x, v0-v4-v2-v6
            p1, p3, p4, p7,

            # -y, v0-v1-v4-v5
            p1, p2, p3, p5,

            # -z, v1-v0-v3-v2
            p2, p1, p6, p4,

            # x, v5-v1-v7-v3
            p5, p2, p8, p6,

            # y, v3-v2-v7-v6
            p6, p4, p8, p7,

            # z, v4-v5-v6-v7
            p3, p5, p7, p8,
        ], dtype=np.float32)

        normals = np.array([
            # -x, v0-v4-v2-v6
            [-1, 0, 0],
            [-1, 0, 0],
            [-1, 0, 0],
            [-1, 0, 0],

            # -y, v0-v1-v4-v5
            [0, -1, 0],
            [0, -1, 0],
            [0, -1, 0],
            [0, -1, 0],

            # -z, v1-v0-v3-v2
            [0, 0, -1],
            [0, 0, -1],
            [0, 0, -1],
            [0, 0, -1],

            # x, v5-v1-v7-v3
            [1, 0, 0],
            [1, 0, 0],
            [1, 0, 0],
            [1, 0, 0],

            # y, v3-v2-v7-v6
            [0, 1, 0],
            [0, 1, 0],
            [0, 1, 0],
            [0, 1, 0],

            # z, v4-v5-v6-v7
            [0, 0, 1],
            [0, 0, 1],
            [0, 0, 1],
            [0, 0, 1],
        ], dtype=np.float32)
        
        triangles = np.array([
            [0, 1, 2], [2, 1, 3],           # -x
            [4, 5, 6], [6, 5, 7],           # -y
            [8, 9, 10], [10, 9, 11],        # -z
            [12, 13, 14], [14, 13, 15],     # x
            [16, 17, 18], [18, 17, 19],     # y
            [20, 21, 22], [22, 21, 23],     # z
        ], dtype=np.int32)
        
        return vertices, normals, triangles
    
    def show_grid(self) -> None:
        """
        Callback for the "Show Grid" button
        - Get minimum and maximum coordinates of the KVFinder-web 3D-grid, dependent on selected parameters.
        :return: Call draw_grid function with minimum and maximum coordinates or return Error.
        """

        global xg, yg, zg

        if self.ui.input.count() > 0:
            # Get minimum and maximum dimensions of target PDB
            pdb = self.ui.input.currentText()
            
            model = self._get_model(pdb)
            
            run(self.session, f"sel {model.atomspec}")
            sel_atoms = selected_atoms(self.session) 
            run(self.session, "sel clear")          
            
            if len(sel_atoms) > 0 :
                minCoords = sel_atoms.coords.min(axis=0)
                maxCoords = sel_atoms.coords.max(axis=0)
                
                min_x, min_y, min_z = minCoords[0], minCoords[1], minCoords[2]
                max_x, max_y, max_z = maxCoords[0], maxCoords[1], maxCoords[2]
            else:
                print(f"I can't find the model {self.input}. Check if the model is open")
                

            print(f"Min coords (Grid): {min_x, min_y, min_z}\nMax coords (Grid): {max_x, max_y, max_z}")

            # Get Probe Out value
            probe_out = self.ui.probe_out.value()
            probe_out = round(probe_out - round(probe_out, 4) % round(0.6, 4), 1)

            # Prepare dimensions
            min_x = round(min_x - (min_x % 0.6), 1) - probe_out
            min_y = round(min_y - (min_y % 0.6), 1) - probe_out
            min_z = round(min_z - (min_z % 0.6), 1) - probe_out
            max_x = round(max_x - (max_x % 0.6) + 0.6, 1) + probe_out
            max_y = round(max_y - (max_y % 0.6) + 0.6, 1) + probe_out
            max_z = round(max_z - (max_z % 0.6) + 0.6, 1) + probe_out

            # Get center of each dimension (x, y, z)
            xg = (min_x + max_x) / 2
            yg = (min_y + max_y) / 2
            zg = (min_z + max_z) / 2

            # Draw Grid
            self.draw_grid(min_x, max_x, min_y, max_y, min_z, max_z)
        else:
            from PyQt5.QtWidgets import QMessageBox

            QMessageBox.critical(self.tool_window, "Error", "Select an input PDB!")
            return

    def draw_grid(self, min_x, max_x, min_y, max_y, min_z, max_z) -> None:
        """
        Draw Grid in ChimeraX.
        :param min_x: minimum X coordinate.
        :param max_x: maximum X coordinate.
        :param min_y: minimum Y coordinate.
        :param max_y: maximum Y coordinate.
        :param min_z: minimum Z coordinate.
        :param max_z: maximum Z coordinate.
        :return: grid object in PyMOL.
        """

        from math import sin, cos
        
        models = all_objects(self.session).models

        for model in models:
            if model.name == "grid":
                model.delete()   

        # Prepare dimensions
        angle1 = 0.0
        angle2 = 0.0
        min_x = xg - min_x
        max_x = max_x - xg
        min_y = yg - min_y
        max_y = max_y - yg
        min_z = zg - min_z
        max_z = max_z - zg

        # Get positions of grid vertices
        # P1
        x1 = (
            -min_x * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + xg
        )

        y1 = -min_y * cos(angle1) + (-min_z) * sin(angle1) + yg

        z1 = (
            min_x * sin(angle2)
            + min_y * sin(angle1) * cos(angle2)
            - min_z * cos(angle1) * cos(angle2)
            + zg
        )

        # P2
        x2 = (
            max_x * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + xg
        )

        y2 = (-min_y) * cos(angle1) + (-min_z) * sin(angle1) + yg

        z2 = (
            (-max_x) * sin(angle2)
            - (-min_y) * sin(angle1) * cos(angle2)
            + (-min_z) * cos(angle1) * cos(angle2)
            + zg
        )

        # P3
        x3 = (
            (-min_x) * cos(angle2)
            - max_y * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + xg
        )

        y3 = max_y * cos(angle1) + (-min_z) * sin(angle1) + yg

        z3 = (
            -(-min_x) * sin(angle2)
            - max_y * sin(angle1) * cos(angle2)
            + (-min_z) * cos(angle1) * cos(angle2)
            + zg
        )

        # P4
        x4 = (
            (-min_x) * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + max_z * cos(angle1) * sin(angle2)
            + xg
        )

        y4 = (-min_y) * cos(angle1) + max_z * sin(angle1) + yg

        z4 = (
            -(-min_x) * sin(angle2)
            - (-min_y) * sin(angle1) * cos(angle2)
            + max_z * cos(angle1) * cos(angle2)
            + zg
        )

        # P5
        x5 = (
            max_x * cos(angle2)
            - max_y * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + xg
        )

        y5 = max_y * cos(angle1) + (-min_z) * sin(angle1) + yg

        z5 = (
            (-max_x) * sin(angle2)
            - max_y * sin(angle1) * cos(angle2)
            + (-min_z) * cos(angle1) * cos(angle2)
            + zg
        )

        # P6
        x6 = (
            max_x * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + max_z * cos(angle1) * sin(angle2)
            + xg
        )

        y6 = (-min_y) * cos(angle1) + max_z * sin(angle1) + yg

        z6 = (
            (-max_x) * sin(angle2)
            - (-min_y) * sin(angle1) * cos(angle2)
            + max_z * cos(angle1) * cos(angle2)
            + zg
        )

        # P7
        x7 = (
            (-min_x) * cos(angle2)
            - max_y * sin(angle1) * sin(angle2)
            + max_z * cos(angle1) * sin(angle2)
            + xg
        )

        y7 = max_y * cos(angle1) + max_z * sin(angle1) + yg

        z7 = (
            -(-min_x) * sin(angle2)
            - max_y * sin(angle1) * cos(angle2)
            + max_z * cos(angle1) * cos(angle2)
            + zg
        )

        # P8
        x8 = (
            max_x * cos(angle2)
            - max_y * sin(angle1) * sin(angle2)
            + max_z * cos(angle1) * sin(angle2)
            + xg
        )

        y8 = max_y * cos(angle1) + max_z * sin(angle1) + yg

        z8 = (
            (-max_x) * sin(angle2)
            - max_y * sin(angle1) * cos(angle2)
            + max_z * cos(angle1) * cos(angle2)
            + zg
        )

        P1 = (x1, y1, z1)
        P2 = (x2, y2, z2)
        P3 = (x3, y3, z3)
        P4 = (x4, y4, z4)
        P5 = (x5, y5, z5)
        P6 = (x6, y6, z6)
        P7 = (x7, y7, z7)
        P8 = (x8, y8, z8)
        
        from chimerax.shape.shape import _show_surface
        
        points = [P1, P2, P3, P4, P5, P6, P7, P8]
        
        from chimerax.markers.markers import MarkerSet, create_link
        
        markers = MarkerSet(self.session, name="grid")
        

        links = [(1, 2), (1, 3), (1, 4), (3, 5), (3, 7), (5, 2), (5, 8),
 (2, 6), (4, 7), (7, 8), (8, 6), (6, 4)]
        
        #       v2 ---- v3
        #        |\      |\
        #        | v6 ---- v7 = urf
        #        |  |    | |
        #        |  |    | |
        # llb = v0 -|---v1 |
        #         \ |     \|
        #          v4 ---- v5
        
        # p1: v0
        # p3: v4
        # p4: v2
        # p7: v6
        # p2: v1
        # p6: v3
        # p7: v6
        # p8: v7
        # p5: v5
        
        """
        links = {
            1:  [2, 3, 4],
            8: [5, 6, 7],
            4: [6, 7],
            5: [2, 3],
            2: [6],
            3: [7]
        }
        """
        
        i = 1
        for point in points:
            markers.create_marker(xyz= point, rgba=(255, 200, 0, 40), radius=1.5, id= i )
            i += 1
            
        self.session.models.add([markers])
        
        model = self._get_model("grid")
        
        if model:
            atoms = model.atoms
            
            for link in links:
                a, b = link[0], link[1]
                create_link(atoms[a-1], atoms[b-1])
              
        
        
        # varray, normals, tarray = self.box_geometry(P1, P2, P3, P4, P5, P6, P7, P8)

        # # Create box object
        # self.s = _show_surface(self.session, varray=varray, tarray=tarray, color = (255, 255, 220, 30), mesh=False,
        #             center=None, rotation=None, qrotation=None, coordinate_system=None, 
        #             slab=None, model_id= None, shape_name="grid")
        
    def draw_box(self) -> None:
        """
        Callback for the "Draw box" button.

        This method calculates each vertice of the custom box. Then, it draws and connects them on the PyMOL viewer as a object named 'box'.
        """

        from math import cos, pi, sin

        # Convert angle
        angle1 = (self.ui.angle1.value() / 180.0) * pi
        angle2 = (self.ui.angle2.value() / 180.0) * pi

        # Get positions of box vertices
        # P1
        x1 = (
            -self.ui.min_x.value() * cos(angle2)
            - (-self.ui.min_y.value()) * sin(angle1) * sin(angle2)
            + (-self.ui.min_z.value()) * cos(angle1) * sin(angle2)
            + self.x
        )

        y1 = (
            -self.ui.min_y.value() * cos(angle1)
            + (-self.ui.min_z.value()) * sin(angle1)
            + self.y
        )

        z1 = (
            self.ui.min_x.value() * sin(angle2)
            + self.ui.min_y.value() * sin(angle1) * cos(angle2)
            - self.ui.min_z.value() * cos(angle1) * cos(angle2)
            + self.z
        )

        # P2
        x2 = (
            self.ui.max_x.value() * cos(angle2)
            - (-self.ui.min_y.value()) * sin(angle1) * sin(angle2)
            + (-self.ui.min_z.value()) * cos(angle1) * sin(angle2)
            + self.x
        )

        y2 = (
            (-self.ui.min_y.value()) * cos(angle1)
            + (-self.ui.min_z.value()) * sin(angle1)
            + self.y
        )

        z2 = (
            (-self.ui.max_x.value()) * sin(angle2)
            - (-self.ui.min_y.value()) * sin(angle1) * cos(angle2)
            + (-self.ui.min_z.value()) * cos(angle1) * cos(angle2)
            + self.z
        )

        # P3
        x3 = (
            (-self.ui.min_x.value()) * cos(angle2)
            - self.ui.max_y.value() * sin(angle1) * sin(angle2)
            + (-self.ui.min_z.value()) * cos(angle1) * sin(angle2)
            + self.x
        )

        y3 = (
            self.ui.max_y.value() * cos(angle1)
            + (-self.ui.min_z.value()) * sin(angle1)
            + self.y
        )

        z3 = (
            -(-self.ui.min_x.value()) * sin(angle2)
            - self.ui.max_y.value() * sin(angle1) * cos(angle2)
            + (-self.ui.min_z.value()) * cos(angle1) * cos(angle2)
            + self.z
        )

        # P4
        x4 = (
            (-self.ui.min_x.value()) * cos(angle2)
            - (-self.ui.min_y.value()) * sin(angle1) * sin(angle2)
            + self.ui.max_z.value() * cos(angle1) * sin(angle2)
            + self.x
        )

        y4 = (
            (-self.ui.min_y.value()) * cos(angle1)
            + self.ui.max_z.value() * sin(angle1)
            + self.y
        )

        z4 = (
            -(-self.ui.min_x.value()) * sin(angle2)
            - (-self.ui.min_y.value()) * sin(angle1) * cos(angle2)
            + self.ui.max_z.value() * cos(angle1) * cos(angle2)
            + self.z
        )

        # P5
        x5 = (
            self.ui.max_x.value() * cos(angle2)
            - self.ui.max_y.value() * sin(angle1) * sin(angle2)
            + (-self.ui.min_z.value()) * cos(angle1) * sin(angle2)
            + self.x
        )

        y5 = (
            self.ui.max_y.value() * cos(angle1)
            + (-self.ui.min_z.value()) * sin(angle1)
            + self.y
        )

        z5 = (
            (-self.ui.max_x.value()) * sin(angle2)
            - self.ui.max_y.value() * sin(angle1) * cos(angle2)
            + (-self.ui.min_z.value()) * cos(angle1) * cos(angle2)
            + self.z
        )

        # P6
        x6 = (
            self.ui.max_x.value() * cos(angle2)
            - (-self.ui.min_y.value()) * sin(angle1) * sin(angle2)
            + self.ui.max_z.value() * cos(angle1) * sin(angle2)
            + self.x
        )

        y6 = (
            (-self.ui.min_y.value()) * cos(angle1)
            + self.ui.max_z.value() * sin(angle1)
            + self.y
        )

        z6 = (
            (-self.ui.max_x.value()) * sin(angle2)
            - (-self.ui.min_y.value()) * sin(angle1) * cos(angle2)
            + self.ui.max_z.value() * cos(angle1) * cos(angle2)
            + self.z
        )

        # P7
        x7 = (
            (-self.ui.min_x.value()) * cos(angle2)
            - self.ui.max_y.value() * sin(angle1) * sin(angle2)
            + self.ui.max_z.value() * cos(angle1) * sin(angle2)
            + self.x
        )

        y7 = (
            self.ui.max_y.value() * cos(angle1) + self.ui.max_z.value() * sin(angle1) + self.y
        )

        z7 = (
            -(-self.ui.min_x.value()) * sin(angle2)
            - self.ui.max_y.value() * sin(angle1) * cos(angle2)
            + self.ui.max_z.value() * cos(angle1) * cos(angle2)
            + self.z
        )

        # P8
        x8 = (
            self.ui.max_x.value() * cos(angle2)
            - self.ui.max_y.value() * sin(angle1) * sin(angle2)
            + self.ui.max_z.value() * cos(angle1) * sin(angle2)
            + self.x
        )

        y8 = (
            self.ui.max_y.value() * cos(angle1) + self.ui.max_z.value() * sin(angle1) + self.y
        )

        z8 = (
            (-self.ui.max_x.value()) * sin(angle2)
            - self.ui.max_y.value() * sin(angle1) * cos(angle2)
            + self.ui.max_z.value() * cos(angle1) * cos(angle2)
            + self.z
        )
        
        P1 = [x1, y1, z1]
        P2 = [x2, y2, z2]
        P3 = [x3, y3, z3]
        P4 = [x4, y4, z4]
        P5 = [x5, y5, z5]
        P6 = [x6, y6, z6]
        P7 = [x7, y7, z7]
        P8 = [x8, y8, z8]
        
        from chimerax.shape.shape import _show_surface
        
        varray, normals, tarray = self.box_geometry(P1, P2, P3, P4, P5, P6, P7, P8)

        self.s = _show_surface(self.session, varray=varray, tarray=tarray, color = (255, 0, 255, 100), mesh=False,
                    center=None, rotation=None, qrotation=None, coordinate_system=None, 
                    slab=None, model_id= None, shape_name="box")
        
    def create_box_parameters(
        self, is_internal_box=False
    ):
        from math import pi, cos, sin

        # Get box parameters
        if self.ui.box_adjustment.isChecked():
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
            min_x += self.ui.probe_out.value()
            max_x += self.ui.probe_out.value()
            min_y += self.ui.probe_out.value()
            max_y += self.ui.probe_out.value()
            min_z += self.ui.probe_out.value()
            max_z += self.ui.probe_out.value()

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
                self.ui.output_dir_path.setText(fname)

        return

    def cprint(self, text):
        return self.session.logger.info(text)

    def get_number_of_cavities(self):
        # Read results (Ubuntu/macOS)
        results = toml.load(
            f"{os.path.join(self.ui.output_dir_path.text(), 'KV_Files', self.ui.base_name.text(), f'{self.ui.base_name.text()}.KVFinder.results.toml')}"
        )

        return len(results["RESULTS"]["VOLUME"].keys())

    def extract_pdb_session(self, list_models, selected=True):
        """Extract the PDB of a specific model.
        By default, the `selected` option is True, so it will extract the PDB
        of selected atoms in the model {name}. When `selected` is False, this function extracts
        the PDB of all atoms inside the model.

        Parameters
        ----------
        name : str
            Name of the model.
        selected : bool, optional
            Controls whether the function will return the entire model or only selected atoms. 
            Defaults to True.

        Raises
        ------
        AssertionError
            If the function cannot find the model {name}, it raises an exception error.

        Returns
        -------
        numpy.ndarray
            An array containing the atomic information.
        """

        if selected:
            sel_atoms = selected_atoms(self.session)
        else:
            sel_atoms = list_models[0].atoms
            for model in list_models[1:]:
                sel_atoms = sel_atoms.merge(model.atoms)

            
            if len(sel_atoms) < 1:
                raise AssertionError(f"WARNING: I didn't find any structure with the name {name}")


        atomNP = np.zeros(shape=(len(sel_atoms), 8), dtype='<U32')
        if self.ui.dictionary.text() != "":            
            path = self.ui.dictionary.text()
            vdw = pyKVFinder.read_vdw(path)
        else:  
            vdw = pyKVFinder.read_vdw()
        for i in range(0, len(sel_atoms)):
            atom = sel_atoms[i]
            residue_name, atom_name, atom_element = str(atom.residue.name).upper(), str(atom.name).upper(), str(atom.element).upper()
            if residue_name in vdw.keys() and atom_name in vdw[residue_name].keys():
                radius = vdw[residue_name][atom_name]
            else:
                radius = vdw["GEN"][atom_element]
                self.cprint(f"Warning: Atom {atom_name} of residue {residue_name} \  not found in dictionary.")
                self.cprint(f"Warning: Using generic atom {atom_element} \radius: {radius} \u00c5.")

            try:
                atomNP[i] = [str(atom.residue.number), str(atom.residue.chain)[1:], residue_name, atom_name, atom.coord[0], atom.coord[1], atom.coord[2], radius ]
            except:
                self.cprint(f"Problem to modify line {str(i)}: {str(atom)}")

        # print(f"atomNP Final: {atomNP}")
        
        return atomNP       
   
    def extract_pdb_session_old(self, name, selected=True):
        """Extract the PDB of a specific model.
        By default, the `selected` option is True, so it will extract the PDB
        of selected atoms in the model {name}. When `selected` is False, this function extracts
        the PDB of all atoms inside the model.

        Parameters
        ----------
        name : str
            Name of the model.
        selected : bool, optional
            Controls whether the function will return the entire model or only selected atoms. 
            Defaults to True.

        Raises
        ------
        AssertionError
            If the function cannot find the model {name}, it raises an exception error.

        Returns
        -------
        numpy.ndarray
            An array containing the atomic information.
        """

        if selected:
            sel_atoms = selected_atoms(self.session)
        else:
            structures = all_atomic_structures(self.session)
            for structure in structures:
                if structure.name == name:
                    sel_atoms = structure.atoms
                    break
            else:
                raise AssertionError(f"WARNING: I didn't find any structure with the name {name}")


        atomNP = np.zeros(shape=(len(sel_atoms), 8), dtype='<U32')
        if self.ui.dictionary.text() != "":            
            path = self.ui.dictionary.text()
            vdw = pyKVFinder.read_vdw(path)
        else:  
            vdw = pyKVFinder.read_vdw()
        for i in range(0, len(sel_atoms)):
            atom = sel_atoms[i]
            residue_name, atom_name, atom_element = str(atom.residue.name).upper(), str(atom.name).upper(), str(atom.element).upper()
            if residue_name in vdw.keys() and atom_name in vdw[residue_name].keys():
                radius = vdw[residue_name][atom_name]
            else:
                radius = vdw["GEN"][atom_element]
                self.cprint(f"Warning: Atom {atom_name} of residue {residue_name} \  not found in dictionary.")
                self.cprint(f"Warning: Using generic atom {atom_element} \radius: {radius} \u00c5.")

            try:
                atomNP[i] = [str(atom.residue.number), str(atom.residue.chain)[1:], residue_name, atom_name, atom.coord[0], atom.coord[1], atom.coord[2], radius ]
            except:
                self.cprint(f"Problem to modify line {str(i)}: {str(atom)}")

        return atomNP       

    def check_resolution(self):
        if self.ui.resolution_label.isChecked():
            self.ui.resolution.setEnabled(True)
            self.ui.resolution.setCurrentText(self._default.resolution)
            self.ui.step_size_label.setChecked(False)
            self.ui.step_size.setEnabled(False)
            self.ui.step_size.setValue(self._default.step)
        else:
            self.ui.resolution.setEnabled(False)
            self.ui.resolution.setCurrentText("Off")
            self.ui.step_size_label.setChecked(True)
            self.ui.step_size.setEnabled(True)
            self.ui.step_size.setValue(0.6)

    def check_step_size(self):
        if self.ui.step_size_label.isChecked():
            self.ui.resolution_label.setChecked(False)
            self.ui.resolution.setEnabled(False)
            self.ui.resolution.setCurrentText("Off")
            self.ui.step_size.setEnabled(True)
            self.ui.step_size.setValue(0.6)
        else:
            self.ui.resolution_label.setChecked(True)
            self.ui.resolution.setEnabled(True)
            self.ui.resolution.setCurrentText(self._default.resolution)
            self.ui.step_size.setEnabled(False)
            self.ui.step_size.setValue(self._default.step)    

    def fill_context_menu(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically filled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.) 
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        from Qt.QtGui import QAction
        clear_action = QAction("Clear", menu)
        clear_action.triggered.connect(lambda *args: self.line_edit.clear())
        menu.addAction(clear_action)

        # We will be adding an item to the tool's context menu, so override
        # the default MainToolWindow fill_context_menu method
        self.tool_window.fill_context_menu = self.fill_context_menu

    def select_file(self, caption, entry, filters) -> None:
        """
        Callback for the "Browse ..." button
        Open a QFileDialog to select a file.
        """

        from PyQt5.QtWidgets import QFileDialog
        from PyQt5.QtCore import QDir

        # Get results file
        fname, _ = QFileDialog.getOpenFileName(
            self.tool_window, caption=caption, directory=os.getcwd(), filter=filters
        )

        if fname:
            fname = QDir.toNativeSeparators(fname)
            if os.path.exists(fname):
                entry.setText(fname)

        return

    def refresh_area(self) -> None:
        # Get cavity indexes
        indexes = sorted(results["RESULTS"]["AREA"].keys())
        # Include Area
        for index in indexes:
            item = f"{index}: {results['RESULTS']['AREA'][index]}"
            self.ui.area_list.addItem(item)
        return

    def refresh_volume(self) -> None:
        # Get cavity indexes
        indexes = sorted(results["RESULTS"]["VOLUME"].keys())
        # Include Volume
        for index in indexes:
            item = f"{index}: {results['RESULTS']['VOLUME'][index]}"
            self.ui.volume_list.addItem(item)
        return
    
    def refresh_avg_depth(self) -> None:
        # Get cavity indexes
        indexes = sorted(results["RESULTS"]["AVG_DEPTH"].keys())
        # Include Average Depth
        for index in indexes:
            item = f"{index}: {results['RESULTS']['AVG_DEPTH'][index]}"
            self.ui.avg_depth_list.addItem(item)
        return

    def refresh_max_depth(self) -> None:
        # Get cavity indexes
        indexes = sorted(results["RESULTS"]["MAX_DEPTH"].keys())
        # Include Maximum Depth
        for index in indexes:
            item = f"{index}: {results['RESULTS']['MAX_DEPTH'][index]}"
            self.ui.max_depth_list.addItem(item)
        return

    def refresh_avg_hydropathy(self) -> None:
        # Get cavity indexes
        indexes = sorted(results["RESULTS"]["AVG_HYDROPATHY"].keys())
        # Include Average Hydropathy
        for index in indexes:
            if index != "EisenbergWeiss":
                item = f"{index}: {results['RESULTS']['AVG_HYDROPATHY'][index]}"
                self.ui.avg_hydropathy_list.addItem(item)
        return

    def refresh_residues(self) -> None:
        # Get cavity indexes
        indexes = sorted(results["RESULTS"]["RESIDUES"].keys())
        # Include Interface Residues
        for index in indexes:
            self.ui.residues_list.addItem(index)
        return
    
    def _get_model(self, name) -> None | AtomicStructure | list:
        
        if isinstance(name, str):
            objects = all_objects(self.session)
            models = objects.models
            
            for model in models:
                if model.name == name:
                    return model
            else:       
                return None
        elif isinstance(name, list):
            objects = all_objects(self.session)
            models_obj = objects.models
            
            models = []
            print(f"name: {name}")
            print(f"models_obj: {models_obj}")
        
            
            for model in models_obj:
                model_name = model.name + " " + model.id_string
                print(f"model_name: {model_name}")
                
                if model_name in name:
                    models.append(model)
            
            return models if models else None
                
                
    def _deselect_all_items(self, list_widget):
        for index in range(list_widget.count()):
            item = list_widget.item(index)
            if item.isSelected():
                item.setSelected(False)

    def _reset_areas(self):
        listWidgets = [self.ui.volume_list, self.ui.area_list, self.ui.avg_depth_list, self.ui.max_depth_list, self.ui.avg_hydropathy_list, self.ui.residues_list]
        for widget in listWidgets:
            self._deselect_all_items(widget)

    def show_residues(self) -> None:
        """Show the residues in structure
        """

        # Select items of list2
        list1 = self.ui.residues_list
        cavs = []
        residues = []
        deselect = []
        deselects_res = []
        number_of_items = list1.count()
        for index in range(number_of_items):
            item1 = list1.item(index)

            if item1.text()[0:3] not in self.res_selected and item1.isSelected():
                cavs.append(item1.text()[0:3])
                self.res_selected.append(item1.text()[0:3])


            elif item1.text()[0:3] in self.res_selected and item1.isSelected() == False:
                self.res_selected.remove(item1.text()[0:3])
                deselect.append(item1.text()[0:3])
            
            elif item1.text()[0:3] in self.res_selected and item1.isSelected():
                cavs.append(item1.text()[0:3])

        for cav in cavs:
            for residue in results["RESULTS"]["RESIDUES"][cav]:
                if residue not in residues:
                    residues.append(residue)

        model = self._get_model(self.input_pdb)

        if model:
            spec = model.atomspec
            if deselect:

                for residue in results["RESULTS"]["RESIDUES"][deselect[0]]:
                    if residue not in deselects_res:
                        deselects_res.append(residue)
                command = spec

                while len(deselects_res) > 0:
                    res, chain, _ = deselects_res.pop(0)
                    command = f"{command}/{chain}:{res}"
                command_hide = f"hide {command}"
                run(self.session, command_hide)
                return 
            # Select residues

            command = spec
            while len(residues) > 0:
                res, chain, _ = residues.pop(0)
                command = f"{command}/{chain}:{res}"
            command_show = f"show {command}"

            run(self.session, command_show)

            command_style = f"style {command} stick"
            run(self.session, command_style)

        else:
            print(f"Didn't find the model {self.input_pdb}")

    def show_default_view(self) -> None:
        """    Color all cavities with white.

        This method is responsible for coloring all cavities within the structure with white.

        Returns
        -------
        None
            This method does not return anything; it only modifies the color of the cavities.
        """

        model = self._get_model(self.cavity_pdb)


        if model:
            spec = model.atomspec
            command = f"color {spec} white"
            run(self.session, command)
            self._reset_areas()
        else:
            print(f"WARNING: Didn't find the model {self.cavity_pdb}")

    def show_depth_view(self) -> None:
        """
        Color all cavities using the depth attribute.

        This method is responsible for coloring all cavities within the structure based on their depth attribute.

        Returns
        -------
        None
            This method does not return anything; it only modifies the color of the cavities.
        """

        model = self._get_model(self.cavity_pdb)

        if model:
            spec = model.atomspec
            command = f"color byattribute bfactor {spec} palette rainbow"
            run(self.session, command)
            self._reset_areas()
        else:
            print(f"WARNING: I Didn't find the model {self.cavity_pdb}")

    def show_hydropathy_view(self) -> None:
        """
        Color all cavities using the hydropath attribute.

        This method is responsible for coloring all cavities within the structure based on their hydropath attribute.

        Returns
        -------
        None
            This method does not return anything; it only modifies the color of the cavities.
        """

        model = self._get_model(self.cavity_pdb)
        inpModel = self._get_model(self.input_pdb)

        if model:
            spec = model.atomspec + "@HA"

            command = f"color byattribute occupancy {spec} palette yellow:white:blue"
            run(self.session, command)
            self._reset_areas()
        else:
            print(f"Didn't find the model {self.cavity_pdb}")

    def show_cavities(self, list1, list2) -> None:
        """   
        Color the atoms of selected cavities.

        This method colors the atoms of selected cavities. It uses blue for pseudoatoms with HA bonds
        and red for pseudoatoms that do not form bonds.

        Parameters
        ----------
        list1 : list
            A list of QtListItem representing pseudoatoms.
        list2 : list
            A list of QtListItem representing pseudoatoms.

        Returns
        -------
        None
        """

        # Get items from list1
        cavs = []
        deselect = []

        # Select items of list2
        number_of_items = list1.count()
        for index in range(number_of_items):
            item1 = list1.item(index)
            item2 = list2.item(index)

            if item1.text()[0:3] not in self.vs_selected and item1.isSelected():
                cavs.append(item1.text()[0:3])
                self.vs_selected.append(item1.text()[0:3])
                item2.setSelected(True)

            elif item1.text()[0:3] in self.vs_selected and item1.isSelected() == False:
                self.vs_selected.remove(item1.text()[0:3])
                deselect.append(item1.text()[0:3])
                item2.setSelected(False)

        model = self._get_model(self.cavity_pdb)

        if model:
            spec = model.atomspec
            # Return if no cavity is selected
            if len(cavs) > 0:
                command = f"color {spec}:{cavs[0]} blue; color {spec}:{cavs[0]}@HA red"
            elif len(deselect) > 0:
                command = f"color {spec}:{deselect[0]} white"
            else:
                return
            
            run(self.session, command)
            # style(self.session, spec, atom_style='ball', dashes=5)
        else:
            print(f"Didn't find the model {self.cavity_pdb}")

    def show_depth(self, list1, list2) -> None:
        """
        Color atoms by depth attribute using a rainbow palette.

        This method colors atoms based on their depth attribute using a rainbow color palette.

        Parameters
        ----------
        list1 : list
            A list of QtListItem representing average depth.
        list2 : list
            A list of QtListItem representing maximum depth.

        Returns
        -------
        None
        """

        # Get items from list1
        cavs = []
        deselect = []

        # Select items of list2
        number_of_items = list1.count()
        for index in range(number_of_items):
            item1 = list1.item(index)
            item2 = list2.item(index)

            if item1.text()[0:3] not in self.am_selected and item1.isSelected():
                cavs.append(item1.text()[0:3])
                self.am_selected.append(item1.text()[0:3])
                item2.setSelected(True)

            elif item1.text()[0:3] in self.am_selected and item1.isSelected() == False:
                self.am_selected.remove(item1.text()[0:3])
                deselect.append(item1.text()[0:3])
                item2.setSelected(False)
            elif item1.text()[0:3] in self.am_selected and item1.isSelected():
                cavs.append(item1.text()[0:3])

        model = self._get_model(self.cavity_pdb)

        if model:
            spec = model.atomspec
        
            # Return if no cavity is selected
            if len(cavs) > 0:
                cavs = ":".join(cavs)
                command = f"color byattribute bfactor {spec}:{cavs} palette rainbow"
                cavs = []
                # while len(cavs) > 0:
                #     command = f"{command}{cavs.pop(0)}:"
                # command = command[:-1] + " palette rainbow"
                run(self.session, command)

            if len(deselect) > 0:
                command = f"color {spec}:{deselect[0]} white"
                run(self.session, command)
            else:
                return
        else:
            print(f"Didn't find the model {self.cavity_pdb}")

    def show_hydropathy(self, list1) -> None:
        """
        Color atoms in selected cavities by hydropathy factor using a yellow-white-blue color palette.

        This method colors atoms within selected cavities based on their hydropathy factor using a yellow-white-blue color palette.

        Parameters
        ----------
        list1 : list
            A list of QtListItem representing atoms within selected cavities.

        Returns
        -------
        None
        """

        # Get items from list1
        cavs = []
        deselect = []

        # Select items of list2
        number_of_items = list1.count()
        for index in range(number_of_items):
            item1 = list1.item(index)

            if item1.text()[0:3] not in self.hyd_selected and item1.isSelected():
                cavs.append(item1.text()[0:3])
                self.hyd_selected.append(item1.text()[0:3])

            elif item1.text()[0:3] in self.hyd_selected and item1.isSelected() == False:
                self.hyd_selected.remove(item1.text()[0:3])
                deselect.append(item1.text()[0:3])
            elif item1.text()[0:3] in self.hyd_selected and item1.isSelected():
                cavs.append(item1.text()[0:3])

        model = self._get_model(self.cavity_pdb)
        inpModel = self._get_model(self.input_pdb)

        if model:
            spec = model.atomspec
        
            # Return if no cavity is selected
            if len(cavs) > 0:
                cavs = "@HA:".join(cavs)
                
                if inpModel:
                    command = f"color byattribute occupancy {inpModel.atomspec}:{cavs}@HA palette yellow:white:blue"
                else:
                    command = f"color byattribute occupancy {spec}:{cavs}@HA palette yellow:white:blue"

                run(self.session, command)

            if len(deselect) > 0:
                command = f"color {spec}:{deselect[0]} white"
                run(self.session, command)
            else:
                return
        else:
            print(f"Didn't find the model {self.cavity_pdb}")

class Ui_pyKVFinder(object):
    def setupUi(self, pyKVFinder):
        pyKVFinder.setObjectName("pyKVFinder")
        self.gui = QtWidgets.QWidget()
        self.gui.setObjectName("gui")
        self.gridLayout = QtWidgets.QGridLayout(self.gui)
        self.gridLayout.setContentsMargins(10, 10, 10, 10)
        self.gridLayout.setVerticalSpacing(10)
        self.gridLayout.setObjectName("gridLayout")
        self.dialog_separator = QtWidgets.QFrame(self.gui)
        self.dialog_separator.setFrameShape(QtWidgets.QFrame.HLine)
        self.dialog_separator.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.dialog_separator.setObjectName("dialog_separator")
        self.gridLayout.addWidget(self.dialog_separator, 2, 0, 1, 1)
        self.dialog_buttons = QtWidgets.QHBoxLayout()
        self.dialog_buttons.setContentsMargins(20, -1, 20, -1)
        self.dialog_buttons.setSpacing(6)
        self.dialog_buttons.setObjectName("dialog_buttons")
        self.button_run = QtWidgets.QPushButton(self.gui)

        #sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        #sizePolicy.setHorizontalStretch(0)
        #sizePolicy.setVerticalStretch(0)
        #sizePolicy.setHeightForWidth(self.button_run.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.button_run)

        self.button_run.setSizePolicy(sizePolicy)
        self.button_run.setText("Run pyKVFinder")
        self.button_run.setObjectName("button_run")
        self.dialog_buttons.addWidget(self.button_run)
        self.button_grid = QtWidgets.QPushButton(self.gui)

        #sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        #sizePolicy.setHorizontalStretch(0)
        #sizePolicy.setVerticalStretch(0)
        #sizePolicy.setHeightForWidth(self.button_grid.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.button_grid)

        self.button_grid.setSizePolicy(sizePolicy)
        self.button_grid.setObjectName("button_grid")
        self.dialog_buttons.addWidget(self.button_grid)
        self.button_save_parameters = QtWidgets.QPushButton(self.gui)

        #sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        #sizePolicy.setHorizontalStretch(0)
        #sizePolicy.setVerticalStretch(0)
        #sizePolicy.setHeightForWidth(self.button_save_parameters.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.button_save_parameters)

        self.button_save_parameters.setSizePolicy(sizePolicy)
        self.button_save_parameters.setText("Save Parameters")
        self.button_save_parameters.setObjectName("button_save_parameters")
        self.dialog_buttons.addWidget(self.button_save_parameters)
        self.button_restore = QtWidgets.QPushButton(self.gui)

        #sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        #sizePolicy.setHorizontalStretch(0)
        #sizePolicy.setVerticalStretch(0)
        #sizePolicy.setHeightForWidth(self.button_restore.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.button_restore)

        self.button_restore.setSizePolicy(sizePolicy)
        self.button_restore.setObjectName("button_restore")
        self.dialog_buttons.addWidget(self.button_restore)
        self.button_exit = QtWidgets.QPushButton(self.gui)

        #sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        #sizePolicy.setHorizontalStretch(0)
        #sizePolicy.setVerticalStretch(0)
        #sizePolicy.setHeightForWidth(self.button_exit.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.button_exit)

        self.button_exit.setSizePolicy(sizePolicy)
        self.button_exit.setText("Exit")
        self.button_exit.setObjectName("button_exit")
        self.dialog_buttons.addWidget(self.button_exit)
        self.gridLayout.addLayout(self.dialog_buttons, 3, 0, 1, 1)
        self.tabs = QtWidgets.QTabWidget(self.gui)

        #sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        #sizePolicy.setHorizontalStretch(0)
        #sizePolicy.setVerticalStretch(0)
        #sizePolicy.setHeightForWidth(self.tabs.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.tabs)

        self.tabs.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.tabs.setFont(font)
        self.tabs.setObjectName("tabs")
        self.main = QtWidgets.QWidget()
        self.main.setObjectName("main")
        self.verticalLayout_8 = QtWidgets.QVBoxLayout(self.main)
        self.verticalLayout_8.setObjectName("verticalLayout_8")

        self.parameters = QtWidgets.QGroupBox(self.main)
        sizePolicy = self._setPolicy(self.parameters)
        self.parameters.setSizePolicy(sizePolicy)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)

        self.parameters.setFont(font)
        self.parameters.setObjectName("parameters")

        self.verticalLayout = QtWidgets.QVBoxLayout(self.parameters)
        self.verticalLayout.setSpacing(3)
        self.verticalLayout.setObjectName("verticalLayout")

        self.hframe1 = QtWidgets.QFrame(self.parameters)
        self.hframe1.setObjectName("hframe1")

        self.horizontalLayout_14 = QtWidgets.QHBoxLayout(self.hframe1)
        self.horizontalLayout_14.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_14.setObjectName("horizontalLayout_14")

        self.input_label = QtWidgets.QLabel(self.hframe1)
        sizePolicy = self._setPolicy(self.input_label)
        self.input_label.setSizePolicy(sizePolicy)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        
        self.input_label.setFont(font)
        self.input_label.setMouseTracking(False)
        self.input_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.input_label.setTextFormat(QtCore.Qt.PlainText)
        self.input_label.setObjectName("input_label")

        self.horizontalLayout_14.addWidget(self.input_label)

        self.input = CheckableComboBox(self.hframe1)
        sizePolicy = self._setPolicy(self.input)
        self.input.setSizePolicy(sizePolicy)
        self.input.setObjectName("input")

        self.horizontalLayout_14.addWidget(self.input)

        self.refresh_input = QtWidgets.QPushButton(self.hframe1)


        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)

        self.refresh_input.setFont(font)
        self.refresh_input.setObjectName("refresh_input")

        sizePolicy = self._setPolicy(self.refresh_input)
        self.refresh_input.setSizePolicy(sizePolicy)
        
        self.horizontalLayout_14.addWidget(self.refresh_input)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_14.addItem(spacerItem)

        self.verticalLayout.addWidget(self.hframe1)

        self.hframe1_5 = QtWidgets.QHBoxLayout()
        self.hframe1_5.setObjectName("hframe1_5")

        self.regionOption_frame = QtWidgets.QFrame(self.parameters)
        self.regionOption_frame.setObjectName("regionoption_frame")

        self.hL_Option = QtWidgets.QHBoxLayout(self.regionOption_frame)  
        self.hL_Option.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.hL_Option.setObjectName("hL_Option")
        
        self.regionOption_label = QtWidgets.QLabel("Structure: ")
        sizePolicy = self._setPolicy(self.regionOption_label)
        self.regionOption_label.setSizePolicy(sizePolicy)      

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        
        self.regionOption_label.setFont(font)
        self.regionOption_label.setMouseTracking(False)
        self.regionOption_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.regionOption_label.setTextFormat(QtCore.Qt.PlainText)
        self.regionOption_label.setObjectName("regionOption_label")

        self.hL_Option.addWidget(self.regionOption_label)

        self.regionOption_rbtn1 = QtWidgets.QRadioButton("Default")
        self.regionOption_rbtn1.setChecked(True)
        self.regionOption_rbtn1.setAutoExclusive(True)
        self.regionOption_rbtn2  = QtWidgets.QRadioButton("Selected")
        self.regionOption_rbtn3 = QtWidgets.QRadioButton("Protein")
        self.regionOption_rbtn4  = QtWidgets.QRadioButton("All ligands without solvent")

        self.groupButton = QtWidgets.QButtonGroup(self.regionOption_frame)

        self.groupButton.addButton(self.regionOption_rbtn1)
        self.groupButton.addButton(self.regionOption_rbtn2)
        self.groupButton.addButton(self.regionOption_rbtn3)
        self.groupButton.addButton(self.regionOption_rbtn4)

        self.hL_Option.addWidget(self.regionOption_rbtn1)
        self.hL_Option.addWidget(self.regionOption_rbtn2)
        self.hL_Option.addWidget(self.regionOption_rbtn3)
        self.hL_Option.addWidget(self.regionOption_rbtn4)

        self.ignore_backbone_checkbox = QtWidgets.QCheckBox('Ignore Backbone')
        self.ignore_backbone_checkbox.setChecked(False)

        sizePolicy = self._setPolicy(self.ignore_backbone_checkbox)
        self.ignore_backbone_checkbox.setSizePolicy(sizePolicy)

        self.hL_Option.addStretch(1)

        self.hL_Option.addWidget(self.ignore_backbone_checkbox)

        self.hL_Option.addStretch(3)

        self.hframe1_5.addWidget(self.regionOption_frame)

        self.verticalLayout.addLayout(self.hframe1_5)

        self.hframe2 = QtWidgets.QHBoxLayout()
        self.hframe2.setObjectName("hframe2")

        self.step_size_frame = QtWidgets.QFrame(self.parameters)
        self.step_size_frame.setObjectName("step_size_frame")

        self.horizontalLayout_20 = QtWidgets.QHBoxLayout(self.step_size_frame)
        self.horizontalLayout_20.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_20.setObjectName("horizontalLayout_20")

        self.step_size_label = QtWidgets.QLabel(self.step_size_frame)
        self.step_size_label.setObjectName("step_size_label")
        self.step_size_label.setToolTip("<ul><li>Low: 0.6 Å</li><li>Medium: 0.5 Å</li><li>High: 0.25 Å</li></ul>")
        self.horizontalLayout_20.addWidget(self.step_size_label)

        self.step_size = QtWidgets.QDoubleSpinBox(self.step_size_frame)
        sizePolicy = self._setPolicy(self.step_size)
        self.step_size.setSizePolicy(sizePolicy)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)

        self.step_size.setFont(font)
        self.step_size.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.step_size.setDecimals(1)
        self.step_size.setMaximum(20.0)
        self.step_size.setSingleStep(0.1)
        self.step_size.setProperty("value", 0.0)
        self.step_size.setObjectName("step_size")

        self.horizontalLayout_20.addWidget(self.step_size)

        self.hframe2.addWidget(self.step_size_frame)

        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.hframe2.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.hframe2)
        self.hframe3 = QtWidgets.QHBoxLayout()
        self.hframe3.setObjectName("hframe3")
        self.probe_in_frame = QtWidgets.QFrame(self.parameters)
        self.probe_in_frame.setObjectName("probe_in_frame")
        self.horizontalLayout_13 = QtWidgets.QHBoxLayout(self.probe_in_frame)
        self.horizontalLayout_13.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_13.setObjectName("horizontalLayout_13")

        self.probe_in_label = QtWidgets.QLabel(self.probe_in_frame)
        sizePolicy = self._setPolicy(self.probe_in_label)
        self.probe_in_label.setSizePolicy(sizePolicy)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)

        self.probe_in_label.setFont(font)
        self.probe_in_label.setMouseTracking(True)
        self.probe_in_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.probe_in_label.setTextFormat(QtCore.Qt.RichText)
        self.probe_in_label.setObjectName("probe_in_label")
        self.horizontalLayout_13.addWidget(self.probe_in_label)

        self.probe_in = QtWidgets.QDoubleSpinBox(self.probe_in_frame)
        sizePolicy = self._setPolicy(self.probe_in)
        self.probe_in.setSizePolicy(sizePolicy)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.probe_in.setFont(font)
        self.probe_in.setDecimals(1)
        self.probe_in.setMaximum(5.0)
        self.probe_in.setSingleStep(0.1)
        self.probe_in.setProperty("value", 1.4)
        self.probe_in.setObjectName("probe_in")
        self.horizontalLayout_13.addWidget(self.probe_in)
        self.hframe3.addWidget(self.probe_in_frame)
        self.probe_out_frame = QtWidgets.QFrame(self.parameters)
        self.probe_out_frame.setObjectName("probe_out_frame")
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout(self.probe_out_frame)
        self.horizontalLayout_10.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")

        self.probe_out_label = QtWidgets.QLabel(self.probe_out_frame)
        sizePolicy = self._setPolicy(self.probe_out_label)
        self.probe_out_label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.probe_out_label.setFont(font)
        self.probe_out_label.setMouseTracking(True)
        self.probe_out_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.probe_out_label.setTextFormat(QtCore.Qt.RichText)
        self.probe_out_label.setObjectName("probe_out_label")
        self.horizontalLayout_10.addWidget(self.probe_out_label)
        self.probe_out = QtWidgets.QDoubleSpinBox(self.probe_out_frame)

        # sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        # sizePolicy.setHorizontalStretch(0)
        # sizePolicy.setVerticalStretch(0)
        # sizePolicy.setHeightForWidth(self.probe_out.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.probe_out)

        self.probe_out.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.probe_out.setFont(font)
        self.probe_out.setDecimals(1)
        self.probe_out.setMaximum(50.0)
        self.probe_out.setSingleStep(0.1)
        self.probe_out.setProperty("value", 4.0)
        self.probe_out.setObjectName("probe_out")
        self.horizontalLayout_10.addWidget(self.probe_out)
        self.hframe3.addWidget(self.probe_out_frame)
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.hframe3.addItem(spacerItem2)
        self.verticalLayout.addLayout(self.hframe3)
        self.hframe4_2 = QtWidgets.QHBoxLayout()
        self.hframe4_2.setObjectName("hframe4_2")
        self.removal_distance_frame = QtWidgets.QFrame(self.parameters)
        self.removal_distance_frame.setObjectName("removal_distance_frame")
        self.horizontalLayout_16 = QtWidgets.QHBoxLayout(self.removal_distance_frame)
        self.horizontalLayout_16.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_16.setObjectName("horizontalLayout_16")
        self.removal_distance_label = QtWidgets.QLabel(self.removal_distance_frame)

        # sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        # sizePolicy.setHorizontalStretch(0)
        # sizePolicy.setVerticalStretch(0)
        # sizePolicy.setHeightForWidth(self.removal_distance_label.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.removal_distance_label)

        self.removal_distance_label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.removal_distance_label.setFont(font)
        self.removal_distance_label.setMouseTracking(True)
        self.removal_distance_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.removal_distance_label.setTextFormat(QtCore.Qt.RichText)
        self.removal_distance_label.setObjectName("removal_distance_label")
        self.horizontalLayout_16.addWidget(self.removal_distance_label)
        self.removal_distance = QtWidgets.QDoubleSpinBox(self.removal_distance_frame)

        # sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        # sizePolicy.setHorizontalStretch(0)
        # sizePolicy.setVerticalStretch(0)
        # sizePolicy.setHeightForWidth(self.removal_distance.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.removal_distance)

        self.removal_distance.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.removal_distance.setFont(font)
        self.removal_distance.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.removal_distance.setDecimals(1)
        self.removal_distance.setMaximum(10.0)
        self.removal_distance.setSingleStep(0.1)
        self.removal_distance.setProperty("value", 2.4)
        self.removal_distance.setObjectName("removal_distance")
        self.horizontalLayout_16.addWidget(self.removal_distance)
        self.hframe4_2.addWidget(self.removal_distance_frame)
        self.volume_cutoff_frame = QtWidgets.QFrame(self.parameters)
        self.volume_cutoff_frame.setObjectName("volume_cutoff_frame")
        self.horizontalLayout_17 = QtWidgets.QHBoxLayout(self.volume_cutoff_frame)
        self.horizontalLayout_17.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_17.setObjectName("horizontalLayout_17")
        self.volume_cutoff_label = QtWidgets.QLabel(self.volume_cutoff_frame)
        
        # sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        # sizePolicy.setHorizontalStretch(0)
        # sizePolicy.setVerticalStretch(0)
        # sizePolicy.setHeightForWidth(self.volume_cutoff_label.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.volume_cutoff_label)

        self.volume_cutoff_label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.volume_cutoff_label.setFont(font)
        self.volume_cutoff_label.setMouseTracking(True)
        self.volume_cutoff_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.volume_cutoff_label.setTextFormat(QtCore.Qt.RichText)
        self.volume_cutoff_label.setObjectName("volume_cutoff_label")
        self.horizontalLayout_17.addWidget(self.volume_cutoff_label)
        self.volume_cutoff = QtWidgets.QDoubleSpinBox(self.volume_cutoff_frame)

        # sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        # sizePolicy.setHorizontalStretch(0)
        # sizePolicy.setVerticalStretch(0)
        # sizePolicy.setHeightForWidth(self.volume_cutoff.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.volume_cutoff)

        self.volume_cutoff.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.volume_cutoff.setFont(font)
        self.volume_cutoff.setDecimals(1)
        self.volume_cutoff.setMaximum(1000000000.0)
        self.volume_cutoff.setSingleStep(1.0)
        self.volume_cutoff.setProperty("value", 5.0)
        self.volume_cutoff.setObjectName("volume_cutoff")
        self.horizontalLayout_17.addWidget(self.volume_cutoff)
        self.hframe4_2.addWidget(self.volume_cutoff_frame)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.hframe4_2.addItem(spacerItem3)
        self.verticalLayout.addLayout(self.hframe4_2)

        self.hframe_5 = QtWidgets.QHBoxLayout()
        self.hframe_5.setObjectName("hframe_5")

        self.surface_representation_frame = QtWidgets.QFrame(self.parameters)
        self.surface_representation_frame.setObjectName("surface_representation_frame")
        
        self.horizontalLayout_26 = QtWidgets.QHBoxLayout(self.surface_representation_frame)
        self.horizontalLayout_26.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_26.setObjectName("horizontalLayout_26")

        self.surface_label = QtWidgets.QLabel(self.surface_representation_frame)
        sizePolicy = self._setPolicy(self.surface_label)
        self.surface_label.setSizePolicy(sizePolicy)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.surface_label.setFont(font)
        self.surface_label.setMouseTracking(True)
        self.surface_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.surface_label.setTextFormat(QtCore.Qt.RichText)
        self.surface_label.setObjectName("surface_label")
        self.horizontalLayout_26.addWidget(self.surface_label)
        self.surface = QtWidgets.QComboBox(self.surface_representation_frame)
        self.surface.setObjectName("surface")
        self.surface.addItem("")
        self.surface.addItem("")
        self.horizontalLayout_26.addWidget(self.surface)
        self.hframe_5.addWidget(self.surface_representation_frame)

        # self.cavity_representation_frame = QtWidgets.QFrame(self.parameters)
        # self.cavity_representation_frame.setObjectName("cavity_representation_frame")
        # self.horizontalLayout_11 = QtWidgets.QHBoxLayout(self.cavity_representation_frame)
        # self.horizontalLayout_11.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        # self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        # self.cavity_representation_label = QtWidgets.QLabel(self.cavity_representation_frame)

        # # sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        # # sizePolicy.setHorizontalStretch(0)
        # # sizePolicy.setVerticalStretch(0)
        # # sizePolicy.setHeightForWidth(self.cavity_representation_label.sizePolicy().hasHeightForWidth())

        # sizePolicy = self._setPolicy(self.cavity_representation_label)

        # self.cavity_representation_label.setSizePolicy(sizePolicy)
        # font = QtGui.QFont()
        # font.setPointSize(10)
        # font.setBold(False)
        # font.setItalic(False)
        # font.setWeight(50)
        # font.setKerning(True)
        # self.cavity_representation_label.setFont(font)
        # self.cavity_representation_label.setMouseTracking(True)
        # self.cavity_representation_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        # self.cavity_representation_label.setTextFormat(QtCore.Qt.RichText)
        # self.cavity_representation_label.setObjectName("cavity_representation_label")
        # self.horizontalLayout_11.addWidget(self.cavity_representation_label)
        # self.cavity_representation = QtWidgets.QComboBox(self.cavity_representation_frame)
        # self.cavity_representation.setObjectName("cavity_representation")
        # self.cavity_representation.addItem("")
        # self.cavity_representation.addItem("")
        # self.horizontalLayout_11.addWidget(self.cavity_representation)
        # self.hframe_5.addWidget(self.cavity_representation_frame)
        spacerItem4 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.hframe_5.addItem(spacerItem4)
        self.verticalLayout.addLayout(self.hframe_5)

        self.hframe6_2 = QtWidgets.QFrame(self.parameters)
        self.hframe6_2.setObjectName("hframe6_2")
        self.horizontalLayout_15 = QtWidgets.QHBoxLayout(self.hframe6_2)
        self.horizontalLayout_15.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_15.setObjectName("horizontalLayout_15")

        self.output_base_name_label = QtWidgets.QLabel(self.hframe6_2)
        sizePolicy = self._setPolicy(self.output_base_name_label)
        self.output_base_name_label.setSizePolicy(sizePolicy)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)

        self.output_base_name_label.setFont(font)
        self.output_base_name_label.setMouseTracking(False)
        self.output_base_name_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.output_base_name_label.setTextFormat(QtCore.Qt.PlainText)
        self.output_base_name_label.setObjectName("output_base_name_label")
        self.horizontalLayout_15.addWidget(self.output_base_name_label)
        self.base_name = QtWidgets.QLineEdit(self.hframe6_2)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)

        self.base_name.setFont(font)
        self.base_name.setText("output")
        self.base_name.setCursorMoveStyle(QtCore.Qt.VisualMoveStyle)
        self.base_name.setClearButtonEnabled(True)
        self.base_name.setObjectName("base_name")
        self.horizontalLayout_15.addWidget(self.base_name)
        spacerItem5 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_15.addItem(spacerItem5)
        self.verticalLayout.addWidget(self.hframe6_2)

        self.hframe7_2 = QtWidgets.QFrame(self.parameters)
        self.hframe7_2.setObjectName("hframe7_2")
        self.horizontalLayout_12 = QtWidgets.QHBoxLayout(self.hframe7_2)
        self.horizontalLayout_12.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_12.setObjectName("horizontalLayout_12")

        self.output_dir_label = QtWidgets.QLabel(self.hframe7_2)
        sizePolicy = self._setPolicy(self.output_dir_label)
        self.output_dir_label.setSizePolicy(sizePolicy)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)

        self.output_dir_label.setFont(font)
        self.output_dir_label.setMouseTracking(False)
        self.output_dir_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.output_dir_label.setTextFormat(QtCore.Qt.PlainText)
        self.output_dir_label.setObjectName("output_dir_label")
        self.horizontalLayout_12.addWidget(self.output_dir_label)
        self.output_dir_path = QtWidgets.QLineEdit(self.hframe7_2)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)

        self.output_dir_path.setFont(font)
        self.output_dir_path.setText("")
        self.output_dir_path.setEchoMode(QtWidgets.QLineEdit.Normal)
        self.output_dir_path.setReadOnly(True)
        self.output_dir_path.setClearButtonEnabled(False)
        self.output_dir_path.setObjectName("output_dir_path")

        self.horizontalLayout_12.addWidget(self.output_dir_path)
        self.button_browse = QtWidgets.QPushButton(self.hframe7_2)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)

        self.button_browse.setFont(font)
        self.button_browse.setText("Browse...")
        self.button_browse.setObjectName("button_browse")

        self.horizontalLayout_12.addWidget(self.button_browse)

        self.verticalLayout.addWidget(self.hframe7_2)


        # self.file_locations = QtWidgets.QGroupBox(self.main)
        # self.file_locations.setObjectName("file_locations")

        # self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.file_locations)
        # self.verticalLayout_3.setSpacing(3)
        # self.verticalLayout_3.setObjectName("verticalLayout_3")

        # self.hframe5_2 = QtWidgets.QFrame(self.parameters)
        # self.hframe5_2.setObjectName("hframe5_2")

        # self.horizontalLayout_24 = QtWidgets.QHBoxLayout(self.hframe5_2)
        # self.horizontalLayout_24.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        # self.horizontalLayout_24.setObjectName("horizontalLayout_24")
        
        # self.verticalLayout_3.addWidget(self.hframe5_2)

        self.hframe5_3 = QtWidgets.QFrame(self.parameters)
        self.hframe5_3.setObjectName("hframe5_3")

        self.horizontalLayout_25 = QtWidgets.QHBoxLayout(self.hframe5_3)
        self.horizontalLayout_25.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_25.setObjectName("horizontalLayout_25")

        self.dictionary_label = QtWidgets.QLabel(self.hframe5_3)
        sizePolicy = self._setPolicy(self.dictionary_label)
        self.dictionary_label.setSizePolicy(sizePolicy)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)

        self.dictionary_label.setFont(font)
        self.dictionary_label.setMouseTracking(False)
        self.dictionary_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.dictionary_label.setTextFormat(QtCore.Qt.PlainText)
        self.dictionary_label.setObjectName("dictionary_label")
        self.horizontalLayout_25.addWidget(self.dictionary_label)
        self.dictionary = QtWidgets.QLineEdit(self.hframe5_3)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)

        self.dictionary.setFont(font)
        self.dictionary.setText("")
        self.dictionary.setEchoMode(QtWidgets.QLineEdit.Normal)
        self.dictionary.setReadOnly(True)
        self.dictionary.setClearButtonEnabled(False)
        self.dictionary.setObjectName("dictionary")
        self.horizontalLayout_25.addWidget(self.dictionary)
        self.button_browse3 = QtWidgets.QPushButton(self.hframe5_3)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)

        self.button_browse3.setFont(font)
        self.button_browse3.setText("Browse...")
        self.button_browse3.setObjectName("button_browse3")
        self.horizontalLayout_25.addWidget(self.button_browse3)
        self.verticalLayout.addWidget(self.hframe5_3)
        self.verticalLayout_8.addWidget(self.parameters)


        spacerItem6 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_8.addItem(spacerItem6)
        self.tabs.addTab(self.main, "")

        self.search_space = QtWidgets.QWidget()
        self.search_space.setObjectName("search_space")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.search_space)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.box_adjustment = QtWidgets.QGroupBox(self.search_space)
        self.box_adjustment.setCheckable(True)
        self.box_adjustment.setChecked(True)
        self.box_adjustment.setObjectName("box_adjustment")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.box_adjustment)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.hframe12 = QtWidgets.QHBoxLayout()
        self.hframe12.setObjectName("hframe12")
        self.min_y_label = QtWidgets.QLabel(self.box_adjustment)
        self.min_y_label.setTextFormat(QtCore.Qt.RichText)
        self.min_y_label.setAlignment(QtCore.Qt.AlignCenter)
        self.min_y_label.setObjectName("min_y_label")
        self.hframe12.addWidget(self.min_y_label)
        self.min_y = QtWidgets.QDoubleSpinBox(self.box_adjustment)
        self.min_y.setEnabled(False)
        self.min_y.setDecimals(1)
        self.min_y.setMaximum(50.0)
        self.min_y.setSingleStep(0.1)
        self.min_y.setObjectName("min_y")
        self.hframe12.addWidget(self.min_y)
        self.gridLayout_3.addLayout(self.hframe12, 6, 0, 1, 1)
        self.hframe13 = QtWidgets.QHBoxLayout()
        self.hframe13.setObjectName("hframe13")
        self.min_z_label = QtWidgets.QLabel(self.box_adjustment)
        self.min_z_label.setTextFormat(QtCore.Qt.RichText)
        self.min_z_label.setAlignment(QtCore.Qt.AlignCenter)
        self.min_z_label.setObjectName("min_z_label")
        self.hframe13.addWidget(self.min_z_label)
        self.min_z = QtWidgets.QDoubleSpinBox(self.box_adjustment)
        self.min_z.setEnabled(False)
        self.min_z.setDecimals(1)
        self.min_z.setMaximum(50.0)
        self.min_z.setSingleStep(0.1)
        self.min_z.setObjectName("min_z")
        self.hframe13.addWidget(self.min_z)
        self.gridLayout_3.addLayout(self.hframe13, 8, 0, 1, 1)
        self.hframe15 = QtWidgets.QHBoxLayout()
        self.hframe15.setObjectName("hframe15")
        self.angle1_label = QtWidgets.QLabel(self.box_adjustment)
        self.angle1_label.setTextFormat(QtCore.Qt.RichText)
        self.angle1_label.setAlignment(QtCore.Qt.AlignCenter)
        self.angle1_label.setObjectName("angle1_label")
        self.hframe15.addWidget(self.angle1_label)
        self.angle1 = QtWidgets.QDoubleSpinBox(self.box_adjustment)
        self.angle1.setEnabled(False)
        self.angle1.setDecimals(0)
        self.angle1.setMaximum(180.0)
        self.angle1.setSingleStep(1.0)
        self.angle1.setObjectName("angle1")
        self.hframe15.addWidget(self.angle1)
        self.gridLayout_3.addLayout(self.hframe15, 10, 0, 1, 1)
        self.hframe7 = QtWidgets.QHBoxLayout()
        self.hframe7.setObjectName("hframe7")
        self.button_draw_box = QtWidgets.QPushButton(self.box_adjustment)
        self.button_draw_box.setObjectName("button_draw_box")
        self.hframe7.addWidget(self.button_draw_box)
        self.button_delete_box = QtWidgets.QPushButton(self.box_adjustment)
        self.button_delete_box.setObjectName("button_delete_box")
        self.hframe7.addWidget(self.button_delete_box)
        self.button_redraw_box = QtWidgets.QPushButton(self.box_adjustment)
        self.button_redraw_box.setEnabled(False)
        self.button_redraw_box.setObjectName("button_redraw_box")
        self.hframe7.addWidget(self.button_redraw_box)
        self.gridLayout_3.addLayout(self.hframe7, 2, 0, 1, 1)
        self.hframe9 = QtWidgets.QHBoxLayout()
        self.hframe9.setObjectName("hframe9")
        self.min_x_label = QtWidgets.QLabel(self.box_adjustment)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.min_x_label.setFont(font)
        self.min_x_label.setStyleSheet("")
        self.min_x_label.setTextFormat(QtCore.Qt.RichText)
        self.min_x_label.setAlignment(QtCore.Qt.AlignCenter)
        self.min_x_label.setObjectName("min_x_label")
        self.hframe9.addWidget(self.min_x_label)
        self.min_x = QtWidgets.QDoubleSpinBox(self.box_adjustment)
        self.min_x.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.min_x.setFont(font)
        self.min_x.setStyleSheet("")
        self.min_x.setDecimals(1)
        self.min_x.setMaximum(50.0)
        self.min_x.setSingleStep(0.1)
        self.min_x.setObjectName("min_x")
        self.hframe9.addWidget(self.min_x)
        self.gridLayout_3.addLayout(self.hframe9, 4, 0, 1, 1)
        self.hframe6 = QtWidgets.QHBoxLayout()
        self.hframe6.setObjectName("hframe6")
        self.box_adjustment_label = QtWidgets.QLabel(self.box_adjustment)
        self.box_adjustment_label.setEnabled(True)

        sizePolicy = self._setPolicy(self.box_adjustment_label)

        self.box_adjustment_label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.box_adjustment_label.setFont(font)
        self.box_adjustment_label.setMouseTracking(False)
        self.box_adjustment_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.box_adjustment_label.setTextFormat(QtCore.Qt.PlainText)
        self.box_adjustment_label.setAlignment(QtCore.Qt.AlignCenter)
        self.box_adjustment_label.setObjectName("box_adjustment_label")
        self.hframe6.addWidget(self.box_adjustment_label)
        self.button_box_adjustment_help = QtWidgets.QToolButton(self.box_adjustment)
        self.button_box_adjustment_help.setMinimumSize(QtCore.QSize(30, 26))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        font.setKerning(True)
        self.button_box_adjustment_help.setFont(font)
        self.button_box_adjustment_help.setCursor(QtGui.QCursor(QtCore.Qt.WhatsThisCursor))
        self.button_box_adjustment_help.setFocusPolicy(QtCore.Qt.NoFocus)
        self.button_box_adjustment_help.setObjectName("button_box_adjustment_help")
        self.hframe6.addWidget(self.button_box_adjustment_help)
        self.gridLayout_3.addLayout(self.hframe6, 1, 0, 1, 1)
        self.hframe16 = QtWidgets.QHBoxLayout()
        self.hframe16.setObjectName("hframe16")
        self.angle2_label = QtWidgets.QLabel(self.box_adjustment)
        self.angle2_label.setTextFormat(QtCore.Qt.RichText)
        self.angle2_label.setAlignment(QtCore.Qt.AlignCenter)
        self.angle2_label.setObjectName("angle2_label")
        self.hframe16.addWidget(self.angle2_label)
        self.angle2 = QtWidgets.QDoubleSpinBox(self.box_adjustment)
        self.angle2.setEnabled(False)
        self.angle2.setDecimals(0)
        self.angle2.setMaximum(180.0)
        self.angle2.setSingleStep(1.0)
        self.angle2.setObjectName("angle2")
        self.hframe16.addWidget(self.angle2)
        self.gridLayout_3.addLayout(self.hframe16, 11, 0, 1, 1)
        self.hframe10 = QtWidgets.QHBoxLayout()
        self.hframe10.setObjectName("hframe10")
        self.max_x_label = QtWidgets.QLabel(self.box_adjustment)
        self.max_x_label.setTextFormat(QtCore.Qt.RichText)
        self.max_x_label.setAlignment(QtCore.Qt.AlignCenter)
        self.max_x_label.setObjectName("max_x_label")
        self.hframe10.addWidget(self.max_x_label)
        self.max_x = QtWidgets.QDoubleSpinBox(self.box_adjustment)
        self.max_x.setEnabled(False)
        self.max_x.setDecimals(1)
        self.max_x.setMaximum(50.0)
        self.max_x.setSingleStep(0.1)
        self.max_x.setObjectName("max_x")
        self.hframe10.addWidget(self.max_x)
        self.gridLayout_3.addLayout(self.hframe10, 5, 0, 1, 1)
        self.hframe14 = QtWidgets.QHBoxLayout()
        self.hframe14.setObjectName("hframe14")
        self.max_z_label = QtWidgets.QLabel(self.box_adjustment)
        self.max_z_label.setTextFormat(QtCore.Qt.RichText)
        self.max_z_label.setAlignment(QtCore.Qt.AlignCenter)
        self.max_z_label.setObjectName("max_z_label")
        self.hframe14.addWidget(self.max_z_label)
        self.max_z = QtWidgets.QDoubleSpinBox(self.box_adjustment)
        self.max_z.setEnabled(False)
        self.max_z.setDecimals(1)
        self.max_z.setMaximum(50.0)
        self.max_z.setSingleStep(0.1)
        self.max_z.setObjectName("max_z")
        self.hframe14.addWidget(self.max_z)
        self.gridLayout_3.addLayout(self.hframe14, 9, 0, 1, 1)
        self.hframe11 = QtWidgets.QHBoxLayout()
        self.hframe11.setObjectName("hframe11")
        self.max_y_label = QtWidgets.QLabel(self.box_adjustment)
        self.max_y_label.setTextFormat(QtCore.Qt.RichText)
        self.max_y_label.setAlignment(QtCore.Qt.AlignCenter)
        self.max_y_label.setObjectName("max_y_label")
        self.hframe11.addWidget(self.max_y_label)
        self.max_y = QtWidgets.QDoubleSpinBox(self.box_adjustment)
        self.max_y.setEnabled(False)
        self.max_y.setDecimals(1)
        self.max_y.setMaximum(50.0)
        self.max_y.setSingleStep(0.1)
        self.max_y.setObjectName("max_y")
        self.hframe11.addWidget(self.max_y)
        self.gridLayout_3.addLayout(self.hframe11, 7, 0, 1, 1)
        self.hframe8 = QtWidgets.QHBoxLayout()
        self.hframe8.setObjectName("hframe8")
        self.padding_label = QtWidgets.QLabel(self.box_adjustment)
        self.padding_label.setTextFormat(QtCore.Qt.RichText)
        self.padding_label.setAlignment(QtCore.Qt.AlignCenter)
        self.padding_label.setObjectName("padding_label")
        self.hframe8.addWidget(self.padding_label)
        self.padding = QtWidgets.QDoubleSpinBox(self.box_adjustment)
        self.padding.setDecimals(1)
        self.padding.setMaximum(10.0)
        self.padding.setSingleStep(0.1)
        self.padding.setProperty("value", 3.5)
        self.padding.setObjectName("padding")
        self.hframe8.addWidget(self.padding)
        self.gridLayout_3.addLayout(self.hframe8, 3, 0, 1, 1)
        spacerItem7 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_3.addItem(spacerItem7, 12, 0, 1, 1)
        self.gridLayout_4.addWidget(self.box_adjustment, 0, 0, 1, 1)
        self.ligand_adjustment = QtWidgets.QGroupBox(self.search_space)
        self.ligand_adjustment.setEnabled(True)
        self.ligand_adjustment.setCheckable(True)
        self.ligand_adjustment.setChecked(True)
        self.ligand_adjustment.setObjectName("ligand_adjustment")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.ligand_adjustment)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.hframe17 = QtWidgets.QFrame(self.ligand_adjustment)
        self.hframe17.setObjectName("hframe17")
        self.horizontalLayout_18 = QtWidgets.QHBoxLayout(self.hframe17)
        self.horizontalLayout_18.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_18.setObjectName("horizontalLayout_18")
        self.ligand_label = QtWidgets.QLabel(self.hframe17)

        # sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        # sizePolicy.setHorizontalStretch(0)
        # sizePolicy.setVerticalStretch(0)
        # sizePolicy.setHeightForWidth(self.ligand_label.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.ligand_label)

        self.ligand_label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.ligand_label.setFont(font)
        self.ligand_label.setMouseTracking(False)
        self.ligand_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.ligand_label.setTextFormat(QtCore.Qt.PlainText)
        self.ligand_label.setObjectName("ligand_label")
        self.horizontalLayout_18.addWidget(self.ligand_label)
        
        self.ligand = CheckableComboBox(self.hframe17)
        sizePolicy = self._setPolicy(self.ligand)
        self.ligand.setSizePolicy(sizePolicy)
        self.ligand.setObjectName("ligand")
        
        self.horizontalLayout_18.addWidget(self.ligand)
        self.extract = QtWidgets.QPushButton()
        self.refresh_ligand = QtWidgets.QPushButton(self.hframe17)
        
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.refresh_ligand.sizePolicy().hasHeightForWidth())
        self.refresh_ligand.setSizePolicy(sizePolicy)
        
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.refresh_ligand.setFont(font)
        self.refresh_ligand.setObjectName("refresh_ligand")
        self.horizontalLayout_18.addWidget(self.refresh_ligand)
        self.verticalLayout_2.addWidget(self.hframe17)
        self.hframe18 = QtWidgets.QFrame(self.ligand_adjustment)
        self.hframe18.setObjectName("hframe18")
        self.horizontalLayout_19 = QtWidgets.QHBoxLayout(self.hframe18)
        self.horizontalLayout_19.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_19.setObjectName("horizontalLayout_19")
        self.ligand_cutoff_label = QtWidgets.QLabel(self.hframe18)
        self.ligand_cutoff_label.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ligand_cutoff_label.sizePolicy().hasHeightForWidth())
        self.ligand_cutoff_label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.ligand_cutoff_label.setFont(font)
        self.ligand_cutoff_label.setMouseTracking(True)
        self.ligand_cutoff_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.ligand_cutoff_label.setTextFormat(QtCore.Qt.RichText)
        self.ligand_cutoff_label.setObjectName("ligand_cutoff_label")
        self.horizontalLayout_19.addWidget(self.ligand_cutoff_label)
        self.ligand_cutoff = QtWidgets.QDoubleSpinBox(self.hframe18)
        self.ligand_cutoff.setEnabled(True)

        sizePolicy = self._setPolicy(self.ligand_cutoff)

        self.ligand_cutoff.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.ligand_cutoff.setFont(font)
        self.ligand_cutoff.setDecimals(1)
        self.ligand_cutoff.setMaximum(1000000000.0)
        self.ligand_cutoff.setSingleStep(0.1)
        self.ligand_cutoff.setProperty("value", 5.0)
        self.ligand_cutoff.setObjectName("ligand_cutoff")
        self.horizontalLayout_19.addWidget(self.ligand_cutoff)
        spacerItem8 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_19.addItem(spacerItem8)
        self.verticalLayout_2.addWidget(self.hframe18)
        spacerItem9 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem9)
        self.gridLayout_4.addWidget(self.ligand_adjustment, 0, 1, 1, 1)
        self.tabs.addTab(self.search_space, "")
        self.results = QtWidgets.QWidget()
        self.results.setObjectName("results")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.results)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.show_descriptors = QtWidgets.QHBoxLayout()
        self.show_descriptors.setObjectName("show_descriptors")
        self.show_descriptors_label = QtWidgets.QLabel(self.results)

        sizePolicy = self._setPolicy(self.show_descriptors_label)

        self.show_descriptors_label.setSizePolicy(sizePolicy)
        self.show_descriptors_label.setObjectName("show_descriptors_label")
        self.show_descriptors.addWidget(self.show_descriptors_label)
        self.default_view = QtWidgets.QRadioButton(self.results)

        sizePolicy = self._setPolicy(self.default_view)

        self.default_view.setSizePolicy(sizePolicy)
        self.default_view.setChecked(True)
        self.default_view.setAutoExclusive(True)
        self.default_view.setObjectName("default_view")
        self.show_descriptors.addWidget(self.default_view)
        self.depth_view = QtWidgets.QRadioButton(self.results)

        sizePolicy = self._setPolicy(self.depth_view)

        self.depth_view.setSizePolicy(sizePolicy)
        self.depth_view.setAcceptDrops(False)
        self.depth_view.setObjectName("depth_view")
        self.show_descriptors.addWidget(self.depth_view)
        self.hydropathy_view = QtWidgets.QRadioButton(self.results)
           
        sizePolicy = self._setPolicy(self.hydropathy_view)

        self.hydropathy_view.setSizePolicy(sizePolicy)
        self.hydropathy_view.setObjectName("hydropathy_view")
        self.show_descriptors.addWidget(self.hydropathy_view)
        spacerItem10 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.show_descriptors.addItem(spacerItem10)
        self.gridLayout_5.addLayout(self.show_descriptors, 1, 0, 1, 1)
        self.results_information = QtWidgets.QGroupBox(self.results)
        self.results_information.setObjectName("results_information")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout(self.results_information)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.hframe26 = QtWidgets.QHBoxLayout()
        self.hframe26.setObjectName("hframe26")
        self.results_file_label = QtWidgets.QLabel(self.results_information)

        sizePolicy = self._setPolicy(self.results_file_label)

        self.results_file_label.setSizePolicy(sizePolicy)
        self.results_file_label.setObjectName("results_file_label")
        self.hframe26.addWidget(self.results_file_label)
        self.results_file_entry = QtWidgets.QLineEdit(self.results_information)
        self.results_file_entry.setReadOnly(True)
        self.results_file_entry.setObjectName("results_file_entry")
        self.hframe26.addWidget(self.results_file_entry)
        self.button_browse4 = QtWidgets.QPushButton(self.results_information)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.button_browse4.setFont(font)
        self.button_browse4.setText("Browse...")
        self.button_browse4.setObjectName("button_browse4")
        self.hframe26.addWidget(self.button_browse4)
        self.button_load_results = QtWidgets.QPushButton(self.results_information)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        font.setKerning(True)
        self.button_load_results.setFont(font)
        self.button_load_results.setText("Load")
        self.button_load_results.setObjectName("button_load_results")
        self.hframe26.addWidget(self.button_load_results)
        self.verticalLayout_7.addLayout(self.hframe26)
        self.hframe27 = QtWidgets.QHBoxLayout()
        self.hframe27.setObjectName("hframe27")
        self.input_file_label = QtWidgets.QLabel(self.results_information)
        self.input_file_label.setObjectName("input_file_label")
        self.hframe27.addWidget(self.input_file_label)
        self.input_file_entry = QtWidgets.QLineEdit(self.results_information)
        self.input_file_entry.setReadOnly(True)
        self.input_file_entry.setObjectName("input_file_entry")
        self.hframe27.addWidget(self.input_file_entry)
        self.verticalLayout_7.addLayout(self.hframe27)
        self.hframe28 = QtWidgets.QHBoxLayout()
        self.hframe28.setObjectName("hframe28")
        self.ligand_file_label = QtWidgets.QLabel(self.results_information)
        self.ligand_file_label.setObjectName("ligand_file_label")
        self.hframe28.addWidget(self.ligand_file_label)
        self.ligand_file_entry = QtWidgets.QLineEdit(self.results_information)
        self.ligand_file_entry.setReadOnly(True)
        self.ligand_file_entry.setObjectName("ligand_file_entry")
        self.hframe28.addWidget(self.ligand_file_entry)
        self.verticalLayout_7.addLayout(self.hframe28)
        self.hframe29 = QtWidgets.QHBoxLayout()
        self.hframe29.setObjectName("hframe29")
        self.cavities_file_label = QtWidgets.QLabel(self.results_information)
        self.cavities_file_label.setObjectName("cavities_file_label")
        self.hframe29.addWidget(self.cavities_file_label)
        self.cavities_file_entry = QtWidgets.QLineEdit(self.results_information)
        self.cavities_file_entry.setReadOnly(True)
        self.cavities_file_entry.setObjectName("cavities_file_entry")
        self.hframe29.addWidget(self.cavities_file_entry)
        self.verticalLayout_7.addLayout(self.hframe29)
        self.hframe30 = QtWidgets.QHBoxLayout()
        self.hframe30.setObjectName("hframe30")
        self.step_size_label_2 = QtWidgets.QLabel(self.results_information)
        self.step_size_label_2.setObjectName("step_size_label_2")
        self.hframe30.addWidget(self.step_size_label_2)
        self.step_size_entry = QtWidgets.QLineEdit(self.results_information)

        sizePolicy = self._setPolicy(self.step_size_entry)

        self.step_size_entry.setSizePolicy(sizePolicy)
        self.step_size_entry.setMaximumSize(QtCore.QSize(50, 16777215))
        self.step_size_entry.setText("")
        self.step_size_entry.setMaxLength(10)
        self.step_size_entry.setAlignment(QtCore.Qt.AlignCenter)
        self.step_size_entry.setReadOnly(True)
        self.step_size_entry.setObjectName("step_size_entry")
        self.hframe30.addWidget(self.step_size_entry)
        spacerItem11 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.hframe30.addItem(spacerItem11)
        self.verticalLayout_7.addLayout(self.hframe30)
        self.gridLayout_5.addWidget(self.results_information, 0, 0, 1, 1)
        self.descriptors = QtWidgets.QGroupBox(self.results)

        sizePolicy = self._setPolicy(self.descriptors)

        self.descriptors.setSizePolicy(sizePolicy)
        self.descriptors.setObjectName("descriptors")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.descriptors)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.vframe1 = QtWidgets.QVBoxLayout()
        self.vframe1.setObjectName("vframe1")
        self.volume_label = QtWidgets.QLabel(self.descriptors)

        sizePolicy = self._setPolicy(self.volume_label)

        self.volume_label.setSizePolicy(sizePolicy)
        self.volume_label.setAlignment(QtCore.Qt.AlignCenter)
        self.volume_label.setObjectName("volume_label")
        self.vframe1.addWidget(self.volume_label)
        self.volume_list = QtWidgets.QListWidget(self.descriptors)

        sizePolicy = self._setPolicy(self.volume_list)

        self.volume_list.setSizePolicy(sizePolicy)
        self.volume_list.setMinimumSize(QtCore.QSize(153, 0))
        self.volume_list.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.volume_list.setObjectName("volume_list")
        self.vframe1.addWidget(self.volume_list)
        self.horizontalLayout.addLayout(self.vframe1)
        self.vframe2 = QtWidgets.QVBoxLayout()
        self.vframe2.setObjectName("vframe2")
        self.area_label = QtWidgets.QLabel(self.descriptors)

        sizePolicy = self._setPolicy(self.area_label)

        self.area_label.setSizePolicy(sizePolicy)
        self.area_label.setAlignment(QtCore.Qt.AlignCenter)
        self.area_label.setObjectName("area_label")
        self.vframe2.addWidget(self.area_label)
        self.area_list = QtWidgets.QListWidget(self.descriptors)

        sizePolicy = self._setPolicy(self.area_list)

        self.area_list.setSizePolicy(sizePolicy)
        self.area_list.setMinimumSize(QtCore.QSize(153, 0))
        self.area_list.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.area_list.setObjectName("area_list")
        self.vframe2.addWidget(self.area_list)
        self.horizontalLayout.addLayout(self.vframe2)
        self.vframe4 = QtWidgets.QVBoxLayout()
        self.vframe4.setObjectName("vframe4")
        self.avg_depth_label = QtWidgets.QLabel(self.descriptors)
        self.avg_depth_label.setAlignment(QtCore.Qt.AlignCenter)
        self.avg_depth_label.setObjectName("avg_depth_label")
        self.vframe4.addWidget(self.avg_depth_label)
        self.avg_depth_list = QtWidgets.QListWidget(self.descriptors)

        sizePolicy = self._setPolicy(self.avg_depth_list)

        self.avg_depth_list.setSizePolicy(sizePolicy)
        self.avg_depth_list.setMinimumSize(QtCore.QSize(153, 0))
        self.avg_depth_list.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.avg_depth_list.setObjectName("avg_depth_list")
        self.vframe4.addWidget(self.avg_depth_list)
        self.horizontalLayout.addLayout(self.vframe4)
        self.vframe5 = QtWidgets.QVBoxLayout()
        self.vframe5.setObjectName("vframe5")
        self.max_depth_label = QtWidgets.QLabel(self.descriptors)
        self.max_depth_label.setAlignment(QtCore.Qt.AlignCenter)
        self.max_depth_label.setObjectName("max_depth_label")
        self.vframe5.addWidget(self.max_depth_label)
        self.max_depth_list = QtWidgets.QListWidget(self.descriptors)

        # sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        # sizePolicy.setHorizontalStretch(0)
        # sizePolicy.setVerticalStretch(0)
        # sizePolicy.setHeightForWidth(self.max_depth_list.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.max_depth_list)

        self.max_depth_list.setSizePolicy(sizePolicy)
        self.max_depth_list.setMinimumSize(QtCore.QSize(153, 0))
        self.max_depth_list.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.max_depth_list.setObjectName("max_depth_list")
        self.vframe5.addWidget(self.max_depth_list)
        self.horizontalLayout.addLayout(self.vframe5)
        self.verticalLayout_6 = QtWidgets.QVBoxLayout()
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.avg_hydropathy_label = QtWidgets.QLabel(self.descriptors)
        self.avg_hydropathy_label.setAlignment(QtCore.Qt.AlignCenter)
        self.avg_hydropathy_label.setObjectName("avg_hydropathy_label")
        self.verticalLayout_6.addWidget(self.avg_hydropathy_label)
        self.avg_hydropathy_list = QtWidgets.QListWidget(self.descriptors)

        # sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        # sizePolicy.setHorizontalStretch(0)
        # sizePolicy.setVerticalStretch(0)
        # sizePolicy.setHeightForWidth(self.avg_hydropathy_list.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.avg_hydropathy_list)

        self.avg_hydropathy_list.setSizePolicy(sizePolicy)
        self.avg_hydropathy_list.setMinimumSize(QtCore.QSize(153, 0))
        self.avg_hydropathy_list.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.avg_hydropathy_list.setObjectName("avg_hydropathy_list")
        self.verticalLayout_6.addWidget(self.avg_hydropathy_list)
        self.horizontalLayout.addLayout(self.verticalLayout_6)
        self.vframe3 = QtWidgets.QVBoxLayout()
        self.vframe3.setObjectName("vframe3")
        self.residues_label = QtWidgets.QLabel(self.descriptors)

        # sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        # sizePolicy.setHorizontalStretch(0)
        # sizePolicy.setVerticalStretch(0)
        # sizePolicy.setHeightForWidth(self.residues_label.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.residues_label)

        self.residues_label.setSizePolicy(sizePolicy)
        self.residues_label.setAlignment(QtCore.Qt.AlignCenter)
        self.residues_label.setObjectName("residues_label")
        self.vframe3.addWidget(self.residues_label)
        self.residues_list = QtWidgets.QListWidget(self.descriptors)

        # sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        # sizePolicy.setHorizontalStretch(0)
        # sizePolicy.setVerticalStretch(0)
        # sizePolicy.setHeightForWidth(self.residues_list.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.residues_list)

        self.residues_list.setSizePolicy(sizePolicy)
        self.residues_list.setMinimumSize(QtCore.QSize(153, 0))
        self.residues_list.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.residues_list.setObjectName("residues_list")
        self.vframe3.addWidget(self.residues_list)
        self.horizontalLayout.addLayout(self.vframe3)
        self.gridLayout_5.addWidget(self.descriptors, 2, 0, 1, 1)
        self.tabs.addTab(self.results, "")

        self.about = QtWidgets.QWidget()
        self.about.setObjectName("about")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.about)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.about_text = QtWidgets.QTextBrowser(self.about)
        self.about_text.setEnabled(True)
        self.about_text.viewport().setProperty("cursor", QtGui.QCursor(QtCore.Qt.IBeamCursor))
        self.about_text.setStyleSheet("background-color: #d3d3d3;color:black; padding: 20px; font: 10pt \"Sans Serif\";")
        self.about_text.setFrameShape(QtWidgets.QFrame.Box)
        self.about_text.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustIgnored)
        self.about_text.setAutoFormatting(QtWidgets.QTextEdit.AutoBulletList)
        self.about_text.setLineWrapMode(QtWidgets.QTextEdit.WidgetWidth)
        self.about_text.setAcceptRichText(True)
        self.about_text.setOpenExternalLinks(True)
        self.about_text.setObjectName("about_text")
        self.gridLayout_2.addWidget(self.about_text, 0, 0, 1, 1)
        self.tabs.addTab(self.about, "")
        self.gridLayout.addWidget(self.tabs, 1, 0, 1, 1)
        self.main_description = QtWidgets.QLabel(self.gui)
        self.main_description.setEnabled(True)

        # sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        # sizePolicy.setHorizontalStretch(0)
        # sizePolicy.setVerticalStretch(0)
        # sizePolicy.setHeightForWidth(self.main_description.sizePolicy().hasHeightForWidth())

        sizePolicy = self._setPolicy(self.main_description)

        self.main_description.setSizePolicy(sizePolicy)
        self.main_description.setMinimumSize(QtCore.QSize(0, 0))
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.main_description.setFont(font)
        self.main_description.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.main_description.setStyleSheet("background-color: #d3d3d3;color:black; padding: 10px")
        self.main_description.setFrameShape(QtWidgets.QFrame.Box)
        self.main_description.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.main_description.setText("parKVFinder software identifies and describes cavities in a target biomolecular structure using a dual probe system.\n"
"\n"
"The description includes spatial, depth, constitutional and hydropathy characterization. The spatial description includes shape, volume, and area. The depth description defines depths for each cavity point, shown in the B-factor, and calculates the average and maximum depth per cavity. The constitutional description includes amino acids that form the identified cavities. The hydropathy description maps Eisenberg & Weiss hydrophobicity scale at surface points, shown in the Q-factor, and estimates average hydropathy per cavity.")
        self.main_description.setTextFormat(QtCore.Qt.PlainText)
        self.main_description.setScaledContents(False)
        self.main_description.setAlignment(QtCore.Qt.AlignJustify|QtCore.Qt.AlignVCenter)
        self.main_description.setWordWrap(True)
        self.main_description.setObjectName("main_description")
        self.gridLayout.addWidget(self.main_description, 0, 0, 1, 1)
        pyKVFinder.setCentralWidget(self.gui)

        self.retranslateUi(pyKVFinder)
        self.tabs.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(pyKVFinder)

        pyKVFinder.setTabOrder(self.tabs, self.button_run)
        pyKVFinder.setTabOrder(self.button_run, self.button_grid)
        pyKVFinder.setTabOrder(self.button_grid, self.button_restore)
        pyKVFinder.setTabOrder(self.button_restore, self.button_exit)
        pyKVFinder.setTabOrder(self.button_exit, self.input)
        pyKVFinder.setTabOrder(self.input, self.refresh_input)
        pyKVFinder.setTabOrder(self.refresh_input, self.regionOption_rbtn4)
        pyKVFinder.setTabOrder(self.regionOption_rbtn4, self.regionOption_rbtn3)
        pyKVFinder.setTabOrder(self.regionOption_rbtn3, self.regionOption_rbtn2)
        pyKVFinder.setTabOrder(self.regionOption_rbtn2, self.regionOption_rbtn1)
        pyKVFinder.setTabOrder(self.regionOption_rbtn1, self.regionOption_label)
        pyKVFinder.setTabOrder(self.regionOption_label, self.probe_out)
        pyKVFinder.setTabOrder(self.probe_out, self.probe_in)
        pyKVFinder.setTabOrder(self.probe_in, self.volume_cutoff)
        pyKVFinder.setTabOrder(self.volume_cutoff, self.removal_distance)
        pyKVFinder.setTabOrder(self.removal_distance, self.base_name)
        pyKVFinder.setTabOrder(self.base_name, self.output_dir_path)
        pyKVFinder.setTabOrder(self.output_dir_path, self.button_browse)
        pyKVFinder.setTabOrder(self.button_browse, self.box_adjustment)
        pyKVFinder.setTabOrder(self.box_adjustment, self.button_draw_box)
        pyKVFinder.setTabOrder(self.button_draw_box, self.button_delete_box)
        pyKVFinder.setTabOrder(self.button_delete_box, self.button_redraw_box)
        pyKVFinder.setTabOrder(self.button_redraw_box, self.padding)
        pyKVFinder.setTabOrder(self.padding, self.min_x)
        pyKVFinder.setTabOrder(self.min_x, self.max_x)
        pyKVFinder.setTabOrder(self.max_x, self.min_y)
        pyKVFinder.setTabOrder(self.min_y, self.max_y)
        pyKVFinder.setTabOrder(self.max_y, self.min_z)
        pyKVFinder.setTabOrder(self.min_z, self.max_z)
        pyKVFinder.setTabOrder(self.max_z, self.angle1)
        pyKVFinder.setTabOrder(self.angle1, self.angle2)
        pyKVFinder.setTabOrder(self.angle2, self.ligand_adjustment)
        pyKVFinder.setTabOrder(self.ligand_adjustment, self.ligand)
        pyKVFinder.setTabOrder(self.ligand, self.refresh_ligand)
        pyKVFinder.setTabOrder(self.refresh_ligand, self.ligand_cutoff)
        pyKVFinder.setTabOrder(self.ligand_cutoff, self.results_file_entry)
        pyKVFinder.setTabOrder(self.results_file_entry, self.button_browse4)
        pyKVFinder.setTabOrder(self.button_browse4, self.button_load_results)
        pyKVFinder.setTabOrder(self.button_load_results, self.volume_list)
        pyKVFinder.setTabOrder(self.volume_list, self.area_list)
        pyKVFinder.setTabOrder(self.area_list, self.residues_list)
        pyKVFinder.setTabOrder(self.residues_list, self.about_text)

        self.gui.adjustSize()
    def _setPolicy(self, element):
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(element.sizePolicy().hasHeightForWidth())
        return sizePolicy

    def retranslateUi(self, pyKVFinder):
        _translate = QtCore.QCoreApplication.translate
        self.button_grid.setText(_translate("pyKVFinder", "Show Grid"))
        self.button_restore.setText(_translate("pyKVFinder", "Restore Default Values"))
        self.parameters.setTitle(_translate("pyKVFinder", "Parameters"))
        self.input_label.setText(_translate("pyKVFinder", "Input PDB:"))
        self.refresh_input.setText(_translate("pyKVFinder", "Refresh"))
        # self.resolution_label.setText(_translate("pyKVFinder", "Resolution:"))
        # self.resolution.setItemText(0, _translate("pyKVFinder", "Low"))
        # self.resolution.setItemText(1, _translate("pyKVFinder", "Medium"))
        # self.resolution.setItemText(2, _translate("pyKVFinder", "High"))
        # self.resolution.setItemText(3, _translate("pyKVFinder", "Off"))
        self.step_size_label.setText(_translate("pyKVFinder", "Step Size (Å):"))
        self.probe_in_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Probe In (Å):</p></body></html>"))
        self.probe_out_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Probe Out (Å):</p></body></html>"))
        self.removal_distance_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Removal Distance (Å):</p></body></html>"))
        self.volume_cutoff_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Volume Cutoff (Å³):</p></body></html>"))
        self.surface_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Surface Representation:</p></body></html>"))
        self.surface.setItemText(0, _translate("pyKVFinder", "Solvent Excluded Surface (SES)"))
        self.surface.setItemText(1, _translate("pyKVFinder", "Solvent Accesible Surface (SAS)"))
        # self.cavity_representation_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Cavity Representation:</p></body></html>"))
        # self.cavity_representation.setItemText(0, _translate("pyKVFinder", "Filtered"))
        # self.cavity_representation.setItemText(1, _translate("pyKVFinder", "Full"))
        self.output_base_name_label.setText(_translate("pyKVFinder", "Output Base Name:"))
        self.output_dir_label.setText(_translate("pyKVFinder", "Output Directory:"))
        # self.file_locations.setTitle(_translate("pyKVFinder", "File Locations"))
        # self.parKVFinder_label.setText(_translate("pyKVFinder", "parKVFinder:"))
        self.dictionary_label.setText(_translate("pyKVFinder", "vdW dictionary:"))
        self.tabs.setTabText(self.tabs.indexOf(self.main), _translate("pyKVFinder", "Main"))
        self.box_adjustment.setTitle(_translate("pyKVFinder", "Box Adjustment"))
        self.min_y_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Minimum Y (Å):</p></body></html>"))
        self.min_z_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Minimum Z (Å):</p></body></html>"))
        self.angle1_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Angle 1 (°):</p></body></html>"))
        self.button_draw_box.setText(_translate("pyKVFinder", "Draw Box"))
        self.button_delete_box.setText(_translate("pyKVFinder", "Delete Box"))
        self.button_redraw_box.setText(_translate("pyKVFinder", "Redraw Box"))
        self.min_x_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Minimum X (Å):</p></body></html>"))
        self.box_adjustment_label.setText(_translate("pyKVFinder", "Select residues and press Draw Box:"))
        self.button_box_adjustment_help.setText(_translate("pyKVFinder", "?"))
        self.angle2_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Angle 2 (°):</p></body></html>"))
        self.max_x_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Maximum X (Å):</p></body></html>"))
        self.max_z_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Maximum Z (Å):</p></body></html>"))
        self.max_y_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Maximum Y (Å):</p></body></html>"))
        self.padding_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Padding (Å):</p></body></html>"))
        self.ligand_adjustment.setTitle(_translate("pyKVFinder", "Ligand Adjustment"))
        self.ligand_label.setText(_translate("pyKVFinder", "Ligand PDB:"))
        self.refresh_ligand.setText(_translate("pyKVFinder", "Refresh"))
        self.ligand_cutoff_label.setText(_translate("pyKVFinder", "<html><head/><body><p>Ligand Cutoff (Å):</p></body></html>"))
        self.tabs.setTabText(self.tabs.indexOf(self.search_space), _translate("pyKVFinder", "Search Space"))
        self.show_descriptors_label.setText(_translate("pyKVFinder", "Show descriptors:"))
        self.default_view.setText(_translate("pyKVFinder", "Default"))
        self.depth_view.setText(_translate("pyKVFinder", "Depth"))
        self.hydropathy_view.setText(_translate("pyKVFinder", "Hydropathy"))
        self.results_information.setTitle(_translate("pyKVFinder", "Information"))
        self.results_file_label.setText(_translate("pyKVFinder", "Results File:"))
        self.input_file_label.setText(_translate("pyKVFinder", "Input File:"))
        self.ligand_file_label.setText(_translate("pyKVFinder", "Ligand File:"))
        self.cavities_file_label.setText(_translate("pyKVFinder", "Cavities File:"))
        self.step_size_label_2.setText(_translate("pyKVFinder", "Step Size (Å):"))
        self.descriptors.setTitle(_translate("pyKVFinder", "Descriptors"))
        self.volume_label.setText(_translate("pyKVFinder", "Volume (Å³)"))
        self.area_label.setText(_translate("pyKVFinder", "<html>Surface Area (Å&#178;)<\\html>"))
        self.avg_depth_label.setText(_translate("pyKVFinder", "Average Depth (Å)"))
        self.max_depth_label.setText(_translate("pyKVFinder", "Maximum Depth (Å)"))
        self.avg_hydropathy_label.setText(_translate("pyKVFinder", "Average Hydropathy"))
        self.residues_label.setText(_translate("pyKVFinder", "Interface Residues"))
        self.tabs.setTabText(self.tabs.indexOf(self.results), _translate("pyKVFinder", "Results"))
        
        # TODO: Change the about text
        self.about_text.setHtml(_translate("pyKVFinder", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans Serif\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">PyMOL2 parKVFinder Tools integrates PyMOL v2.x (http://PyMOL.org/) with parKVFinder (https://github.com/LBC-LNBio/parKVFinder).</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">In the simplest case to run parKVFinder:</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">1) Load a target biomolecular structure into PyMOL v2.x.</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">2) Start PyMOL2 parKVFinder Tools plugin.</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">3) Select an input PDB on \'Main\' tab.</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">4) Ensure that parKVFinder executable path is correct on the &quot;Program Locations&quot; tab.</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">5) Click the &quot;Run parKVFinder&quot; button.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">Completed runs are available on \'Results\' tab, where users can check run information (i.e., input file, ligand file, output directory, step size) and spatial properties (i.e., volume, surface area and interface residues). In addition, the results can be loaded directly from a results file (.KVFinder.results.toml) by selecting a \'Results File\' by clicking on \'Browse ...\' and then clicking on \'Load\'. </span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">In addition to whole structure cavity detection, there are two search space adjustments: Box and Ligand adjustments.</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">- The \'Box adjustment\' mode creates a custom search box around a selection of interest by clicking on \'Draw Box\' button, which can be adapted by changing one box parameter (minimum and maximum XYZ, padding and angles) at a time by clicking on \'Redraw Box\'. For more information, there is a help button in \'Box adjustment\' group.</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">- The \'Ligand adjustment\' keeps cavity points around a target ligand PDB within a radius defined by the \'Ligand Cutoff\' parameter.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">parKVFinder and PyMOL2 parKVFinder Tools was developed by:</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">- João Victor da Silva Guerra</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">- Helder Veras Filho</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">- Leandro Oliveira Bortot</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">- Rodrigo Vargas Honorato</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">- José Geraldo de Carvalho Pereira</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">- Paulo Sergio Lopes de Oliveira (paulo.oliveira@lnbio.cnpem.br)</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">Brazilian Center for Research in Energy and Materials - CNPEM</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">Brazilian Biosciences National Laboratory - LNBio</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\">Please refer and cite our papers if you use it in a publication.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; font-weight:600; text-decoration: underline; color:#000000;\">Citations</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Droid Sans Mono\',\'monospace\',\'monospace\',\'Droid Sans Fallback\'; color:#000000;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">If you use <span style=\" text-decoration: underline;\">parKVFinder</span> or <span style=\" text-decoration: underline;\">PyMOL2 parKVFinder Tools</span>, please cite:</p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">João Victor da Silva Guerra, Helder Veras Ribeiro Filho, Leandro Oliveira Bortot, Rodrigo Vargas Honorato, José Geraldo de Carvalho Pereira, Paulo Sergio Lopes de Oliveira. ParKVFinder: A thread-level parallel approach in biomolecular cavity detection. SoftwareX (2020). <a href=\"https://doi.org/10.1016/j.softx.2020.100606\"><span style=\" text-decoration: underline; color:#0000ff;\">https://doi.org/10.1016/j.softx.2020.100606</span></a>.</p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">If you use <span style=\" text-decoration: underline;\">depth and hydropathy characterization</span>, please also cite:</p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Guerra, J.V.d., Ribeiro-Filho, H.V., Jara, G.E. et al. pyKVFinder: an efficient and integrable Python package for biomolecular cavity detection and characterization in data science. BMC Bioinformatics 22, 607 (2021). <a href=\"https://doi.org/10.1186/s12859-021-04519-4\"><span style=\" text-decoration: underline; color:#0000ff;\">https://doi.org/10.1186/s12859-021-04519-4</span></a>.</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">PyMOL citation may be found here:</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"http://pymol.sourceforge.net/faq.html#CITE\"><span style=\" text-decoration: underline; color:#0000ff;\">https://pymol.org/2/support.html?</span></a></p></body></html>"))
        self.tabs.setTabText(self.tabs.indexOf(self.about), _translate("pyKVFinder", "About"))
        
"""
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    pyKVFinder = QtWidgets.QMainWindow()
    ui = Ui_pyKVFinder()
    ui.setupUi(pyKVFinder)
    pyKVFinder.show()
    sys.exit(app.exec_())
"""

from PyQt5.QtWidgets import QApplication, QComboBox, QMainWindow
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QStyledItemDelegate, qApp
from PyQt5.QtGui import QStandardItemModel, QPalette, QStandardItem, QFontMetrics
from PyQt5.QtCore import Qt, QEvent
import typing
import sys
from io import StringIO
class CheckableComboBox(QComboBox):

    # Subclass Delegate to increase item height
    class Delegate(QStyledItemDelegate):
        def sizeHint(self, option, index):
            size = super().sizeHint(option, index)
            size.setHeight(20)
            return size

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Make the combo editable to set a custom text, but readonly
        self.setEditable(True)
        self.lineEdit().setReadOnly(True)
        # Make the lineedit the same color as QPushButton
        palette = qApp.palette()
        palette.setBrush(QPalette.Base, palette.button())
        self.lineEdit().setPalette(palette)

        # Use custom delegate
        self.setItemDelegate(CheckableComboBox.Delegate())

        # Update the text when an item is toggled
        self.model().dataChanged.connect(self.updateText)

        # Hide and show popup when clicking the line edit
        self.lineEdit().installEventFilter(self)
        self.closeOnLineEditClick = False

        # Prevent popup from closing when clicking on an item
        self.view().viewport().installEventFilter(self)

    def resizeEvent(self, event):
        # Recompute text to elide as needed
        self.updateText()
        super().resizeEvent(event)

    def eventFilter(self, object, event):

        if object == self.lineEdit():
            if event.type() == QEvent.MouseButtonRelease:
                if self.closeOnLineEditClick:
                    self.hidePopup()
                else:
                    self.showPopup()
                return True
            return False

        if object == self.view().viewport():
            if event.type() == QEvent.MouseButtonRelease:
                index = self.view().indexAt(event.pos())
                item = self.model().item(index.row())

                if item.checkState() == Qt.Checked:
                    item.setCheckState(Qt.Unchecked)
                else:
                    item.setCheckState(Qt.Checked)
                return True
        return False

    def showPopup(self):
        super().showPopup()
        # When the popup is displayed, a click on the lineedit should close it
        self.closeOnLineEditClick = True

    def hidePopup(self):
        super().hidePopup()
        # Used to prevent immediate reopening when clicking on the lineEdit
        self.startTimer(100)
        # Refresh the display text when closing
        self.updateText()

    def timerEvent(self, event):
        # After timeout, kill timer, and reenable click on line edit
        self.killTimer(event.timerId())
        self.closeOnLineEditClick = False

    def updateText(self):
        texts = []
        for i in range(self.model().rowCount()):
            if self.model().item(i).checkState() == Qt.Checked:
                texts.append(self.model().item(i).text())
        text = ", ".join(texts)

        # Compute elided text (with "...")
        metrics = QFontMetrics(self.lineEdit().font())
        elidedText = metrics.elidedText(text, Qt.ElideRight, self.lineEdit().width())
        self.lineEdit().setText(elidedText)

    def addItem(self, text, data=None):
        item = QStandardItem()
        item.setText(text)
        if data is None:
            item.setData(text)
        else:
            item.setData(data)
        item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsUserCheckable)
        item.setData(Qt.Unchecked, Qt.CheckStateRole)
        self.model().appendRow(item)

    def addItems(self, texts, datalist=None):
        for i, text in enumerate(texts):
            try:
                data = datalist[i]
            except (TypeError, IndexError):
                data = None
            self.addItem(text, data)

    def currentData(self):
        # Return the list of selected items data
        res = []
        for i in range(self.model().rowCount()):
            if self.model().item(i).checkState() == Qt.Checked:
                res.append(self.model().item(i).data())
        return res



class InputGUI():
    def __init__(self, parentWidget):
        self.parentWidget = parentWidget

    def readline(self):
        text, ok = QtWidgets.QInputDialog.getText(self.parentWidget, 'Introduce value', 'Value:')
        if ok: 
            return str(text)
        else:
            return ''
        
class SampleGUI(QtWidgets.QWidget):
    def __init__(self, tool):
        super(SampleGUI, self).__init__()
        self.tool = tool
        self.u = self.tool.ui
        self.s = self.tool.session
        self.all_objects = all_objects
        self.cx = cx
        self.log = ""
        
        self.initGUI()

    def initGUI(self):
        self.code = QtWidgets.QTextEdit()
        self.result = QtWidgets.QTextEdit()

        evalBtn = QtWidgets.QPushButton('Evaluate')
        evalBtn.clicked.connect(self.evaluate)

        clearBtn = QtWidgets.QPushButton('Clear')
        clearBtn.clicked.connect(self.clearT)

        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(evalBtn)
        hbox.addWidget(clearBtn)

        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.code)
        vbox.addLayout(hbox)
        vbox.addWidget(self.result)

        self.setLayout(vbox)
        #self.show()

    def evaluate(self):
        source_code = str(self.code.toPlainText())
        streams = sys.stdin, sys.stdout
        sys.stdin = InputGUI(self)                                                                                                                                                              
        redirected_output = sys.stdout = StringIO()
        exec(source_code)
        sys.stdin, sys.stdout = streams
        self.log += '>>> '+redirected_output.getvalue()+'\n'
        self.result.setText(self.log)

    def clearT(self):
        self.log = ""
        self.code.clear()
        self.result.clear()
