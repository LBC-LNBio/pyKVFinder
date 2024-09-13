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

from PyQt5 import QtWidgets
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QMainWindow
from PyQt5.uic import loadUi


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

    def _initialize_gui(self) -> None:
        """
        Initialize the graphical user interface of the plugin.
        """
        # populate the QMainWindow from our *.ui file
        uifile = os.path.join(os.path.dirname(__file__), "PyMOL-pyKVFinder-Tools.ui")
        loadUi(uifile, self)

    
