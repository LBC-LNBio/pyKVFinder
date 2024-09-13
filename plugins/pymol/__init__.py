# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This is the source code for PyMOL pyKVFinder Tools. It was developed using the
PyMOL API, and is intended to be used as a plugin for PyMOL. The code is
written in Python and uses the PyQt5 library for the graphical user interface.
Changes in this file are not advised, as it controls all interactions with
pyKVFinder.
"""

from .gui import PyMOLpyKVFinderTools
from .io import *

# global reference to avoid garbage collection of our dialog
dialog = None


def __init_plugin__(app=None):
    """
    Add an entry to the PyMOL "Plugin" menu. This function is called when the
    plugin is loaded by PyMOL.
    """
    from pymol.plugins import addmenuitemqt

    addmenuitemqt("PyMOL pyKVFinder Tools", run_plugin_gui)


def run_plugin_gui():
    """
    Open our custom dialog. If it is already open, bring it to the front.
    """
    global dialog

    if dialog is None:
        dialog = PyMOLpyKVFinderTools()

    dialog.show()
