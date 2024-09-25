# Installation

First, you need to install PyMOL v2.x or v3.x on your computer. If you do not have it installed, please refer to PyMOL [website](https://pymol.org/).

_PyMOL pyKVFinder Tools_ is available to use with PyMOL and has been been developed using Qt interface and Python3.

After installing PyMOL, plug-in requires the installation of [pyKVFinder](https://pypi.org/project/pyKVFinder/), [toml](https://pypi.org/project/toml/), and [pyqt5](https://pypi.org/project/PyQt5/) modules and PyMOL's native Python does not have it installed. So you need to install it:

```bash
# PyMOL v3 (python 3)
$ pip3 install pyKVFinder toml pyqt5
```

_Note:_ [pip](https://pypi.org/project/pip/) package management system
installation is required.

Then, to install the _PyMOL pyKVFinder Tools_ on PyMOL, download the latest `PyMOL-pyKVFinder-Tools.zip` from [here](https://github.com/LBC-LNBio/pyKVFinder/releases/latest/download/PyMOL-pyKVFinder-Tools.zip) and follow these steps:

1. Open PyMOL.
2. Go to **Plugin** menu, click on **Plugin Manager**.
3. The **Plugin Manager** window will open, go to the **Install New Plugin** tab.
4. Under **Install from local file** panel, click on **Choose file...**.
5. The **Install Plugin** window will open, select the `PyMOL-pyKVFinder-Tools.zip` that you downloaded earliar.
6. The **Select plugin directory** window will open, select
   **/home/\<user\>/.pymol/startup** and click **OK**.
7. The **Confirm** window will open, click **OK**.
8. The **Sucess** window will appear, confirming that the plug-in has
   been installed.
9. Restart PyMOL.
10. _PyMOL pyKVFinder Tools_ is ready to use.

Or, instead of selecting `PyMOL-pyKVFinder-Tools.zip` (Step 5), you can select `__init__.py` file on `pyKVFinder/plugins/PyMOL-pyKVFinder-Tools` directory tree.
