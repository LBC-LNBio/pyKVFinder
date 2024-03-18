# Installation

Before installing the _Chimerax pyKVFinder Tools_ plugin for ChimeraX, ensure that ChimeraX is installed on your system. If you haven't installed ChimeraX yet, please visit the [ChimeraX website](https://www.cgl.ucsf.edu/chimerax/download.html) for instructions on installation.

The _ChimeraX pyKVFinder Tools_ plugin is compatible with ChimeraX and has been developed using the Qt interface and Python 3.

Once ChimeraX is installed, the plugin requires the installation of the [toml](https://pypi.org/project/toml/) and [PyQt5](https://pypi.org/project/PyQt5/) modules. Ensure that you have these modules installed using your preferred package manager.

```bash
# Install required modules for ChimeraX (Python 3)
$ pip3 install toml PyQt5 pyKVFinder
# Alternatively, in a conda environment via Anaconda Cloud
$ conda install -c conda-forge toml
$ conda install -c anaconda pyqt 
```

_Note:_ ChimeraX has an internal Python interpreter. The _ChimeraX pyKVFinder Tools_ will attempt to install the required packages automatically. Restarting ChimeraX is necessary for the changes to take effect. If this automated installation process fails, please run the following code in the **ChimeraX CLI**:

```bash
# Execute this command in ChimeraX CLI
$ pip3 install pyKVFinder PyQt5 toml
```

_Note:_ Installation via [pip](https://pypi.org/project/pip/) or [Anaconda](https://www.anaconda.com/) package management system is required.

To install the _ChimeraX pyKVFinder Tools_ plugin for ChimeraX, you have two options. You can either use a precompiled **wheel** file: **`ChimeraX_pyKVFinder-0.1-py3-none-any.whl`**, or you can compile it directly within ChimeraX:

## Using precompiled wheel file

1. Download the file **`ChimeraX_pyKVFinder-0.1-py3-none-any.whl`**
2. Launc ChimeraX.
3. Run the following command in ChimeraX CLI:
```bash
$ toolshed install [path_to_binary]
```
4. Relaunch ChimeraX after installed.
5. The _ChimeraX pyKVFinder Tools_ will install the dependencies. Relaunch ChimeraX if asked.
6. The _ChimeraX pyKVFinder Tools_ plugin will now be ready for use.

## Building using ChimeraX

1. Download our repository from GitHub, or via the following command in the terminal: **`git clone [repo_link]`**
2. Launch ChimeraX.
3. Execute the following commands in the ChimeraX CLI:
```bash
$ devel build [repository_path]
$ devel install [repository_path]
```
4. Relaunch ChimeraX after installed.
5. The _ChimeraX pyKVFinder Tools_ will install the dependencies. Relaunch ChimeraX if asked.
6. The _ChimeraX pyKVFinder Tools_ plugin will now be ready for use.