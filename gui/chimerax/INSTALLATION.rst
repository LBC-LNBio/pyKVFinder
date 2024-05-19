Installation
============

Before installing the *Chimerax pyKVFinder Tools* plugin for ChimeraX, ensure that ChimeraX is installed on your system. If you haven't installed ChimeraX yet, please visit the `ChimeraX website <https://www.cgl.ucsf.edu/chimerax/download.html>`_ for instructions on installation.

The *ChimeraX pyKVFinder Tools* plugin is compatible with ChimeraX and has been developed using the Qt interface and Python 3.

Once ChimeraX is installed, the plugin requires the installation of the `toml <https://pypi.org/project/toml/>`_ and `PyQt5 <https://pypi.org/project/PyQt5/>`_ modules. Ensure that you have these modules installed using your preferred package manager.

.. code-block:: bash

    # Install required modules for ChimeraX (Python 3)
    $ pip3 install toml PyQt5 pyKVFinder
    # Alternatively, in a conda environment via Anaconda Cloud
    $ conda install -c conda-forge toml
    $ conda install -c anaconda pyqt

.. note::

    ChimeraX has an internal Python interpreter. The *ChimeraX pyKVFinder Tools* will attempt to install the required packages automatically. Restarting ChimeraX is necessary for the changes to take effect. If this automated installation process fails, please run the following code in the **ChimeraX CLI**:

.. code-block:: bash

    # Execute this command in ChimeraX CLI
    $ pip3 install pyKVFinder PyQt5 toml

.. note::

    Installation via `pip <https://pypi.org/project/pip/>`_ or `Anaconda <https://www.anaconda.com/>`_ package management system is required.

To install the *ChimeraX pyKVFinder Tools* plugin for ChimeraX, you have two options. You can either use a precompiled **wheel** file: **``ChimeraX_pyKVFinder-0.1-py3-none-any.whl``**, or you can compile it directly within ChimeraX:

Using precompiled wheel file
----------------------------

1. Download the file **``ChimeraX_pyKVFinder-0.1-py3-none-any.whl``**
2. Launch ChimeraX.
3. Run the following command in ChimeraX CLI:

.. code-block:: bash

    $ toolshed install [path_to_binary]

4. Relaunch ChimeraX after installed.
5. The *ChimeraX pyKVFinder Tools* will install the dependencies. Relaunch ChimeraX if asked.
6. The *ChimeraX pyKVFinder Tools* plugin will now be ready for use.

Building using ChimeraX
-----------------------

1. Download our repository from GitHub, or via the following command in the terminal: **``git clone [repo_link]``**
2. Launch ChimeraX.
3. Execute the following commands in the ChimeraX CLI:

.. code-block:: bash

    $ devel build [repository_path]
    $ devel install [repository_path]

4. Relaunch ChimeraX after installed.
5. The *ChimeraX pyKVFinder Tools* will install the dependencies. Relaunch ChimeraX if asked.
6. The *ChimeraX pyKVFinder Tools* plugin will now be ready for use.