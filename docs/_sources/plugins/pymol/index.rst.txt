**********************
PyMOL pyKVFinder Tools
**********************

pyKVFinder package is integrated into PyMOL through a plugin, called PyMOL pyKVFinder Tools. This plugin allows the user to run pyKVFinder on PyMOL and visualize the results directly on the PyMOL viewer.

Installation
============

First, it is required to install Anaconda or Miniconda (Python 3.10 or later). If you do not have it installed, please refer to `Anaconda <https://www.anaconda.com/products/distribution>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ websites.

In the `Anaconda Prompt`, create a new environment and install PyMOL:

.. code-block:: bash

    conda create -n pymol python=3.10
    conda activate pymol
    # PyMOL v3.0
    conda install -c conda-forge -c schrodinger pymol-bundle
    # PyMOL v2.6
    conda install -c conda-forge -c schrodinger pymol-bundle=2.6 
    # PyMOL Open-Source
    conda install -c conda-forge pymol-open-source


After installing one of the PyMOL versions (v3.0, v2.6 or Open-Source), you need to install the `pyKVFinder <https://pypi.org/project/pyKVFinder/>`_ package, if you have not done it yet:

.. code-block:: bash

    pip install pyKVFinder

Then, to install the `PyMOL pyKVFinder Tools` on PyMOL, download the latest `PyMOL-pyKVFinder-Tools.zip` from `here <https://github.com/LBC-LNBio/pyKVFinder/releases/latest/download/PyMOL-pyKVFinder-Tools.zip>`_ and follow these steps:

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
10. `PyMOL pyKVFinder Tools` is ready to use.

Or, instead of selecting `PyMOL-pyKVFinder-Tools.zip` (Step 5), you can select `__init__.py` file on `pyKVFinder/plugins/PyMOL-pyKVFinder-Tools` directory tree.

.. admonition:: Summary

    1. Download and install `Anaconda <https://www.anaconda.com/products/distribution>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_.

    2. Open the `Anaconda Prompt` and execute the following commands:

    .. code-block:: bash

        conda create -n pymol python=3.10
        conda activate pymol
        conda install -c conda-forge -c schrodinger pymol-bundle
        pip install pyKVFinder
        pymol

    3. Download the latest `PyMOL-pyKVFinder-Tools.zip <https://github.com/LBC-LNBio/pyKVFinder/releases/latest/download/PyMOL-pyKVFinder-Tools.zip>`_ and install it on PyMOL: `Plugin` > `Plugin Manager` > `Install New Plugin` > `Choose file...` > `PyMOL-pyKVFinder-Tools.zip` > `OK` > `OK`.


Tutorial
========

This section offers step-by-step tutorials on how to use the PyMOL pyKVFinder Tools to detect and analyze cavities in biomolecular and supramolecular cage structures. The tutorials are organized into two categories:

.. toctree::
    :maxdepth: 1

    Biomolecular cavity detection <tutorial/biomolecular-cavity-detection>
    Supramolecular cage cavity detection <tutorial/cage-cavity-detection>
