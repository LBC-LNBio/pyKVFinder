Requirements
============

If you don’t have Python v3 and/or SWIG yet, the installation procedure differs depending on the operating system.

Python and SWIG can be installed with conda, or with a package manager on Linux and macOS.

Package managers
----------------

On Linux:

.. code-block:: bash
    
  sudo apt install python3 swig

On macOS:

.. code-block:: bash
    
  brew install python3 swig

.. note:: 

  Users can use their preferred package manager to install SWIG and Python v3.

Conda
-----

If you use conda, you can install Python v3 and SWIG from the defaults channel:

.. code-block:: bash
    
  # Use an environment rather than install in base environement (recommended)
  conda create -n myenv python=3
  conda activate myenv
  # The SWIG install command
  conda install swig

Installation
============

The prerequisites for installing pyKVFinder is Python v3 and SWIG. If you don’t have Python and/or SWIG yet, please refer to this `section <index.html#requirements>`_.

To install the latest release on `PyPI <https://pypi.org/project/pyKVFinder>`_, 
run:

.. code-block:: bash

  pip install pyKVFinder

Or to install the latest developmental version, run:

.. code-block:: bash

  git clone https://github.com/LBC-LNBio/pyKVFinder.git
  pip install pyKVFinder
