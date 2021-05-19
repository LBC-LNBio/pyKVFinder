"""
Python-C parallel KVFinder

pyKVFinder detects and characterizes cavities in biomolecular structures. 

The characterization includes shape, volume, area, depth, hydropathy and interface residues and their frequencies. 

In addition to the set of function that can be imported on Python scritps, it contains a command line interface (CLI).

Python package
--------------
import pyKVFinder

Command Line Interface
----------------------
Usage: pyKVFinder [-h] [-v] [--version] [-b <str>] [-O <path>] [--nthreads <int>] [-d <file>] [-s <float>] [-i <float>] [-o <float>] [-V <float>] [-R <float>] [-S <enum>]
                [--ignore_backbone] [-B <.toml>] [-L <.pdb>] [--ligand_cutoff <float>]
                <.pdb>

GitHub repository: https://github.com/LBC-LNBio/pyKVFinder
Documentation: https://lbc-lnbio.github.io/pyKVFinder
"""

__name__ = "pyKVFinder"
__version__ = "0.1.0"
license = "GNU GPL-3.0 License"

try:
    from .utils import *
    from .grid import *
    from .main import *
except SyntaxError:
    pass
except ModuleNotFoundError:
    pass
