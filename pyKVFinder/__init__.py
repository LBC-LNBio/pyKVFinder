"""
Python-C parallel KVFinder

pyKVFinder detects and characterizes cavities in biomolecular structures. The characterization includes shape, volume, area and interface residues. In addition to the set of function that can be imported on Python scritps, it contains a command line interface (CLI).

Command Line Interface
----------------------
Usage: pyKVFinder [-h] [-v] [--version] [-b <str>] [-O <path>] [--nthreads <int>] [-d <file>] [-s <float>] [-i <float>] [-o <float>] [-V <float>] [-R <float>] [-S <enum>]
                [--ignore_backbone] [-B <.toml>] [-L <.pdb>] [--ligand_cutoff <float>]
                <.pdb>

For more information: https://github.com/LBC-LNBio/pyKVFinder
"""

__name__ = "pyKVFinder"
__version__ = "0.1"
license = "GNU GPL-3.0 License"

try:
    from .utils import *
    from .grid import *
    from .main import *
except SyntaxError:
    pass
