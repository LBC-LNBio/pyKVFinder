"""Python-C parallel KVFinder.

pyKVFinder detects and characterizes cavities in biomolecular structures.

The characterization includes shape, volume, area, depth, hydropathy and
interface residues and their frequencies.

In addition to the set of functions that can be imported into Python scripts,
it contains a command line interface (CLI).

Python package
--------------
>>> import pyKVFinder

Command Line Interface
----------------------
  Usage: pyKVFinder [-h] [-v] [--version] [-b <str>] [-O <str>]\
                    [-m <int>] [--nthreads <int>] [-d <str>] [-s <float>]\
                    [-i <float>] [-o <float>] [-V <float>] [-R <float>]\
                    [-S <str>] [--ignore_backbone] [-D] [--plot_frequencies]\
                    [--hydropathy [{EisenbergWeiss, HessaHeijne, KyteDoolittle,\
                    MoonFleming, RadzickaWolfenden, WimleyWhite, ZhaoLondon, <.toml>}]]\
                    [-B <.toml>] [-L (<.pdb> | <.xyz>)] [--ligand_cutoff <float>]\
                    (<.pdb> | <.xyz>)

See also
--------
* GitHub repository: https://github.com/LBC-LNBio/pyKVFinder

* Documentation: https://lbc-lnbio.github.io/pyKVFinder
"""

__name__ = "pyKVFinder"
__version__ = "0.6.15"
license = "GNU GPL-3.0 License"

from .utils import *
from .grid import *
from .main import *
