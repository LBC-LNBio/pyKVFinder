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

__version__ = "0.9.0"
__license__ = "GPL-3.0-or-later"

from .utils import (
    _write_parameters as _write_parameters,
    calculate_frequencies as calculate_frequencies,
    plot_frequencies as plot_frequencies,
    read_cavity as read_cavity,
    read_pdb as read_pdb,
    read_vdw as read_vdw,
    read_xyz as read_xyz,
    write_results as write_results,
)

from .grid import (
    _get_dimensions as _get_dimensions,
    _get_sincos as _get_sincos,
    constitutional as constitutional,
    depth as depth,
    detect as detect,
    export as export,
    export_openings as export_openings,
    get_vertices as get_vertices,
    get_vertices_from_file as get_vertices_from_file,
    hydropathy as hydropathy,
    openings as openings,
    spatial as spatial,
)

from .main import (
    run_workflow as run_workflow,
    pyKVFinderResults as pyKVFinderResults,
    Molecule as Molecule,
)

__all__ = [
    # utils
    "_write_parameters",
    "calculate_frequencies",
    "plot_frequencies",
    "read_cavity",
    "read_pdb",
    "read_vdw",
    "read_xyz",
    "write_results",
    # grid
    "_get_dimensions",
    "_get_sincos",
    "constitutional",
    "depth",
    "detect",
    "export",
    "export_openings",
    "get_vertices",
    "get_vertices_from_file",
    "hydropathy",
    "openings",
    "spatial",
    # main
    "run_workflow",
    "pyKVFinderResults",
    "Molecule",
]
