import argparse
import os
from typing import Any

from . import __name__, __version__

__all__ = ["argparser"]


def _check_pdb_xyz(x: str) -> str:
    """Checks if x is a .pdb or .xyz extension and returns its abspath.

    Parameters
    ----------
    x : Any
        Option to be checked.

    Returns
    -------
    x : str
        Valid option.

    Raises
    ------
    argparse.ArgumentTypeError
        `x` must have .pdb or .xyz extension.
    """
    if x.endswith(".pdb"):
        x = os.path.abspath(x)
    elif x.endswith(".xyz"):
        x = os.path.abspath(x)
    else:
        raise (argparse.ArgumentTypeError("%r must have .pdb or .xyz extension." % x))

    return x


def _check_hydropathy_options(x: Any) -> str:
    """Checks if x is a acceptable option for hydropathy argument.

    Parameters
    ----------
    x : Any
        Option to be checked.

    Returns
    -------
    x : str
        Valid option.

    Raises
    ------
    argparse.ArgumentTypeError
        Invalid TOML-formatted file: {x} (define a TOML-formatted
        file with hydrophobicity values for standard residues).
    argparse.ArgumentTypeError
        Invalid choice: '{x}' (choose from 'EisenbergWeiss', 'HessaHeijne',
        'KyteDoolittle', 'MoonFleming', 'RadzickaWolfenden', 'WimleyWhite',
        'ZhaoLondon' or define a TOML-formatted file with hydrophobicity
        values for standard residues).
    """
    if x in [
        "EisenbergWeiss",
        "HessaHeijne",
        "KyteDoolittle",
        "RadzickaWolfenden",
        "MoonFleming",
        "WimleyWhite",
        "ZhaoLondon",
    ]:
        return x
    elif os.path.exists(x):
        import toml

        try:
            toml.load(x)
        except toml.decoder.TomlDecodeError:
            msg = f"Invalid TOML-formatted file: {x} (define a TOML-formatted \
file with hydrophobicity values for standard residues)."
            raise (argparse.ArgumentTypeError(msg))
        return x
    else:
        msg = f"Invalid choice: '{x}' (choose from 'EisenbergWeiss', 'HessaHeijne', \
'KyteDoolittle', 'MoonFleming', 'RadzickaWolfenden', 'WimleyWhite', 'ZhaoLondon' or define a \
TOML-formatted file with hydrophobicity values for standard residues)."
        raise (argparse.ArgumentTypeError(msg))


def _positive_float(x: Any) -> float:
    """Checks if x is a float.

    Parameters
    ----------
    x : Any
        Variable to be checked.

    Returns
    -------
    x : float
        A valid float.

    Raises
    ------
    argparse.ArgumentTypeError
        `x` is not a floating-point literal.
    argparse.ArgumentTypeError
        `x` is not a positive floating-point.
    """
    try:
        x = float(x)
    except ValueError:
        msg = "%r is not a floating-point literal." % (x,)
        raise (argparse.ArgumentTypeError(msg))

    if x < 0.0:
        msg = "%r is not a positive floating-point." % (x,)
        raise (argparse.ArgumentTypeError(msg))
    return x


def _restricted_step_size(x: Any) -> float:
    """Checks if x is in range (0.0, 2.0].

    Parameters
    ----------
    x : Any
        Grid spacing (A).

    Returns
    -------
    x : float
        Grid spacing (A) between (0.0, 2.0].

    Raises
    ------
    argparse.ArgumentTypeError
        `x` is not a floating-point literal.
    argparse.ArgumentTypeError
        `x` is not in range (0.0, 2.0].
    """
    try:
        x = float(x)
    except ValueError:
        msg = "%r is not a floating-point literal." % (x,)
        raise (argparse.ArgumentTypeError(msg))

    if not 0.0 < x <= 2.0:
        msg = "%r is not in range (0.0, 2.0]." % (x,)
        raise (argparse.ArgumentTypeError(msg))
    return x


def argparser() -> argparse.ArgumentParser:
    """Defines pyKVFinder argument parser.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        pyKVFinder argument parser
    """

    # Overrides method in HelpFormatter
    class CapitalisedHelpFormatter(argparse.HelpFormatter):
        def add_usage(self, usage, actions, groups, prefix=None):
            if prefix is None:
                prefix = "Usage: "
            return super(CapitalisedHelpFormatter, self).add_usage(
                usage, actions, groups, prefix
            )

    parser = argparse.ArgumentParser(
        prog="pyKVFinder",
        description="pyKVFinder software detects and characterizes cavities.",
        formatter_class=CapitalisedHelpFormatter,
        add_help=True,
    )

    # Change parser titles
    parser._positionals.title = "Positional arguments"
    parser._optionals.title = "Optional arguments"

    # Positional arguments
    parser.add_argument(
        "input",
        metavar="(<.pdb> | <.xyz>)",
        type=_check_pdb_xyz,
        help="Path to a target PDB or XYZ file.",
    )

    # Optional arguments
    parser.add_argument(
        "-v",
        "--verbose",
        help="Print extra information to standard output.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"{__name__} v{__version__}",
        help="Show pyKVFinder version number and exit.",
    )
    parser.add_argument(
        "-b",
        "--base_name",
        metavar="<str>",
        type=str,
        help="Prefix to output files.",
    )
    parser.add_argument(
        "-O",
        "--output_directory",
        metavar="<str>",
        type=os.path.abspath,
        default=os.path.abspath("."),
        help="Output directory. (default: %(default)s)",
    )
    parser.add_argument(
        "-m",
        "--model",
        metavar="<int>",
        type=int,
        default=None,
        help="Define a specific model for a multi-model PDB file. \
(default: %(default)s (all models))",
    )
    parser.add_argument(
        "--nthreads",
        metavar="<int>",
        type=int,
        default=os.cpu_count() - 1,
        help="Number of threads to apply in parallel routines \
(default: %(default)s)",
    )

    # Create argument group
    parameters = parser.add_argument_group("Parameters")

    parameters.add_argument(
        "-d",
        "--dictionary",
        default=os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "data/vdw.dat"
        ),
        metavar="<str>",
        type=str,
        help="Path to van der Waals radii file. (default: %(default)s)",
    )
    parameters.add_argument(
        "-s",
        "--step",
        metavar="<float>",
        default=0.6,
        type=_restricted_step_size,
        help="Step size (grid spacing). (default: %(default).1lf)",
    )
    parameters.add_argument(
        "-i",
        "--probe_in",
        metavar="<float>",
        default=1.4,
        type=_positive_float,
        help="Probe In size (\u212B). (default: %(default).1lf)",
    )
    parameters.add_argument(
        "-o",
        "--probe_out",
        metavar="<float>",
        default=4.0,
        type=_positive_float,
        help="Probe Out size (\u212B). (default: %(default).1lf)",
    )
    parameters.add_argument(
        "-V",
        "--volume_cutoff",
        metavar="<float>",
        default=5.0,
        type=_positive_float,
        help="Cavities volume filter (\u212B\u00b3). (default: %(default).1f)",
    )
    parameters.add_argument(
        "-R",
        "--removal_distance",
        metavar="<float>",
        default=2.4,
        type=_positive_float,
        help="Length to be removed from the cavity-bulk frontier (\u212B). \
(default: %(default).1lf)",
    )
    parameters.add_argument(
        "-S",
        "--surface",
        metavar="<str>",
        type=str,
        default="SES",
        choices=["SES", "SAS"],
        help="Surface representation. Options: %(choices)s. \
SES specifies solvent excluded surface. SAS specifies solvent accessible \
surface. (default: %(default)s)",
    )
    parameters.add_argument(
        "--ignore_backbone",
        action="store_true",
        default=False,
        help="Ignore backbone contacts to cavity when defining interface \
residues. (default: %(default)s)",
    )

    # Extra modes
    extra_modes = parser.add_argument_group("Additional characterization")
    extra_modes.add_argument(
        "-D",
        "--depth",
        action="store_true",
        default=False,
        help="Depth characterization of the detected cavities. \
Map depth in B-factor in the cavity PDB file and maximum and average depth of \
the detected cavities. (default: %(default)s)",
    )
    extra_modes.add_argument(
        "--plot_frequencies",
        action="store_true",
        default=False,
        help="Plot bar charts of calculated frequencies of the detected \
cavities in a PDF file. (default: %(default)s)",
    )
    extra_modes.add_argument(
        "--hydropathy",
        type=_check_hydropathy_options,
        metavar="{EisenbergWeiss, HessaHeijne, KyteDoolittle, MoonFleming, \
RadzickaWolfenden, WimleyWhite, ZhaoLondon, <.toml>}",
        nargs="?",
        const="EisenbergWeiss",
        default=False,
        help="Hydropathy characterization of the detected cavities. Map \
hydrophobicity scale values as B-factor at surface points of detected \
cavities. Hydrophobicity scales options: %(metavar)s. A custom hydrophobicity \
scale can be defined via a TOML-formatted file. (default: %(default)s) \
(constant: %(const)s)",
    )

    # Create argument group
    box_adjustment = parser.add_argument_group("Box adjustment arguments")

    # Box adjustment arguments
    box_adjustment.add_argument(
        "-B",
        "--box",
        metavar="<.toml>",
        type=os.path.abspath,
        help='A path to TOML-formatted file with box parameters. Adjust the \
3D grid based on a list of residues (["resnum", "chain", "resname"]) and a \
padding or a set of four vertices (p1: origin, p2: X-axis max, p3: Y-axis \
max, p4: Z-axis max) with xyz coordinates ([x, y, z]). Define the custom 3D \
grid with only one of this box parameters.',
    )

    # Create argument group
    ligand_adjustment = parser.add_argument_group("Ligand adjustment arguments")

    # Ligand adjustment arguments
    ligand_adjustment.add_argument(
        "-L",
        "--ligand",
        metavar="(<.pdb> | <.xyz>)",
        type=_check_pdb_xyz,
        help="Path to a ligand PDB or XYZ file to limit the cavities within a \
radius (ligand_cutoff) around it.",
    )
    ligand_adjustment.add_argument(
        "--ligand_cutoff",
        metavar="<float>",
        default=5.0,
        type=_positive_float,
        help="Radius value to limit a space around the defined ligand. \
(default: %(default).1lf)",
    )

    return parser
