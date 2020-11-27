import os
import argparse
from . import __name__, __version__

__all__ = ["argparser"]


def _positive_float(x: float) -> float:
    """
    Checks if x is a float

    Parameters
    ----------
        x (float): variable to be checked

    Returns
    -------
        x (float): float variable
    """
    try:
        x = float(x)
    except ValueError:
        msg = "%r not a floating-point literal" % (x,)
        raise argparse.ArgumentTypeError(msg)

    if x < 0.0:
        msg = "%r not a positive floating-point" % (x,)
        raise argparse.ArgumentTypeError(msg)
    return x


def _restricted_step_size(x: float) -> float:
    """
    Checks: 0.0 < x <= 2.0

    Parameters
    ----------
        x (float): grid spacing (A)

    Returns
    -------
        x (float): grid spacing (A) between 0.0 and 2.0
    """
    try:
        x = float(x)
    except ValueError:
        msg = "%r not a floating-point literal" % (x,)
        raise argparse.ArgumentTypeError(msg)

    if not 0.0 < x <= 2.0:
        msg = "%r not in range (0.0, 2.0]" % (x,)
        raise argparse.ArgumentTypeError(msg)
    return x


def argparser() -> argparse.ArgumentParser:
    """
    Defines pyKVFinder argument parser

    Parameters
    ----------
        None

    Returns
    -------
        parser (argparse.ArgumentParser): pyKVFinder argument parser
    """
    # Overrides method in HelpFormatter
    class CapitalisedHelpFormatter(argparse.HelpFormatter):

        def add_usage(self, usage, actions, groups, prefix=None):
            if prefix is None:
                prefix = 'Usage: '
            return super(CapitalisedHelpFormatter, self).add_usage(
                usage, actions, groups, prefix)

    parser = argparse.ArgumentParser(
        prog="pyKVFinder",
        description="pyKVFinder software detects and characterizes cavities.",
        formatter_class=CapitalisedHelpFormatter,
        add_help=True)

    # Change parser titles
    parser._positionals.title = 'Positional arguments'
    parser._optionals.title = 'Optional arguments'

    # Positional arguments
    parser.add_argument(
        "pdb",
        metavar="<.pdb>",
        type=os.path.abspath,
        help="Path to a target PDB file."
    )

    # Optional arguments
    parser.add_argument("-v", "--verbose",
                        help="Print extra information to standard output.",
                        action="store_true")
    parser.add_argument('--version',
                        action='version',
                        version=f'{__name__} v{__version__}',
                        help="Show pyKVFinder version number and exit.")
    parser.add_argument("-b",
                        "--base_name",
                        metavar="<str>",
                        type=str,
                        help="Prefix to output files.")
    parser.add_argument("-O",
                        "--output_directory",
                        metavar="<path>",
                        type=os.path.abspath,
                        default=os.path.abspath("."),
                        help="Output directory. (default: %(default)s)")
    parser.add_argument("--nthreads",
                        metavar="<int>",
                        type=int,
                        default=os.cpu_count() - 1,
                        help="Number of threads to apply in parallel routines (default: %(default)s)")

    # Create argument group
    parameters = parser.add_argument_group("Parameters")

    parameters.add_argument("-d",
                            "--dictionary",
                            default=os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/vdw.dat"),
                            metavar="<file>",
                            type=str,
                            help="Path to van der Waals radii file. (default: %(default)s)")
    parameters.add_argument("-s",
                            "--step",
                            metavar="<float>",
                            default=0.6,
                            type=_restricted_step_size,
                            help="Step size (grid spacing). (default: %(default).1lf)")
    parameters.add_argument("-i",
                            "--probe_in",
                            metavar="<float>",
                            default=1.4,
                            type=_positive_float,
                            help=u"Probe In size (\u212B). (default: %(default).1lf)")
    parameters.add_argument("-o",
                            "--probe_out",
                            metavar="<float>",
                            default=4.0,
                            type=_positive_float,
                            help=u"Probe Out size (\u212B). (default: %(default).1lf)")
    parameters.add_argument("-V",
                            "--volume_cutoff",
                            metavar="<float>",
                            default=5.0,
                            type=_positive_float,
                            help=u"Cavities volume filter (\u212B\u00b3). (default: %(default).1lf)")
    parameters.add_argument("-R",
                            "--removal_distance",
                            metavar="<float>",
                            default=2.4,
                            type=_positive_float,
                            help=u"Length to be removed from the cavity-bulk frontier (\u212B). (default: %(default).1lf)")
    parameters.add_argument("-S",
                            "--surface",
                            metavar="<enum>",
                            type=str,
                            default="SES",
                            choices=["SES", "SAS"],
                            help="Surface representation. Options: %(choices)s. SAS specifies solvent accessible surface. SES specifies solvent excluded surface. (default: %(default)s)")
    parameters.add_argument("--ignore_backbone",
                            action='store_true',
                            default=False,
                            help="Ignore backbone contacts to cavity when defining interface residues. (default: %(default)s)")

    # Create argument group
    box_adjustment = parser.add_argument_group("Box adjustment arguments")

    # Box adjustment arguments
    box_adjustment.add_argument("-B",
                                "--box",
                                metavar="<.toml>",
                                type=os.path.abspath,
                                help="Adjust the 3D grid based on a list of residues ([\"resnum\",\"chain\"]) and a padding (default: 3.5) or a set of four points (p1: origin, p2: x-axis max, p3: y-axis max, p4: z-axis max) with xyz coordinates ([x, y, z]). Define one of them using a toml file.")

    # Create argument group
    ligand_adjustment = parser.add_argument_group("Ligand adjustment arguments")

    # Ligand adjustment arguments
    ligand_adjustment.add_argument("-L",
                                   "--ligand",
                                   metavar="<.pdb>",
                                   type=os.path.abspath,
                                   help="Path to a ligand PDB file to limit the cavities within a radius (ligand_cutoff) around it.")
    ligand_adjustment.add_argument("--ligand_cutoff",
                                   metavar="<float>",
                                   default=5.0,
                                   type=_positive_float,
                                   help="Radius value to limit a space around the defined ligand. (default: %(default).1lf)")

    return parser
