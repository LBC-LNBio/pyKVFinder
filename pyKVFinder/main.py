import os
import time
import logging
import numpy
from datetime import datetime
from .argparser import argparser
from .utils import (
    read_vdw,
    read_pdb,
    calculate_frequencies,
    plot_frequencies,
    write_results,
    _write_parameters,
)
from .grid import (
    get_vertices,
    get_grid_from_file,
    get_dimensions,
    get_sincos,
    detect,
    spatial,
    depth,
    constitutional,
    hydropathy,
    export,
)

__all__ = ["pyKVFinder", "pyKVFinderResults"]

here = os.path.abspath(os.path.dirname(__file__))
_dictionary = os.path.join(here, "data/vdw.dat")


def cli():
    """
    pyKVFinder Command Line Interface (CLI)

    Parameters
    ----------
        None

    Returns
    -------
        None

    Example
    -------
    Usage: pyKVFinder [-h] [-v] [--version] [-b <str>] [-O <str>] [--nthreads <int>] [-d <str>] [-s <float>] [-i <float>] [-o <float>] [-V <float>] [-R <float>] [-S <str>] [--ignore_backbone]
                    [-D] [--plot_frequencies] [-B <.toml>] [-L <.pdb>] [--ligand_cutoff <float>]
                    <.pdb>
    """
    # Start time
    start_time = time.time()

    # Load pyKVFinder argument parser
    parser = argparser()

    # Parse command-line arguments
    args = parser.parse_args()

    # Get base name from pdb file if not defined by user
    if not args.base_name:
        args.base_name = os.path.basename(args.pdb.replace(".pdb", ""))

    # Create output directory
    os.makedirs(args.output_directory, exist_ok=True)

    # Print message to stdout
    print(f"[PID {os.getpid()}] Running pyKVFinder for: {args.pdb}")

    # Start logging
    logging.basicConfig(
        filename=f"{os.path.join(args.output_directory, 'KVFinder.log')}",
        level=logging.INFO,
        format="%(message)s",
    )
    logging.info("=" * 80)
    logging.info(
        f"Date: {datetime.now().strftime('%a %d %B, %Y')}\nTime: {datetime.now().strftime('%H:%M:%S')}\n"
    )
    logging.info(f"[ Running pyKVFinder for: {args.pdb} ]")
    logging.info(f"> vdW radii file: {args.dictionary}")

    if args.verbose:
        print("> Loading atomic dictionary file")
    vdw = read_vdw(args.dictionary)

    if args.verbose:
        print("> Reading PDB coordinates")
    atominfo, xyzr = read_pdb(args.pdb, vdw)

    if args.ligand:
        if args.verbose:
            print("> Reading ligand coordinates")
        _, lxyzr = read_pdb(args.ligand, vdw)
    else:
        lxyzr = None

    if args.verbose:
        print("> Calculating 3D grid dimensions")
    if args.box:
        # Get vertices from file
        args.vertices, atominfo, xyzr, args.sincos, nx, ny, nz = get_grid_from_file(
            args.box,
            atominfo,
            xyzr,
            args.step,
            args.probe_in,
            args.probe_out,
            args.nthreads,
        )

        # Set flag to boolean
        args.box = True
    else:
        # Get vertices from pdb
        args.vertices = get_vertices(xyzr, args.probe_out, args.step)

        # Calculate distance between points
        nx, ny, nz = get_dimensions(args.vertices, args.step)
        if args.verbose:
            print(f"Dimensions: (nx:{nx}, ny:{ny}, nz:{nz})")

        # Calculate sin and cos of angles a and b
        args.sincos = get_sincos(args.vertices)
        if args.verbose:
            print(f"sina: {args.sincos[0]:.2f}\tsinb: {args.sincos[2]:.2f}")
            print(f"cosa: {args.sincos[1]:.2f}\tcosb: {args.sincos[3]:.2f}")

        # Set flag to boolean
        args.box = False

    # Logging parameters
    logging.info(f"> Step: {args.step} \u00c5")
    logging.info(f"> Probe In: {args.probe_in} \u00c5")
    logging.info(f"> Probe Out: {args.probe_out} \u00c5")
    logging.info(f"> Voxel volume: {args.step * args.step * args.step} \u00c5\u00b3")
    logging.info(f"> Dimensions: (nx:{nx}, ny:{ny}, nz:{nz})")
    logging.info(f"> sina: {args.sincos[0]:.2f}\tcosa: {args.sincos[1]:.2f}")
    logging.info(f"> sinb: {args.sincos[2]:.2f}\tcosb: {args.sincos[3]:.2f}")

    # Cavity detection
    ncav, cavities = detect(
        nx,
        ny,
        nz,
        xyzr,
        args.vertices,
        args.sincos,
        args.step,
        args.probe_in,
        args.probe_out,
        args.removal_distance,
        args.volume_cutoff,
        lxyzr,
        args.ligand_cutoff,
        args.box,
        args.surface,
        args.nthreads,
        args.verbose,
    )

    # Cavities were found
    if ncav > 0:
        # Spatial characterization
        surface, volume, area = spatial(
            cavities, ncav, args.step, args.nthreads, args.verbose
        )

        # Depth characterization
        if args.depth:
            depths, max_depth, avg_depth = depth(
                cavities, ncav, args.step, args.nthreads, args.verbose
            )
        else:
            depths, max_depth, avg_depth = None, None, None

        # Constitutional characterization
        residues = constitutional(
            cavities,
            atominfo,
            xyzr,
            args.vertices,
            args.sincos,
            ncav,
            args.step,
            args.probe_in,
            args.ignore_backbone,
            args.nthreads,
            args.verbose,
        )
        frequencies = calculate_frequencies(residues)

        # Plot histograms of frequencies
        if args.plot_frequencies:
            output_plot = os.path.join(
                args.output_directory, f"{args.base_name}.histograms.pdf"
            )
            plot_frequencies(frequencies, output_plot)

        # Hydropathy characterization
        if args.hydropathy:
            # Map hydrophobicity scales
            scales, avg_hydropathy = hydropathy(
                surface,
                atominfo,
                xyzr,
                args.vertices,
                args.sincos,
                ncav,
                args.step,
                args.probe_in,
                args.hydropathy,
                args.ignore_backbone,
                args.nthreads,
                args.verbose,
            )
            output_hydropathy = os.path.join(
                args.output_directory,
                f"{args.base_name}.{list(avg_hydropathy.keys())[-1]}.pdb",
            )
        else:
            scales, avg_hydropathy, output_hydropathy = None, None, None

        # Export cavities
        output_cavity = os.path.join(
            args.output_directory, f"{args.base_name}.KVFinder.output.pdb"
        )
        export(
            output_cavity,
            cavities,
            surface,
            args.vertices,
            args.sincos,
            ncav,
            args.step,
            depths,
            output_hydropathy,
            scales,
            args.nthreads,
        )

        # Write results
        output_results = os.path.join(
            args.output_directory, f"{args.base_name}.KVFinder.results.toml"
        )
        write_results(
            output_results,
            args.pdb,
            args.ligand,
            output_cavity,
            output_hydropathy,
            volume,
            area,
            max_depth,
            avg_depth,
            avg_hydropathy,
            residues,
            frequencies,
            args.step,
        )

        # Write parameters
        _write_parameters(args)
    else:
        print("> No cavities detected!")

    # Elapsed time
    elapsed_time = time.time() - start_time
    print(f"[ \033[1mElapsed time:\033[0m {elapsed_time:.4f} ]")
    logging.info(f"[ Elapsed time (s): {elapsed_time:.4f} ]\n")

    return True


class pyKVFinderResults(object):
    f"""
    A class containing pyKVFinder detection and characterization results

    Attributes
    ----------
        cavities (numpy.ndarray): cavity points in the 3D grid (cavities[nx][ny][nz])
        surface (numpy.ndarray): surface points in the 3D grid (surface[nx][ny][nz])
        depths (numpy.ndarray): a numpy array with depth of cavity points (depth[nx][ny][nz])
        scales (numpy.ndarray): a numpy array with hydrophobicity scale value mapped at surface points (scales[nx][ny][nz])
        volume (dict): a dictionary with volume of each detected cavity
        area (dict): a dictionary with area of each detected cavity
        max_depth (dict): a dictionary with maximum depth of each detected cavity
        avg_depth (dict): a dictionary with average depth of each detected cavity
        avg_hydropathy (dict): a dictionary with average hydropathy of each detected cavity and the range of the hydrophobicity scale (min, max)
        residues (dict): a dictionary of list of interface residues pairs of each detected cavity
        frequencies (dict): a dictionary with frequencies of interface residues and classes of residues of each detected cavity
        _vertices (numpy.ndarray): a numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)
        _step (float): grid spacing (A)
        _ncav (int): number of cavities
        _pdb (str): a path to input PDB file
        _ligand (str): a path to ligand PDB file

    Methods
    -------
        export(output = 'cavity.pdb', output_hydropathy = 'hydropathy.pdb', nthreads = {os.cpu_count() - 1}):
            Exports cavities and characterizations to PDB-formatted files
        write(fn = 'results.toml', output = None, output_hydropathy = None):
            Writes file paths and characterizations to a TOML-formatted results file
        plot_frequencies(pdf = 'histograms.pdf')
            Plot histograms of frequencies (residues and classes of residues) in a PDF file
        export_all(fn = 'results.toml', output = 'cavity.pdb', output_hydropathy = 'hydropathy.pdb', include_frequencies_pdf = False, pdf = 'histogtrams.pdf', nthreads = {os.cpu_count() - 1}):
            Exports cavities and characterizations to PDB-formatted files, writes file paths and characterizations to a TOML-formatted results file, and optionally plot histograms of frequencies in a PDF file
    """

    def __init__(
        self,
        cavities: numpy.ndarray,
        surface: numpy.ndarray,
        depths: numpy.ndarray,
        scales: numpy.ndarray,
        volume: dict,
        area: dict,
        max_depth: dict,
        avg_depth: dict,
        avg_hydropathy: dict,
        residues: dict,
        frequencies: dict,
        _vertices: numpy.ndarray,
        _step: float,
        _ncav: int,
        _pdb: str = None,
        _ligand: str = None,
    ):
        """
        Constructs attributes for pyKVFinderResults object

        Parameters
        ----------
            cavities (numpy.ndarray): cavity points in the 3D grid (cavities[nx][ny][nz])
            surface (numpy.ndarray): surface points in the 3D grid (surface[nx][ny][nz])
            depths (numpy.ndarray): a numpy array with depth of cavity points (depth[nx][ny][nz])
            scales (numpy.ndarray): a numpy array with hydrophobicity scale value mapped at surface points (scales[nx][ny][nz])
            volume (dict): a dictionary with volume of each detected cavity
            area (dict): a dictionary with area of each detected cavity
            max_depth (dict): a dictionary with maximum depth of each detected cavity
            avg_depth (dict): a dictionary with average depth of each detected cavity
            avg_hydropathy (dict): a dictionary with average hydropathy of each detected cavity and the range of the hydrophobicity scale (min, max)
            residues (dict): a dictionary of list of interface residues pairs of each detected cavity
            frequencies (dict): a dictionary with frequencies of interface residues and classes of residues of each detected cavity
            _vertices (numpy.ndarray): a numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)
            _step (float): grid spacing (A)
            _ncav (int): number of cavities
            _pdb (str): a path to input PDB file
            _ligand (str): a path to ligand PDB file
        """
        self.cavities = cavities
        self.surface = surface
        self.depths = depths
        self.scales = scales
        self.volume = volume
        self.area = area
        self.max_depth = max_depth
        self.avg_depth = avg_depth
        self.avg_hydropathy = avg_hydropathy
        self.residues = residues
        self.frequencies = frequencies
        self._vertices = _vertices
        self._step = _step
        self._ncav = _ncav
        self._pdb = os.path.abspath(_pdb)
        self._ligand = os.path.abspath(_ligand) if _ligand else None

    def __repr__(self):
        return "<pyKVFinderResults object>"

    def export(
        self,
        output: str = "cavity.pdb",
        output_hydropathy: str = "hydropathy.pdb",
        nthreads: int = os.cpu_count() - 1,
    ) -> None:
        """
        Exports cavitiy (H) and surface (HA) points to PDB-formatted file with a variable (B; optional) in B-factor column, and hydropathy to PDB-formatted file in B-factor column at surface points (HA)

        Parameters
        ----------
            output (str): a path to PDB file for writing cavities
            output_hydropathy (str): a path to PDB file for writing hydropathy at surface points
            nthreads (int): number of threads

        Returns
        -------
            None

        Note
        ----
            The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.
        """
        sincos = get_sincos(self._vertices)
        export(
            output,
            self.cavities,
            self.surface,
            self._vertices,
            sincos,
            self._ncav,
            self._step,
            self.depths,
            output_hydropathy,
            self.scales,
            nthreads,
        )

    def write(
        self,
        fn: str = "results.toml",
        output: str = None,
        output_hydropathy: str = None,
    ) -> None:
        """
        Writes file paths and cavity characterization to TOML-formatted file

         Parameters
         ----------
             fn (str): a path to TOML-formatted file for writing file paths and cavity characterization (volume, area, depth and interface residues) per cavity detected
             output (str): a path to a cavity PDB file
             output_hydropathy (str): a path to PDB file for writing hydropathy at surface points

         Returns
         -------
             None

         Note
         ----
             The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.
        """
        output = os.path.abspath(output) if output else None
        output_hydropathy = (
            os.path.abspath(output_hydropathy) if output_hydropathy else None
        )
        write_results(
            fn,
            self._pdb,
            self._ligand,
            output,
            output_hydropathy,
            self.volume,
            self.area,
            self.max_depth,
            self.avg_depth,
            self.avg_hydropathy,
            self.residues,
            self.frequencies,
            self._step,
        )

    def plot_frequencies(self, pdf: str = "histograms.pdf"):
        """
        Plot histograms of frequencies (residues and classes of residues) in a PDF file

        Parameters
        ----------
            pdf (str): a path to a PDF file

        Returns
        -------
            None

        Note
        ----
            The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

        Classes
        -------
            Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine
            Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine
            Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine
            Negatively charged (R4): Aspartate, Glutamate
            Positively charged (R5): Arginine, Histidine, Lysine
            Non-standard (RX): Non-standard residues
        """
        plot_frequencies(self.frequencies, pdf)

    def export_all(
        self,
        fn: str = "results.toml",
        output: str = "cavity.pdb",
        output_hydropathy: str = "hydropathy.pdb",
        include_frequencies_pdf: bool = False,
        pdf: str = "histograms.pdf",
        nthreads: int = os.cpu_count() - 1,
    ) -> None:
        """
        Exports cavities and characterization to PDB-formatted files, writes file paths and characterization to a TOML-formatted file, and optionally plot histograms of frequencies (residues and classes of residues) in a PDF file

        Parameters
        ----------
            fn (str): a path to TOML-formatted file for writing file paths and cavity characterization (volume, area and interface residues) per cavity detected
            output (str): a path to PDB file for writing cavities
            output_hydropathy (str): a path to PDB file for writing hydropathy at surface points
            include_frequencies_pdf (bool): whether to plot frequencies (residues and classes of residues) to PDF file
            pdf (str): a path to a PDF file
            nthreads (int): number of threads

        Returns
        -------
            None

        Note
        ----
            The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

        Classes
        -------
            Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine
            Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine
            Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine
            Negatively charged (R4): Aspartate, Glutamate
            Positively charged (R5): Arginine, Histidine, Lysine
            Non-standard (RX): Non-standard residues
        """
        # Export cavity PDB file
        self.export(output, output_hydropathy, nthreads)
        # Write KVFinder results TOML
        self.write(fn, output, output_hydropathy)
        # Plot histograms of frequencies
        if include_frequencies_pdf:
            self.plot_frequencies(pdf)


def pyKVFinder(
    pdb: str,
    ligand: str = None,
    dictionary: str = _dictionary,
    box: str = None,
    step: float = 0.6,
    probe_in: float = 1.4,
    probe_out: float = 4.0,
    removal_distance: float = 2.4,
    volume_cutoff: float = 5.0,
    ligand_cutoff: float = 5.0,
    include_depth: bool = False,
    include_hydropathy: bool = False,
    hydrophobicity_scale: str = "EisenbergWeiss",
    surface: str = "SES",
    ignore_backbone: bool = False,
    nthreads: int = os.cpu_count() - 1,
    verbose: bool = False,
) -> pyKVFinderResults:
    """
    Detects and characterizes cavities (volume, area, depth [optional], hydropathy [optional] and interface residues)

    Parameters
    ----------
        pdb (str): a path to input PDB file
        ligand (str): a path to ligand PDB file
        dictionary (str): a path to van der Waals radii file
        box (str): a path to box configuration file (TOML-formatted)
        step (float): grid spacing (A)
        probe_in (float): Probe In size (A)
        probe_out (float): Probe Out size (A)
        removal_distance (float): length to be removed from the cavity-bulk frontier (A)
        volume_cutoff (float): cavities volume filter (A3)
        ligand_cutoff (float): radius value to limit a space around a ligand (A)
        include_depth (bool): whether to characterize the depth of the detected cavities
        include_hydropathy (bool): whether to characterize the hydropathy of the detected cavities
        hydrophobicity_scale (str): name of a native hydrophobicity scale (EisenbergWeiss, HessaHeijne, KyteDoolitte, MoonFleming, WimleyWhite, ZhaoLondon) or a path to a TOML-formatted file with a custom hydrophobicity scale.
        surface (str): keywords options are SES (Solvent Excluded Surface) or SAS (Solvent Accessible Surface)
        ignore_backbone (bool): whether to ignore backbone atoms (C, CA, N, O) when defining interface residues
        nthreads (int): number of threads
        verbose: print extra information to standard output

    Returns
    -------
        results (pyKVFinderResults): a class with the following attributes defined:
            * cavities (numpy.ndarray): cavity points in the 3D grid (cavities[nx][ny][nz])
            * surface (numpy.ndarray): surface points in the 3D grid (surface[nx][ny][nz])
            * depths (numpy.ndarray): a numpy array with depth of cavity points (depth[nx][ny][nz])
            * scales (numpy.ndarray): a numpy array with hydrophobicity scale value mapped at surface points (scales[nx][ny][nz])
            * volume (dict): a dictionary with volume of each detected cavity
            * area (dict): a dictionary with area of each detected cavity
            * max_depth (dict): a dictionary with maximum depth of each detected cavity
            * avg_depth (dict): a dictionary with average depth of each detected cavity
            * avg_hydropathy (dict): a dictionary with average hydropathy of each detected cavity and the range of the hydrophobicity scale (min, max)
            * residues (dict): a dictionary of list of interface residues pairs of each detected cavity
            * frequencies (dict): a dictionary with frequencies of interface residues and classes of residues of each detected cavity
            * _vertices (numpy.ndarray): a numpy array with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis)
            * _step (float): grid spacing (A)
            * _ncav (int): number of cavities
            * _pdb (str): a path to input PDB file
            * _ligand (str): a path to ligand PDB file

    Notes
    -----
        The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on

        Cavities array has integer labels in each position, that are:
            * -1: bulk points
            * 0: biomolecule points
            * 1: empty space points
            * >=2: cavity points

        Surface array has integer labels in each position, that are:
            * -1: bulk points
            * 0: biomolecule or empty space points
            * >=2: cavity points

    Classes
    -------
        Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine
        Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine
        Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine
        Negatively charged (R4): Aspartate, Glutamate
        Positively charged (R5): Arginine, Histidine, Lysine
        Non-standard (RX): Non-standard residues
    """
    if verbose:
        print("> Loading atomic dictionary file")
    vdw = read_vdw(dictionary)

    if verbose:
        print("> Reading PDB coordinates")
    atominfo, xyzr = read_pdb(pdb, vdw)

    if ligand:
        if verbose:
            print("> Reading ligand coordinates")
        _, lxyzr = read_pdb(ligand, vdw)
    else:
        lxyzr = None

    if verbose:
        print("> Calculating 3D grid dimensions")
    if box:
        # Get vertices from file
        vertices, atominfo, xyzr, sincos, nx, ny, nz = get_grid_from_file(
            box, atominfo, xyzr, step, probe_in, probe_out, nthreads
        )

        # Set flag to boolean
        box = True
    else:
        # Get vertices from pdb
        vertices = get_vertices(xyzr, probe_out, step)
        # Calculate distance between points
        nx, ny, nz = get_dimensions(vertices, step)
        if verbose:
            print(f"Dimensions: (nx:{nx}, ny:{ny}, nz:{nz})")

        # Calculate sin and cos of angles a and b
        sincos = get_sincos(vertices)
        if verbose:
            print(f"sina: {sincos[0]:.2f}\tsinb: {sincos[2]:.2f}")
            print(f"cosa: {sincos[1]:.2f}\tcosb: {sincos[3]:.2f}")

        # Set flag to boolean
        box = False

    # Cavity detection
    ncav, cavities = detect(
        nx,
        ny,
        nz,
        xyzr,
        vertices,
        sincos,
        step,
        probe_in,
        probe_out,
        removal_distance,
        volume_cutoff,
        lxyzr,
        ligand_cutoff,
        box,
        surface,
        nthreads,
        verbose,
    )

    if ncav > 0:
        # Spatial characterization
        surface, volume, area = spatial(cavities, ncav, step, nthreads, verbose)

        # Depth characterization
        if include_depth:
            depths, max_depth, avg_depth = depth(
                cavities, ncav, step, nthreads, verbose
            )
        else:
            depths, max_depth, avg_depth = None, None, None

        # Constitutional characterization
        residues = constitutional(
            cavities,
            atominfo,
            xyzr,
            vertices,
            sincos,
            ncav,
            step,
            probe_in,
            ignore_backbone,
            nthreads,
            verbose,
        )
        frequencies = calculate_frequencies(residues)

        # Hydropathy hydrophobicity scales
        if include_hydropathy:
            scales, avg_hydropathy = hydropathy(
                surface,
                atominfo,
                xyzr,
                vertices,
                sincos,
                ncav,
                step,
                probe_in,
                hydrophobicity_scale,
                ignore_backbone,
                nthreads,
                verbose,
            )
        else:
            scales, avg_hydropathy = None, None
    else:
        print("Warning: No cavities detected, returning None!")
        return None

    # Return dict
    results = pyKVFinderResults(
        cavities,
        surface,
        depths,
        scales,
        volume,
        area,
        max_depth,
        avg_depth,
        avg_hydropathy,
        residues,
        frequencies,
        vertices,
        step,
        ncav,
        pdb,
        ligand,
    )

    return results
