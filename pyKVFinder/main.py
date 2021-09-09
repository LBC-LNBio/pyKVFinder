import os
import time
import logging
import numpy
import pathlib
from datetime import datetime
from typing import Union, Optional, Dict, List
from .argparser import argparser
from .utils import (
    read_vdw,
    read_pdb,
    read_xyz,
    calculate_frequencies,
    plot_frequencies,
    write_results,
    _write_parameters,
)
from .grid import (
    get_vertices,
    get_vertices_from_file,
    _get_dimensions,
    _get_sincos,
    detect,
    spatial,
    depth,
    constitutional,
    hydropathy,
    export,
)

__all__ = ["run_workflow", "pyKVFinderResults"]

VDW = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/vdw.dat")


def cli() -> None:
    """pyKVFinder Command Line Interface (CLI).

    Parameters
    ----------
    None

    Returns
    -------
    None

    Example
    -------
    Usage: pyKVFinder [-h] [-v] [--version] [-b <str>] [-O <str>]
                      [--nthreads <int>] [-d <str>] [-s <float>] [-i <float>]
                      [-o <float>] [-V <float>] [-R <float>] [-S <str>]
                      [--ignore_backbone] [-D] [--plot_frequencies]
                      [-B <.toml>] [-L (<.pdb> | <.xyz>)] [--ligand_cutoff <float>]
                      (<.pdb> | <.xyz>)
    """
    # Start time
    start_time = time.time()

    # Load pyKVFinder argument parser
    parser = argparser()

    # Parse command-line arguments
    args = parser.parse_args()

    # Get base name from input file if not defined by user
    if args.base_name is None:
        args.base_name = os.path.basename(
            args.input.replace(".pdb", "").replace(".xyz", "")
        )

    # Create output directory
    os.makedirs(args.output_directory, exist_ok=True)

    # Print message to stdout
    print(f"[PID {os.getpid()}] Running pyKVFinder for: {args.input}")

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
    logging.info(f"[ Running pyKVFinder for: {args.input} ]")
    logging.info(f"> vdW radii file: {args.dictionary}")

    if args.verbose:
        print("> Loading atomic dictionary file")
    vdw = read_vdw(args.dictionary)

    if args.verbose:
        print("> Reading PDB coordinates")
    if args.input.endswith(".pdb"):
        atomic = read_pdb(args.input, vdw)
    elif args.input.endswith(".xyz"):
        atomic = read_xyz(args.input, vdw)

    if args.ligand:
        if args.verbose:
            print("> Reading ligand coordinates")
        if args.ligand.endswith('.pdb'):
            latomic = read_pdb(args.ligand, vdw)
        elif args.ligand.endswith('.xyz'):
            latomic = read_xyz(args.ligand, vdw)
    else:
        latomic = None

    if args.verbose:
        print("> Calculating 3D grid dimensions")
    if args.box:
        # Get vertices from file
        args.vertices, atomic = get_vertices_from_file(
            args.box,
            atomic,
            args.step,
            args.probe_in,
            args.probe_out,
            args.nthreads,
        )

        # Set flag to boolean
        args.box = True
    else:
        # Get vertices from input
        args.vertices = get_vertices(atomic, args.probe_out, args.step)

        # Set flag to boolean
        args.box = False

    # Calculate distance between points
    nx, ny, nz = _get_dimensions(args.vertices, args.step)

    # Calculate sin and cos of angles a and b
    args.sincos = _get_sincos(args.vertices)

    if args.verbose:
        print(f"p1: {args.vertices[0]}")
        print(f"p2: {args.vertices[1]}")
        print(f"p3: {args.vertices[2]}")
        print(f"p4: {args.vertices[3]}")
        print(f"Dimensions: (nx:{nx}, ny:{ny}, nz:{nz})")
        print(f"sina: {args.sincos[0]:.2f}\tsinb: {args.sincos[2]:.2f}")
        print(f"cosa: {args.sincos[1]:.2f}\tcosb: {args.sincos[3]:.2f}")

    # Logging parameters
    logging.info(f"> Step: {args.step} \u00c5")
    logging.info(f"> Probe In: {args.probe_in} \u00c5")
    logging.info(f"> Probe Out: {args.probe_out} \u00c5")
    logging.info(f"> Voxel volume: {args.step * args.step * args.step} \u00c5\u00b3")
    logging.info(f"> p1: {args.vertices[0]}")
    logging.info(f"> p2: {args.vertices[1]}")
    logging.info(f"> p3: {args.vertices[2]}")
    logging.info(f"> p4: {args.vertices[3]}")
    logging.info(f"> Dimensions: (nx:{nx}, ny:{ny}, nz:{nz})")
    logging.info(f"> sina: {args.sincos[0]:.2f}\tcosa: {args.sincos[1]:.2f}")
    logging.info(f"> sinb: {args.sincos[2]:.2f}\tcosb: {args.sincos[3]:.2f}")

    # Cavity detection
    ncav, cavities = detect(
        atomic,
        args.vertices,
        args.step,
        args.probe_in,
        args.probe_out,
        args.removal_distance,
        args.volume_cutoff,
        latomic,
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
            cavities, args.step, None, args.nthreads, args.verbose
        )

        # Constitutional characterization
        residues = constitutional(
            cavities,
            atomic,
            args.vertices,
            args.step,
            args.probe_in,
            args.ignore_backbone,
            None,
            args.nthreads,
            args.verbose,
        )
        frequencies = calculate_frequencies(residues)

        # Depth characterization
        if args.depth:
            depths, max_depth, avg_depth = depth(
                cavities, args.step, None, args.nthreads, args.verbose
            )
        else:
            depths, max_depth, avg_depth = None, None, None

        # Plot bar charts of frequencies
        if args.plot_frequencies:
            output_plot = os.path.join(
                args.output_directory, f"{args.base_name}.barplot.pdf"
            )
            plot_frequencies(frequencies, output_plot)

        # Hydropathy characterization
        if args.hydropathy:
            # Map hydrophobicity scales
            scales, avg_hydropathy = hydropathy(
                surface,
                atomic,
                args.vertices,
                args.step,
                args.probe_in,
                args.hydropathy,
                args.ignore_backbone,
                None,
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
            args.step,
            depths,
            output_hydropathy,
            scales,
            None,
            args.nthreads,
        )

        # Write results
        output_results = os.path.join(
            args.output_directory, f"{args.base_name}.KVFinder.results.toml"
        )
        write_results(
            output_results,
            args.input,
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

    return 0


class pyKVFinderResults(object):
    """A class containing pyKVFinder detection and characterization results.

    Parameters
    ----------
    cavities : numpy.ndarray
        Cavity points in the 3D grid (cavities[nx][ny][nz]).
        Cavities array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule points;

            * 1: empty space points;

            * >=2: cavity points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.
    surface : numpy.ndarray
        Surface points in the 3D grid (surface[nx][ny][nz]).
        Surface array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule or empty space points;

            * >=2: surface points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.
    depths : numpy.ndarray, optional
        A numpy.ndarray with depth of cavity points (depth[nx][ny][nz]).
    scales : numpy.ndarray, optional
        A numpy.ndarray with hydrophobicity scale value mapped at surface
        points (scales[nx][ny][nz]).
    volume : Dict[str, float]
        A dictionary with volume of each detected cavity.
    area : Dict[str, float]
        A dictionary with area of each detected cavity.
    max_depth : Dict[str, float], optional
        A dictionary with maximum depth of each detected cavity.
    avg_depth : Dict[str, float], optional
        A dictionary with average depth of each detected cavity.
    avg_hydropathy : Dict[str, float], optional
        A dictionary with average hydropathy for each detected cavity and
        the range of the hydrophobicity scale (min, max).
    residues: Dict[str, List[List[str]]]
        A dictionary with a list of interface residues for each detected
        cavity.
    frequencies : Dict[str, Dict[str, Dict[str, int]]], optional
        A dictionary with frequencies of residues and class for
        residues of each detected cavity.
    _vertices : numpy.ndarray
        A numpy.ndarray or a list with xyz vertices coordinates (origin,
        X-axis, Y-axis, Z-axis).
    _step : float
        Grid spacing (A).
    _input : Union[str, pathlib.Path], optional
        A path to input PDB or XYZ file, by default None.
    _ligand : Union[str, pathlib.Path], optional
        A path to ligand PDB or XYZ file, by default None.

    Attributes
    ----------
    cavities : numpy.ndarray
        Cavity points in the 3D grid (cavities[nx][ny][nz]).
        Cavities array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule points;

            * 1: empty space points;

            * >=2: cavity points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.
    surface : numpy.ndarray
        Surface points in the 3D grid (surface[nx][ny][nz]).
        Surface array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule or empty space points;

            * >=2: surface points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.
    depths : numpy.ndarray, optional
        A numpy.ndarray with depth of cavity points (depth[nx][ny][nz]).
    scales : numpy.ndarray, optional
        A numpy.ndarray with hydrophobicity scale value mapped at surface
        points (scales[nx][ny][nz]).
    ncav : int
        Number of cavities.
    volume : Dict[str, float]
        A dictionary with volume of each detected cavity.
    area : Dict[str, float]
        A dictionary with area of each detected cavity.
    max_depth : Dict[str, float], optional
        A dictionary with maximum depth of each detected cavity.
    avg_depth : Dict[str, float], optional
        A dictionary with average depth of each detected cavity.
    avg_hydropathy : Dict[str, float], optional
        A dictionary with average hydropathy for each detected cavity and
        the range of the hydrophobicity scale (min, max).
    residues: Dict[str, List[List[str]]]
        A dictionary with a list of interface residues for each detected
        cavity.
    frequencies : Dict[str, Dict[str, Dict[str, int]]], optional
        A dictionary with frequencies of residues and class for
        residues of each detected cavity.
    _vertices : numpy.ndarray
        A numpy.ndarray or a list with xyz vertices coordinates (origin,
        X-axis, Y-axis, Z-axis).
    _step : float
        Grid spacing (A).
    _input : Union[str, pathlib.Path], optional
        A path to input PDB or XYZ file, by default None.
    _ligand : Union[str, pathlib.Path], optional
        A path to ligand PDB or XYZ file, by default None.
    """

    def __init__(
        self,
        cavities: numpy.ndarray,
        surface: numpy.ndarray,
        depths: Optional[numpy.ndarray],
        scales: Optional[numpy.ndarray],
        volume: Dict[str, float],
        area: Dict[str, float],
        max_depth: Optional[Dict[str, float]],
        avg_depth: Optional[Dict[str, float]],
        avg_hydropathy: Optional[Dict[str, float]],
        residues: Dict[str, List[List[str]]],
        frequencies: Optional[Dict[str, Dict[str, Dict[str, int]]]],
        _vertices: numpy.ndarray,
        _step: Union[float, int],
        _input: Optional[Union[str, pathlib.Path]] = None,
        _ligand: Optional[Union[str, pathlib.Path]] = None,
    ):
        self.cavities = cavities
        self.surface = surface
        self.depths = depths
        self.scales = scales
        self.volume = volume
        self.ncav = cavities.max() - 1
        self.area = area
        self.max_depth = max_depth
        self.avg_depth = avg_depth
        self.avg_hydropathy = avg_hydropathy
        self.residues = residues
        self.frequencies = frequencies
        self._vertices = _vertices
        self._step = _step
        self._input = os.path.abspath(_input)
        self._ligand = os.path.abspath(_ligand) if _ligand else None

    def __repr__(self):
        return "<pyKVFinderResults object>"

    def export(
        self,
        output: Union[str, pathlib.Path] = "cavity.pdb",
        output_hydropathy: Union[str, pathlib.Path] = "hydropathy.pdb",
        nthreads: Optional[int] = None,
    ) -> None:
        """Exports cavitiy (H) and surface (HA) points to PDB-formatted file
        with a variable (B; optional) in B-factor column, and hydropathy to
        PDB-formatted file in B-factor column at surface points (HA).

        Parameters
        ----------
        output : Union[str, pathlib.Path]), optional
            A path to PDB file for writing cavities, by default `cavity.pdb`.
        output_hydropathy : Union[str, pathlib.Path], optional
            A path to PDB file for writing hydropathy at surface points, by
            default `hydropathy.pdb`.
        nthreads : int, optional
            Number of threads, by default None. If None, the number of threads is
            `os.cpu_count() - 1`.

        Note
        ----
        The cavity nomenclature is based on the integer label. The cavity
        marked with 2, the first integer corresponding to a cavity, is KAA,
        the cavity marked with 3 is KAB, the cavity marked with 4 is KAC
        and so on.
        """
        export(
            output,
            self.cavities,
            self.surface,
            self._vertices,
            self._step,
            self.depths,
            output_hydropathy,
            self.scales,
            None,
            nthreads,
        )

    def write(
        self,
        fn: Union[str, pathlib.Path] = "results.toml",
        output: Optional[Union[str, pathlib.Path]] = None,
        output_hydropathy: Optional[Union[str, pathlib.Path]] = None,
    ) -> None:
        """
        Writes file paths and cavity characterization to TOML-formatted file

        Parameters
        ----------
        fn : Union[str, pathlib.Path], optional
            A path to TOML-formatted file for writing file paths and cavity
            characterization (volume, area, depth and interface residues)
            per cavity detected, by default `results.toml`.
        output : Union[str, pathlib.Path], optional
            A path to a cavity PDB file, by default None.
        output_hydropathy : Union[str, pathlib.Path], optional
            A path to PDB file for writing hydropathy at surface points, by
            default None.

        Note
        ----
        The cavity nomenclature is based on the integer label. The cavity
        marked with 2, the first integer corresponding to a cavity, is KAA,
        the cavity marked with 3 is KAB, the cavity marked with 4 is KAC
        and so on.
        """
        write_results(
            fn,
            self._input,
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

    def plot_frequencies(self, pdf: Union[str, pathlib.Path] = "barplots.pdf"):
        """Plot bar charts of frequencies (residues and classes of residues) in
        a PDF file.

        Parameters
        ----------
        pdf : Union[str, pathlib.Path], optional
            A path to a PDF file, by default `barplots.pdf`.

        Note
        ----
        The cavity nomenclature is based on the integer label. The cavity
        marked with 2, the first integer corresponding to a cavity, is KAA,
        the cavity marked with 3 is KAB, the cavity marked with 4 is KAC
        and so on.

        Note
        ----
        The classes of residues are:

        * Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine.

        * Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine.

        * Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine.

        * Negatively charged (R4): Aspartate, Glutamate.

        * Positively charged (R5): Arginine, Histidine, Lysine.

        * Non-standard (RX): Non-standard residues.
        """
        plot_frequencies(self.frequencies, pdf)

    def export_all(
        self,
        fn: Union[str, pathlib.Path] = "results.toml",
        output: Union[str, pathlib.Path] = "cavity.pdb",
        output_hydropathy: Union[str, pathlib.Path] = "hydropathy.pdb",
        include_frequencies_pdf: bool = False,
        pdf: Union[str, pathlib.Path] = "barplots.pdf",
        nthreads: Optional[int] = None,
    ) -> None:
        """Exports cavities and characterization to PDB-formatted files,
        writes file paths and characterization to a TOML-formatted file, and
        optionally plot bar charts of frequencies (residues and classes of
        residues) in a PDF file.

        Parameters
        ----------
        fn : Union[str, pathlib.Path], optional
            A path to TOML-formatted file for writing file paths and
            cavity characterization (volume, area and interface residues)
            per cavity detected, by default `results.toml`.
        output : Union[str, pathlib.Path], optional
            A path to PDB file for writing cavities, by default `cavity.pdb`.
        output_hydropathy : Union[str, pathlib.Path], optional
            A path to PDB file for writing hydropathy at surface points,
            by default `hydropathy.pdb`.
        include_frequencies_pdf : bool, optional
            Whether to plot frequencies (residues and classes of residues)
            to PDF file, by default False.
        pdf : Union[str, pathlib.Path], optional
            A path to a PDF file, by default `barplots.pdf`.
        nthreads : int, optional
            Number of threads, by default None. If None, the number of threads is
            `os.cpu_count() - 1`.

        Note
        ----
        The cavity nomenclature is based on the integer label. The cavity
        marked with 2, the first integer corresponding to a cavity, is KAA,
        the cavity marked with 3 is KAB, the cavity marked with 4 is KAC
        and so on.

        Note
        ----
        The classes of residues are:

        * Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine.

        * Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine.

        * Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine.

        * Negatively charged (R4): Aspartate, Glutamate.

        * Positively charged (R5): Arginine, Histidine, Lysine.

        * Non-standard (RX): Non-standard residues.
        """
        # Export cavity PDB file
        self.export(output, output_hydropathy, nthreads)
        # Write KVFinder results TOML
        self.write(fn, output, output_hydropathy)
        # Plot bar charts of frequencies
        if include_frequencies_pdf:
            self.plot_frequencies(pdf)


def run_workflow(
    input: Union[str, pathlib.Path],
    ligand: Optional[Union[str, pathlib.Path]] = None,
    vdw: Optional[Union[str, pathlib.Path]] = None,
    box: Optional[Union[str, pathlib.Path]] = None,
    step: Union[float, int] = 0.6,
    probe_in: Union[float, int] = 1.4,
    probe_out: Union[float, int] = 4.0,
    removal_distance: Union[float, int] = 2.4,
    volume_cutoff: Union[float, int] = 5.0,
    ligand_cutoff: Union[float, int] = 5.0,
    include_depth: bool = False,
    include_hydropathy: bool = False,
    hydrophobicity_scale: Union[str, pathlib.Path] = "EisenbergWeiss",
    surface: str = "SES",
    ignore_backbone: bool = False,
    nthreads: Optional[int] = None,
    verbose: bool = False,
) -> pyKVFinderResults:
    """Detects and characterizes cavities (volume, area, depth [optional],
    hydropathy [optional] and interface residues).

    Parameters
    ----------
    input : Union[str, pathlib.Path]
        A path to a target structure file, in PDB or XYZ format, to detect and characterize cavities. 
    ligand : Union[str, pathlib.Path], optional
        A path to ligand file, in PDB or XYZ format, by default None.
    vdw : Union[str, pathlib.Path], optional
        A path to a van der Waals radii file, by default None. If None, apply the built-in van der
        Waals radii file: `vdw.dat`.
    box : Union[str, pathlib.Path], optional
        A path to box configuration file (TOML-formatted), by default None.
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6.
    probe_in : Union[float, int], optional
        Probe In size (A), by default 1.4.
    probe_out : Union[float, int], optional
        Probe Out size (A), by default 4.0.
    removal_distance : Union[float, int], optional
        Length to be removed from the cavity-bulk frontier (A), by default 2.4.
    volume_cutoff : Union[float, int], optional
        Cavities volume filter (A3), by default 5.0.
    ligand_cutoff : Union[float, int], optional
        Radius value to limit a space around a ligand (A), by default 5.0.
    include_depth : bool, optional
        Whether to characterize the depth of the detected cavities, by
        default False.
    include_hydropathy : bool, optional
        Whether to characterize the hydropathy of the detected cavities, by
        default False.
    hydrophobicity_scale : Union[str, pathlib.Path], optional
        Name of a built-in hydrophobicity scale (EisenbergWeiss, HessaHeijne,
        KyteDoolitte, MoonFleming, WimleyWhite, ZhaoLondon) or a path to a
        TOML-formatted file with a custom hydrophobicity scale, by default
        `EisenbergWeiss`.
    surface : str, optional
        Keywords options are SES (Solvent Excluded Surface) or SAS (Solvent
        Accessible Surface), by default SES.
    ignore_backbone : bool, optional
        Whether to ignore backbone atoms (C, CA, N, O) when defining interface
        residues, by default False.
    nthreads : int, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    results : pyKVFinderResults
        A class with the following attributes defined:

        * cavities : numpy.ndarray
            Cavity points in the 3D grid (cavities[nx][ny][nz]).
            Cavities array has integer labels in each position, that are:

                * -1: bulk points;

                * 0: biomolecule points;

                * 1: empty space points;

                * >=2: cavity points.

            The empty space points are regions that do not meet the chosen
            volume cutoff to be considered a cavity.
        * surface : numpy.ndarray
            Surface points in the 3D grid (surface[nx][ny][nz]).
            Surface array has integer labels in each position, that are:

                * -1: bulk points;

                * 0: biomolecule or empty space points;

                * >=2: surface points.

            The empty space points are regions that do not meet the chosen
            volume cutoff to be considered a cavity.
        * depths : numpy.ndarray, optional
            A numpy.ndarray with depth of cavity points (depth[nx][ny][nz]).
        * scales : numpy.ndarray, optional
            A numpy.ndarray with hydrophobicity scale value mapped at surface
            points (scales[nx][ny][nz]).
        * ncav : int
            Number of cavities.
        * volume : Dict[str, float]
            A dictionary with volume of each detected cavity.
        * area : Dict[str, float]
            A dictionary with area of each detected cavity.
        * max_depth : Dict[str, float], optional
            A dictionary with maximum depth of each detected cavity.
        * avg_depth : Dict[str, float], optional
            A dictionary with average depth of each detected cavity.
        * avg_hydropathy : Dict[str, float], optional
            A dictionary with average hydropathy for each detected cavity and
            the range of the hydrophobicity scale (min, max).
        * residues: Dict[str, List[List[str]]]
            A dictionary with a list of interface residues for each detected
            cavity.
        * frequencies : Dict[str, Dict[str, Dict[str, int]]], optional
            A dictionary with frequencies of residues and class for
            residues of each detected cavity.
        * _vertices : numpy.ndarray
            A numpy.ndarray or a list with xyz vertices coordinates (origin,
            X-axis, Y-axis, Z-axis).
        * _step : float
            Grid spacing (A).
        * _input : Union[str, pathlib.Path], optional
            A path to input PDB or XYZ file.
        * _ligand : Union[str, pathlib.Path], optional
            A path to ligand PDB or XYZ file.

    Raises
    ------
    TypeError
        `input` must have .pdb or .xyz extension.
    TypeError
        `ligand` must have .pdb or .xyz extension.

    Note
    ----
    The cavity nomenclature is based on the integer label. The cavity marked
    with 2, the first integer corresponding to a cavity, is KAA, the cavity
    marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

    Note
    ----
    The classes of residues are:

    * Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine.

    * Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine.

    * Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine.

    * Negatively charged (R4): Aspartate, Glutamate.

    * Positively charged (R5): Arginine, Histidine, Lysine.

    * Non-standard (RX): Non-standard residues.
    """
    if verbose:
        print("> Loading atomic dictionary file")
    if vdw is not None:
        vdw = read_vdw(vdw)
    else:
        vdw = read_vdw(VDW)

    if verbose:
        print("> Reading PDB coordinates")
    if input.endswith('.pdb'):
        atomic = read_pdb(input, vdw)
    elif input.endswith('.xyz'):
        atomic = read_xyz(input, vdw)
    else:
        raise TypeError("`target` must have .pdb or .xyz extension.")
    
    if ligand:
        if verbose:
            print("> Reading ligand coordinates")
        if ligand.endswith('.pdb'):
            latomic = read_pdb(ligand, vdw)
        elif ligand.endswith('.xyz'):
            latomic = read_xyz(ligand, vdw)
        else:
            raise TypeError("`ligand` must have .pdb or .xyz extension.")
    else:
        latomic = None

    if verbose:
        print("> Calculating 3D grid dimensions")
    if box:
        # Get vertices from file
        vertices, atomic = get_vertices_from_file(
            box, atomic, step, probe_in, probe_out, nthreads
        )

        # Set flag to boolean
        box = True
    else:
        # Get vertices from input
        vertices = get_vertices(atomic, probe_out, step)

        # Set flag to boolean
        box = False

    # Calculate distance between points
    nx, ny, nz = _get_dimensions(vertices, step)
    if verbose:
        print(f"Dimensions: (nx:{nx}, ny:{ny}, nz:{nz})")

    # Calculate sin and cos of angles a and b
    sincos = _get_sincos(vertices)
    if verbose:
        print(f"sina: {sincos[0]:.2f}\tsinb: {sincos[2]:.2f}")
        print(f"cosa: {sincos[1]:.2f}\tcosb: {sincos[3]:.2f}")

    # Cavity detection
    ncav, cavities = detect(
        atomic,
        vertices,
        step,
        probe_in,
        probe_out,
        removal_distance,
        volume_cutoff,
        latomic,
        ligand_cutoff,
        box,
        surface,
        nthreads,
        verbose,
    )

    if ncav > 0:
        # Spatial characterization
        surface, volume, area = spatial(cavities, step, None, nthreads, verbose)

        # Depth characterization
        if include_depth:
            depths, max_depth, avg_depth = depth(
                cavities, step, None, nthreads, verbose
            )
        else:
            depths, max_depth, avg_depth = None, None, None

        # Constitutional characterization
        residues = constitutional(
            cavities,
            atomic,
            vertices,
            step,
            probe_in,
            ignore_backbone,
            None,
            nthreads,
            verbose,
        )
        frequencies = calculate_frequencies(residues)

        # Hydropathy hydrophobicity scales
        if include_hydropathy:
            scales, avg_hydropathy = hydropathy(
                surface,
                atomic,
                vertices,
                step,
                probe_in,
                hydrophobicity_scale,
                ignore_backbone,
                None,
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
        input,
        ligand,
    )

    return results
