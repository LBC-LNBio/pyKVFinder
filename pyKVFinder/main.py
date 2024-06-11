import logging
import os
import pathlib
import time
from datetime import datetime
from typing import Any, Dict, List, Tuple, Optional, Union

import numpy

from .argparser import argparser
from .grid import (
    _get_dimensions,
    _get_sincos,
    constitutional,
    depth,
    detect,
    export,
    get_vertices,
    get_vertices_from_file,
    hydropathy,
    spatial,
)
from .utils import (
    _write_parameters,
    calculate_frequencies,
    plot_frequencies,
    read_pdb,
    read_vdw,
    read_xyz,
    write_results,
)

__all__ = ["run_workflow", "pyKVFinderResults", "Molecule"]

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
        atomic = read_pdb(args.input, vdw, args.model)
    elif args.input.endswith(".xyz"):
        atomic = read_xyz(args.input, vdw)

    if args.ligand:
        if args.verbose:
            print("> Reading ligand coordinates")
        if args.ligand.endswith(".pdb"):
            latomic = read_pdb(args.ligand, vdw)
        elif args.ligand.endswith(".xyz"):
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
        else:
            scales, avg_hydropathy = None, None

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
    print(f"[ \033[1mElapsed time:\033[0m {elapsed_time:.4f}s ]")
    logging.info(f"[ Elapsed time (s): {elapsed_time:.4f}s ]\n")

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
        nthreads: Optional[int] = None,
    ) -> Optional[str]:
        """Exports cavitiy (H) and surface (HA) points to PDB-formatted file
        with a variable (B; optional) in B-factor column, and hydropathy to
        PDB-formatted file in B-factor column at surface points (HA).

        Parameters
        ----------
        output : Union[str, pathlib.Path]), optional
            A path to PDB file for writing cavities, by default `cavity.pdb`.
        nthreads : int, optional
            Number of threads, by default None. If None, the number of threads is
            `os.cpu_count() - 1`.

        Returns
        -------
        Optional[str]
            A raw string with the PDB-formatted file.

        Note
        ----
        The cavity nomenclature is based on the integer label. The cavity
        marked with 2, the first integer corresponding to a cavity, is KAA,
        the cavity marked with 3 is KAB, the cavity marked with 4 is KAC
        and so on.

        Example
        -------
        >>> import os
        >>> import pyKVFinder
        >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
        >>> results = pyKVFinder.run_workflow(pdb)
        >>> results.export()
        """
        string = export(
            output,
            self.cavities,
            self.surface,
            self._vertices,
            self._step,
            self.depths,
            self.scales,
            None,
            nthreads,
        )
        if output is None:
            return string

    def write(
        self,
        fn: Union[str, pathlib.Path] = "results.toml",
        output: Optional[Union[str, pathlib.Path]] = None,
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

        Note
        ----
        The cavity nomenclature is based on the integer label. The cavity
        marked with 2, the first integer corresponding to a cavity, is KAA,
        the cavity marked with 3 is KAB, the cavity marked with 4 is KAC
        and so on.

        Example
        -------
        >>> import os
        >>> import pyKVFinder
        >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
        >>> results = pyKVFinder.run_workflow(pdb)
        >>> results.write()
        """
        write_results(
            fn,
            self._input,
            self._ligand,
            output,
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

        * Non-standard (RX): Non-standard residues

        Example
        -------
        >>> import os
        >>> import pyKVFinder
        >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
        >>> results = pyKVFinder.run_workflow(pdb)
        >>> results.plot_frequencies()
        """
        plot_frequencies(self.frequencies, pdf)

    def export_all(
        self,
        fn: Union[str, pathlib.Path] = "results.toml",
        output: Union[str, pathlib.Path] = "cavity.pdb",
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

        Example
        -------
        >>> import os
        >>> import pyKVFinder
        >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
        >>> results = pyKVFinder.run_workflow(pdb)
        >>> results.export_all()

        Yet, we can set a ``include_frequencies_pdf`` flag to True to plot the bar charts of the frequencies in a PDF file.

        >>> results.export_all(include_frequencies_pdf=True)
        """
        # Export cavity PDB file
        string = self.export(output, nthreads)
        # Write KVFinder results TOML
        self.write(fn, output)
        # Plot bar charts of frequencies
        if include_frequencies_pdf:
            self.plot_frequencies(pdf)
        # Return PDB-formatted file as a string
        if output is None:
            return string


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
    model: Optional[int] = None,
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
        KyteDoolitte, MoonFleming, RadzickaWolfenden, WimleyWhite, ZhaoLondon)
        or a path to a TOML-formatted file with a custom hydrophobicity scale,
        by default `EisenbergWeiss`.
    surface : str, optional
        Keywords options are SES (Solvent Excluded Surface) or SAS (Solvent
        Accessible Surface), by default SES.
    ignore_backbone : bool, optional
        Whether to ignore backbone atoms (C, CA, N, O) when defining interface
        residues, by default False.
    model : int, optional
        Model number, by default None. If None, keep atoms from all models.
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

    See Also
    --------
    pyKVFinderResults

    Example
    -------
    The **standard workflow** for cavity detection with spatial (surface points, volume, area) and constitutional (interface residues and their frequencies) characterization  can be run at once with one command:

    >>> import os
    >>> import pyKVFinder
    >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
    >>> results = pyKVFinder.run_workflow(pdb)
    >>> results
    <pyKVFinderResults object>
    >>> results.cavities
    array([[[-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            ...,
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1]],
           ...,
           [[-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            ...,
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)
    >>> results.surface
    array([[[-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            ...,
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1]],
           ...,
           [[-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            ...,
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)
    >>> results.ncav
    >>> 18
    >>> results.volume
    {'KAA': 137.16, 'KAB': 47.52, 'KAC': 66.96, 'KAD': 8.21, 'KAE': 43.63, 'KAF': 12.53, 'KAG': 6.26, 'KAH': 520.13, 'KAI': 12.31, 'KAJ': 26.57, 'KAK': 12.31, 'KAL': 33.91, 'KAM': 23.11, 'KAN': 102.82, 'KAO': 6.05, 'KAP': 15.55, 'KAQ': 7.99, 'KAR': 7.78}
    >>> results.area
    {'KAA': 126.41, 'KAB': 62.37, 'KAC': 74.57, 'KAD': 19.06, 'KAE': 57.08, 'KAF': 22.77, 'KAG': 15.38, 'KAH': 496.97, 'KAI': 30.58, 'KAJ': 45.64, 'KAK': 30.58, 'KAL': 45.58, 'KAM': 45.25, 'KAN': 129.77, 'KAO': 12.28, 'KAP': 25.04, 'KAQ': 13.46, 'KAR': 16.6}
    >>> results.residues
    {'KAA': [['14', 'E', 'SER'], ['15', 'E', 'VAL'], ['18', 'E', 'PHE'], ['19', 'E', 'LEU'], ['100', 'E', 'PHE'], ['152', 'E', 'LEU'], ['155', 'E', 'GLU'], ['156', 'E', 'TYR'], ['292', 'E', 'LYS'], ['302', 'E', 'TRP'], ['303', 'E', 'ILE'], ['306', 'E', 'TYR']], 'KAB': [['18', 'E', 'PHE'], ['22', 'E', 'ALA'], ['25', 'E', 'ASP'], ['26', 'E', 'PHE'], ['29', 'E', 'LYS'], ['97', 'E', 'ALA'], ['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['156', 'E', 'TYR']], 'KAC': [['141', 'E', 'PRO'], ['142', 'E', 'HIS'], ['144', 'E', 'ARG'], ['145', 'E', 'PHE'], ['148', 'E', 'ALA'], ['299', 'E', 'THR'], ['300', 'E', 'THR'], ['305', 'E', 'ILE'], ['310', 'E', 'VAL'], ['311', 'E', 'GLU'], ['313', 'E', 'PRO']], 'KAD': [['122', 'E', 'TYR'], ['124', 'E', 'ALA'], ['176', 'E', 'GLN'], ['318', 'E', 'PHE'], ['320', 'E', 'GLY'], ['321', 'E', 'PRO'], ['322', 'E', 'GLY'], ['323', 'E', 'ASP']], 'KAE': [['95', 'E', 'LEU'], ['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['100', 'E', 'PHE'], ['103', 'E', 'LEU'], ['104', 'E', 'VAL'], ['105', 'E', 'LYS'], ['106', 'E', 'LEU']], 'KAF': [['123', 'E', 'VAL'], ['124', 'E', 'ALA'], ['175', 'E', 'ASP'], ['176', 'E', 'GLN'], ['181', 'E', 'GLN']], 'KAG': [['34', 'E', 'SER'], ['37', 'E', 'THR'], ['96', 'E', 'GLN'], ['106', 'E', 'LEU'], ['107', 'E', 'GLU'], ['108', 'E', 'PHE'], ['109', 'E', 'SER']], 'KAH': [['49', 'E', 'LEU'], ['50', 'E', 'GLY'], ['51', 'E', 'THR'], ['52', 'E', 'GLY'], ['53', 'E', 'SER'], ['54', 'E', 'PHE'], ['55', 'E', 'GLY'], ['56', 'E', 'ARG'], ['57', 'E', 'VAL'], ['70', 'E', 'ALA'], ['72', 'E', 'LYS'], ['74', 'E', 'LEU'], ['84', 'E', 'GLN'], ['87', 'E', 'HIS'], ['88', 'E', 'THR'], ['91', 'E', 'GLU'], ['104', 'E', 'VAL'], ['120', 'E', 'MET'], ['121', 'E', 'GLU'], ['122', 'E', 'TYR'], ['123', 'E', 'VAL'], ['127', 'E', 'GLU'], ['166', 'E', 'ASP'], ['168', 'E', 'LYS'], ['170', 'E', 'GLU'], ['171', 'E', 'ASN'], ['173', 'E', 'LEU'], ['183', 'E', 'THR'], ['184', 'E', 'ASP'], ['186', 'E', 'GLY'], ['187', 'E', 'PHE'], ['201', 'E', 'THR'], ['327', 'E', 'PHE']], 'KAI': [['131', 'E', 'HIS'], ['138', 'E', 'PHE'], ['142', 'E', 'HIS'], ['146', 'E', 'TYR'], ['174', 'E', 'ILE'], ['314', 'E', 'PHE']], 'KAJ': [['33', 'E', 'PRO'], ['89', 'E', 'LEU'], ['92', 'E', 'LYS'], ['93', 'E', 'ARG'], ['96', 'E', 'GLN'], ['349', 'E', 'GLU'], ['350', 'E', 'PHE']], 'KAK': [['157', 'E', 'LEU'], ['162', 'E', 'LEU'], ['163', 'E', 'ILE'], ['164', 'E', 'TYR'], ['185', 'E', 'PHE'], ['188', 'E', 'ALA']], 'KAL': [['49', 'E', 'LEU'], ['127', 'E', 'GLU'], ['129', 'E', 'PHE'], ['130', 'E', 'SER'], ['326', 'E', 'ASN'], ['327', 'E', 'PHE'], ['328', 'E', 'ASP'], ['330', 'E', 'TYR']], 'KAM': [['51', 'E', 'THR'], ['55', 'E', 'GLY'], ['56', 'E', 'ARG'], ['73', 'E', 'ILE'], ['74', 'E', 'LEU'], ['75', 'E', 'ASP'], ['115', 'E', 'ASN'], ['335', 'E', 'ILE'], ['336', 'E', 'ARG']], 'KAN': [['165', 'E', 'ARG'], ['166', 'E', 'ASP'], ['167', 'E', 'LEU'], ['199', 'E', 'CYS'], ['200', 'E', 'GLY'], ['201', 'E', 'THR'], ['204', 'E', 'TYR'], ['205', 'E', 'LEU'], ['206', 'E', 'ALA'], ['209', 'E', 'ILE'], ['219', 'E', 'VAL'], ['220', 'E', 'ASP'], ['223', 'E', 'ALA']], 'KAO': [['48', 'E', 'THR'], ['51', 'E', 'THR'], ['56', 'E', 'ARG'], ['330', 'E', 'TYR'], ['331', 'E', 'GLU']], 'KAP': [['222', 'E', 'TRP'], ['238', 'E', 'PHE'], ['253', 'E', 'GLY'], ['254', 'E', 'LYS'], ['255', 'E', 'VAL'], ['273', 'E', 'LEU']], 'KAQ': [['207', 'E', 'PRO'], ['208', 'E', 'GLU'], ['211', 'E', 'LEU'], ['213', 'E', 'LYS'], ['275', 'E', 'VAL'], ['277', 'E', 'LEU']], 'KAR': [['237', 'E', 'PRO'], ['238', 'E', 'PHE'], ['249', 'E', 'LYS'], ['254', 'E', 'LYS'], ['255', 'E', 'VAL'], ['256', 'E', 'ARG']]}
    >>> results.frequencies
    {'KAA': {'RESIDUES': {'GLU': 1, 'ILE': 1, 'LEU': 2, 'LYS': 1, 'PHE': 2, 'SER': 1, 'TRP': 1, 'TYR': 2, 'VAL': 1}, 'CLASS': {'R1': 4, 'R2': 5, 'R3': 1, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAB': {'RESIDUES': {'ALA': 2, 'ASN': 1, 'ASP': 1, 'LYS': 1, 'PHE': 2, 'TYR': 1, 'VAL': 1}, 'CLASS': {'R1': 3, 'R2': 3, 'R3': 1, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAC': {'RESIDUES': {'ALA': 1, 'ARG': 1, 'GLU': 1, 'HIS': 1, 'ILE': 1, 'PHE': 1, 'PRO': 2, 'THR': 2, 'VAL': 1}, 'CLASS': {'R1': 5, 'R2': 1, 'R3': 2, 'R4': 1, 'R5': 2, 'RX': 0}}, 'KAD': {'RESIDUES': {'ALA': 1, 'ASP': 1, 'GLN': 1, 'GLY': 2, 'PHE': 1, 'PRO': 1, 'TYR': 1}, 'CLASS': {'R1': 4, 'R2': 2, 'R3': 1, 'R4': 1, 'R5': 0, 'RX': 0}}, 'KAE': {'RESIDUES': {'ASN': 1, 'LEU': 3, 'LYS': 1, 'PHE': 1, 'VAL': 2}, 'CLASS': {'R1': 5, 'R2': 1, 'R3': 1, 'R4': 0, 'R5': 1, 'RX': 0}}, 'KAF': {'RESIDUES': {'ALA': 1, 'ASP': 1, 'GLN': 2, 'VAL': 1}, 'CLASS': {'R1': 2, 'R2': 0, 'R3': 2, 'R4': 1, 'R5': 0, 'RX': 0}}, 'KAG': {'RESIDUES': {'GLN': 1, 'GLU': 1, 'LEU': 1, 'PHE': 1, 'SER': 2, 'THR': 1}, 'CLASS': {'R1': 1, 'R2': 1, 'R3': 4, 'R4': 1, 'R5': 0, 'RX': 0}}, 'KAH': {'RESIDUES': {'ALA': 1, 'ARG': 1, 'ASN': 1, 'ASP': 2, 'GLN': 1, 'GLU': 4, 'GLY': 4, 'HIS': 1, 'LEU': 3, 'LYS': 2, 'MET': 1, 'PHE': 3, 'SER': 1, 'THR': 4, 'TYR': 1, 'VAL': 3}, 'CLASS': {'R1': 11, 'R2': 4, 'R3': 8, 'R4': 6, 'R5': 4, 'RX': 0}}, 'KAI': {'RESIDUES': {'HIS': 2, 'ILE': 1, 'PHE': 2, 'TYR': 1}, 'CLASS': {'R1': 1, 'R2': 3, 'R3': 0, 'R4': 0, 'R5': 2, 'RX': 0}}, 'KAJ': {'RESIDUES': {'ARG': 1, 'GLN': 1, 'GLU': 1, 'LEU': 1, 'LYS': 1, 'PHE': 1, 'PRO': 1}, 'CLASS': {'R1': 2, 'R2': 1, 'R3': 1, 'R4': 1, 'R5': 2, 'RX': 0}}, 'KAK': {'RESIDUES': {'ALA': 1, 'ILE': 1, 'LEU': 2, 'PHE': 1, 'TYR': 1}, 'CLASS': {'R1': 4, 'R2': 2, 'R3': 0, 'R4': 0, 'R5': 0, 'RX': 0}}, 'KAL': {'RESIDUES': {'ASN': 1, 'ASP': 1, 'GLU': 1, 'LEU': 1, 'PHE': 2, 'SER': 1, 'TYR': 1}, 'CLASS': {'R1': 1, 'R2': 3, 'R3': 2, 'R4': 2, 'R5': 0, 'RX': 0}}, 'KAM': {'RESIDUES': {'ARG': 2, 'ASN': 1, 'ASP': 1, 'GLY': 1, 'ILE': 2, 'LEU': 1, 'THR': 1}, 'CLASS': {'R1': 4, 'R2': 0, 'R3': 2, 'R4': 1, 'R5': 2, 'RX': 0}}, 'KAN': {'RESIDUES': {'ALA': 2, 'ARG': 1, 'ASP': 2, 'CYS': 1, 'GLY': 1, 'ILE': 1, 'LEU': 2, 'THR': 1, 'TYR': 1, 'VAL': 1}, 'CLASS': {'R1': 7, 'R2': 1, 'R3': 2, 'R4': 2, 'R5': 1, 'RX': 0}}, 'KAO': {'RESIDUES': {'ARG': 1, 'GLU': 1, 'THR': 2, 'TYR': 1}, 'CLASS': {'R1': 0, 'R2': 1, 'R3': 2, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAP': {'RESIDUES': {'GLY': 1, 'LEU': 1, 'LYS': 1, 'PHE': 1, 'TRP': 1, 'VAL': 1}, 'CLASS': {'R1': 3, 'R2': 2, 'R3': 0, 'R4': 0, 'R5': 1, 'RX': 0}}, 'KAQ': {'RESIDUES': {'GLU': 1, 'LEU': 2, 'LYS': 1, 'PRO': 1, 'VAL': 1}, 'CLASS': {'R1': 4, 'R2': 0, 'R3': 0, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAR': {'RESIDUES': {'ARG': 1, 'LYS': 2, 'PHE': 1, 'PRO': 1, 'VAL': 1}, 'CLASS': {'R1': 2, 'R2': 1, 'R3': 0, 'R4': 0, 'R5': 3, 'RX': 0}}}

    However, users may opt to perform cavity detection in a segmented space through ligand adjustment and/or box adjustment modes.

    The cavity detection can be limited around the target ligand(s), which will be passed to pyKVFinder through a *.pdb* or a *.xyz* files. Thus, the detected cavities are limited within a radius (``ligand_cutoff``) of the target ligand(s).

    >>> ligand = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', 'ADN.pdb')
    >>> results = pyKVFinder.run_workflow(pdb, ligand)
    >>> results
    <pyKVFinderResults object>
    >>> results.cavities
    array([[[-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            ...,
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1]],
           ...,
           [[-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            ...,
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)
    >>> results.surface
    array([[[-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            ...,
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1]],
           ...,
           [[-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            ...,
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)
    >>> results.ncav
    >>> 18
    >>> results.volume
    {'KAA': 365.04, 'KAB': 16.85}
    >>> results.area
    {'KAA': 328.79, 'KAB': 23.15}
    >>> results.residues
    {'KAA': [['49', 'E', 'LEU'], ['50', 'E', 'GLY'], ['51', 'E', 'THR'], ['52', 'E', 'GLY'], ['53', 'E', 'SER'], ['55', 'E', 'GLY'], ['56', 'E', 'ARG'], ['57', 'E', 'VAL'], ['70', 'E', 'ALA'], ['72', 'E', 'LYS'], ['104', 'E', 'VAL'], ['120', 'E', 'MET'], ['121', 'E', 'GLU'], ['122', 'E', 'TYR'], ['123', 'E', 'VAL'], ['127', 'E', 'GLU'], ['166', 'E', 'ASP'], ['168', 'E', 'LYS'], ['170', 'E', 'GLU'], ['171', 'E', 'ASN'], ['173', 'E', 'LEU'], ['183', 'E', 'THR'], ['184', 'E', 'ASP'], ['327', 'E', 'PHE']], 'KAB': [['49', 'E', 'LEU'], ['127', 'E', 'GLU'], ['130', 'E', 'SER'], ['326', 'E', 'ASN'], ['327', 'E', 'PHE'], ['328', 'E', 'ASP'], ['330', 'E', 'TYR']]}
    >>> results.frequencies
    {'KAA': {'RESIDUES': {'ALA': 1, 'ARG': 1, 'ASN': 1, 'ASP': 2, 'GLU': 3, 'GLY': 3, 'LEU': 2, 'LYS': 2, 'MET': 1, 'PHE': 1, 'SER': 1, 'THR': 2, 'TYR': 1, 'VAL': 3}, 'CLASS': {'R1': 9, 'R2': 2, 'R3': 5, 'R4': 5, 'R5': 3, 'RX': 0}}, 'KAB': {'RESIDUES': {'ASN': 1, 'ASP': 1, 'GLU': 1, 'LEU': 1, 'PHE': 1, 'SER': 1, 'TYR': 1}, 'CLASS': {'R1': 1, 'R2': 2, 'R3': 2, 'R4': 2, 'R5': 0, 'RX': 0}}}

    Further, we can also perform cavity detection on a custom 3D grid, where we can explore closed regions with a custom box, which can be defined by a *.toml* file (see `Box configuration file template`).

    >>> fn = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', 'custom-box.toml')
    >>> with open(fn, 'r') as f:
    ...     print(f.read())
    [box]
    p1 = [3.11, 7.34, 1.59]
    p2 = [11.51, 7.34, 1.59]
    p3 = [3.11, 10.74, 1.59]
    p4 = [3.11, 7.34, 6.19]
    >>> results = pyKVFinder.run_workflow(pdb, box=fn)
    >>> results
    <pyKVFinderResults object>
    >>> results.cavities
    array([[[-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            ...,
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1]],
           ...,
           [[-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            ...,
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)
    >>> results.surface
    array([[[-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            ...,
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1]],
           ...,
           [[-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            ...,
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1],
            [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)
    >>> results.ncav
    >>> 1
    >>> results.volume
    {'KAA': 115.78}
    >>> results.area
    {'KAA': 33.91}
    >>> results.residues
    {'KAA': [['49', 'E', 'LEU'], ['50', 'E', 'GLY'], ['51', 'E', 'THR'], ['57', 'E', 'VAL'], ['70', 'E', 'ALA'], ['104', 'E', 'VAL'], ['121', 'E', 'GLU'], ['122', 'E', 'TYR'], ['123', 'E', 'VAL'], ['127', 'E', 'GLU'], ['170', 'E', 'GLU'], ['171', 'E', 'ASN'], ['173', 'E', 'LEU'], ['183', 'E', 'THR'], ['327', 'E', 'PHE']]}
    >>> results.frequencies
    {'KAA': {'RESIDUES': {'ALA': 1, 'ASN': 1, 'GLU': 3, 'GLY': 1, 'LEU': 2, 'PHE': 1, 'THR': 2, 'TYR': 1, 'VAL': 3}, 'CLASS': {'R1': 7, 'R2': 2, 'R3': 3, 'R4': 3, 'R5': 0, 'RX': 0}}}

    However, users may opt to perform the **full workflow** for cavity detection with spatial (surface points, volume and area), constitutional (interface residues and their frequencies), hydropathy and depth characterization. This full workflow can be run with one command by setting some parameters of ``run_workflow`` function:

    >>> results = pyKVFinder.run_workflow(pdb, include_depth=True, include_hydropathy=True, hydrophobicity_scale='EisenbergWeiss')
    >>> results.depths
    array([[[0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            ...,
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.]],
           ...,
           [[0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            ...,
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.]]])
    >>> results.scales
    array([[[0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            ...,
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.]],
           ...,
           [[0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            ...,
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.],
            [0., 0., 0., ..., 0., 0., 0.]]])
    >>> results.avg_depth
    {'KAA': 1.35, 'KAB': 0.91, 'KAC': 0.68, 'KAD': 0.32, 'KAE': 0.99, 'KAF': 0.24, 'KAG': 0.1, 'KAH': 3.91, 'KAI': 0.0, 'KAJ': 0.96, 'KAK': 0.0, 'KAL': 1.07, 'KAM': 0.24, 'KAN': 0.0, 'KAO': 0.29, 'KAP': 0.7, 'KAQ': 0.22, 'KAR': 0.12}
    >>> results.max_depth
    {'KAA': 3.79, 'KAB': 2.68, 'KAC': 2.62, 'KAD': 0.85, 'KAE': 3.0, 'KAF': 0.85, 'KAG': 0.6, 'KAH': 10.73, 'KAI': 0.0, 'KAJ': 2.24, 'KAK': 0.0, 'KAL': 3.0, 'KAM': 1.2, 'KAN': 0.0, 'KAO': 1.04, 'KAP': 2.08, 'KAQ': 0.85, 'KAR': 0.6}
    >>> results.avg_hydropathy
    {'KAA': -0.73, 'KAB': -0.05, 'KAC': -0.07, 'KAD': -0.62, 'KAE': -0.81, 'KAF': -0.14, 'KAG': -0.33, 'KAH': -0.17, 'KAI': -0.4, 'KAJ': 0.62, 'KAK': -0.99, 'KAL': 0.36, 'KAM': -0.33, 'KAN': 0.18, 'KAO': 0.88, 'KAP': -0.96, 'KAQ': 0.48, 'KAR': 0.24, 'EisenbergWeiss': [-1.42, 2.6]}
    """
    if verbose:
        print("> Loading atomic dictionary file")
    if vdw is not None:
        vdw = read_vdw(vdw)
    else:
        vdw = read_vdw(VDW)

    if verbose:
        print("> Reading PDB coordinates")
    if input.endswith(".pdb"):
        atomic = read_pdb(input, vdw, model)
    elif input.endswith(".xyz"):
        atomic = read_xyz(input, vdw)
    else:
        raise TypeError("`target` must have .pdb or .xyz extension.")

    if ligand:
        if verbose:
            print("> Reading ligand coordinates")
        if ligand.endswith(".pdb"):
            latomic = read_pdb(ligand, vdw)
        elif ligand.endswith(".xyz"):
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

        # Depth characterization
        if include_depth:
            depths, max_depth, avg_depth = depth(
                cavities, step, None, nthreads, verbose
            )
        else:
            depths, max_depth, avg_depth = None, None, None

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


class Molecule(object):
    """A class for representing molecular structures.

    Parameters
    ----------
    molecule : Union[str, pathlib.Path]
        A file path to the molecule in either PDB or XYZ format
    radii : Union[str, pathlib.Path, Dict[str, Any]], optional
        A file path to a van der Waals radii file or a dictionary of VDW radii, by default None. If None, apply the built-in van der Waals radii file: `vdw.dat`.
    model : int, optional
        The model number of a multi-model PDB file, by default None. If None, keep atoms from all models.
    nthreads : int, optional
        Number of threads, by default None. If None, the number of threads is `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Attributes
    ----------
    _atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates and radius) for each atom.
    _dim : tuple
        Grid dimensions.
    _grid : numpy.ndarray
        Molecule points in the 3D grid (grid[nx][ny][nz]).
        Grid array has integer labels in each position, that are:

            * 0: molecule points;

            * 1: solvent points.
    _molecule : Union[str, pathlib.Path]
        A file path to the molecule in either PDB or XYZ format.
    _padding : float
        The length to add to each direction of the 3D grid.
    _probe : float
        Spherical probe size to define the molecular surface based on a molecular representation.
    _radii : Dict[str, Any]
        A dictionary containing radii values, by default None.
    _representation : str, optional
        Molecular surface representation. Keywords options are vdW (van der Waals surface), SES (Solvent Excluded Surface) or SAS (Solvent Accessible Surface), by default SES.
    _rotation : numpy.ndarray
        A numpy.ndarray with sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb).
    _step : float
        Grid spacing (A).
    _vertices : numpy.ndarray
        A numpy.ndarray or a list with xyz vertices coordinates (origin, X-axis, Y-axis, Z-axis).
    nthreads : int
        Number of threads for parallel processing.
    verbose : bool
        Whether to print extra information to standard output.

    Note
    ----
    The van der Waals radii file defines the radius values for each atom by residue and when not defined, it uses a generic value based on the atom type. The function by default loads the built-in van der Waals radii file: ``vdw.dat``.

    See Also
    --------
    read_vdw

    Example
    -------
    The ``Molecule`` class loads the target molecular structure (ClO4) into pyKVFinder. class.

    >>> import os
    >>> from pyKVFinder import Molecule
    >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', 'ClO4.pdb')
    >>> molecule = Molecule(pdb)
    >>> molecule
    >>> <pyKVFinder.main.Molecule object at 0x7f5ddacf2230>

    The van der Waals radii can be define by:

        * creating a Python dictionary:

        >>> # PyMOL (v2.5.0) vdW radii values
        >>> vdw = {'GEN': {'CL': 1.75, 'O': 1.52}}
        >>> molecule = Molecule(pdb, radii=vdw)
        >>> molecule.radii
        {'GEN': {'CL': 1.75, 'O': 1.52}}

        * specifying a *.dat* file following template of `van der Waals radii file`.

        >>> from pyKVFinder import read_vdw
        >>> # ChimeraX vdW radii values
        >>> with open('vdw.dat', 'w') as f:
        ...     f.write('>GEN\\nCL\\t\\t1.98\\nO\\t\\t1.46\\n')
        >>> vdw = read_vdw('vdw.dat')
        >>> molecule = Molecule(pdb, radii=vdw)
        >>> molecule.radii
        {'GEN': {'CL': 1.98, 'O': 1.46}}
    """

    def __init__(
        self,
        molecule: Union[str, pathlib.Path],
        radii: Union[str, pathlib.Path, Dict[str, Any]] = None,
        model: Optional[int] = None,
        nthreads: Optional[int] = None,
        verbose: bool = False,
    ):
        """Initialize the Molecule object with molecule, radii, model, nthreads and verbose.

        Parameters
        ----------
        molecule : Union[str, pathlib.Path]
            A file path to the molecule in either PDB or XYZ format.
        radii : Union[str, pathlib.Path, Dict[str, Any]], optional
            A file path to a van der Waals radii file or a dictionary of VDW radii, by default None. If None, apply the built-in van der Waals radii file: `vdw.dat`.
        model : int, optional
            The model number of a multi-model PDB file, by default None. If None, keep atoms from all models.
        nthreads : int, optional
            Number of threads, by default None. If None, the number of threads is `os.cpu_count() - 1`.
        verbose : bool, optional
            Print extra information to standard output, by default False.

        Raises
        ------
        TypeError
            `molecule` must be a string or a pathlib.Path.
        TypeError
            `molecule` must have .pdb or .xyz extension.
        TypeError
            `nthreads` must be a positive integer.
        ValueError
            `nthreads` must be a positive integer.
        """

        # Attributes
        self._grid = None
        self._step = None
        self._padding = None
        self._probe = None
        self._representation = None
        self._vertices = None
        self._dim = None
        self._rotation = None
        self.verbose = verbose

        # Molecule
        if type(molecule) not in [str, pathlib.Path]:
            raise TypeError("`molecule` must be a string or a pathlib.Path.")
        self._molecule = os.path.realpath(molecule)

        # van der Waals radii
        if self.verbose:
            print("> Loading van der Waals radii")
        if radii is None:
            # default
            self._radii = read_vdw(VDW)
        elif type(radii) in [str, pathlib.Path]:
            # vdw file
            self._radii = read_vdw(radii)
        elif type(radii) in [dict]:
            # Processed dictionary
            self._radii = radii

        # Atomic information
        if self.verbose:
            print("> Reading molecule coordinates")
        if molecule.endswith(".pdb"):
            self._atomic = read_pdb(molecule, self.radii, model)
        elif molecule.endswith(".xyz"):
            self._atomic = read_xyz(molecule, self.radii)
        else:
            raise TypeError("`molecule` must have .pdb or .xyz extension.")

        # Number of threads
        if nthreads is not None:
            if type(nthreads) not in [int]:
                raise TypeError("`nthreads` must be a positive integer.")
            elif nthreads <= 0:
                raise ValueError("`nthreads` must be a positive integer.")
            else:
                self.nthreads = nthreads
        else:
            self.nthreads = os.cpu_count() - 1

    @property
    def atomic(self) -> numpy.ndarray:
        """Get _atomic attribute."""
        return self._atomic

    @property
    def dim(self) -> Tuple[int, int, int]:
        """Get _dim attribute"""
        return self._dim

    @property
    def grid(self) -> numpy.ndarray:
        """Get _grid attribute."""
        return self._grid

    @property
    def molecule(self) -> Union[str, pathlib.Path]:
        """Get _molecule attribute."""
        return self._molecule

    @property
    def nx(self) -> int:
        """Get grid units in X-axis."""
        if self._dim is not None:
            return self._dim[0]

    @property
    def ny(self) -> int:
        """Get grid units in Y-axis."""
        if self._dim is not None:
            return self._dim[1]

    @property
    def nz(self) -> int:
        """Get grid units in Z-axis."""
        if self._dim is not None:
            return self._dim[2]

    @property
    def p1(self) -> numpy.ndarray:
        """Get origin of the 3D grid."""
        if self._vertices is not None:
            return self._vertices[0]

    @property
    def p2(self) -> numpy.ndarray:
        """Get X-axis max of the 3D grid."""
        if self._vertices is not None:
            return self._vertices[1]

    @property
    def p3(self) -> numpy.ndarray:
        """Get Y-axis max of the 3D grid."""
        if self._vertices is not None:
            return self._vertices[2]

    @property
    def p4(self) -> numpy.ndarray:
        """Get Z-axis max of the 3D grid."""
        if self._vertices is not None:
            return self._vertices[3]

    @property
    def padding(self) -> float:
        """Get _padding attribute."""
        return self._padding

    @property
    def probe(self) -> float:
        """Get _probe attribute."""
        return self._probe

    @property
    def radii(self) -> Dict[str, Any]:
        """Get _radii attribute."""
        return self._radii

    @property
    def representation(self) -> str:
        """Get _representation attribute."""
        return self._representation

    @property
    def rotation(self) -> numpy.ndarray:
        """Get _rotation attribute."""
        return self._rotation

    @property
    def step(self) -> float:
        """Get _step attribute."""
        if self._step is not None:
            return self._step

    @property
    def vertices(self) -> numpy.ndarray:
        """Get _vertices attribute."""
        return self._vertices

    @property
    def xyzr(self) -> numpy.ndarray:
        """Get xyz coordinates and radius of molecule atoms."""
        return self._atomic[:, 4:].astype(numpy.float64)

    def _set_grid(self, padding: Optional[float] = None) -> None:
        """Define the 3D grid for the target molecule.

        Parameters
        ----------
        padding : float, optional
            The length to add to each direction of the 3D grid, by default None. If None, automatically define the length based on molecule coordinates, probe size, grid spacing and atom radii.

        Raises
        ------
        TypeError
            `padding` must be a non-negative real number.
        ValueError
            `padding` must be a non-negative real number.
        """
        # Padding
        if padding is not None:
            if type(padding) not in [int, float]:
                raise TypeError("`padding` must be a non-negative real number.")
            elif padding < 0.0:
                raise ValueError("`padding` must be a non-negative real number.")
            else:
                self._padding = padding
        else:
            self._padding = self._get_padding()

        # 3D grid
        if self.verbose:
            print("> Calculating 3D grid")
        self._vertices = get_vertices(self.atomic, self.padding, self.step)
        self._dim = _get_dimensions(self.vertices, self.step)
        self._rotation = _get_sincos(self.vertices)
        if self.verbose:
            print(f"p1: {self.vertices[0]}")
            print(f"p2: {self.vertices[1]}")
            print(f"p3: {self.vertices[2]}")
            print(f"p4: {self.vertices[3]}")
            print("nx: {}, ny: {}, nz: {}".format(*self.dim))
            print("sina: {}, sinb: {}, cosa: {}, cosb: {}".format(*self.rotation))

    def _get_padding(self) -> float:
        """Automatically define the padding based on molecule coordinates, probe size, grid spacing and atom radii.

        Returns
        -------
        padding : float
            The length to add to each direction of the 3D grid.
        """
        padding = 1.1 * self.xyzr[:, 3].max()
        if self.representation in ["SES", "SAS"]:
            padding += self._probe
        return float(padding.round(decimals=1))

    def vdw(self, step: float = 0.6, padding: Optional[float] = None) -> None:
        """Fill the 3D grid with the molecule as the van der Waals surface representation.

        Parameters
        ----------
        step : float, optional
            Grid spacing (A), by default 0.6.
        padding : float, optional
            The length to add to each direction of the 3D grid, by default None. If None, automatically define the length based on molecule coordinates, probe size, grid spacing and atom radii.

        Raises
        ------
        TypeError
            `step` must be a positive real number.
        ValueError
            `step` must be a positive real number.

        Example
        -------
        The ``Molecule.vdw()`` method takes a grid spacing and returns a NumPy array with the molecule points representing the vdW surface in the 3D grid.

        >>> # Grid Spacing (step): 0.1
        >>> step = 0.1
        >>> molecule.vdw(step=step)
        >>> molecule.grid
        array([[[1, 1, 1, ..., 1, 1, 1],
                [1, 1, 1, ..., 1, 1, 1],
                [1, 1, 1, ..., 1, 1, 1],
                ...,
                [1, 1, 1, ..., 1, 1, 1],
                [1, 1, 1, ..., 1, 1, 1],
                [1, 1, 1, ..., 1, 1, 1]],
               ...,
               [[1, 1, 1, ..., 1, 1, 1],
                [1, 1, 1, ..., 1, 1, 1],
                [1, 1, 1, ..., 1, 1, 1],
                ...,
                [1, 1, 1, ..., 1, 1, 1],
                [1, 1, 1, ..., 1, 1, 1],
                [1, 1, 1, ..., 1, 1, 1]]], dtype=int32)
        """
        from _pyKVFinder import _fill_receptor

        # Check arguments
        if type(step) not in [int, float]:
            raise TypeError("`step` must be a postive real number.")
        elif step <= 0.0:
            raise ValueError("`step` must be a positive real number.")
        else:
            self._step = step

        # Attributes
        self._representation = "vdW"
        self._probe = None

        # Define 3D grid
        self._set_grid(padding)

        # van der Waals atoms (hard sphere model) to grid
        if self.verbose:
            print("> Inserting atoms with van der Waals radii into 3D grid")
        self._grid = _fill_receptor(
            self.nx * self.ny * self.nz,
            self.nx,
            self.ny,
            self.nz,
            self.xyzr,
            self.p1,
            self.rotation,
            self.step,
            0.0,
            False,
            self.nthreads,
            self.verbose,
        ).reshape(self.nx, self.ny, self.nz)

    def surface(
        self,
        step: float = 0.6,
        probe: float = 1.4,
        surface: str = "SES",
        padding: Optional[float] = None,
    ) -> None:
        """Fill the 3D grid with the molecule as the van der Waals surface representation.

        Parameters
        ----------
        step : float, optional
            Grid spacing (A), by default 0.6.
        probe : float, optional
            Spherical probe size to define the molecular surface based on a molecular representation, by default 1.4.
        surface : str, optional
            Molecular surface representation. Keywords options are vdW (van der Waals surface), SES (Solvent Excluded Surface) or SAS (Solvent Accessible Surface), by default "SES".
        padding : float, optional
            The length to add to each direction of the 3D grid, by default None. If None, automatically define the length based on molecule coordinates, probe size, grid spacing and atom radii.

        Raises
        ------
        TypeError
            `step` must be a positive real number.
        ValueError
            `step` must be a positive real number.
        TypeError
            `probe_out` must be a positive real number.
        ValueError
            `probe_out` must be a positive real number.

        Example
        -------
        The ``Molecule.surface()`` method takes the grid spacing, the spherical probe size to model the surface, the surface representation and returns a NumPy array with the molecule points representing the SES in the 3D grid.

        The molecular surface can be modelled as:

            * Solvent Excluded Surface (SES):

            >>> # Surface Representation: SES
            >>> surface = 'SES'
            >>> # Grid Spacing (step): 0.1
            >>> step = 0.1
            >>> # Spherical Probe (probe): 1.4
            >>> probe = 1.4
            >>> molecule.surface(step=step, probe=probe, surface=surface)
            >>> molecule.grid
            array([[[1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1],
                    ...,
                    [1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1]],
                   ...,
                   [[1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1],
                    ...,
                    [1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1]]], dtype=int32)

        The molecular surface can be modelled as:

            * Solvent Accessible Surface (SAS):

            >>> # Surface Representation: SAS
            >>> surface = 'SAS'
            >>> # Grid Spacing (step): 0.1
            >>> step = 0.1
            >>> # Spherical Probe (probe): 1.4
            >>> probe = 1.4
            >>> molecule.surface(step=step, probe=probe, surface=surface)
            >>> molecule.grid
            array([[[1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1],
                    ...,
                    [1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1]],
                   ...,
                   [[1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1],
                    ...,
                    [1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1],
                    [1, 1, 1, ..., 1, 1, 1]]], dtype=int32)
        """
        from _pyKVFinder import _fill_receptor

        # Check arguments
        if type(step) not in [int, float]:
            raise TypeError("`step` must be a postive real number.")
        elif step <= 0.0:
            raise ValueError("`step` must be a positive real number.")
        else:
            self._step = step

        # Probe
        if type(probe) not in [int, float, numpy.float64]:
            raise TypeError("`probe_out` must be a non-negative real number.")
        elif probe <= 0.0:
            raise ValueError("`probe_out` must be a non-negative real number.")
        self._probe = probe

        # Surface
        if surface == "SES":
            if self.verbose:
                print("> Surface representation: Solvent Excluded Surface (SES).")
            self._representation = surface
            surface = True
        elif surface == "SAS":
            if self.verbose:
                print("> Surface representation: Solvent Accessible Surface (SAS).")
            self._representation = surface
            surface = False
        else:
            raise ValueError(f"`surface` must be SAS or SES, not {surface}.")

        # Define 3D grid
        self._set_grid(padding)

        # Molecular surface (SES or SAS) to grid
        self._grid = _fill_receptor(
            self.nx * self.ny * self.nz,
            self.nx,
            self.ny,
            self.nz,
            self.xyzr,
            self.p1,
            self.rotation,
            self.step,
            self.probe,
            surface,
            self.nthreads,
            self.verbose,
        ).reshape(self.nx, self.ny, self.nz)

    def volume(self) -> float:
        """Estimate the volume of the molecule based on the molecular surface representation, ie, vdW, SES or SAS representations.

        Returns
        -------
        volume : float
            Molecular volume (A³).

        Example
        -------
        With the molecular surface modelled by ``Molecule.vdw()`` or ``Molecule.surface()``, the volume can be estimated by running:

        >>> molecule.volume()
        90.8
        """
        from _pyKVFinder import _volume

        if self.grid is not None:
            volume = _volume(
                (self.grid == 0).astype(numpy.int32) * 2, self.step, 1, self.nthreads
            )
            return volume.round(decimals=2).item()

    def preview(self, **kwargs) -> None:
        """Preview the molecular surface in the 3D grid.

        Example
        -------
        With the molecular surface modelled by ``Molecule.vdw()`` or ``Molecule.surface()``, the modelled molecule in the 3D grid can be previewed by running:

        >>> molecule.preview()
        """
        if self.grid is not None:
            from plotly.express import scatter_3d

            x, y, z = numpy.nonzero(self.grid == 0)
            fig = scatter_3d(x=x, y=y, z=z, **kwargs)
            fig.update_layout(
                scene_xaxis_showticklabels=False,
                scene_yaxis_showticklabels=False,
                scene_zaxis_showticklabels=False,
            )
            fig.show()

    def export(
        self,
        fn: Union[str, pathlib.Path] = "molecule.pdb",
    ) -> None:
        """Export molecule points (H) to a PDB-formatted file.

        Parameters
        ----------
        fn : Union[str, pathlib.Path], optional
            A file path to the molecular volume in the grid-based rerpesentation in PDB format, by default "molecule.pdb".

        Raises
        ------
        TypeError
            `fn` must be a string or a pathlib.Path.

        Example
        -------
        With the molecular surface modelled by ``Molecule.vdw()`` or ``Molecule.surface()``, the modelled molecule in the 3D grid can be exported to a PDB-formatted file by running:

        >>> molecule.export('model.pdb')
        """
        # Filename (fn)
        if type(fn) not in [str, pathlib.Path]:
            raise TypeError("`fn` must be a string or a pathlib.Path.")
        os.makedirs(os.path.abspath(os.path.dirname(fn)), exist_ok=True)

        # Save grid to PDB file
        export(
            fn, (self.grid == 0).astype(numpy.int32) * 2, None, self.vertices, self.step
        )
