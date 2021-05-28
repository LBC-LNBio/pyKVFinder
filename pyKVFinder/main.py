import os
import time
import logging
import numpy
import pathlib
from datetime import datetime
from typing import Union, Dict, List
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

_vdw = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/vdw.dat")


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
                      [-B <.toml>] [-L <.pdb>] [--ligand_cutoff <float>]
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

        # Calculate sin and cos of angles a and b
        args.sincos = get_sincos(args.vertices)

        # Set flag to boolean
        args.box = False

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
    f"""A class containing pyKVFinder detection and characterization results.

    Attributes
    ----------
    cavities : numpy.ndarray
        Cavity points in the 3D grid (cavities[nx][ny][nz]).
    surface : numpy.ndarray
        Surface points in the 3D grid (surface[nx][ny][nz])
    depths : Union[numpy.ndarray, None]
        A numpy.ndarray with depth of cavity points (depth[nx][ny][nz]).
    scales : Union[numpy.ndarray, None]
        A numpy.ndarray with hydrophobicity scale value mapped at surface
        points (scales[nx][ny][nz]).
    volume : Dict[str, float]
        A dictionary with volume of each detected cavity.
    area : Dict[str, float]
        A dictionary with area of each detected cavity.
    max_depth : Union[Dict[str, float], None]
        A dictionary with maximum depth of each detected cavity.
    avg_depth : Union[Dict[str, float], None]
        A dictionary with average depth of each detected cavity.
    avg_hydropathy : Union[Dict[str, float], None]
        A dictionary with average hydropathy for each detected cavity and
        the range of the hydrophobicity scale (min, max).
    residues: Dict[str, List[List[str]]]
        A dictionary with a list of interface residues for each detected
        cavity.
    frequencies : Union[Dict[str, Dict[str, Dict[str, int]]], None]
        A dictionary with frequencies of residues and class for
        residues of each detected cavity.
    _vertices : numpy.ndarray
        A numpy.ndarray or a list with xyz vertices coordinates (origin,
        X-axis, Y-axis, Z-axis).
    _step : float
        Grid spacing (A).
    _ncav : int
        Number of cavities.
    _pdb : Union[str, pathlib.Path, None], optional
        A path to input PDB file, by default None.
    _ligand : Union[str, pathlib.Path, None], optional
        A path to ligand PDB file, by default None.

    Methods
    -------
    export(output = 'cavity.pdb', output_hydropathy = 'hydropathy.pdb',
    nthreads = {os.cpu_count() - 1})
        Exports cavities and characterizations to PDB-formatted files.
    write(fn = 'results.toml', output = None, output_hydropathy = None)
        Writes file paths and characterizations to a TOML-formatted
        results file.
    plot_frequencies(pdf = 'histograms.pdf')
        Plot histograms of frequencies (residues and classes of residues) in a
        PDF file.
    export_all(fn = 'results.toml', output = 'cavity.pdb',
    output_hydropathy = 'hydropathy.pdb', include_frequencies_pdf = False,
    pdf = 'histogtrams.pdf', nthreads = {os.cpu_count() - 1})
        Exports cavities and characterizations to PDB-formatted files, writes
        file paths and characterizations to a TOML-formatted results file,
        and optionally plot histograms of frequencies in a PDF file.
    """

    def __init__(
        self,
        cavities: numpy.ndarray,
        surface: numpy.ndarray,
        depths: Union[numpy.ndarray, None],
        scales: Union[numpy.ndarray, None],
        volume: Dict[str, float],
        area: Dict[str, float],
        max_depth: Union[Dict[str, float], None],
        avg_depth: Union[Dict[str, float], None],
        avg_hydropathy: Union[Dict[str, float], None],
        residues: Dict[str, List[List[str]]],
        frequencies: Union[Dict[str, Dict[str, Dict[str, int]]], None],
        _vertices: numpy.ndarray,
        _step: Union[float, int],
        _ncav: int,
        _pdb: Union[str, pathlib.Path, None] = None,
        _ligand: Union[str, pathlib.Path, None] = None,
    ):
        """Constructs attributes for pyKVFinderResults object.

        Parameters
        ----------
        cavities : numpy.ndarray
            Cavity points in the 3D grid (cavities[nx][ny][nz]).
        surface : numpy.ndarray
            Surface points in the 3D grid (surface[nx][ny][nz])
        depths : Union[numpy.ndarray, None]
            A numpy.ndarray with depth of cavity points (depth[nx][ny][nz]).
        scales : Union[numpy.ndarray, None]
            A numpy.ndarray with hydrophobicity scale value mapped at surface
            points (scales[nx][ny][nz]).
        volume : Dict[str, float]
            A dictionary with volume of each detected cavity.
        area : Dict[str, float]
            A dictionary with area of each detected cavity.
        max_depth : Union[Dict[str, float], None]
            A dictionary with maximum depth of each detected cavity.
        avg_depth : Union[Dict[str, float], None]
            A dictionary with average depth of each detected cavity.
        avg_hydropathy : Union[Dict[str, float], None]
            A dictionary with average hydropathy for each detected cavity and
            the range of the hydrophobicity scale (min, max).
        residues: Dict[str, List[List[str]]]
            A dictionary with a list of interface residues for each detected
            cavity.
        frequencies : Union[Dict[str, Dict[str, Dict[str, int]]], None]
            A dictionary with frequencies of residues and class for
            residues of each detected cavity.
        _vertices : numpy.ndarray
            A numpy.ndarray or a list with xyz vertices coordinates (origin,
            X-axis, Y-axis, Z-axis).
        _step : float
            Grid spacing (A).
        _ncav : int
            Number of cavities.
        _pdb : Union[str, pathlib.Path, None], optional
            A path to input PDB file, by default None.
        _ligand : Union[str, pathlib.Path, None], optional
            A path to ligand PDB file, by default None.
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
        output: Union[str, pathlib.Path] = "cavity.pdb",
        output_hydropathy: Union[str, pathlib.Path] = "hydropathy.pdb",
        nthreads: int = os.cpu_count() - 1,
    ) -> None:
        f"""Exports cavitiy (H) and surface (HA) points to PDB-formatted file
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
            Number of threads, by default {os.cpu_count() - 1}.

        Returns
        -------
        None

        Note
        ----
        The cavity nomenclature is based on the integer label. The cavity
        marked with 2, the first integer corresponding to a cavity, is KAA,
        the cavity marked with 3 is KAB, the cavity marked with 4 is KAC
        and so on.
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
        fn: Union[str, pathlib.Path] = "results.toml",
        output: Union[str, pathlib.Path, None] = None,
        output_hydropathy: Union[str, pathlib.Path, None] = None,
    ) -> None:
        """
        Writes file paths and cavity characterization to TOML-formatted file

        Parameters
        ----------
        fn : Union[str, pathlib.Path], optional
            A path to TOML-formatted file for writing file paths and cavity
            characterization (volume, area, depth and interface residues)
            per cavity detected, by default `results.toml`.
        output : Union[str, pathlib.Path, None], optional
            A path to a cavity PDB file, by default None.
        output_hydropathy : Union[str, pathlib.Path, None]
            A path to PDB file for writing hydropathy at surface points, by
            default None.

        Returns
        -------
        None

        Note
        ----
        The cavity nomenclature is based on the integer label. The cavity
        marked with 2, the first integer corresponding to a cavity, is KAA,
        the cavity marked with 3 is KAB, the cavity marked with 4 is KAC
        and so on.
        """
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

    def plot_frequencies(self, pdf: Union[str, pathlib.Path] = "histograms.pdf"):
        """Plot histograms of frequencies (residues and classes of residues) in
        a PDF file.

        Parameters
        ----------
        pdf : Union[str, pathlib.Path], optional
            A path to a PDF file, by default `histograms.pdf`.

        Returns
        -------
        None

        Note
        ----
        The cavity nomenclature is based on the integer label. The cavity
        marked with 2, the first integer corresponding to a cavity, is KAA,
        the cavity marked with 3 is KAB, the cavity marked with 4 is KAC
        and so on.

        Classes
        -------
        Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine,
        Methionine, Valine.
        Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine.
        Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline,
        Serine, Threonine.
        Negatively charged (R4): Aspartate, Glutamate.
        Positively charged (R5): Arginine, Histidine, Lysine.
        Non-standard (RX): Non-standard residues.
        """
        plot_frequencies(self.frequencies, pdf)

    def export_all(
        self,
        fn: Union[str, pathlib.Path] = "results.toml",
        output: Union[str, pathlib.Path] = "cavity.pdb",
        output_hydropathy: Union[str, pathlib.Path] = "hydropathy.pdb",
        include_frequencies_pdf: bool = False,
        pdf: Union[str, pathlib.Path] = "histograms.pdf",
        nthreads: int = os.cpu_count() - 1,
    ) -> None:
        f"""Exports cavities and characterization to PDB-formatted files,
        writes file paths and characterization to a TOML-formatted file, and
        optionally plot histograms of frequencies (residues and classes of
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
            A path to a PDF file, by default `histograms.pdf`.
        nthreads : int, optional
            Number of threads, by default {os.cpu_count() - 1}.

        Returns
        -------
        None

        Note
        ----
        The cavity nomenclature is based on the integer label. The cavity
        marked with 2, the first integer corresponding to a cavity, is KAA,
        the cavity marked with 3 is KAB, the cavity marked with 4 is KAC
        and so on.

        Classes
        -------
        Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine,
        Methionine, Valine.
        Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine.
        Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline,
        Serine, Threonine.
        Negatively charged (R4): Aspartate, Glutamate.
        Positively charged (R5): Arginine, Histidine, Lysine.
        Non-standard (RX): Non-standard residues.
        """
        # Export cavity PDB file
        self.export(output, output_hydropathy, nthreads)
        # Write KVFinder results TOML
        self.write(fn, output, output_hydropathy)
        # Plot histograms of frequencies
        if include_frequencies_pdf:
            self.plot_frequencies(pdf)


def pyKVFinder(
    pdb: Union[str, pathlib.Path],
    ligand: Union[str, pathlib.Path, None] = None,
    dictionary: Union[str, pathlib.Path] = _vdw,
    box: Union[str, pathlib.Path, None] = None,
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
    nthreads: int = os.cpu_count() - 1,
    verbose: bool = False,
) -> pyKVFinderResults:
    f"""Detects and characterizes cavities (volume, area, depth [optional],
    hydropathy [optional] and interface residues).

    Parameters
    ----------
    pdb : Union[str, pathlib.Path]
        A path to input PDB file.
    ligand : Union[str, pathlib.Path], optional
        A path to ligand PDB file, by default None.
    dictionary : Union[str, pathlib.Path], optional
        A path to van der Waals radii file, by default {_vdw}.
    box : Union[str, pathlib.Path, None], optional
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
        Number of threads, by default {os.cpu_count()-1}.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    results : pyKVFinderResults
        A class with the following attributes defined:
        * cavities : numpy.ndarray
            Cavity points in the 3D grid (cavities[nx][ny][nz]).
        * surface : numpy.ndarray
            Surface points in the 3D grid (surface[nx][ny][nz])
        * depths : Union[numpy.ndarray, None]
            A numpy.ndarray with depth of cavity points (depth[nx][ny][nz]).
        * scales : Union[numpy.ndarray, None]
            A numpy.ndarray with hydrophobicity scale value mapped at surface
            points (scales[nx][ny][nz]).
        * volume : Dict[str, float]
            A dictionary with volume of each detected cavity.
        * area : Dict[str, float]
            A dictionary with area of each detected cavity.
        * max_depth : Union[Dict[str, float], None]
            A dictionary with maximum depth of each detected cavity.
        * avg_depth : Union[Dict[str, float], None]
            A dictionary with average depth of each detected cavity.
        * avg_hydropathy : Union[Dict[str, float], None]
            A dictionary with average hydropathy for each detected cavity and
            the range of the hydrophobicity scale (min, max).
        * residues: Dict[str, List[List[str]]]
            A dictionary with a list of interface residues for each detected
            cavity.
        * frequencies : Union[Dict[str, Dict[str, Dict[str, int]]], None]
            A dictionary with frequencies of residues and class for
            residues of each detected cavity.
        * _vertices : numpy.ndarray
            A numpy.ndarray or a list with xyz vertices coordinates (origin,
            X-axis, Y-axis, Z-axis).
        * _step : float
            Grid spacing (A).
        * _ncav : int
            Number of cavities.
        * _pdb : Union[str, pathlib.Path]
            A path to input PDB file.
        * _ligand : Union[str, pathlib.Path, None]
            A path to ligand PDB file.

    Notes
    -----
    The cavity nomenclature is based on the integer label. The cavity marked
    with 2, the first integer corresponding to a cavity, is KAA, the cavity
    marked with 3 is KAB, the cavity marked with 4 is KAC and so on.

    Cavities array has integer labels in each position, that are:
        * -1: bulk points;
        * 0: biomolecule points;
        * 1: empty space points;
        * >=2: cavity points.

    Surface array has integer labels in each position, that are:
        * -1: bulk points;
        * 0: biomolecule or empty space points;
        * >=2: cavity points.

    Classes
    -------
    Aliphatic apolar (R1): Alanine, Glycine, Isoleucine, Leucine,
    Methionine, Valine.
    Aromatic (R2): Phenylalanine, Tryptophan, Tyrosine.
    Polar Uncharged (R3): Asparagine, Cysteine, Glutamine, Proline,
    Serine, Threonine.
    Negatively charged (R4): Aspartate, Glutamate.
    Positively charged (R5): Arginine, Histidine, Lysine.
    Non-standard (RX): Non-standard residues.
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
