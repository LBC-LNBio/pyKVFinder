import os
import time
import logging
import numpy as np
from datetime import datetime
from .argparser import argparser
from .utils import read_vdw, read_pdb, write_results, _write_parameters
from .grid import calculate_vertices, calculate_dimensions, calculate_sincos, get_vertices, detect, spatial, constitutional, export

__all__ = ['pyKVFinder']

_here = os.path.abspath(os.path.dirname(__file__))
_dictionary = os.path.join(_here, "data/vdw.dat")


def cli():
    # Start time
    start_time = time.time()

    # Load pyKVFinder argument parser
    parser = argparser()

    # Parse command-line arguments
    args = parser.parse_args()

    # Get base name from pdb file if not defined by user
    if not args.base_name:
        args.base_name = os.path.basename(args.pdb.replace('.pdb', ''))

    # Create output directory
    os.makedirs(args.output_directory, exist_ok=True)

    # Print message to stdout
    print(f"[PID {os.getpid()}] Running pyKVFinder for: {args.pdb}")

    # Start logging
    logging.basicConfig(filename=f"{os.path.join(args.output_directory, 'KVFinder.log')}", level=logging.INFO, format='%(message)s')
    logging.info("=" * 80)
    logging.info(f"Date: {datetime.now().strftime('%a %d %B, %Y')}\nTime: {datetime.now().strftime('%H:%M:%S')}\n")
    logging.info(f"[ Running pyKVFinder for: {args.pdb} ]")
    logging.info(f"> vdW radii file: {args.dictionary}")

    if args.verbose:
        print("> Loading atomic dictionary file")
    vdw = read_vdw(args.dictionary)

    if args.verbose:
        print("> Reading PDB coordinates")
    pdb, xyzr = read_pdb(args.pdb, vdw)

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
        args.vertices, pdb, xyzr, args.sincos, nx, ny, nz = get_vertices(args.box, pdb, xyzr, args.probe_in, args.probe_out, args.step, args.nthreads)

        # Set flag to boolean
        args.box = True
    else:
        # Get vertices from pdb
        args.vertices = calculate_vertices(xyzr, args.probe_out, args.step)

        # Calculate distance between points
        nx, ny, nz = calculate_dimensions(args.vertices, args.step)
        if args.verbose:
            print(f"Dimensions: (nx:{nx}, ny:{ny}, nz:{nz})")

        # Calculate sin and cos of angles a and b
        args.sincos = calculate_sincos(args.vertices)
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

    if args.surface == 'SES':
        args.surface = True
        if args.verbose:
            print("> Surface representation: Solvent Excluded Surface (SES)")
    else:
        args.surface = False
        if args.verbose:
            print("> Surface representation: Solvent Accessible Surface (SAS)")

    # Cavity detection
    ncav, cavities = detect(nx, ny, nz, xyzr, args.vertices, args.sincos, args.step, args.probe_in, args.probe_out, args.removal_distance, args.volume_cutoff, lxyzr, args.ligand_cutoff, args.box, args.surface, args.nthreads, args.verbose)

    # Cavities were found
    if ncav > 0:
        # Spatial characterization
        surface, volume, area = spatial(cavities, nx, ny, nz, ncav, args.step, args.nthreads, args.verbose)

        # Constitutional characterization
        residues = constitutional(cavities, pdb, xyzr, args.vertices, args.sincos, ncav, args.step, args.probe_in, args.ignore_backbone, args.nthreads, args.verbose)

        # Export cavities
        output_cavity = os.path.join(args.output_directory, f"{args.base_name}.KVFinder.output.pdb")
        export(output_cavity, cavities, surface, args.vertices, args.sincos, ncav, args.step, args.nthreads)

        # Write results
        output_results = os.path.join(args.output_directory, f"{args.base_name}.KVFinder.results.toml")
        write_results(output_results, args.pdb, args.ligand, output_cavity, volume, area, residues, args.step)

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

    def __init__(self, cavities: np.ndarray, surface: np.ndarray, volume: dict, area: dict, residues: dict, vertices: np.ndarray, step: float, ncav: int):
        self.cavities = cavities
        self.surface = surface
        self.volume = volume
        self.area = area
        self.residues = residues
        self._vertices = vertices
        self._step = step
        self._ncav = ncav

    def __repr__(self):
        return '<pyKVFinderResults class>'

    def export(self, fn: str = 'cavity.pdb', nthreads: int = os.cpu_count() - 1):
        sincos = calculate_sincos(self.vertices)
        export(fn, self.cavities, self.surface, self._vertices, sincos, self._ncav, self._step, nthreads)


def pyKVFinder(pdb: str, ligand: str = None, dictionary: str = _dictionary, box: dict = None, step: float = 0.6, probe_in: float = 1.4, probe_out: float = 4.0, removal_distance: float = 2.4, volume_cutoff: float = 5.0, ligand_cutoff: float = 5.0, surface: str = 'SES', ignore_backbone: bool = False, nthreads: int = os.cpu_count() - 1, verbose: bool = False) -> pyKVFinderResults:

    if verbose:
        print("> Loading atomic dictionary file")
    vdw = read_vdw(dictionary)

    if verbose:
        print("> Reading PDB coordinates")
    pdb, xyzr = read_pdb(pdb, vdw)

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
        vertices, pdb, xyzr, sincos, nx, ny, nz = get_vertices(box, pdb, xyzr, probe_in, probe_out, step, nthreads)

        # Set flag to boolean
        box = True
    else:
        # Get vertices from pdb
        vertices = calculate_vertices(xyzr, probe_out, step)
        # Calculate distance between points
        nx, ny, nz = calculate_dimensions(vertices, step)
        if verbose:
            print(f"Dimensions: (nx:{nx}, ny:{ny}, nz:{nz})")

        # Calculate sin and cos of angles a and b
        sincos = calculate_sincos(vertices)
        if verbose:
            print(f"sina: {sincos[0]:.2f}\tsinb: {sincos[2]:.2f}")
            print(f"cosa: {sincos[1]:.2f}\tcosb: {sincos[3]:.2f}")

        # Set flag to boolean
        box = False

    if surface == 'SES':
        surface = True
        if verbose:
            print("> Surface representation: Solvent Excluded Surface (SES)")
    elif surface == 'SAS':
        surface = False
        if verbose:
            print("> Surface representation: Solvent Accessible Surface (SAS)")
    else:
        raise Exception(f'surface variable should be SAS or SES, not {surface}')

    # Cavity detection
    ncav, cavities = detect(nx, ny, nz, xyzr, vertices, sincos, step, probe_in, probe_out, removal_distance, volume_cutoff, lxyzr, ligand_cutoff, box, surface, nthreads, verbose)

    if ncav > 0:
        # Spatial characterization
        surface, volume, area = spatial(cavities, nx, ny, nz, ncav, step, nthreads, verbose)

        # Constitutional characterization
        residues = constitutional(cavities, pdb, xyzr, vertices, sincos, ncav, step, probe_in, ignore_backbone, nthreads, verbose)
    else:
        volume, area, residues = None, None, None

    # Return dict
    return pyKVFinderResults(cavities, surface, volume, area, residues, vertices, step, ncav)
