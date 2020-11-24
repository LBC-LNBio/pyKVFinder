import os
import time
import toml
import logging
from datetime import datetime
from .argparser import argparser
from .utils import read_vdw, read_pdb, write_results, _write_parameters
from .grid import calculate_vertices, calculate_dimensions, calculate_sincos, prepare_box_vertices, detect, spatial, constitutional, export


def run():
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
        args.vertices, pdb, xyzr, args.sincos, nx, ny, nz = prepare_box_vertices(args.box, pdb, xyzr, args.probe_in, args.probe_out, args.step, args.nthreads)

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
