import os
import sys
import time

import numpy as np

from argparser import argparser
from modules.utils import read_vdw_dat, read_pdb, write_results
from _gridprocessing import detect, characterize, export

def run(args):
    
    if args.verbose:
        print (f"[PID {os.getpid()}] Running pyKVFinder for: {args.pdb}")

    if args.verbose:
        print("> Loading atomic dictionary file")
    vdw = read_vdw_dat(args.dictionary)

    if args.verbose:
        print("> Reading PDB coordinates")
    pdb, xyzr = read_pdb(args.pdb, vdw)

    if args.verbose: 
        print("> Calculating 3D grid dimensions")
    P1 = np.min(xyzr[:, 0:3], axis=0) - args.probe_out 
    xmax, ymax, zmax = np.max(xyzr[:, 0:3], axis=0) + args.probe_out 
    P2 = np.array([xmax, P1[1], P1[2]])
    P3 = np.array([P1[0], ymax, P1[2]])
    P4 = np.array([P1[0], P1[1], zmax])

    # Calculate distance between points
    norm1 = np.linalg.norm(P2-P1)
    norm2 = np.linalg.norm(P3-P1)
    norm3 = np.linalg.norm(P4-P1)

    # Calculate grid dimensions
    nx = int ( norm1 / args.step ) + 1
    ny = int ( norm2 / args.step ) + 1
    nz = int ( norm3 / args.step ) + 1
    if args.verbose:
        print(f"nx: {nz}\tny: {ny}\tnz: {nz}")

    # Calculate sin and cos of angles a and b
    sincos = np.array([
        ( P4[1] - P1[1] ) / norm3, # sin a
        ( P3[1] - P1[1] ) / norm2, # cos a
        ( P2[2] - P1[2] ) / norm1, # sin b
        ( P2[0] - P1[0] ) / norm1  # cos b
    ])
    if args.verbose:
        print(f"sina: {sincos[0]}\tsinb: {sincos[2]}")
        print(f"cosa: {sincos[1]}\tcosb: {sincos[3]}")

    if args.surface == 'SES':
        args.surface = True
        if args.verbose:
            print ("> Surface representation: Solvent Excluded Surface (SES)") 
    else:
        args.surface = False
        if args.verbose:
            print ("> Surface representation: Solvent Accessible Surface (SAS)")

    nvoxels = nx * ny * nz
    cavities = detect(nvoxels, nx, ny, nz, xyzr, P1, sincos, args.step, args.probe_in, args.probe_out, args.removal_distance, args.volume_cutoff, args.surface, 15, args.verbose).reshape(nx, ny, nz)
    
    # Number of cavities 
    # FIXME: check if it better to return this value from detect
    ncav = int(np.sum((np.unique(cavities) > 1)))

    # No cavities were found
    if (not ncav):
        return True

    # Backbone
    # FIXME: need to be tested
    backbone = False
    if backbone:
        mask = np.where(pdb[:,2] != 'C') and np.where(pdb[:,2] != 'CA') and np.where(pdb[:,2] != 'N') and np.where(pdb[:,2] != 'O')
        pdb = pdb[mask[0], 0]
        xyzr = xyzr[mask[0]]

    # Characterization
    residues, surface, volume, area = characterize(cavities, nvoxels, ncav, ncav, xyzr, P1, sincos, args.step, args.probe_in, args.probe_out, 15, args.verbose)
    surface = surface.reshape(nx, ny, nz)
    print(type(residues), type(surface), type(volume), type(area))
    print(residues)

    # Export cavities
    export("tests/cavity.pdb", cavities, surface, P1, sincos, args.step, ncav, 15)

    # Export results
    # FIXME: bad formatting
    write_results("tests/results.toml", volume, area)

    return True

if __name__ == "__main__":
    # Start time
    start_time = time.time()

    # Load pyKVFinder argument parser
    parser = argparser()

    # Parse command-line arguments
    args = parser.parse_args()
    print(args)

    # Run pyKVFinder
    run(args)

    # Elapsed time
    elapsed_time = time.time() - start_time
    print(f"[ \033[1mElapsed time:\033[0m {elapsed_time:.6f} ]")

