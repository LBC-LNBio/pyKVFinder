import os
import sys
import time

import numpy as np

from argparser import argparser
from modules.utils import read_vdw_dat, read_pdb
from _gridprocessing import igrid, fgrid, dgrid, fill, subtract, export

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
    dx = int ( norm1 / args.step ) + 1
    dy = int ( norm2 / args.step ) + 1
    dz = int ( norm3 / args.step ) + 1
    if args.verbose:
        print(f"dx: {dx}\tdy: {dy}\tdz: {dz}")

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

    if args.verbose: 
        print("> Creating 3D grid")
    A = igrid(dx * dy * dz).reshape(dx, dy, dz)
    S = igrid(dx * dy * dz).reshape(dx, dy, dz)

    # if args.verbose:
    #     print ("> Filling 3D grid with Probe In")
    # fill_grid(A, xyzr, P1, sincos, args.step, args.probe_in, 15)
    # if args.surface == 'SES':
    #     if args.verbose:
    #         print ("> Surface representation: Solvent Excluded Surface (SES)") 
    #         ses(A, args.step, args.probe_in, 15) 

    if args.verbose:
        print ("> Filling 3D grid with Probe In, ", end="")
    if args.surface == 'SES':
        args.surface = True
        if args.verbose:
            print ("using Solvent Excluded Surface (SES) representation")
    else:
        args.surface = False
        if args.verbose:
            print ("using Solvent Accessible Surface (SAS) representation")
    fill (A, xyzr, P1, sincos, args.step, args.probe_in, 15, args.surface)

    if args.verbose:
        print ("> Filling 3D grid with Probe Out")
    fill (S, xyzr, P1, sincos, args.step, args.probe_out, 15, False)

    if args.verbose:
        print ("> Defining biomolecular cavities")
        subtract(A, S, args.step, args.removal_distance, 15)

    # print(np.max(A), np.min(A), np.sum(A)*0.6*0.6*0.6)
    export(A, S, args.step, "tests/cavity.pdb", P1, sincos, args.filled)
    # from sklearn.cluster import DBSCAN
    # from sklearn.cluster import AgglomerativeClustering
    # model = DBSCAN(eps=args.step, min_samples=round(args.volume_cutoff/0.6), n_jobs=15)
    # pred = model.fit_predict()
    
    # print(A.shape)

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

