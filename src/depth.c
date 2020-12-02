#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "depth.h"


/*
 * Function: dgrid
 * ---------------
 * 
 * Fill double grid with 0.0
 * 
 * grid: empty 3D grid
 * size: number of voxels
 * 
 */
void 
dgrid (double *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = 0.0;
               
}

/*
 * Function: define_frontier_points
 * --------------------------------
 * 
 * Identify cavity-bulk frontier points based on neighboring points
 * 
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * i: x coordinate of cavity point
 * j: y coordinate of cavity point
 * k: z coordinate of cavity point
 * 
 * returns: cavity identifier (tag) or cavity-bulk frontier identifier (-tag)
 */
int 
define_frontier_points (int *cavities, int nx, int ny, int nz, int i, int j, int k)
{
    if (i-1>=0 && i+1<nx && j-1>=0 && j+1<ny && k-1>=0 && k+1<nz)
    {
        if (cavities[ k + nz * (j + ( ny * (i - 1) ) ) ] == -1)
            return -cavities[k + nz * (j + ( ny * i ) )];
        if (cavities[ k + nz * (j + ( ny * (i + 1) ) ) ] == -1)
            return -cavities[k + nz * (j + ( ny * i ) )];
        if (cavities[ k + nz * ( (j - 1) + ( ny * i ) ) ] == -1)
            return -cavities[k + nz * (j + ( ny * i ) )];
        if (cavities[ k + nz * ( (j + 1) + ( ny * i ) ) ] == -1)
            return -cavities[k + nz * (j + ( ny * i ) )];
        if (cavities[ (k - 1) + nz * (j + ( ny * i ) ) ] == -1)
            return -cavities[k + nz * (j + ( ny * i ) )];
        if (cavities[ (k + 1) + nz * (j + ( ny * i ) ) ] == -1)
            return -cavities[k + nz * (j + ( ny * i ) )];
    }
	return cavities[k + nz * (j + ( ny * i ) )];
} 

/*
 * Function: filter_frontier
 * -------------------------
 * 
 * Inspect cavities 3D grid and mark detected cavity-bulk frontier points
 * 
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * nthreads: number of threads for OpenMP
 * 
 */
void 
filter_frontier (int *cavities, int nx, int ny, int nz, int nthreads)
{

}

void 
_depth (int *cavities, int nx, int ny, int nz, double *depth, int size, double *max_depth, int n, double *avg_depth, int nn, double step, int nthreads, int verbose)
{
    dgrid(depth, size);
}

