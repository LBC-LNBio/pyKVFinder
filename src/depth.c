#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

# define min(x,y) ( ((x) < (y)) ? (x) : (y))
# define max(x,y) ( ((x) > (y)) ? (x) : (y))


/*
 * Struct: POINTS
 * --------------
 * 
 * Two 3D grid points with xyz coordinates (P1 and P2)
 * 
 * X1: x coordinate of P1
 * Y1: y coordinate of P1
 * Z1: z coordinate of P1
 * X2: x coordinate of P2
 * Y2: y coordinate of P2
 * Z2: z coordinate of P2
 *  
 */
typedef struct POINTS 
{
    double X1;
    double Y1;
    double Z1;
    double X2;
    double Y2;
    double Z2;
} pts;

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
filter_frontier (int *cavities, int nx, int ny, int nz, pts *cavs, pts *frontiers, int nthreads)
{
    int i, j, k, tag;

    // Set number of threads in OpenMP
    omp_set_num_threads(nthreads);
    omp_set_nested(1);

    #pragma omp parallel default(none), shared(cavities, nx, ny, nz, cavs, frontiers), private(i, j, k, tag)
    {
        #pragma omp for collapse(3) schedule(static)
        for (i=0; i<nx; i++)
            for (j=0; j<ny; j++)
                for (k=0; k<nz; k++)
                    if (cavities[k + nz * (j + ( ny * i ) )] > 1)
                    {
                        // Get cavity identifier
                        tag = cavities[k + nz * (j + ( ny * i ) )] - 2;
                        
                        // Get min and max coordinates of each cavity
                        cavs[tag].X1 = min(cavs[tag].X1, i);
                        cavs[tag].Y1 = min(cavs[tag].Y1, i);
                        cavs[tag].Z1 = min(cavs[tag].Z1, i);
                        cavs[tag].X2 = max(cavs[tag].X2, i);
                        cavs[tag].Y2 = max(cavs[tag].Y2, i);
                        cavs[tag].Z2 = max(cavs[tag].Z2, i);

                        // Define cavity-bulk frontier points
                        cavities[k + nz * (j + ( ny * i ) )] = define_frontier_points (cavities, nx, ny, nz, i, j, k);
                        
                        // Get min and max coordinates of each cavity-bulk frontier
                        if (cavities[k + nz * (j + ( ny * i ) )] < -1)
                        {
                            frontiers[tag].X1 = min(frontiers[tag].X1, i);
                            frontiers[tag].Y1 = min(frontiers[tag].Y1, i);
                            frontiers[tag].Z1 = min(frontiers[tag].Z1, i);
                            frontiers[tag].X2 = max(frontiers[tag].X2, i);
                            frontiers[tag].Y2 = max(frontiers[tag].Y2, i);
                            frontiers[tag].Z2 = max(frontiers[tag].Z2, i);
                        }
                    }
    }
}

/*
 * Function: _depth
 * ----------------
 * 
 * Description
 * 
 * cavities: 
 * nx:
 * ny:
 * nz:
 * depth:
 * size:
 * max_depth:
 * n:
 * avg_depth:
 * nn:
 * step:
 * nthreads:
 * verbose:
 * 
 * returns:
 * 
 */
void 
_depth (int *cavities, int nx, int ny, int nz, double *depth, int size, double *max_depth, int n, double *avg_depth, int nn, double step, int nthreads, int verbose)
{
    int i, j, k, i2, j2, k2, count, tag;
    double distance, tmp;
    pts *cavs, *frontiers;

    // Set number of threads in OpenMP
    omp_set_num_threads(nthreads);
    omp_set_nested(1);

    // Fill depth 3D grid with 0.0
    dgrid(depth, size);

    // Allocate memory
    cavs = (pts*) calloc (n, sizeof(pts));
    frontiers = (pts*) calloc (n, sizeof(pts));

    // Fill frontiers and cavs points
    for (i=0; i<n; i++)
    {
        frontiers[i].X1 = cavs[i].X1 = nx; 
        frontiers[i].Y1 = cavs[i].Y1 = ny;
        frontiers[i].Z1 = cavs[i].Z1 = nz;
        frontiers[i].X2 = cavs[i].X2 = 0.0;
        frontiers[i].Y2 = cavs[i].Y2 = 0.0;
        frontiers[i].Z2 = cavs[i].Z2 = 0.0;
    }

    // Define cavity-bulk frontier points
    filter_frontier(cavities, nx, ny, nz, cavs, frontiers, nthreads);

    #pragma omp parallel default(none), shared(cavities, depth, max_depth, avg_depth, cavs, frontiers, n, nx, ny, nz, step), private(tmp, tag, i, j, k, i2, j2, k2, distance, count)
    {
        #pragma omp for schedule(dynamic)
        for (tag=0; tag < n; tag++)
        {
            max_depth[tag] = 0.0;
            avg_depth[tag] = 0.0;
            count = 0;
            
            for (i=cavs[tag].X1; i<=cavs[tag].X2; i++)
                for (j=cavs[tag].Y1; j<=cavs[tag].Y2; j++)
                    for (k=cavs[tag].Z1; i<=cavs[tag].Z2; k++)
                        if (abs(cavities[k + nz * (j + ( ny * i ) )]) + 2 == tag)
                        {   
                            tmp = sqrt( pow(nx, 2) + pow(ny, 2) + pow(nz, 2));
                            count++;

                            if (frontiers[tag].X1 == nx && frontiers[tag].Y1 == ny && frontiers[tag].Z1 == nz && frontiers[tag].X2 == 0 && frontiers[tag].Y2 == 0 && frontiers[tag].Z2 == 0)
                            {
                                tmp = 0.0;
                            }
                            else
                            {
                                for (i2=frontiers[tag].X1; i2<=frontiers[tag].X2; i2++)
                                    for (j2=frontiers[tag].Y1; j2<=frontiers[tag].Y2; j2++)
                                        for (k2=frontiers[tag].Z1; i2<=frontiers[tag].Z2; k2++)
                                            if (cavities[k2 + nz * (j2 + ( ny * i2 ) )] == -(tag + 2))
                                            {
                                                distance = sqrt( pow(i2-i, 2) + pow(j2-j, 2) + pow(k2-k, 2) );
                                                if (distance < tmp)
                                                    tmp = distance;
                                            }
                            }
                            // Save depth for cavity point
                            depth[k + nz * (j + ( ny * i ) )] = tmp;
                            if (max_depth[tag] < tmp)
                                max_depth[tag] = tmp;
                            avg_depth[tag] += tmp;
                        }
            avg_depth[tag] /= count;
        }
    }
}
