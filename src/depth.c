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
 * Function: define_boundary_points
 * --------------------------------
 * 
 * Identify cavity-bulk boundary points based on neighboring points
 * 
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * i: x coordinate of cavity point
 * j: y coordinate of cavity point
 * k: z coordinate of cavity point
 * 
 * returns: cavity identifier (tag) or cavity-bulk boundary identifier (-tag)
 */
int 
define_boundary_points (int *cavities, int nx, int ny, int nz, int i, int j, int k)
{
    if (i-1>=0)
        if (cavities[ k + nz * (j + ( ny * (i - 1) ) ) ] == -1)
            return -(cavities[k + nz * (j + ( ny * i ) )]);
    if (i+1<nx)
        if (cavities[ k + nz * (j + ( ny * (i + 1) ) ) ] == -1)
            return -(cavities[k + nz * (j + ( ny * i ) )]);
    if (j-1>=0)
        if (cavities[ k + nz * ( (j - 1) + ( ny * i ) ) ] == -1)
            return -(cavities[k + nz * (j + ( ny * i ) )]);
    if (j+1<ny)
        if (cavities[ k + nz * ( (j + 1) + ( ny * i ) ) ] == -1)
            return -(cavities[k + nz * (j + ( ny * i ) )]);
    if (k-1>=0)
        if (cavities[ (k - 1) + nz * (j + ( ny * i ) ) ] == -1)
            return -(cavities[k + nz * (j + ( ny * i ) )]);
    if (k+1<nz)
        if (cavities[ (k + 1) + nz * (j + ( ny * i ) ) ] == -1)
            return -(cavities[k + nz * (j + ( ny * i ) )]);

	return cavities[k + nz * (j + ( ny * i ) )];
} 

/*
 * Function: filter_boundary
 * -------------------------
 * 
 * Inspect cavities 3D grid and mark detected cavity-bulk boundary points
 * 
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * nthreads: number of threads for OpenMP
 * 
 */
void 
filter_boundary (int *cavities, int nx, int ny, int nz, pts *cavs, pts *boundaries, int nthreads)
{
    int i, j, k, tag;

    // Set number of threads in OpenMP
    omp_set_num_threads(nthreads);
    omp_set_nested(1);

    #pragma omp parallel default(none), shared(cavities, nx, ny, nz, cavs, boundaries), private(i, j, k, tag)
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
                        cavs[tag].Y1 = min(cavs[tag].Y1, j);
                        cavs[tag].Z1 = min(cavs[tag].Z1, k);
                        cavs[tag].X2 = max(cavs[tag].X2, i);
                        cavs[tag].Y2 = max(cavs[tag].Y2, j);
                        cavs[tag].Z2 = max(cavs[tag].Z2, k);

                        // Define cavity-bulk boundary points
                        cavities[k + nz * (j + ( ny * i ) )] = define_boundary_points (cavities, nx, ny, nz, i, j, k);
                        
                        // Get min and max coordinates of each cavity-bulk boundary
                        if (cavities[k + nz * (j + ( ny * i ) )] < -1)
                        {
                            boundaries[tag].X1 = min(boundaries[tag].X1, i);
                            boundaries[tag].Y1 = min(boundaries[tag].Y1, j);
                            boundaries[tag].Z1 = min(boundaries[tag].Z1, k);
                            boundaries[tag].X2 = max(boundaries[tag].X2, i);
                            boundaries[tag].Y2 = max(boundaries[tag].Y2, j);
                            boundaries[tag].Z2 = max(boundaries[tag].Z2, k);
                        }
                    }
    }
}

/*
 * Function: define_depth
 * ----------------------
 * 
 * Description
 * 
 * cavities: 
 * depth:
 * nx:
 * ny:
 * nz:
 * max_depth:
 * avg_depth:
 * ncav:
 * step:
 * cavs:
 * boundaries:
 * nthreads:
 * 
 */
void
define_depth(int *cavities, double *depth, int nx, int ny, int nz, double *max_depth, double *avg_depth, int ncav, pts *cavs, pts *boundaries, int nthreads)
{
    int i, j, k, i2, j2, k2, count, tag;
    double distance, tmp;

    // Set number of threads in OpenMP
    omp_set_num_threads(nthreads);
    omp_set_nested(1);

    #pragma omp parallel default(none), shared(cavities, depth, max_depth, avg_depth, cavs, boundaries, ncav, nx, ny, nz), private(tmp, tag, i, j, k, i2, j2, k2, distance, count)
    {
        #pragma omp for schedule(dynamic)
            for (tag=0; tag < ncav; tag++)
            {
                max_depth[tag] = 0.0;
                avg_depth[tag] = 0.0;
                count = 0;
                
                for (i=cavs[tag].X1; i<=cavs[tag].X2; i++)
                    for (j=cavs[tag].Y1; j<=cavs[tag].Y2; j++)
                        for (k=cavs[tag].Z1; k<=cavs[tag].Z2; k++)
                            if ( abs(cavities[k + nz * (j + ( ny * i ) )]) == (tag + 2) )
                            {   
                                tmp = sqrt( pow(nx, 2) + pow(ny, 2) + pow(nz, 2));
                                count++;

                                if (boundaries[tag].X1 == nx && boundaries[tag].Y1 == ny && boundaries[tag].Z1 == nz && boundaries[tag].X2 == 0.0 && boundaries[tag].Y2 == 0.0 && boundaries[tag].Z2 == 0.0)
                                {
                                    tmp = 0.0;
                                }
                                else
                                {
                                    for (i2=boundaries[tag].X1; i2<=boundaries[tag].X2; i2++)
                                        for (j2=boundaries[tag].Y1; j2<=boundaries[tag].Y2; j2++)
                                            for (k2=boundaries[tag].Z1; k2<=boundaries[tag].Z2; k2++)
                                                if ( cavities[k2 + nz * (j2 + ( ny * i2 ) )] == -(tag + 2) )
                                                {
                                                    distance = sqrt( pow(i2-i, 2) + pow(j2-j, 2) + pow(k2-k, 2) );
                                                    if (distance < tmp)
                                                        tmp = distance;
                                                }
                                }
                                // Save depth for cavity point
                                depth[k + nz * (j + ( ny * i ) )] = tmp;
                                if (tmp > max_depth[tag])
                                    max_depth[tag] = tmp;
                                avg_depth[tag] += tmp;
                            }
                avg_depth[tag] /= count;
            }
            
    }

}

void
remove_boundary (int *cavities, int nx, int ny, int nz, int ncav, pts *boundaries, int nthreads)
{
    int i, j, k, tag;

    // Set number of threads in OpenMP
    omp_set_num_threads(nthreads);
    omp_set_nested(1);

    #pragma omp parallel default(none), shared(cavities, boundaries, ncav, nx, ny, nz), private(tag, i, j, k)
        #pragma omp for schedule(dynamic)
            for (tag=0; tag < ncav; tag++)
                for (i=boundaries[tag].X1; i<=boundaries[tag].X2; i++)
                    for (j=boundaries[tag].Y1; j<=boundaries[tag].Y2; j++)
                        for (k=boundaries[tag].Z1; k<=boundaries[tag].Z2; k++)
                            if ( cavities[k + nz * (j + ( ny * i ) )] < -1 )
                                cavities[k + nz * (j + ( ny * i ) )] = abs(cavities[k + nz * (j + ( ny * i ) )]);
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
 * nn:calculate_depth
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
    int i;
    pts *cavs, *boundaries;

    // Fill depth 3D grid with 0.0
    dgrid (depth, size);

    // Allocate memory
    cavs = (pts*) calloc (n, sizeof(pts));
    boundaries = (pts*) calloc (n, sizeof(pts));

    // Fill boundaries and cavs points
    for (i=0; i<n; i++)
    {
        boundaries[i].X1 = cavs[i].X1 = nx; 
        boundaries[i].Y1 = cavs[i].Y1 = ny;
        boundaries[i].Z1 = cavs[i].Z1 = nz;
        boundaries[i].X2 = cavs[i].X2 = 0.0;
        boundaries[i].Y2 = cavs[i].Y2 = 0.0;
        boundaries[i].Z2 = cavs[i].Z2 = 0.0;
    }

    if (verbose)
        fprintf (stdout, "> Defining bulk-cavity boundary points\n");
    // Define cavity-bulk boundary points
    filter_boundary (cavities, nx, ny, nz, cavs, boundaries, nthreads);

    if (verbose)
        fprintf (stdout, "> Estimating depth\n");
    define_depth (cavities, depth, nx, ny, nz, max_depth, avg_depth, n, cavs, boundaries, nthreads);

    // Multiply depth values by grid spacing (step)
    for (i=0; i<n; i++)
    {
        max_depth[i] *= step;
        avg_depth[i] *= step;
    }

    remove_boundary (cavities, nx, ny, nz, n, boundaries, nthreads);

    // Free pts
    free (cavs);
    free (boundaries);
}
