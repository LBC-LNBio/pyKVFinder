#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void igrid (int *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = 1;              

}

void fgrid (float *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = 0.0;
               
}

void dgrid (double *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = 0.0;
               
}

void cgrid (int *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = '\0';
               
}

void filter (int *grid, int dx, int dy, int dz)
{
    int i, j, k;

    for (i=0; i<dx; i++)
        for (j=0; j<dy; j++)
            for (k=0; k<dz; k++)
                printf("(%d, %d, %d): %d\n", i, j, k, grid[k + dz * (j + ( dy * i ) ) ]);

}


/******* sincos ******
* sincos[0] = sin a  *
* sincos[1] = cos a  *
* sincos[2] = sin b  *
* sincos[3] = cos b  *
*********************/

void fill_grid (int *grid, int dx, int dy, int dz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int ncores)
{
    int i, j, k, imax, jmax, kmax, atom;
    double distance, H, x, y, z, xaux, yaux, zaux;


    /* Set number of processes in OpenMP */
	omp_set_num_threads (ncores);
    omp_set_nested (1);

    #pragma omp parallel default(none), shared(grid, reference, step, probe, natoms, dx, dy, dz, sincos, atoms, ncores), private(atom, i, j, k, imax, jmax, kmax, distance, H, x, y, z, xaux, yaux, zaux)
    {
    #pragma omp for schedule(dynamic) nowait
        for (atom=0; atom<natoms; atom++) {

            // Convert atom coordinates in 3D grid coordinates
            x = atoms[atom * 4] / step; 
            y = atoms[1 + (atom * 4)] / step; 
            z = atoms[2 + (atom * 4)] / step;

            xaux = x * sincos[3] + z * sincos[2];
            yaux = y;
            zaux = -x * sincos[2] + z * sincos[3];

            x = xaux;
            y = yaux * sincos[1] - zaux * sincos[0];
            z = yaux * sincos[0] + zaux * sincos[1];

            // Create a radius (H) for space occupied by probe and atom
            H = ( probe + atoms[4 + (atom * 4)] ) / step;

            imax = ceil(x + H)+1;
            jmax = ceil(y + H)+1;
            kmax = ceil(z + H)+1;

            // #pragma omp simd
            // #pragma preomp for collapse(3) schedule(dynamic) ordered nowait
                // Loop around radius from atom center
                for (i=floor(x - H); i<imax; i++)
                    for (j=floor(y - H); j<jmax; j++)
                        for (k=floor(z - H); k<kmax; k++) {

                            // Get distance between atom center and point inspected
                            distance = sqrt( pow(i - x, 2) + pow(j - y, 2) + pow(k - z, 2));
                            if (distance < H)
                                if (i >= 0 && i < dx && j >= 0 && j < dy && k >= 0 && i < dz)
                                    grid[ k + dz * (j + ( dy * i ) ) ] = 0;

                        }
        }
    }
}