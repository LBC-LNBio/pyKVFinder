#include <stdio.h>
#include <stdlib.h>

void igrid (int *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = i;              

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

void insert_atom (int *grid, int dx, int dy, int dz, double *atom, int xyz, double *reference, int ndims, int npoints, double step, double probe_in)
{
    return 0;
}