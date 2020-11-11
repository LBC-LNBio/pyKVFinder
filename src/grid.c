#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define DEBUG 1

/******* sincos ******
* sincos[0] = sin a  *
* sincos[1] = cos a  *
* sincos[2] = sin b  *
* sincos[3] = cos b  *
*********************/

void 
detect (
    int *PI, int size, 
    int nx, int ny, int nz,
    double *atoms, int natoms, int xyzr, 
    double *reference, int ndims, 
    double *sincos, int nvalues, 
    double step, 
    double probe_in,
    double probe_out,
    double removal_threshold,
    int is_ses,
    int ncores)
{
    int *PO;

    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {
            igrid(PI, size);
            fill(PI, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_in, ncores);
        }

        #pragma omp section
        {
            PO = (int *) calloc (size, sizeof (int));
            igrid(PO, size);
            fill(PO, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_out, ncores);
        }
    }

    if (is_ses)
        ses(PI, nx, ny, nz, step, probe_in, 15);
    // ses(PO, nx, ny, nz, step, probe_in, 15);

    subtract(PI, PO, nx, ny, nz, step, removal_threshold, 15);

    if (DEBUG)
    {
        export ("tests/cavity.pdb", PI, nx, ny, nz, reference, ndims, sincos, nvalues, step, 1);
    }

    // Free PO
    free(PO);
    
}

void 
fill (int *grid, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int ncores)
{
    int i, j, k, atom;
    double x, y, z, xaux, yaux, zaux, distance, H;

    /* Set number of processes in OpenMP */
	omp_set_num_threads (ncores);
    omp_set_nested (1);

    #pragma omp parallel default(none), shared(grid, reference, step, probe, natoms, nx, ny, nz, sincos, atoms, ncores), private(atom, i, j, k, distance, H, x, y, z, xaux, yaux, zaux)
    {
    #pragma omp for schedule(dynamic) nowait
        for (atom=0; atom<natoms; atom++)
        {
            // Convert atom coordinates in 3D grid coordinates
            x = ( atoms[atom * 4] - reference[0] ) / step; 
            y = ( atoms[1 + (atom * 4)] - reference[1] ) / step; 
            z = ( atoms[2 + (atom * 4)] - reference[2] ) / step;

            xaux = x * sincos[3] + z * sincos[2];
            yaux = y;
            zaux = (-x) * sincos[2] + z * sincos[3];

            x = xaux;
            y = yaux * sincos[1] - zaux * sincos[0];
            z = yaux * sincos[0] + zaux * sincos[1];

            // Create a radius (H) for space occupied by probe and atom
            H = ( probe + atoms[3 + (atom * 4)] ) / step;
        
            // Loop around radius from atom center
            for (i=floor(x - H); i<ceil(x + H); i++)
                for (j=floor(y - H); j<ceil(y + H); j++)
                    for (k=floor(z - H); k<ceil(z + H); k++) {

                        // Get distance between atom center and point inspected
                        distance = sqrt( pow(i - x, 2) + pow(j - y, 2) + pow(k - z, 2));
                        if (distance < H)
                            if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && i < nz)
                                grid[ k + nz * (j + ( ny * i ) ) ] = 0;

                    }
        }
    }
}

/*
 * Function: check_protein_neighbours
 * ----------------------------------
 * 
 * Checks if a cavity point on the grid is next to a protein point (0 or -2)
 * 
 * grid: 3D grid
 * dx: x grid units
 * dy: y grid units
 * dz: z grid units
 * i: x coordinate of cavity point
 * j: y coordinate of cavity point
 * k: z coordinate of cavity point
 * 
 * returns: true (int 1) or false (int 0)
 */
int
check_protein_neighbours (int *grid, int nx, int ny, int nz, int i, int j, int k)
{
	int x, y, z;

	// Loop around neighboring points
	for (x=i-1; x<=i+1; x++)
        for (y=j-1; y<=j+1; y++)
			for (z=k-1; z<=k+1; z++) 
            {
			    // Check if point is inside 3D grid
				if (x >= 0 || y >= 0 || z >= 0 || x < nx || y < ny || z < nz)
                {
				    if (grid[ z + nz * (y + ( ny * x ) ) ] == 0 || grid[ z + nz * (y + ( ny * x ) ) ] == -2)
                        return 1;
                }
			}
	return 0;
}

void 
ses (int *grid, int nx, int ny, int nz, double step, double probe_in, int ncores)
{
    int i, j, k, i2, j2, k2, aux;
    double distance;    
    
    // Set number of processes in OpenMP
	omp_set_num_threads (ncores);
    omp_set_nested (1);
    
    // Calculate sas limit in 3D grid units
    aux = (int) (probe_in / step) + 1;

    #pragma omp parallel default(none), shared(grid,step,probe_in,aux,nx,ny,nz), private(i,j,k,i2,j2,k2,distance)
    {   
        #pragma omp for schedule(dynamic) collapse(3) nowait
        // Loop around 3D grid
        for (i=0; i<nx; i++)
            for (j=0; j<ny; j++)
                for (k=0; k<nz; k++) 
                {
                    // Check if a cavity point
                    if (grid[ k + nz * (j + ( ny * i ) ) ] == 1) 
                        if ( check_protein_neighbours(grid, nx, ny, nz, i, j, k) )
                        {
                            // Loop around sas limit from cavity point next to protein point
                            for (i2=i-aux; i2<=i+aux; i2++)
                                for (j2=j-aux; j2<=j+aux; j2++)
                                    for (k2=k-aux; k2<=k+aux; k2++)
                                    { 
                                        if (i2>0 && j2>0 && k2>0 && i2<nx-1 && j2<ny-1 && k2<nz-1)
                                        {
                                            // Get distance between point inspected and cavity point
                                            distance = sqrt ( pow(i - i2, 2) + pow(j - j2, 2) + pow(k - k2, 2));
                                            // Check if inspected point is inside sas limit
                                            if ( distance < (probe_in / step) )
                                                if (grid[ k2 + nz * (j2 + ( ny * i2 ) ) ] == 0)
                                                    // Mark cavity point
                                                    grid[ k2 + nz * (j2 + ( ny * i2 ) ) ] = -2;
                                        }
                                    }
                        }
                }

        #pragma omp for collapse(3)
            // Loop around 3D grid
            for (i=0; i<nx; i++)
                for (j=0; j<ny; j++)
                    for (k=0; k<nz; k++) 
                    {
                        // Mark space occupied by sas limit from protein surface
                        if (grid[ k + nz * (j + ( ny * i ) ) ] == -2)
                            grid[ k + nz * (j + ( ny * i ) ) ] = 1;

                    }
    }

}

void
subtract (int *PI, int *PO, int nx, int ny, int nz, double step, double removal_threshold, int ncores)
{
	int i, j, k, i2, j2, k2, rt;
    
    rt = ceil (removal_threshold / step);

    // Set number of processes in OpenMP
	omp_set_num_threads (ncores);
    omp_set_nested (1);

    /* Create a parallel region */
    #pragma omp parallel default(none), shared(PI, PO, nx, ny, nz, i, j, k, step, rt, removal_threshold), private(j2,i2,k2)
    {
        #pragma omp for schedule(dynamic) collapse(3)
            // Loop around the search box
            for (i=0; i<nx; i++)
                for (j=0; j<ny; j++)
                    for (k=0; k<nz; k++) 
                    {
                        // Check if point is a cavity point in grid filled with PO
                        if ( PO[ k + nz * (j + ( ny * i ) ) ] ) 
                        {
                            // Loops around space occupied by probe from atom position
                            for(i2=i-rt; i2<=i+rt; i2++)
                                for(j2=j-rt; j2<=j+rt; j2++)
                                    for(k2=k-rt; k2<=k+rt; k2++)
                                        // Check if inside 3D grid
                                        if(i2>=0 && i2<nx && j2>=0 && j2<ny && k2>=0 && k2<nz)
                                            // Mark points in grid filled with Probe In, where Probe Out reached
                                            PI[ k2 + nz * (j2 + ( ny * i2 ) ) ] = -1;
                        }
                    }
    }
}

void
igrid (int *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = 1;              

}

void
fgrid (float *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = 0.0;
               
}

void 
dgrid (double *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = 0.0;
               
}

void 
cgrid (int *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = '\0';
               
}

void 
filter (int *grid, int dx, int dy, int dz)
{
    int i, j, k;

    for (i=0; i<dx; i++)
        for (j=0; j<dy; j++)
            for (k=0; k<dz; k++)
                printf("(%d, %d, %d): %d\n", i, j, k, grid[k + dz * (j + ( dy * i ) ) ]);

}

void
export (char *fn, int *cavities, int nx, int ny, int nz, double *reference, int ndims, double *sincos, int nvalues, double step, int is_filled)
{
	int i, j, k, count = 1, control = 1, tag = 1;
	double x, y, z, xaux, yaux, zaux;
	FILE *output;

	// Open cavity PDB file
	output = fopen (fn, "w");

	for (i=0; i<nx; i++)
        for (j=0; j<ny; j++)
			for (k=0; k<nz; k++) 
            {
                    // Check if cavity point with value tag
					if ( cavities[k + nz * (j + ( ny * i ) )] == 1 ) 
                    {
						// Convert 3D grid coordinates to real coordinates
						x = i * step; 
                        y = j * step; 
                        z = k * step;
						
                        xaux = (x * sincos[3]) + (y * sincos[0] * sincos[2]) - (z * sincos[1] * sincos[2]) + reference[0];
						yaux =  (y * sincos[1]) + (z * sincos[0]) + reference[1];
						zaux = (x * sincos[2]) - (y * sincos[0] * sincos[3]) + (z * sincos[1] * sincos[3]) + reference[2];

						// Write each cavity point
						fprintf (
                            output, 
                            "ATOM  %5.d  H   K%c%c   259    %8.3lf%8.3lf%8.3lf  1.00%6.2lf\n",
                            count,
                            65+(((cavities[k + nz * (j + ( ny * i ) )]-2)/26)%26),
                            65+((cavities[k + nz * (j + ( ny * i ) )]-2)%26),
                            xaux,
                            yaux,
                            zaux,
                            0.0
                            );

						count++;

						// If count equal to 100,000, restart count
						if (count == 100000)
						    count = 1;
					}
			}
    // Close cavities pdb
	fclose (output);
}

void 
count (int *grid, int nx, int ny, int nz)
{
    int i, j, k, count;

    count = 0;

    for (i=0; i<nx; i++)
        for (j=0; j<ny; j++)
            for (k=0; k<nz; k++)
            {
                if (grid[k + nz * (j + ( ny * i ) )] == 0)
                    count++;
            }
    printf("%d\n", count);

}
