#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

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


/******* sincos ******
* sincos[0] = sin a  *
* sincos[1] = cos a  *
* sincos[2] = sin b  *
* sincos[3] = cos b  *
*********************/
// TODO: Improve parallization in OpenMP
void 
fill_grid (int *grid, int dx, int dy, int dz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int ncores, int is_sas)
{
    int i, j, k, imax, jmax, kmax, atom;
    double distance, H, x, y, z, xaux, yaux, zaux;


    // Set number of processes in OpenMP
	omp_set_num_threads (ncores);
    omp_set_nested (1);

    #pragma omp parallel default(none), shared(grid, reference, step, probe, natoms, dx, dy, dz, sincos, atoms, ncores, is_sas), private(atom, i, j, k, imax, jmax, kmax, distance, H, x, y, z, xaux, yaux, zaux)
    {
    #pragma omp for schedule(dynamic) nowait
        for (atom=0; atom<natoms; atom++) 
        {
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
            if (is_sas)
                H = ( probe + atoms[4 + (atom * 4)] ) / step;
            else 
                H = ( atoms[4 + (atom * 4)] ) / step;

            imax = ceil(x + H)+1;
            jmax = ceil(y + H)+1;
            kmax = ceil(z + H)+1;

            // #pragma omp simd
            // #pragma preomp for collapse(3) schedule(dynamic) ordered nowait
                // Loop around radius from atom center
                for (i=floor(x - H); i<imax; i++)
                    for (j=floor(y - H); j<jmax; j++)
                        for (k=floor(z - H); k<kmax; k++) 
                        {

                            // Get distance between atom center and point inspected
                            distance = sqrt( pow(i - x, 2) + pow(j - y, 2) + pow(k - z, 2));
                            if (distance < H)
                                if (i >= 0 && i < dx && j >= 0 && j < dy && k >= 0 && i < dz)
                                    grid[ k + dz * (j + ( dy * i ) ) ] = 0;

                        }
        }
    }
}


void
subtract (int *grid, int dx, int dy, int dz, int *grid2, int dx2, int dy2, int dz2, double step, double removal_threshold, int ncores)
{
	int i, j, k, i2, j2, k2, removal;
    
    removal = ceil (removal_threshold / step);

    // Set number of processes in OpenMP
	omp_set_num_threads (ncores);
    omp_set_nested (1);

    /* Create a parallel region */
    #pragma omp parallel default(none), shared(grid,grid2,removal,dx,dy,dz,dx2,dy2,dz2,i,j,k), private(j2,i2,k2)
    {
        #pragma omp for schedule(dynamic) collapse(3)
            // Loop around the search box
            for (i=0; i<dx; i++)
                for (j=0; j<dy; j++)
                    for (k=0; k<dz; k++) 
                    {
                        // Check if point is a cavity point
                        if (grid[ k + dz * (j + ( dy * i ) ) ]) 
                        {

                            // Loops around space occupied by probe from atom position
                            for(i2=i-removal; i2<=i+removal; i2++)
                                for(j2=j-removal; j2<=j+removal; j2++)
                                    for(k2=k-removal; k2<=k+removal; k2++)
                                        // Check if inside 3D grid
                                        if(i2 >= 0 && i2 < dx2 && j2 >= 0 && j2 < dy2 && k2 >= 0 && k2 < dz2)
                                            // Mark points where big probe passed in cavities in A
                                            grid2[ k2 + dz2 * (j2 + ( dy2 * i2 ) ) ] = -1;

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
check_protein_neighbours (int *grid, int dx, int dy, int dz, int i, int j, int k)
{
	int x, y, z;

	// Loop around neighboring points
	for (x=i-1; x<=i+1; x++)
        for (y=j-1; y<=j+1; y++)
			for (z=k-1; z<=k+1; z++) 
            {

			    // Check if point is inside 3D grid
				if (x < 0 || y < 0 || z < 0 || x >= dx || y >= dy || z >= dz);
				else
				    // Check if point is a protein protein
				    if (grid[ z + dz * (y + ( dy * x ) ) ] == 0 || grid[ z + dz * (y + ( dy * x ) ) ] == -2)
				        return 1;

			}

	return 0;

}

void
ses (int *grid, int dx, int dy, int dz, double step, double sas, int ncores)
{
    int i, j, k, i2, j2, k2, aux;
    double distance;    
    
    // Set number of processes in OpenMP
	omp_set_num_threads (ncores);
    omp_set_nested (1);
    
    // Calculate sas limit in 3D grid units
    aux = (int) (sas / step) + 1;

    #pragma omp parallel default(none), shared(grid,dx,dy,dz,step,sas,aux), private(i,j,k,i2,j2,k2,distance)
    {   
        #pragma omp for schedule(dynamic) collapse(3)
        // Loop around 3D grid
        for (i=0; i<dx; i++)
            for (j=0; j<dy; j++)
                for (k=0; k<dz; k++) 
                {
                    // Check if a cavity point
                    if (grid[ k + dz * (j + ( dy * i ) ) ] == 1)
                    {
                        if (check_protein_neighbours(grid, dx, dy, dz, i, j, k))
                        {
                            // Loop around sas limit from cavity point next to protein point
                            for (i2=i-aux; i2<=i+aux; i2++)
                                for (j2=j-aux; j2<=j+aux; j2++)
                                    for (k2=k-aux; k2<=k+aux; k2++)
                                    { 
                                        if (i2>0 && j2>0 && k2>0 && i2<=dx && j2<=dy && k2<=dz)
                                        {
                                            // Get distance between point inspected and cavity point
                                            distance = sqrt ( pow(i - i2, 2) + pow(j - j2, 2) + pow(k - k2, 2));
                                            // Check if inspected point is inside sas limit
                                            if ( distance < sas / step )
                                                // Mark cavity point
                                                grid[ k2 + dz * (j2 + ( dy * i2 ) ) ] = -2;
                                        }
                                    }
                        }
                    }
                }

        #pragma omp for collapse(3)
            // Loop around 3D grid
            for (i=0; i<dx; i++)
                for (j=0; j<dy; j++)
                    for (k=0; k<dz; k++) 
                    {
                        // Mark space occupied by sas limit from protein surface
                        if (grid[ k + dz * (j + ( dy * i ) ) ] == -2)
                            grid[ k + dz * (j + ( dy * i ) ) ] = 1;

                    }
    }

}

void
export (int *grid, int dx, int dy, int dz, int *grid2, int dx2, int dy2, int dz2, double step, char *fn, double *reference, int ndims, double *sincos, int nvalues, int cavity_representation)
{
	int i, j, k, count = 1, control = 1, tag = 1;
	double x, y, z, xaux, yaux, zaux;
	FILE *output;

	// Open cavity PDB file
	output = fopen (fn, "w");

	while (control == 1) 
    {
	    control = 0;

		// Loop around grid
		for (i=0; i<dx; i++)
			for (j=0; j<dy; j++)
				for (k=0; k<dz; k++) 
                {
                    // Check if cavity point with value tag
					if (grid[k + dz * (j + ( dy * i ) )] == tag) 
                    {

						control = 1;

						// Convert 3D grid coordinates to real coordinates
						x = i * step; 
                        y = j * step; 
                        z = k * step;
						
                        xaux = (x * sincos[3]) + (y * sincos[0] * sincos[2]) - (z * sincos[1] * sincos[2]) + reference[0];
						yaux =  (y * sincos[1]) + (z * sincos[0]) + reference[1];
						zaux = (x * sincos[2]) - (y * sincos[0] * sincos[3]) + (z * sincos[1] * sincos[3]) + reference[2];

						/* Write each cavity point */
						fprintf (
                            output, 
                            "ATOM  %5.d  HS  K%c%c   259    %8.3lf%8.3lf%8.3lf  1.00%6.2lf\n",
                            count,
                            65+(((grid[k + dz * (j + ( dy * i ) )]-2)/26)%26),
                            65+((grid[k + dz * (j + ( dy * i ) )]-2)%26),
                            xaux,
                            yaux,
                            zaux,
                            0.0
                            );

						count++;

						/* If count equal to 100,000, restart count */
						if (count == 100000)
						    count = 1;

					}
				}

				/* Next cavity tag */
				tag++;

	}

    /* Close output PDB file */
	fclose (output);

}
