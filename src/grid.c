#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define DEBUG 1

/******* sincos ******
* sincos[0] = sin a  *
* sincos[1] = cos a  *
* sincos[2] = sin b  *
* sincos[3] = cos b  *
*********************/

int 
_detect (int *PI, int size, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, double probe_out, double removal_threshold, double volume_cutoff, int is_ses, int nthreads, int verbose)
{
    int *PO, ncav;

    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {   

            if (verbose)
                fprintf(stdout, "> Filling grid with Probe In\n");
            igrid(PI, size);
            fill(PI, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_in, nthreads);

        }
        #pragma omp section
        {
            if (verbose)
                fprintf(stdout, "> Filling grid with Probe Out\n");
            PO = (int *) calloc (size, sizeof (int));
            igrid(PO, size);
            fill(PO, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_out, nthreads);
        }
    }

    if (is_ses)
        ses(PI, nx, ny, nz, step, probe_in, nthreads);
    ses(PO, nx, ny, nz, step, probe_out, nthreads);

    if (verbose)
        fprintf (stdout, "> Defining biomolecular cavities\n");
    subtract(PI, PO, nx, ny, nz, step, removal_threshold, nthreads);
    filter_noise(PI, nx, ny, nz, nthreads);

    if (verbose)
        fprintf (stdout, "> Clustering cavity points\n");
    ncav = cluster(PI, nx, ny, nz, step, volume_cutoff, nthreads);

    // Free PO
    free(PO);

    return ncav;
}

int 
_detect_badj (int *PI, int size, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, double probe_out, double removal_threshold, double volume_cutoff, int box_adjustment, double *P2, int nndims, int is_ses, int nthreads, int verbose)
{
    int *PO, ncav;

    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {   

            if (verbose)
                fprintf(stdout, "> Filling grid with Probe In\n");
            igrid(PI, size);
            fill(PI, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_in, nthreads);

        }
        #pragma omp section
        {
            if (verbose)
                fprintf(stdout, "> Filling grid with Probe Out\n");
            PO = (int *) calloc (size, sizeof (int));
            igrid(PO, size);
            fill(PO, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_out, nthreads);
        }
    }

    if (is_ses)
        ses(PI, nx, ny, nz, step, probe_in, nthreads);
    ses(PO, nx, ny, nz, step, probe_out, nthreads);

    if (verbose)
        fprintf (stdout, "> Defining biomolecular cavities\n");
    subtract(PI, PO, nx, ny, nz, step, removal_threshold, nthreads);

    if (verbose)
        fprintf (stdout, "> Adjusting biomolecular cavities to box\n");
    filter (PI, nx, ny, nz, reference, ndims, P2, nndims, sincos, nvalues, step, probe_out, nthreads);

    filter_noise(PI, nx, ny, nz, nthreads);

    if (verbose)
        fprintf (stdout, "> Clustering cavity points\n");
    ncav = cluster(PI, nx, ny, nz, step, volume_cutoff, nthreads);

    // Free PO
    free(PO);

    return ncav;
}

int 
_detect_ladj (int *PI, int size, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *ligand, int lnatoms, int lxyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, double probe_out, double removal_threshold, double volume_cutoff, int ligand_adjustment, double ligand_cutoff, int box_adjustment, double *P2, int nndims, int is_ses, int nthreads, int verbose)
{
    int *PO, ncav;

    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {   

            if (verbose)
                fprintf(stdout, "> Filling grid with Probe In\n");
            igrid(PI, size);
            fill(PI, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_in, nthreads);

        }
        #pragma omp section
        {
            if (verbose)
                fprintf(stdout, "> Filling grid with Probe Out\n");
            PO = (int *) calloc (size, sizeof (int));
            igrid(PO, size);
            fill(PO, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_out, nthreads);
        }
    }

    if (is_ses)
        ses(PI, nx, ny, nz, step, probe_in, nthreads);
    ses(PO, nx, ny, nz, step, probe_out, nthreads);

    if (verbose)
        fprintf (stdout, "> Defining biomolecular cavities\n");
    subtract(PI, PO, nx, ny, nz, step, removal_threshold, nthreads);

    if (ligand_adjustment)
    {
        if (verbose)
            fprintf (stdout, "> Adjusting biomolecular cavities to ligand\n");
        adjust(PI, nx, ny, nz, ligand, lnatoms, lxyzr, reference, ndims, sincos, nvalues, step, ligand_cutoff, nthreads);
    }

    if (box_adjustment)
    {
        if (verbose)
            fprintf (stdout, "> Adjusting biomolecular cavities to box\n");
        filter (PI, nx, ny, nz, reference, ndims, P2, nndims, sincos, nvalues, step, probe_out, nthreads);
    }

    filter_noise(PI, nx, ny, nz, nthreads);

    if (verbose)
        fprintf (stdout, "> Clustering cavity points\n");
    ncav = cluster(PI, nx, ny, nz, step, volume_cutoff, nthreads);

    // Free PO
    free(PO);

    return ncav;
}

void 
fill (int *grid, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int nthreads)
{
    int i, j, k, atom;
    double x, y, z, xaux, yaux, zaux, distance, H;

    // Set number of processes in OpenMP
	omp_set_num_threads (nthreads);
    omp_set_nested (1);

    #pragma omp parallel default(none), shared(grid, reference, step, probe, natoms, nx, ny, nz, sincos, atoms, nthreads), private(atom, i, j, k, distance, H, x, y, z, xaux, yaux, zaux)
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
                for (i=floor(x - H); i<=ceil(x + H); i++)
                    for (j=floor(y - H); j<=ceil(y + H); j++)
                        for (k=floor(z - H); k<=ceil(z + H); k++) 
                        {
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
				if (x < 0 || y < 0 || z < 0 || x > nx-1 || y > ny-1 || z > nz-1);
                else
    			    if (grid[ z + nz * (y + ( ny * x ) ) ] == 0 || grid[ z + nz * (y + ( ny * x ) ) ] == -2)
                        return 1;
			}
	return 0;
}

void 
ses (int *grid, int nx, int ny, int nz, double step, double probe, int nthreads)
{
    int i, j, k, i2, j2, k2, aux;
    double distance;    

    // Calculate sas limit in 3D grid units
    aux = ceil(probe / step);

    // Set number of processes in OpenMP
	omp_set_num_threads (nthreads);
    omp_set_nested (1);
    
    #pragma omp parallel default(none), shared(grid,step,probe,aux,nx,ny,nz), private(i,j,k,i2,j2,k2,distance)
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
                                        if (i2>0 && j2>0 && k2>0 && i2<nx && j2<ny && k2<nz)
                                        {
                                            // Get distance between point inspected and cavity point
                                            distance = sqrt ( pow(i - i2, 2) + pow(j - j2, 2) + pow(k - k2, 2));
                                            // Check if inspected point is inside sas limit
                                            if ( distance < (probe / step) )
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
subtract (int *PI, int *PO, int nx, int ny, int nz, double step, double removal_threshold, int nthreads)
{
	int i, j, k, i2, j2, k2, rt;

    rt = ceil (removal_threshold / step);

    // Set number of processes in OpenMP
	omp_set_num_threads (nthreads);
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
filter_noise (int *grid, int nx, int ny, int nz, int nthreads)
{
    int i, j, k, contacts;

    for (i=0; i<nx; i++)
        for (j=0; j<ny; j++)
            for (k=0; k<nz; k++)
            {
                if (grid[ k + nz * (j + ( ny * i ) ) ] == 1) 
                {
                    contacts = 0;
                    
                    // Check if a protein point (0) or a medium point (-1) is next to a cavity point (==1)
                    if (i-1>=0 && i+1<nx && j-1>=0 && j+1<ny && k-1>=0 && k+1<nz)
                    {
                        if (grid[ k + nz * (j + ( ny * (i - 1) ) ) ] == 0 || grid[ k + nz * (j + ( ny * (i - 1) ) ) ] == -1)
                            contacts++;
                        if (grid[ k + nz * (j + ( ny * (i + 1) ) ) ] == 0 || grid[ k + nz * (j + ( ny * (i + 1) ) ) ] == -1)
                            contacts++;
                        if (grid[ k + nz * ( (j - 1) + ( ny * i ) ) ] == 0 || grid[ k + nz * ( (j - 1) + ( ny * i ) ) ] == -1)
                            contacts++;
                        if (grid[ k + nz * ( (j + 1) + ( ny * i ) ) ] == 0 || grid[ k + nz * ( (j + 1) + ( ny * i ) ) ] == -1)
                            contacts++;
                        if (grid[ (k - 1) + nz * (j + ( ny * i ) ) ] == 0 || grid[ (k - 1) + nz * (j + ( ny * i ) ) ] == -1)
                            contacts++;
                        if (grid[ (k + 1) + nz * (j + ( ny * i ) ) ] == 0 || grid[ (k + 1) + nz * (j + ( ny * i ) ) ] == -1)
                            contacts++;

                        // Cavity point is a medium point
                        if (contacts == 6)
                            grid[ k + nz * (j + ( ny * i ) ) ] = -1;
                    }
                }
            }
}

void
adjust (int *grid, int nx, int ny, int nz, double *ligand, int lnatoms, int lxyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double ligand_cutoff, int nthreads)
{
    int i, j, k, atom, inside;
    double x, y, z, xaux, yaux, zaux, distance;

    // Set number of processes in OpenMP
	omp_set_num_threads (nthreads);
    omp_set_nested (1);

    #pragma omp parallel default(none), shared(grid, nx, ny, nz, step, sincos, reference, ligand, ligand_cutoff, lnatoms), private(inside, i, j, k, x, y, z, xaux, yaux, zaux, atom, distance)
    {
        #pragma omp for collapse(3) schedule(static) nowait
        for (i=0; i<nx; i++)
            for (j=0; j<ny; j++)
                for (k=0; k<nz; k++)
                {
                    inside = 0;
                    for (atom=0; atom<lnatoms; atom++)
                    {
                        // Get 3D grid point coordinate
                        x = i * step; 
                        y = j * step; 
                        z = k * step;

                        xaux = (x * sincos[3]) + (y * sincos[0] * sincos[2]) - (z * sincos[1] * sincos[2]) + reference[0];
                        yaux =  (y * sincos[1]) + (z * sincos[0]) + reference[1];
                        zaux = (x * sincos[2]) - (y * sincos[0] * sincos[3]) + (z * sincos[1] * sincos[3]) + reference[2];

                        // Get distance between ligand and 3D grid point evaluated
                        distance = sqrt ( pow(xaux - ligand[atom * 4], 2) + pow(yaux - ligand[1 + (atom * 4)], 2) + pow(zaux - ligand[2 + (atom * 4)], 2) );

                        if (distance < ligand_cutoff)
                            inside = 1;

                    }
                    if (inside == 0 && grid[ k + nz * (j + ( ny * i ) ) ])
                        grid[ k + nz * (j + ( ny * i ) ) ] = -1;
                }             
    }
}

int
filter (int *grid, int nx, int ny, int nz, double *P1, int ndims, double *P2, int nndims, double *sincos, int nvalues, double step, double probe_out, int nthreads)
{
    int i, j, k;
    double aux, normB, norm1, X1, Y1, Z1, X2, Y2, Z2;

    // Set number of processes in OpenMP
	omp_set_num_threads (nthreads);
    omp_set_nested (1);

    // Get norm between X1 and X2 in 3D grid
    X1 = P1[0]; Y1 = P1[1]; Z1 = P1[2];
    X2 = P2[0]; Y2 = P2[1]; Z2 = P2[2];
    norm1 = sqrt (pow(X2 - X1, 2) + pow(Y2 - Y1, 2) + pow(Z2 - Z1, 2));

    // Remove Probe out from grid to recreate box
    X1 -= (- (probe_out * sincos[3]) - (probe_out * sincos[0] * sincos[2]) + (probe_out * sincos[1] * sincos[2]));
    Y1 -= (- (probe_out * sincos[1]) - (probe_out * sincos[0]));
    Z1 -= (- (probe_out * sincos[2]) + (probe_out * sincos[0] * sincos[3]) - (probe_out * sincos[1] * sincos[3]));
    X2 -= ((probe_out * sincos[3]) - (probe_out * sincos[0] * sincos[2]) + (probe_out * sincos[1] * sincos[2]));
    Y2 -= (- (probe_out * sincos[1]) - (probe_out * sincos[0]));
    Z2 -= ((probe_out * sincos[2]) + (probe_out * sincos[0] * sincos[3]) - (probe_out * sincos[1] * sincos[3]));

    // Get norm between X1 and X2 in box
	normB = sqrt (pow(X2 - X1, 2) + pow(Y2 - Y1, 2) + pow(Z2 - Z1, 2));

    // Prepare grid units
    aux = (int) (norm1 - normB) / (2 * step);

    #pragma omp parallel default(shared), private(i, j, k)
    {
        for (i=0; i<=aux; i++)
            #pragma omp for collapse(2) nowait
            for (j=0; j<ny; j++)
                for (k=0; k<nz; k++)
                    grid[ k + nz * (j + ( ny * i ) ) ] = -1;

        for (i=nx-1; i>=nx-aux-1; i--)
            #pragma omp for collapse(2) nowait
            for (j=0; j<ny; j++)
                for (k=0; k<nz; k++)
                    grid[ k + nz * (j + ( ny * i ) ) ] = -1;    

        for (j=0; j<=aux; j++)
            #pragma omp for collapse(2) nowait
            for (i=0; i<nx; i++)
                for (k=0; k<nz; k++)
                    grid[ k + nz * (j + ( ny * i ) ) ] = -1;

        for (j=ny-1; j>=ny-aux-1; j--)
            #pragma omp for collapse(2) nowait
            for (i=0; i<nx; i++)
                for (k=0; k<nz; k++)
                    grid[ k + nz * (j + ( ny * i ) ) ] = -1;    

        for (k=0; k<=aux; k++)
            #pragma omp for collapse(2) nowait
            for (i=0; i<nx; i++)
                for (k=0; k<nz; k++)
                    grid[ k + nz * (j + ( ny * i ) ) ] = -1;

        for (k=nz-1; k>=nz-aux-1; k--)
            #pragma omp for collapse(2) nowait
            for (i=0; i<nx; i++)
                for (j=0; j<nz; j++)
                    grid[ k + nz * (j + ( ny * i ) ) ] = -1;  
    }
}

int vol;

int
cluster (int *grid, int nx, int ny, int nz, double step, double volume_cutoff, int nthreads)
{
    int i, j, k, tag;

    tag = 1;

    for (i=0; i<nx; i++)
        for (j=0; j<ny; j++)
            for (k=0; k<nz; k++)
                if (grid[ k + nz * (j + ( ny * i ) ) ] == 1)
                {
                    tag++;
                    vol = 0;
                    
                    // Clustering procedure
                    DFS(grid, nx, ny, nz, i, j, k, tag);
                    
                    // if (DEBUG)
                    //     printf("Volume (%d): %lf\n", tag-2, (double) vol * step * step * step);

                    // Check if cavity reached cutoff
                    if ( (double) vol * pow(step, 3) < volume_cutoff )
                    {
                        remove_cavity(grid, nx, ny, nz, tag, nthreads);
                        tag--;
                    }
                }
    return tag-1;
}

int 
DFS (int *grid, int nx, int ny, int nz, int i, int j, int k, int tag)
{
    int x, y, z;

    if (i == 0 || i == nx-1 || j == 0 || j == ny-1 || k == 0 || k == nz-1)
        return;

    if (grid[ k + nz * (j + ( ny * i ) ) ] == 1)
    {
        grid[ k + nz * (j + ( ny * i ) ) ] = tag;
        vol++;
        #pragma omp taskloop shared(i, j, k, nx, ny, nz, tag, grid), private(x, y, z)
        for (x=i-1 ; x<=i+1 ; x++)
            for (y=j-1 ; y<=j+1 ; y++)
                for (z = k-1; z<=k+1; z++)
                    DFS(grid, nx, ny, nz, x, y, z, tag);
    }
}

void
remove_cavity (int *grid, int nx, int ny, int nz, int tag, int nthreads)
{
    int i, j, k;

    // Set number of threads in OpenMP
	omp_set_num_threads (nthreads);
    omp_set_nested(1);

    #pragma omp parallel default(shared)
        #pragma omp for schedule(static) collapse(3) nowait
            for (i = 0; i < nx; i++)
                for (j = 0; j < ny; j++)
                    for (k = 0; k < nz; k++)
                        // Remove points based on tag
                        if (grid[ k + nz * (j + ( ny * i ) ) ] == tag)
                            grid[ k + nz * (j + ( ny * i ) ) ] = 0;
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
_spatial (int *cavities, int nx, int ny, int nz, int *surface, int size, double *volumes, int nvol, double *areas, int narea, double step, int nthreads, int verbose)
{
    if (verbose)
        fprintf (stdout, "> Defining surface points\n");
    filter_surface (cavities, surface, nx, ny, nz, nthreads);

    #pragma omp sections
    {
        #pragma omp section
        {
            if (verbose)
                fprintf (stdout, "> Estimating volume\n");
            volume (cavities, nx, ny, nz, nvol, step, volumes, nthreads);
            // if (DEBUG) 
            //     for (int i=0; i<nvol; i++) 
            //         printf("%d: %lf\n", i, volumes[i]);
        }

        #pragma omp section
        {
            if (verbose)
                fprintf (stdout, "> Estimating area\n");
            area (surface, nx, ny, nz, narea, step, areas, nthreads);
            // if (DEBUG) 
            //     for (int i=0; i<narea; i++) 
            //         printf("%d: %lf\n", i, areas[i]);
        }
    }    
}

int
define_surface_points (int *grid, int nx, int ny, int nz, int i, int j, int k)
{
    if (i-1>=0 && i+1<nx && j-1>=0 && j+1<ny && k-1>=0 && k+1<nz)
    {
        if (grid[ k + nz * (j + ( ny * (i - 1) ) ) ] == 0)
            return grid[k + nz * (j + ( ny * i ) )];
        if (grid[ k + nz * (j + ( ny * (i + 1) ) ) ] == 0)
            return grid[k + nz * (j + ( ny * i ) )];
        if (grid[ k + nz * ( (j - 1) + ( ny * i ) ) ] == 0)
            return grid[k + nz * (j + ( ny * i ) )];
        if (grid[ k + nz * ( (j + 1) + ( ny * i ) ) ] == 0)
            return grid[k + nz * (j + ( ny * i ) )];
        if (grid[ (k - 1) + nz * (j + ( ny * i ) ) ] == 0)
            return grid[k + nz * (j + ( ny * i ) )];
        if (grid[ (k + 1) + nz * (j + ( ny * i ) ) ] == 0)
            return grid[k + nz * (j + ( ny * i ) )];
    }
	return -1;
}

void 
filter_surface (int *cavities, int *surface, int nx, int ny, int nz, int nthreads)
{
    int i, j, k;

    // Set number of threads in OpenMP
    omp_set_num_threads(nthreads);
    omp_set_nested(1);

    #pragma omp parallel default(none), shared(cavities, surface, nx, ny, nz), private(i, j, k)
    {
        #pragma omp for collapse(3) schedule(static)
        for (i=0; i<nx; i++)
            for (j=0; j<ny; j++)
                for (k=0; k<nz; k++)
                    if (cavities[k + nz * (j + ( ny * i ) )] > 1)
                    {
                        // Define surface cavity points
                        surface[k + nz * (j + ( ny * i ) )] = define_surface_points (cavities, nx, ny, nz, i, j, k);
                    }
                    else
                    {
                        if (cavities[k + nz * (j + ( ny * i ) )] == 0)
                            surface[k + nz * (j + ( ny * i ) )] = 0;
                        else
                            surface[k + nz * (j + ( ny * i ) )] = -1;
                    }
    }

}

void
volume (int *cavities, int nx, int ny, int nz, int ncav, double step, double *volumes, int nthreads) 
{
    int i, j, k;

    // Set number of threads in OpenMP
    omp_set_num_threads(nthreads);
    omp_set_nested(1);

    for (i=0; i<ncav; i++)
        volumes[i] = 0.0;

    #pragma omp parallel default(none), shared(volumes, cavities, ncav, step, nx, ny, nz), private(i, j, k)
    {
        #pragma omp for collapse(3) reduction(+:volumes[:ncav])
            for (i=0; i<nx; i++)
                for (j=0; j<ny; j++)
                    for (k=0; k<nz; k++)
                        if (cavities[k + nz * (j + ( ny * i ) )] > 1)
                            volumes[cavities[k + nz * (j + ( ny * i ) )] - 2] += pow (step, 3);
    }
}

double 
check_voxel_class (int *grid, int nx, int ny, int nz, int i, int j, int k)
{
    int contacts = 0;
    double weight = 1.0;

    // Count face contacts
    if (grid[ k + nz * (j + ( ny * (i - 1) ) ) ] == 0)
        contacts++;
    if (grid[ k + nz * (j + ( ny * (i + 1) ) ) ] == 0)
        contacts++;
    if (grid[ k + nz * ( (j - 1) + ( ny * i ) ) ] == 0)
        contacts++;
    if (grid[ k + nz * ( (j + 1) + ( ny * i ) ) ] == 0)
        contacts++;
    if (grid[ (k - 1) + nz * (j + ( ny * i ) ) ] == 0)
        contacts++;
    if (grid[ (k + 1) + nz * (j + ( ny * i ) ) ] == 0)
        contacts++;

    // Get weight
    switch (contacts)
    {
        // One face in contact with biomolecule
        case 1:
			weight = 0.894;
			break;
        // Two faces in contact with biomolecule
		case 2:
			weight = 1.3409;
			break;
		// Three faces in contact with biomolecule
		case 3:
			if ( ( grid[ k + nz * (j + ( ny * (i + 1) ) ) ] == 0 && grid[ k + nz * (j + ( ny * (i - 1) ) ) ] ) || ( grid[ k + nz * ( (j - 1) + ( ny * i ) ) ] == 0 && grid[ k + nz * ( (j - 1) + ( ny * i ) ) ] == 0 ) || ( grid[ (k + 1) + nz * (j + ( ny * i ) ) ] == 0 && grid[ (k + 1) + nz * (j + ( ny * i ) ) ] == 0 ) )
			    weight = 2;
			else
			    weight = 1.5879;
			break;
	    // Four faces in contact with biomolecule
		case 4:
			weight = 2.6667;
			break;
        // Five in contact with biomolecule
		case 5:
			weight = 3.3333;
			break;
    }
    return weight;
}

void
area (int *surface, int nx, int ny, int nz, int ncav, double step, double *areas, int nthreads)
{
    int i, j, k;

    // Set number of threads in OpenMP
    omp_set_num_threads (nthreads);
    omp_set_nested (1);

    for (i=0; i<ncav; i++)
        areas[i] = 0.0;

    for (i=0; i<nx; i++)
        for (j=0; j<ny; j++)
            for (k=0; k<nz; k++)
                if (surface[k + nz * (j + ( ny * i ) )] > 1)
                    areas[surface[k + nz * (j + ( ny * i ) )] - 2] += check_voxel_class(surface, nx, ny, nz, i, j, k) * pow (step, 2);
}

typedef struct node {
    int pos;
    struct node* next;
} res;

res* 
create (int pos)
{
    res* new = (res*) malloc (sizeof(res));

    new->pos = pos;
    new->next = NULL;

    return new; 
}

void 
insert (res** head, res* new)
{
    res* current;

    if (*head == NULL || (*head)->pos >= new->pos)
    {
        new->next = *head;
        *head = new;
    } 
    else
    {
        current = *head;
        while (current->next != NULL && current->next->pos < new->pos) 
        {
            current = current->next;
        }
        new->next = current->next;
        current->next = new;
    }
}

char
**interface (int *grid, int nx, int ny, int nz, char **pdb, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, int ncav, int nthreads)
{
    int i, j, k, atom, tag, count=0, old_atom=-1, old_tag=-1;
    double x, y, z, xaux, yaux, zaux, distance, H;
    char **residues;

    // Allocate memory for reslist structure    
    res *reslist[ncav], *new;
    
    // Initialize linked list
    for (i=0; i<ncav; i++)
        reslist[i] = NULL;

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
            H = ( probe_in + atoms[3 + (atom * 4)] ) / step;

            // Loop around radius from atom center        
            for (i=floor(x - H); i<=ceil(x + H); i++)
                for (j=floor(y - H); j<=ceil(y + H); j++)
                    for (k=floor(z - H); k<=ceil(z + H); k++) 
                    {
                        if (i < nx && i > 0 && j < ny && j > 0 && k < nz && k > 0)
                            if (abs(grid[ k + nz * (j + ( ny * i ) ) ]) > 1)
                            {
                                tag = grid[ k + nz * (j + ( ny * i ) ) ] - 2;
                                distance = sqrt( pow(i - x, 2) + pow(j - y, 2) + pow(k - z, 2));
                                if (distance <= H)
                                {   
                                    if (old_atom != atom || old_tag != tag)
                                    {
                                        new = create(atom);
                                        insert(&reslist[tag], new);
                                        count++;
                                    }
                                    old_atom = atom;
                                    old_tag = tag;
                                }
                            }
                    }
    }
    
    // Pass res information to char **
    residues = calloc(count + ncav + 1, sizeof(char*));
    for (i=0, j=0; i<ncav; i++)
    {   
        new = reslist[i];
        while (new != NULL) { 
            residues[j++] = pdb[new->pos];
            new = new->next;
        }
        free(reslist[i]);
        residues[j++] = "-1";
    }
    residues[j] = NULL;
    return residues;
}

char
** _constitutional (int *cavities, int nx, int ny, int nz, char **pdb, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, int ncav, int nthreads, int verbose)
{
    char **residues;

    if (verbose)
        fprintf (stdout, "> Retrieving interface residues\n");
    residues = interface(cavities, nx, ny, nz, pdb, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_in, ncav, nthreads);

    return residues;
}

void
_export (char *fn, int *cavities, int nx, int ny, int nz, int *surf, int nxx, int nyy, int nzz, double *reference, int ndims, double *sincos, int nvalues, double step, int ncav, int nthreads)
{
	int i, j, k, count, tag;
	double x, y, z, xaux, yaux, zaux;
	FILE *output;

    // Set number of threads in OpenMP
    omp_set_num_threads(nthreads);
    omp_set_nested(1);

	// Open cavity PDB file
	output = fopen (fn, "w");

    for (count=1, tag=2; tag<=ncav+2; tag++)
        #pragma omp parallel default(none) shared(cavities, surf, reference, sincos, step, ncav, tag, count, nx, ny, nz, output), private(i, j, k, x, y, z, xaux, yaux, zaux)
        {
            #pragma omp for schedule(static) collapse(3) ordered nowait
            for (i=0; i<nx; i++)
                for (j=0; j<ny; j++)
                    for (k=0; k<nz; k++) 
                    {
                        // Check if cavity point with value tag
                        if ( cavities[k + nz * (j + ( ny * i ) )] == tag ) 
                        {
                            // Convert 3D grid coordinates to real coordinates
                            x = i * step; 
                            y = j * step; 
                            z = k * step;
                            
                            xaux = (x * sincos[3]) + (y * sincos[0] * sincos[2]) - (z * sincos[1] * sincos[2]) + reference[0];
                            yaux =  (y * sincos[1]) + (z * sincos[0]) + reference[1];
                            zaux = (x * sincos[2]) - (y * sincos[0] * sincos[3]) + (z * sincos[1] * sincos[3]) + reference[2];

                            // Write cavity point
                            #pragma omp critical
                            if ( surf[k + nz * (j + ( ny * i ) )] == tag )
                                fprintf (output, "ATOM  %5.d  HA  K%c%c   259    %8.3lf%8.3lf%8.3lf  1.00%6.2lf\n", count % 100000, 65 + (((surf[k + nz * (j + ( ny * i ) )]-2) / 26) % 26), 65 + ((surf[k + nz * (j + ( ny * i ) )]-2) % 26), xaux, yaux, zaux, 0.0);
                            else
                                fprintf (output, "ATOM  %5.d  H   K%c%c   259    %8.3lf%8.3lf%8.3lf  1.00%6.2lf\n", count % 100000, 65 + (((cavities[k + nz * (j + ( ny * i ) )]-2) / 26) % 26), 65 + ((cavities[k + nz * (j + ( ny * i ) )]-2) % 26), xaux, yaux, zaux, 0.0);
                            count++;
                        }
                    }
        }
        fclose (output);
}
