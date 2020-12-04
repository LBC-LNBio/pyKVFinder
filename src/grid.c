#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

# define min(x,y) ( ((x) < (y)) ? (x) : (y))
# define max(x,y) ( ((x) > (y)) ? (x) : (y))

/******* sincos ******
* sincos[0] = sin a  *
* sincos[1] = cos a  *
* sincos[2] = sin b  *
* sincos[3] = cos b  *
*********************/

/* Grid initialization */

/*
 * Function: igrid
 * ---------------
 * 
 * Fill integer grid with 1
 * 
 * grid: empty 3D grid
 * size: number of voxels
 * 
 */
void
igrid (int *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = 1;              

}

/*
 * Function: fgrid
 * ---------------
 * 
 * Fill float grid with 0.0
 * 
 * grid: empty 3D grid
 * size: number of voxels
 * 
 */
void
fgrid (float *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = 0.0;
               
}

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
 * Function: cgrid
 * ---------------
 * 
 * Fill char grid with '\0'
 * 
 * grid: empty 3D grid
 * size: number of voxels
 * 
 */
void 
cgrid (int *grid, int size)
{
    int i;

    for (i=0; i<size; i++)
        grid[i] = '\0';
               
}

/* Grid filling */

/*
 * Function: fill
 * --------------
 * 
 * Insert atoms with a probe addition inside a 3D grid
 * 
 * grid: 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * atoms: xyz coordinates and radii of input pdb
 * natoms: number of atoms
 * xyzr: number of data per atom (4: xyzr)
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * probe: Probe size (A)
 * nthreads: number of threads for OpenMP
 * 
 */
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
                                if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz)
                                    grid[ k + nz * (j + ( ny * i ) ) ] = 0;
                        }
            }
        }
}

/* Biomolecular surface representation */

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

/*
 * Function: ses
 * --------------
 * 
 * Adjust surface representation to Solvent Excluded Surface (SES)
 * 
 * grid: 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * step: 3D grid spacing (A)
 * probe: Probe size (A)
 * nthreads: number of threads for OpenMP
 * 
 */
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

/* Grid subtract (Probe In - Probe Out) */

/*
 * Function: subtract
 * ------------------
 * 
 * Compare Probe In and Probe Out 3D grids to define biomolecular cavities
 * 
 * PI: Probe In 3D grid
 * PO: Probe Out 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * step: 3D grid spacing (A)
 * removal_distance: Length to be removed from the cavity-bulk frontier (A)
 * nthreads: number of threads for OpenMP
 * 
 */
void
subtract (int *PI, int *PO, int nx, int ny, int nz, double step, double removal_distance, int nthreads)
{
	int i, j, k, i2, j2, k2, rd;

    rd = ceil (removal_distance / step);

    // Set number of processes in OpenMP
	omp_set_num_threads (nthreads);
    omp_set_nested (1);

    /* Create a parallel region */
    #pragma omp parallel default(none), shared(PI, PO, nx, ny, nz, i, j, k, step, rd, removal_distance), private(j2,i2,k2)
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
                            for(i2=i-rd; i2<=i+rd; i2++)
                                for(j2=j-rd; j2<=j+rd; j2++)
                                    for(k2=k-rd; k2<=k+rd; k2++)
                                        // Check if inside 3D grid
                                        if(i2>=0 && i2<nx && j2>=0 && j2<ny && k2>=0 && k2<nz)
                                            // Mark points in grid filled with Probe In, where Probe Out reached
                                            PI[ k2 + nz * (j2 + ( ny * i2 ) ) ] = -1;
                        }
                    }
    }
}

/* Filter noise from Grid */

/*
 * Function: filter_noise
 * ----------------------
 * 
 * Removes cavities points (1) surrounded by biomolecule (0) or medium (-1) points
 * 
 * grid: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * nthreads: number of threads for OpenMP
 * 
 */
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

/* Ligand adjustment */

/*
 * Function: adjust
 * ----------------
 * 
 * Adjust cavities to a radius around a ligand atoms
 * 
 * grid: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * ligand: xyz coordinates and radii of ligand
 * lnatoms: number of ligand atoms
 * lxyzr: number of data per ligand atom (4: xyzr)
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * ligand_cutoff: Radius value to limit a space around a ligand (A)
 * nthreads: number of threads for OpenMP
 * 
 */
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

/* Box adjustment */

/*
 * Function: _filter_pdb
 * ---------------------
 * 
 * Select pdb atoms inside the 3D grid (box adjustment mode)
 * 
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * atoms: xyz coordinates and radii of input pdb
 * natoms: number of atoms
 * xyzr: number of data per atom (4: xyzr)
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * probe_in: Probe In size (A)
 * nthreads: number of threads for OpenMP
 * 
 */
void
_filter_pdb (int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, int nthreads)
{
    int atom;
    double x, y, z, xaux, yaux, zaux;

    // Set number of processes in OpenMP
	omp_set_num_threads (nthreads);

    #pragma omp parallel default(none), shared(atoms, natoms, reference, sincos, step, probe_in, nx, ny, nz), private(atom, x, y, z, xaux, yaux, zaux)
    {   
        #pragma omp for schedule(static) //nowait
            for (atom=0; atom<natoms; atom++)
            {
                x = ( atoms[atom * 4] - reference[0] ) / step; 
                y = ( atoms[1 + (atom * 4)] - reference[1] ) / step; 
                z = ( atoms[2 + (atom * 4)] - reference[2] ) / step;

                xaux = x * sincos[3] + z * sincos[2];
                yaux = y;
                zaux = (-x) * sincos[2] + z * sincos[3];

                x = xaux;
                y = yaux * sincos[1] - zaux * sincos[0];
                z = yaux * sincos[0] + zaux * sincos[1];

                if (x >  0.0 - (probe_in + atoms[3 + (atom * 4)]) / step && 
                    x < (double)nx + (probe_in + atoms[3 + (atom * 4)]) / step && 
                    y > 0.0 - (probe_in + atoms[3 + (atom * 4)]) / step && 
                    y < (double)ny + (probe_in + atoms[3 + (atom * 4)]) / step && 
                    z >  0.0 - (probe_in + atoms[3 + (atom * 4)]) / step && 
                    z < (double)nz + (probe_in + atoms[3 + (atom * 4)]) / step);
                else
                {
                    atoms[atom * 4] = 0.0;
                    atoms[1 + (atom * 4)] = 0.0;
                    atoms[2 + (atom * 4)] = 0.0;
                    atoms[3 + (atom * 4)] = 0.0;
                }

            }
    }
}

/*
 * Function: filter
 * ----------------
 * 
 * Adjust cavities to a search box
 * 
 * grid: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * P1: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * P2: xyz coordinates of x-axis vertice
 * nndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * probe_out: Radius value to limit a space around a ligand (A)
 * nthreads: number of threads for OpenMP
 * 
 */
void
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
    aux = floor ( (norm1 - normB) / (2 * step) );

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
                for (j=0; j<ny; j++)
                    grid[ k + nz * (j + ( ny * i ) ) ] = -1;

        for (k=nz-1; k>=nz-aux-1; k--)
            #pragma omp for collapse(2) nowait
            for (i=0; i<nx; i++)
                for (j=0; j<ny; j++)
                    grid[ k + nz * (j + ( ny * i ) ) ] = -1;  
    }
}

/* Cavity clustering */

/*
 * Variable: vol
 * -------------
 * 
 * Accumulate volume of a cavity while clustering a cavity
 *
 */
int vol;

/*
 * Function: DFS
 * -------------
 * 
 * Recursive DFS algorithm
 * 
 * grid: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * i: x coordinate of cavity point
 * j: y coordinate of cavity point
 * k: z coordinate of cavity point
 * tag: cavity integer identifier
 * 
 */
void 
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

/*
 * Function: remove_cavity
 * -----------------------
 * 
 * Untag cavity that does not reach volume cutoff
 * 
 * grid: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * tag: cavity integer identifier
 * nthreads: number of threads for OpenMP
 * 
 */
void
remove_cavity (int *grid, int nx, int ny, int nz, int tag, int nthreads)
{
    int i, j, k;

    // Set number of threads in OpenMP
	omp_set_num_threads (nthreads);
    omp_set_nested(1);

    #pragma omp parallel default(shared)
        #pragma omp for schedule(static) collapse(3) //nowait
            for (i = 0; i < nx; i++)
                for (j = 0; j < ny; j++)
                    for (k = 0; k < nz; k++)
                        // Remove points based on tag
                        if (grid[ k + nz * (j + ( ny * i ) ) ] == tag)
                            grid[ k + nz * (j + ( ny * i ) ) ] = 0;
}

/*
 * Function: cluster
 * -----------------
 * 
 * Cluster consecutive cavity points together
 * 
 * grid: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * step: 3D grid spacing (A)
 * volume_cutoff: Cavities volume filter (A3)
 * nthreads: number of threads for OpenMP
 * 
 */
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

                    // Check if cavity reached cutoff
                    if ( (double) vol * pow(step, 3) < volume_cutoff )
                    {
                        remove_cavity(grid, nx, ny, nz, tag, nthreads);
                        tag--;
                    }
                }
    return tag-1;
}

/* Cavity detection */

/*
 * Function: _detect
 * -----------------
 * 
 * Detect and cluster cavities
 * 
 * PI: 3D grid
 * size: number of voxels in 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * atoms: xyz coordinates and radii of input pdb
 * natoms: number of atoms
 * xyzr: number of data per atom (4: xyzr)
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * probe_in: Probe In size (A)
 * probe_out: Probe Out size (A)
 * removal_distance: Length to be removed from the cavity-bulk frontier (A)
 * volume_cutoff: Cavities volume filter (A3)
 * box_adjustment: Box adjustment mode
 * P2: xyz coordinates of x-axis vertice
 * nndims: number of coordinates (3: xyz)
 * is_ses: surface mode (1: SES or 0: SAS)
 * nthreads: number of threads for OpenMP
 * verbose: print extra information to standard output
 * 
 * returns: PI[size] (cavities 3D grid) and ncav (number of cavities)
 */
int 
_detect (int *PI, int size, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, double probe_out, double removal_distance, double volume_cutoff, int box_adjustment, double *P2, int nndims, int is_ses, int nthreads, int verbose)
{
    int *PO, ncav;

    if (verbose)
        fprintf(stdout, "> Filling grid with Probe In\n");
    igrid(PI, size);
    fill(PI, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_in, nthreads);

    if (verbose)
        fprintf(stdout, "> Filling grid with Probe Out\n");
    PO = (int *) calloc (size, sizeof (int));
    igrid(PO, size);
    fill(PO, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_out, nthreads);

    if (is_ses)
        ses(PI, nx, ny, nz, step, probe_in, nthreads);
    ses(PO, nx, ny, nz, step, probe_out, nthreads);

    if (verbose)
        fprintf (stdout, "> Defining biomolecular cavities\n");
    subtract(PI, PO, nx, ny, nz, step, removal_distance, nthreads);

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

/*
 * Function: _detect_ladj
 * ----------------------
 * 
 * Detect and cluster cavities with ligand adjustment
 * 
 * PI: 3D grid
 * size: number of voxels in 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * atoms: xyz coordinates and radii of input pdb
 * natoms: number of atoms
 * xyzr: number of data per atom (4: xyzr)
 * ligand: xyz coordinates and radii of ligand
 * lnatoms: number of ligand atoms
 * lxyzr: number of data per ligand atom (4: xyzr)
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * probe_in: Probe In size (A)
 * probe_out: Probe Out size (A)
 * removal_distance: Length to be removed from the cavity-bulk frontier (A)
 * volume_cutoff: Cavities volume filter (A3)
 * ligand_adjustment: ligand adjustment mode
 * ligand_cutoff: Radius value to limit a space around a ligand (A)
 * box_adjustment: Box adjustment mode
 * P2: xyz coordinates of x-axis vertice
 * nndims: number of coordinates (3: xyz)
 * is_ses: surface mode (1: SES or 0: SAS)
 * nthreads: number of threads for OpenMP
 * verbose: print extra information to standard output
 * 
 * returns: PI[size] (cavities 3D grid) and ncav (number of cavities)
 */
int 
_detect_ladj (int *PI, int size, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *ligand, int lnatoms, int lxyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, double probe_out, double removal_distance, double volume_cutoff, int ligand_adjustment, double ligand_cutoff, int box_adjustment, double *P2, int nndims, int is_ses, int nthreads, int verbose)
{
    int *PO, ncav;

    if (verbose)
        fprintf(stdout, "> Filling grid with Probe In\n");
    igrid(PI, size);
    fill(PI, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_in, nthreads);

    if (verbose)
        fprintf(stdout, "> Filling grid with Probe Out\n");
    PO = (int *) calloc (size, sizeof (int));
    igrid(PO, size);
    fill(PO, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_out, nthreads);

    if (is_ses)
        ses(PI, nx, ny, nz, step, probe_in, nthreads);
    ses(PO, nx, ny, nz, step, probe_out, nthreads);

    if (verbose)
        fprintf (stdout, "> Defining biomolecular cavities\n");
    subtract(PI, PO, nx, ny, nz, step, removal_distance, nthreads);

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

/* Cavity surface points */

/*
 * Function: define_surface_points
 * -------------------------------
 * 
 * Identify surface points based on neighboring points
 * 
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * i: x coordinate of cavity point
 * j: y coordinate of cavity point
 * k: z coordinate of cavity point
 * 
 * returns: cavity identifier (>1) or medium point (-1)
 */
int
define_surface_points (int *cavities, int nx, int ny, int nz, int i, int j, int k)
{
    if (i-1>=0)
        if (cavities[ k + nz * (j + ( ny * (i - 1) ) ) ] == 0)
            return -(cavities[k + nz * (j + ( ny * i ) )]);
    if (i+1<nx)
        if (cavities[ k + nz * (j + ( ny * (i + 1) ) ) ] == 0)
            return -(cavities[k + nz * (j + ( ny * i ) )]);
    if (j-1>=0)
        if (cavities[ k + nz * ( (j - 1) + ( ny * i ) ) ] == 0)
            return -(cavities[k + nz * (j + ( ny * i ) )]);
    if (j+1<ny)
        if (cavities[ k + nz * ( (j + 1) + ( ny * i ) ) ] == 0)
            return -(cavities[k + nz * (j + ( ny * i ) )]);
    if (k-1>=0)
        if (cavities[ (k - 1) + nz * (j + ( ny * i ) ) ] == 0)
            return -(cavities[k + nz * (j + ( ny * i ) )]);
    if (k+1<nz)
        if (cavities[ (k + 1) + nz * (j + ( ny * i ) ) ] == 0)
            return -(cavities[k + nz * (j + ( ny * i ) )]);

	return -1;
}

/*
 * Function: filter_surface
 * ------------------------
 * 
 * Inspect cavities 3D grid and mark detected surface points on a surface 3D grid
 * 
 * cavities: cavities 3D grid
 * surface: surface points 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * nthreads: number of threads for OpenMP
 * 
 */
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

/* Estimate area */

/*
 * Function: check_voxel_class
 * ---------------------------
 * 
 * Identify voxel class of surface voxel and return class weight
 * 
 * surface: surface points 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * i: x coordinate of cavity point
 * j: y coordinate of cavity point
 * k: z coordinate of cavity point
 * 
 * returns: voxel class weight (double)
 */
double 
check_voxel_class (int *surface, int nx, int ny, int nz, int i, int j, int k)
{
    int contacts = 0;
    double weight = 1.0;

    // Count face contacts
    if (surface[ k + nz * (j + ( ny * (i - 1) ) ) ] == 0)
        contacts++;
    if (surface[ k + nz * (j + ( ny * (i + 1) ) ) ] == 0)
        contacts++;
    if (surface[ k + nz * ( (j - 1) + ( ny * i ) ) ] == 0)
        contacts++;
    if (surface[ k + nz * ( (j + 1) + ( ny * i ) ) ] == 0)
        contacts++;
    if (surface[ (k - 1) + nz * (j + ( ny * i ) ) ] == 0)
        contacts++;
    if (surface[ (k + 1) + nz * (j + ( ny * i ) ) ] == 0)
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
			if ( ( surface[ k + nz * (j + ( ny * (i + 1) ) ) ] == 0 && surface[ k + nz * (j + ( ny * (i - 1) ) ) ] ) || ( surface[ k + nz * ( (j - 1) + ( ny * i ) ) ] == 0 && surface[ k + nz * ( (j - 1) + ( ny * i ) ) ] == 0 ) || ( surface[ (k + 1) + nz * (j + ( ny * i ) ) ] == 0 && surface[ (k + 1) + nz * (j + ( ny * i ) ) ] == 0 ) )
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

/*
 * Function: area
 * --------------
 * 
 * Calculate area of cavities
 * 
 * surface: surface points 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * ncav: number of cavities
 * step: 3D grid spacing (A)
 * areas: empty array of areas
 * nthreads: number of threads for OpenMP
 * 
 */
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

/* Estimate volume */

/*
 * Function: volume
 * ----------------
 * 
 * Calculate volume of cavities
 * 
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * ncav: number of cavities
 * step: 3D grid spacing (A)
 * volumes: empty array of volumes
 * nthreads: number of threads for OpenMP
 * 
 */
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

/* Spatial characterization */

/*
 * Function: _spatial
 * -----------------
 * 
 * Spatial characterization (volume and area) of the detected cavities
 * 
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * surface: surface points 3D grid
 * size: number of voxels
 * volumes: empty array of volumes
 * nvol: size of array of volumes
 * areas: empty array of areas
 * narea: size of array of areas
 * step: 3D grid spacing (A)
 * nthreads: number of threads for OpenMP
 * verbose: print extra information to standard output
 * 
 * returns: surface (surface points 3D grid), volumes (array of volumes) and area (array of areas)
 */
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
        }

        #pragma omp section
        {
            if (verbose)
                fprintf (stdout, "> Estimating area\n");
            area (surface, nx, ny, nz, narea, step, areas, nthreads);
        }
    }    
}

/* Bulk-cavity boundary points */

/*
 * Struct: points
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
 * Function: remove_boundary
 * -------------------------
 * 
 * Inspect cavities 3D grid and unmark detected cavity-bulk boundary points
 * 
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * nthreads: number of threads for OpenMP
 * 
 */
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
                                // Untag cavity-bulk boundary points
                                cavities[k + nz * (j + ( ny * i ) )] = abs(cavities[k + nz * (j + ( ny * i ) )]);
}

/* Estimate depth */

/*
 * Function: define_depth
 * ----------------------
 * 
 * Calculate depth of each cavity point and maximum and average depth of cavities
 * 
 * cavities: cavities 3D grid
 * depths: depth of cavities 3D grid points
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * max_depth: empty array of maximum depths
 * avg_depth: empty array of average depths
 * ncav: number of cavities
 * step: 3D grid spacing (A)
 * cavs: Array with minimum and maximum grid coordinates of each cavity
 * boundaries: Array with minimum and maximum grid coordinates of each cavity-bulk boundary
 * step: 3D grid spacing (A)
 * nthreads: number of threads for OpenMP
 * 
 */
void
estimate_depth(int *cavities, double *depths, int nx, int ny, int nz, double *max_depth, double *avg_depth, int ncav, pts *cavs, pts *boundaries, double step, int nthreads)
{
    int i, j, k, i2, j2, k2, count, tag;
    double distance, tmp;

    // Set number of threads in OpenMP
    omp_set_num_threads(nthreads);
    omp_set_nested(1);

    #pragma omp parallel default(none), shared(cavities, depths, max_depth, avg_depth, cavs, boundaries, ncav, nx, ny, nz), private(tmp, tag, i, j, k, i2, j2, k2, distance, count)
    {
        #pragma omp for schedule(dynamic)
            for (tag=0; tag < ncav; tag++)
            {
                // Initialize depths for cavity tag
                max_depth[tag] = 0.0;
                avg_depth[tag] = 0.0;

                // Initialize number of cavity points in cavity
                count = 0;
                
                for (i=cavs[tag].X1; i<=cavs[tag].X2; i++)
                    for (j=cavs[tag].Y1; j<=cavs[tag].Y2; j++)
                        for (k=cavs[tag].Z1; k<=cavs[tag].Z2; k++)
                            if ( abs(cavities[k + nz * (j + ( ny * i ) )]) == (tag + 2) )
                            {   
                                // Initialize tmp depth value
                                tmp = sqrt( pow(nx, 2) + pow(ny, 2) + pow(nz, 2));

                                // Count cavity point
                                count++;

                                if (boundaries[tag].X1 == nx && boundaries[tag].Y1 == ny && boundaries[tag].Z1 == nz && boundaries[tag].X2 == 0.0 && boundaries[tag].Y2 == 0.0 && boundaries[tag].Z2 == 0.0)
                                {
                                    // Cavity without boundary (void)
                                    tmp = 0.0;
                                }
                                else
                                {
                                    // Depth is calculate as the minimum distance between cavity point and bulk-cavity boundary
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
                                depths[k + nz * (j + ( ny * i ) )] = tmp;

                                // Save maximum depth for cavity tag
                                if (tmp > max_depth[tag])
                                    max_depth[tag] = tmp;

                                // Add cavity point depth to average depth for cavity tag
                                avg_depth[tag] += tmp;
                            }
                // Divide sum of depths by number of cavity points for cavity tag
                avg_depth[tag] /= count;
            }
    }
    
    // Multiply depth values by grid spacing (step)
    for (i=0; i<ncav; i++)
    {
        max_depth[i] *= step;
        avg_depth[i] *= step;
    }
}

/* Depth characterization */

/*
 * Function: _depth
 * ----------------
 * 
 * Characterization of the depth of the detected cavities. Calculate depth per cavity point and maximum and average depths of detected cavities.
 * 
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * depth: depth 3D grid
 * size: number of voxels in 3D grid
 * max_depth: empty array of maximum depths
 * nmax: size of array of maximum depths
 * avg_depth: empty array of average depths
 * navg: size of array of average depths
 * step: 3D grid spacing (A)
 * nthreads: number of threads for OpenMP
 * verbose: print extra information to standard output
 * 
 * returns: depths (3D grid with depth of cavity points), max_depth (array of maximum depths) and avg_depth (array of average depths)
 * 
 */
void 
_depth (int *cavities, int nx, int ny, int nz, double *depths, int size, double *max_depth, int nmax, double *avg_depth, int navg, double step, int nthreads, int verbose)
{
    int i, ncav;
    pts *cavs, *boundaries;

    // Get number of cavities
    ncav = nmax;

    // Fill depth 3D grid with 0.0
    dgrid (depths, size);

    // Allocate memory
    cavs = (pts*) calloc (ncav, sizeof(pts));
    boundaries = (pts*) calloc (ncav, sizeof(pts));

    // Initialize boundaries and cavs points
    for (i=0; i<ncav; i++)
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
    filter_boundary (cavities, nx, ny, nz, cavs, boundaries, nthreads);

    if (verbose)
        fprintf (stdout, "> Estimating depth\n");
    estimate_depth (cavities, depths, nx, ny, nz, max_depth, avg_depth, ncav, cavs, boundaries, step, nthreads);

    // Untag bulk-cavity boundary points
    remove_boundary (cavities, nx, ny, nz, ncav, boundaries, nthreads);

    // Free pts
    free (cavs);
    free (boundaries);
}

/* Retrieve interface residues */

/*
 * Struct: node
 * ------------
 * 
 * A linked list node for atom index in xyzr array
 * 
 * pos: atom index in xyzr array (coordinates and radii of pdb)
 * struct node* next: pointer to next linked list node
 *  
 */
typedef struct node 
{ 
    int pos; 
    struct node* next; 
} res;

/*
 * Function: create
 * ----------------
 * 
 * Create a res node
 * 
 * pos: atom index in xyzr array (coordinates and radii of pdb)
 * 
 * returns: res node with atom index
 */
res* 
create (int pos)
{
    res* new = (res*) malloc (sizeof(res));

    new->pos = pos;
    new->next = NULL;

    return new; 
}

/*
 * Function: insert
 * ----------------
 * 
 * Insert res node in linked list
 * 
 * res: pointer to linked list head
 * new: res node
 * 
 */
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

/*
 * Function: interface
 * -------------------
 * 
 * Retrieve interface residues surrounding cavities
 * 
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * pdb: 1D-array of residues information (resnum_chain)
 * atoms: xyz coordinates and radii of input pdb
 * natoms: number of atoms
 * xyzr: number of data per atom (4: xyzr)
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * probe_in: Probe In size (A)
 * ncav: number of cavities
 * nthreads: number of threads for OpenMP
 * 
 * returns: array of strings with interface residues with cavities separated by '-1'
 */
char
**interface (int *cavities, int nx, int ny, int nz, char **pdb, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, int ncav, int nthreads)
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
                            if (abs(cavities[ k + nz * (j + ( ny * i ) ) ]) > 1)
                            {
                                tag = cavities[ k + nz * (j + ( ny * i ) ) ] - 2;
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

/* Constitutional characterization */

/*
 * Function: _constitutional
 * -------------------------
 * 
 * Constitutional characterization (interface residues) of the detected cavities
 * 
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * pdb: 1D-array of residues information (resnum_chain)
 * atoms: xyz coordinates and radii of input pdb
 * natoms: number of atoms
 * xyzr: number of data per atom (4: xyzr)
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * probe_in: Probe In size (A)
 * ncav: number of cavities
 * nthreads: number of threads for OpenMP
 * verbose: print extra information to standard output
 * 
 * returns: array of strings with interface residues with cavities separated by '-1'
 */
char
** _constitutional (int *cavities, int nx, int ny, int nz, char **pdb, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, int ncav, int nthreads, int verbose)
{
    char **residues;

    if (verbose)
        fprintf (stdout, "> Retrieving interface residues\n");
    residues = interface(cavities, nx, ny, nz, pdb, atoms, natoms, xyzr, reference, ndims, sincos, nvalues, step, probe_in, ncav, nthreads);

    return residues;
}

/* Export cavity PDB */

/*
 * Function: _export
 * -----------------
 * 
 * Export cavities to PDB file
 * 
 * fn: cavity PDB filename
 * cavities: cavities 3D grid
 * nx: x grid units (cavities)
 * ny: y grid units (cavities)
 * nz: z grid units (cavities)
 * surf: surface points 3D grid
 * nxx: x grid units (surf)
 * nyy: y grid units (surf)
 * nzz: z grid units (surf)
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * ncav: number of cavities
 * nthreads: number of threads for OpenMP
 * append: append cavities to PDB file
 * 
 */
void
_export (char *fn, int *cavities, int nx, int ny, int nz, int *surf, int nxx, int nyy, int nzz, double *reference, int ndims, double *sincos, int nvalues, double step, int ncav, int nthreads, int append)
{
	int i, j, k, count, tag;
	double x, y, z, xaux, yaux, zaux;
	FILE *output;

    // Set number of threads in OpenMP
    omp_set_num_threads(nthreads);
    omp_set_nested(1);

	// Open cavity PDB file
    if (append)
        output = fopen (fn, "a+");
    else
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

/*
 * Function: _export_b
 * -------------------
 * 
 * Export cavities with a variable as B-factor to PDB file
 * 
 * fn: cavity pdb filename
 * cavities: cavities 3D grid
 * nx: x grid units (cavities)
 * ny: y grid units (cavities)
 * nz: z grid units (cavities)
 * surf: surface points 3D grid
 * nxx: x grid units (surf)
 * nyy: y grid units (surf)
 * nzz: z grid units (surf)
 * B: b-factor 3D grid (depths or hydropathy)
 * nxx: x grid units (B)
 * nyy: y grid units (B)
 * nzz: z grid units (B)
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * ncav: number of cavities
 * nthreads: number of threads for OpenMP
 * append: append cavities to PDB file
 * 
 */
void 
_export_b (char *fn, int *cavities, int nx, int ny, int nz, int *surf, int nxx, int nyy, int nzz, double *B, int nxxx, int nyyy, int nzzz, double *reference, int ndims, double *sincos, int nvalues, double step, int ncav, int nthreads, int append)
{
	int i, j, k, count, tag;
	double x, y, z, xaux, yaux, zaux;
	FILE *output;

    // Set number of threads in OpenMP
    omp_set_num_threads(nthreads);
    omp_set_nested(1);

	// Open cavity PDB file
    if (append)
        output = fopen (fn, "a+");
    else
	    output = fopen (fn, "w");

    for (count=1, tag=2; tag<=ncav+2; tag++)
        #pragma omp parallel default(none) shared(cavities, surf, B, reference, sincos, step, ncav, tag, count, nx, ny, nz, output), private(i, j, k, x, y, z, xaux, yaux, zaux)
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
                                fprintf (output, "ATOM  %5.d  HA  K%c%c   259    %8.3lf%8.3lf%8.3lf  1.00%6.2lf\n", count % 100000, 65 + (((surf[k + nz * (j + ( ny * i ) )]-2) / 26) % 26), 65 + ((surf[k + nz * (j + ( ny * i ) )]-2) % 26), xaux, yaux, zaux, B[k + nz * (j + ( ny * i ) )]);
                            else
                                fprintf (output, "ATOM  %5.d  H   K%c%c   259    %8.3lf%8.3lf%8.3lf  1.00%6.2lf\n", count % 100000, 65 + (((cavities[k + nz * (j + ( ny * i ) )]-2) / 26) % 26), 65 + ((cavities[k + nz * (j + ( ny * i ) )]-2) % 26), xaux, yaux, zaux, B[k + nz * (j + ( ny * i ) )]);
                            count++;
                        }
                    }
        }
        fclose (output);
}
