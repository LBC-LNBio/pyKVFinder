#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pyKVFinder.h"

#define min(x, y) (((x) < (y)) ? (x) : (y))
#define max(x, y) (((x) > (y)) ? (x) : (y))

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
void igrid(int *grid, int size)
{
  int i;

  for (i = 0; i < size; i++)
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
void fgrid(float *grid, int size)
{
  int i;

  for (i = 0; i < size; i++)
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
void dgrid(double *grid, int size)
{
  int i;

  for (i = 0; i < size; i++)
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
void cgrid(int *grid, int size)
{
  int i;

  for (i = 0; i < size; i++)
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
void fill(int *grid, int nx, int ny, int nz, double *atoms, int natoms,
          int xyzr, double *reference, int ndims, double *sincos, int nvalues,
          double step, double probe, int nthreads)
{
  int i, j, k, atom;
  double x, y, z, xaux, yaux, zaux, distance, H;

  // Set number of processes in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

#pragma omp parallel default(none),                                         \
    shared(grid, reference, step, probe, natoms, nx, ny, nz, sincos, atoms, \
           nthreads),                                                       \
    private(atom, i, j, k, distance, H, x, y, z, xaux, yaux, zaux)
  {
#pragma omp for schedule(dynamic) // nowait
    for (atom = 0; atom < natoms; atom++)
    {
      // Convert atom coordinates in 3D grid coordinates
      x = (atoms[atom * 4] - reference[0]) / step;
      y = (atoms[1 + (atom * 4)] - reference[1]) / step;
      z = (atoms[2 + (atom * 4)] - reference[2]) / step;

      xaux = x * sincos[3] + z * sincos[2];
      yaux = y;
      zaux = (-x) * sincos[2] + z * sincos[3];

      x = xaux;
      y = yaux * sincos[1] - zaux * sincos[0];
      z = yaux * sincos[0] + zaux * sincos[1];

      // Create a radius (H) for space occupied by probe and atom
      H = (probe + atoms[3 + (atom * 4)]) / step;

      // Loop around radius from atom center
      for (i = floor(x - H); i <= ceil(x + H); i++)
        for (j = floor(y - H); j <= ceil(y + H); j++)
          for (k = floor(z - H); k <= ceil(z + H); k++)
          {
            // Get distance between atom center and point inspected
            distance = sqrt(pow(i - x, 2) + pow(j - y, 2) + pow(k - z, 2));
            if (distance < H)
              if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz)
                grid[k + nz * (j + (ny * i))] = 0;
          }
    }
  }
}

/*
 * Function: _fill_cavity
 * ----------------------
 *
 * Insert cavities inside a 3D grid
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
void _fill_cavity(int *cavities, int nx, int ny, int nz, double *atoms,
                  int natoms, int xyzr, double *reference, int ndims,
                  double *sincos, int nvalues, double step, int nthreads)
{
  int i, j, k, atom;
  double x, y, z, xaux, yaux, zaux;

  // Set number of processes in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

#pragma omp parallel default(none),                                      \
    shared(cavities, reference, step, natoms, nx, ny, nz, sincos, atoms, \
           nthreads),                                                    \
    private(atom, i, j, k, x, y, z, xaux, yaux, zaux)
  {
#pragma omp for schedule(dynamic) // nowait
    for (atom = 0; atom < natoms; atom++)
    {
      // Convert atom coordinates in 3D grid coordinates
      x = (atoms[atom * 4] - reference[0]) / step;
      y = (atoms[1 + (atom * 4)] - reference[1]) / step;
      z = (atoms[2 + (atom * 4)] - reference[2]) / step;

      xaux = x * sincos[3] + z * sincos[2];
      yaux = y;
      zaux = (-x) * sincos[2] + z * sincos[3];

      i = (int)(xaux);
      j = (int)(yaux * sincos[1] - zaux * sincos[0]);
      k = (int)(yaux * sincos[0] + zaux * sincos[1]);

      cavities[k + nz * (j + (ny * i))] = (int)(atoms[3 + (atom * 4)]);
    }
  }

#pragma omp parallel default(none), shared(cavities, nx, ny, nz), \
    private(i, j, k)
  {
#pragma omp for collapse(3) ordered
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
          if (cavities[k + nz * (j + (ny * i))] == 1)
            cavities[k + nz * (j + (ny * i))] = -1;
  }
}

/*
 * Function: _fill_receptor
 * ------------------------
 *
 * Insert atoms with a SES or SAS representation into a 3D grid
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
 * probe_in: Probe In size (A)
 * is_ses: surface mode (1: SES or 0: SAS)
 * nthreads: number of threads for OpenMP
 * verbose: print extra information to standard output
 *
 */
void _fill_receptor(int *receptor, int size, int nx, int ny, int nz,
                    double *atoms, int natoms, int xyzr, double *reference,
                    int ndims, double *sincos, int nvalues, double step,
                    double probe_in, int is_ses, int nthreads, int verbose)
{

  // Set number of processes in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

  if (verbose)
    fprintf(stdout, "> Creating a SAS representation of the receptor\n");
  igrid(receptor, size);
  fill(receptor, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos,
       nvalues, step, probe_in, nthreads);

  if (is_ses)
  {
    if (verbose)
      fprintf(stdout, "> Adjusting a SES representation of the receptor\n");
    ses(receptor, nx, ny, nz, step, probe_in, nthreads);
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
int check_protein_neighbours(int *grid, int nx, int ny, int nz, int i, int j,
                             int k)
{
  int x, y, z;

  // Loop around neighboring points
  for (x = i - 1; x <= i + 1; x++)
    for (y = j - 1; y <= j + 1; y++)
      for (z = k - 1; z <= k + 1; z++)
      {
        // Check if point is inside 3D grid
        if (x < 0 || y < 0 || z < 0 || x > nx - 1 || y > ny - 1 || z > nz - 1)
          ;
        else if (grid[z + nz * (y + (ny * x))] == 0 ||
                 grid[z + nz * (y + (ny * x))] == -2)
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
void ses(int *grid, int nx, int ny, int nz, double step, double probe,
         int nthreads)
{
  int i, j, k, i2, j2, k2, aux;
  double distance;

  // Calculate sas limit in 3D grid units
  aux = ceil(probe / step);

  // Set number of processes in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

#pragma omp parallel default(none),             \
    shared(grid, step, probe, aux, nx, ny, nz), \
    private(i, j, k, i2, j2, k2, distance)
  {
#pragma omp for schedule(dynamic) collapse(3) // nowait
    // Loop around 3D grid
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
        {
          // Check if a cavity point
          if (grid[k + nz * (j + (ny * i))] == 1)
            if (check_protein_neighbours(grid, nx, ny, nz, i, j, k))
            {
              // Loop around sas limit from cavity point next to protein point
              for (i2 = i - aux; i2 <= i + aux; i2++)
                for (j2 = j - aux; j2 <= j + aux; j2++)
                  for (k2 = k - aux; k2 <= k + aux; k2++)
                  {
                    if (i2 > 0 && j2 > 0 && k2 > 0 && i2 < nx && j2 < ny &&
                        k2 < nz)
                    {
                      // Get distance between point inspected and cavity point
                      distance = sqrt(pow(i - i2, 2) + pow(j - j2, 2) +
                                      pow(k - k2, 2));
                      // Check if inspected point is inside sas limit
                      if (distance < (probe / step))
                        if (grid[k2 + nz * (j2 + (ny * i2))] == 0)
                          // Mark cavity point
                          grid[k2 + nz * (j2 + (ny * i2))] = -2;
                    }
                  }
            }
        }

#pragma omp for collapse(3)
    // Loop around 3D grid
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
        {
          // Mark space occupied by sas limit from protein surface
          if (grid[k + nz * (j + (ny * i))] == -2)
            grid[k + nz * (j + (ny * i))] = 1;
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
void subtract(int *PI, int *PO, int nx, int ny, int nz, double step,
              double removal_distance, int nthreads)
{
  int i, j, k, i2, j2, k2, rd;

  rd = ceil(removal_distance / step);

  // Set number of processes in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

// Create a parallel region */
#pragma omp parallel default(none),                                  \
    shared(PI, PO, nx, ny, nz, i, j, k, step, rd, removal_distance), \
    private(j2, i2, k2)
  {
#pragma omp for schedule(dynamic) collapse(3)
    // Loop around the search box
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
        {
          // Check if point is a cavity point in grid filled with PO
          if (PO[k + nz * (j + (ny * i))])
          {
            // Loops around space occupied by probe from atom position
            for (i2 = i - rd; i2 <= i + rd; i2++)
              for (j2 = j - rd; j2 <= j + rd; j2++)
                for (k2 = k - rd; k2 <= k + rd; k2++)
                  // Check if inside 3D grid
                  if (i2 >= 0 && i2 < nx && j2 >= 0 && j2 < ny && k2 >= 0 &&
                      k2 < nz)
                    // Mark points in grid filled with Probe In, where Probe Out
                    // reached
                    if (PI[k2 + nz * (j2 + (ny * i2))] == 1)
                      PI[k2 + nz * (j2 + (ny * i2))] = -1;
          }
        }
  }
}

/* Filter noise from Grid */

/*
 * Function: filter_noise
 * ----------------------
 *
 * Removes cavities points (1) surrounded by biomolecule (0) or medium (-1)
 * points
 *
 * grid: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * nthreads: number of threads for OpenMP
 *
 */
void filter_noise(int *grid, int nx, int ny, int nz, int nthreads)
{
  int i, j, k, contacts;

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++)
      {
        if (grid[k + nz * (j + (ny * i))] == 1)
        {
          contacts = 0;

          // Check if a protein point (0) or a medium point (-1) is next to a
          // cavity point (==1)
          if (i - 1 >= 0 && i + 1 < nx && j - 1 >= 0 && j + 1 < ny &&
              k - 1 >= 0 && k + 1 < nz)
          {
            if (grid[k + nz * (j + (ny * (i - 1)))] == 0 ||
                grid[k + nz * (j + (ny * (i - 1)))] == -1)
              contacts++;
            if (grid[k + nz * (j + (ny * (i + 1)))] == 0 ||
                grid[k + nz * (j + (ny * (i + 1)))] == -1)
              contacts++;
            if (grid[k + nz * ((j - 1) + (ny * i))] == 0 ||
                grid[k + nz * ((j - 1) + (ny * i))] == -1)
              contacts++;
            if (grid[k + nz * ((j + 1) + (ny * i))] == 0 ||
                grid[k + nz * ((j + 1) + (ny * i))] == -1)
              contacts++;
            if (grid[(k - 1) + nz * (j + (ny * i))] == 0 ||
                grid[(k - 1) + nz * (j + (ny * i))] == -1)
              contacts++;
            if (grid[(k + 1) + nz * (j + (ny * i))] == 0 ||
                grid[(k + 1) + nz * (j + (ny * i))] == -1)
              contacts++;

            // Cavity point is a medium point
            if (contacts == 6)
              grid[k + nz * (j + (ny * i))] = -1;
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
void adjust(int *grid, int nx, int ny, int nz, double *ligand, int lnatoms,
            int lxyzr, double *reference, int ndims, double *sincos,
            int nvalues, double step, double ligand_cutoff, int nthreads)
{
  int i, j, k, atom, inside;
  double x, y, z, xaux, yaux, zaux, distance;

  // Set number of processes in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

#pragma omp parallel default(none),                                          \
    shared(grid, nx, ny, nz, step, sincos, reference, ligand, ligand_cutoff, \
           lnatoms),                                                         \
    private(inside, i, j, k, x, y, z, xaux, yaux, zaux, atom, distance)
  {
#pragma omp for collapse(3) schedule(static) nowait
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
        {
          inside = 0;
          for (atom = 0; atom < lnatoms; atom++)
          {
            // Get 3D grid point coordinate
            x = i * step;
            y = j * step;
            z = k * step;

            xaux = (x * sincos[3]) + (y * sincos[0] * sincos[2]) -
                   (z * sincos[1] * sincos[2]) + reference[0];
            yaux = (y * sincos[1]) + (z * sincos[0]) + reference[1];
            zaux = (x * sincos[2]) - (y * sincos[0] * sincos[3]) +
                   (z * sincos[1] * sincos[3]) + reference[2];

            // Get distance between ligand and 3D grid point evaluated
            distance = sqrt(pow(xaux - ligand[atom * 4], 2) +
                            pow(yaux - ligand[1 + (atom * 4)], 2) +
                            pow(zaux - ligand[2 + (atom * 4)], 2));

            if (distance < ligand_cutoff)
              inside = 1;
          }
          if (inside == 0 && grid[k + nz * (j + (ny * i))])
            grid[k + nz * (j + (ny * i))] = -1;
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
void _filter_pdb(int nx, int ny, int nz, double *atoms, int natoms, int xyzr,
                 double *reference, int ndims, double *sincos, int nvalues,
                 double step, double probe_in, int nthreads)
{
  int atom;
  double x, y, z, xaux, yaux, zaux;

  // Set number of processes in OpenMP
  omp_set_num_threads(nthreads);

#pragma omp parallel default(none),                                       \
    shared(atoms, natoms, reference, sincos, step, probe_in, nx, ny, nz), \
    private(atom, x, y, z, xaux, yaux, zaux)
  {
#pragma omp for schedule(static) // nowait
    for (atom = 0; atom < natoms; atom++)
    {
      x = (atoms[atom * 4] - reference[0]) / step;
      y = (atoms[1 + (atom * 4)] - reference[1]) / step;
      z = (atoms[2 + (atom * 4)] - reference[2]) / step;

      xaux = x * sincos[3] + z * sincos[2];
      yaux = y;
      zaux = (-x) * sincos[2] + z * sincos[3];

      x = xaux;
      y = yaux * sincos[1] - zaux * sincos[0];
      z = yaux * sincos[0] + zaux * sincos[1];

      if (x > 0.0 - (probe_in + atoms[3 + (atom * 4)]) / step &&
          x < (double)nx + (probe_in + atoms[3 + (atom * 4)]) / step &&
          y > 0.0 - (probe_in + atoms[3 + (atom * 4)]) / step &&
          y < (double)ny + (probe_in + atoms[3 + (atom * 4)]) / step &&
          z > 0.0 - (probe_in + atoms[3 + (atom * 4)]) / step &&
          z < (double)nz + (probe_in + atoms[3 + (atom * 4)]) / step)
        ;
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
void filter(int *grid, int nx, int ny, int nz, double *P1, int ndims,
            double *P2, int nndims, double *sincos, int nvalues, double step,
            double probe_out, int nthreads)
{
  int i, j, k;
  double aux, normB, norm1, X1, Y1, Z1, X2, Y2, Z2;

  // Set number of processes in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

  // Get norm between X1 and X2 in 3D grid
  X1 = P1[0];
  Y1 = P1[1];
  Z1 = P1[2];
  X2 = P2[0];
  Y2 = P2[1];
  Z2 = P2[2];
  norm1 = sqrt(pow(X2 - X1, 2) + pow(Y2 - Y1, 2) + pow(Z2 - Z1, 2));

  // Remove Probe out from grid to recreate box
  X1 -= (-(probe_out * sincos[3]) - (probe_out * sincos[0] * sincos[2]) +
         (probe_out * sincos[1] * sincos[2]));
  Y1 -= (-(probe_out * sincos[1]) - (probe_out * sincos[0]));
  Z1 -= (-(probe_out * sincos[2]) + (probe_out * sincos[0] * sincos[3]) -
         (probe_out * sincos[1] * sincos[3]));
  X2 -= ((probe_out * sincos[3]) - (probe_out * sincos[0] * sincos[2]) +
         (probe_out * sincos[1] * sincos[2]));
  Y2 -= (-(probe_out * sincos[1]) - (probe_out * sincos[0]));
  Z2 -= ((probe_out * sincos[2]) + (probe_out * sincos[0] * sincos[3]) -
         (probe_out * sincos[1] * sincos[3]));

  // Get norm between X1 and X2 in box
  normB = sqrt(pow(X2 - X1, 2) + pow(Y2 - Y1, 2) + pow(Z2 - Z1, 2));

  // Prepare grid units
  aux = floor((norm1 - normB) / (2 * step));

#pragma omp parallel default(shared), private(i, j, k)
  {
    for (i = 0; i <= aux; i++)
#pragma omp for collapse(2) nowait
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
          grid[k + nz * (j + (ny * i))] = -1;

    for (i = nx - 1; i >= nx - aux - 1; i--)
#pragma omp for collapse(2) nowait
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
          grid[k + nz * (j + (ny * i))] = -1;

    for (j = 0; j <= aux; j++)
#pragma omp for collapse(2) nowait
      for (i = 0; i < nx; i++)
        for (k = 0; k < nz; k++)
          grid[k + nz * (j + (ny * i))] = -1;

    for (j = ny - 1; j >= ny - aux - 1; j--)
#pragma omp for collapse(2) nowait
      for (i = 0; i < nx; i++)
        for (k = 0; k < nz; k++)
          grid[k + nz * (j + (ny * i))] = -1;

    for (k = 0; k <= aux; k++)
#pragma omp for collapse(2) nowait
      for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
          grid[k + nz * (j + (ny * i))] = -1;

    for (k = nz - 1; k >= nz - aux - 1; k--)
#pragma omp for collapse(2) nowait
      for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
          grid[k + nz * (j + (ny * i))] = -1;
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
 * Variable: big
 *
 * Flag that marks big cavities on clustering
 *
 */
int big;

/*
 * Function: check_unclustered_neighbours
 * --------------------------------------
 *
 * Checks if a cavity point on the grid is next to a unclustered cavity point
 * (1)
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
int check_unclustered_neighbours(int *grid, int nx, int ny, int nz, int i,
                                 int j, int k)
{
  int x, y, z;

  // Loop around neighboring points
  for (x = i - 1; x <= i + 1; x++)
    for (y = j - 1; y <= j + 1; y++)
      for (z = k - 1; z <= k + 1; z++)
      {
        // Check if point is inside 3D grid
        if (x < 0 || y < 0 || z < 0 || x > nx - 1 || y > ny - 1 || z > nz - 1)
          ;
        else if (grid[z + nz * (y + (ny * x))] > 1)
          return grid[z + nz * (y + (ny * x))];
      }
  return 0;
}

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
void DFS(int *grid, int nx, int ny, int nz, int i, int j, int k, int tag)
{
  int x, y, z;

  if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1 || k == 0 || k == nz - 1)
    return;

  if (grid[k + nz * (j + (ny * i))] == 1 && !big)
  {
    grid[k + nz * (j + (ny * i))] = tag;
    vol++;

    if (vol == 10000)
      big = 1;

    if (!big)
    {
      for (x = i - 1; x <= i + 1; x++)
        for (y = j - 1; y <= j + 1; y++)
          for (z = k - 1; z <= k + 1; z++)
            DFS(grid, nx, ny, nz, x, y, z, tag);
    }
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
void remove_cavity(int *grid, int nx, int ny, int nz, int tag, int nthreads)
{
  int i, j, k;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

#pragma omp parallel default(shared)
#pragma omp for schedule(static) collapse(3) // nowait
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++)
        // Remove points based on tag
        if (grid[k + nz * (j + (ny * i))] == tag)
          grid[k + nz * (j + (ny * i))] = -1;
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
int _cluster(int *grid, int nx, int ny, int nz, double step,
             double volume_cutoff, int nthreads)
{
  int i, j, k, i2, j2, k2, tag, vol_aux;

  // Initialize variables
  tag = 1;
  vol_aux = 0;
  big = 0;

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++)
        if (grid[k + nz * (j + (ny * i))] == 1)
        {
          tag++;
          vol = 0;

          // Clustering procedure
          DFS(grid, nx, ny, nz, i, j, k, tag);
          vol_aux = vol;

          // Loop for big cavities
          while (big)
          {
            vol_aux = 0;

            for (i2 = 0; i2 < nx; i2++)
              for (j2 = 0; j2 < ny; j2++)
                for (k2 = 0; k2 < nz; k2++)
                {
                  big = 0;
                  vol_aux += vol;
                  vol = 0;
                  if (grid[k2 + nz * (j2 + (ny * i2))] == 1 &&
                      check_unclustered_neighbours(grid, nx, ny, nz, i2, j2,
                                                   k2) == tag)
                    DFS(grid, nx, ny, nz, i2, j2, k2, tag);
                }
          }
          vol = vol_aux;

          // Check if cavity reached cutoff
          if ((double)vol * pow(step, 3) < volume_cutoff)
          {
            remove_cavity(grid, nx, ny, nz, tag, nthreads);
            tag--;
          }
        }
  return tag - 1;
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
int _detect(int *PI, int size, int nx, int ny, int nz, double *atoms,
            int natoms, int xyzr, double *reference, int ndims, double *sincos,
            int nvalues, double step, double probe_in, double probe_out,
            double removal_distance, double volume_cutoff, int box_adjustment,
            double *P2, int nndims, int is_ses, int nthreads, int verbose)
{
  int *PO, ncav;

  if (verbose)
    fprintf(stdout, "> Filling grid with Probe In\n");
  igrid(PI, size);
  fill(PI, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues,
       step, probe_in, nthreads);

  if (verbose)
    fprintf(stdout, "> Filling grid with Probe Out\n");
  PO = (int *)calloc(size, sizeof(int));
  igrid(PO, size);
  fill(PO, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues,
       step, probe_out, nthreads);

  if (is_ses)
    ses(PI, nx, ny, nz, step, probe_in, nthreads);
  ses(PO, nx, ny, nz, step, probe_out, nthreads);

  if (verbose)
    fprintf(stdout, "> Defining biomolecular cavities\n");
  subtract(PI, PO, nx, ny, nz, step, removal_distance, nthreads);

  if (box_adjustment)
  {
    if (verbose)
      fprintf(stdout, "> Adjusting biomolecular cavities to box\n");
    filter(PI, nx, ny, nz, reference, ndims, P2, nndims, sincos, nvalues, step,
           probe_out, nthreads);
  }

  filter_noise(PI, nx, ny, nz, nthreads);

  if (verbose)
    fprintf(stdout, "> Clustering cavity points\n");
  ncav = _cluster(PI, nx, ny, nz, step, volume_cutoff, nthreads);

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
int _detect_ladj(int *PI, int size, int nx, int ny, int nz, double *atoms,
                 int natoms, int xyzr, double *ligand, int lnatoms, int lxyzr,
                 double *reference, int ndims, double *sincos, int nvalues,
                 double step, double probe_in, double probe_out,
                 double removal_distance, double volume_cutoff,
                 int ligand_adjustment, double ligand_cutoff,
                 int box_adjustment, double *P2, int nndims, int is_ses,
                 int nthreads, int verbose)
{
  int *PO, ncav;

  if (verbose)
    fprintf(stdout, "> Filling grid with Probe In\n");
  igrid(PI, size);
  fill(PI, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues,
       step, probe_in, nthreads);

  if (verbose)
    fprintf(stdout, "> Filling grid with Probe Out\n");
  PO = (int *)calloc(size, sizeof(int));
  igrid(PO, size);
  fill(PO, nx, ny, nz, atoms, natoms, xyzr, reference, ndims, sincos, nvalues,
       step, probe_out, nthreads);

  if (is_ses)
    ses(PI, nx, ny, nz, step, probe_in, nthreads);
  ses(PO, nx, ny, nz, step, probe_out, nthreads);

  if (verbose)
    fprintf(stdout, "> Defining biomolecular cavities\n");
  subtract(PI, PO, nx, ny, nz, step, removal_distance, nthreads);

  if (ligand_adjustment)
  {
    if (verbose)
      fprintf(stdout, "> Adjusting biomolecular cavities to ligand\n");
    adjust(PI, nx, ny, nz, ligand, lnatoms, lxyzr, reference, ndims, sincos,
           nvalues, step, ligand_cutoff, nthreads);
  }

  if (box_adjustment)
  {
    if (verbose)
      fprintf(stdout, "> Adjusting biomolecular cavities to box\n");
    filter(PI, nx, ny, nz, reference, ndims, P2, nndims, sincos, nvalues, step,
           probe_out, nthreads);
  }

  filter_noise(PI, nx, ny, nz, nthreads);

  if (verbose)
    fprintf(stdout, "> Clustering cavity points\n");
  ncav = _cluster(PI, nx, ny, nz, step, volume_cutoff, nthreads);

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
int define_surface_points(int *cavities, int nx, int ny, int nz, int i, int j,
                          int k)
{
  if (i - 1 >= 0)
    if (cavities[k + nz * (j + (ny * (i - 1)))] == 0)
      return cavities[k + nz * (j + (ny * i))];
  if (i + 1 < nx)
    if (cavities[k + nz * (j + (ny * (i + 1)))] == 0)
      return cavities[k + nz * (j + (ny * i))];
  if (j - 1 >= 0)
    if (cavities[k + nz * ((j - 1) + (ny * i))] == 0)
      return cavities[k + nz * (j + (ny * i))];
  if (j + 1 < ny)
    if (cavities[k + nz * ((j + 1) + (ny * i))] == 0)
      return cavities[k + nz * (j + (ny * i))];
  if (k - 1 >= 0)
    if (cavities[(k - 1) + nz * (j + (ny * i))] == 0)
      return cavities[k + nz * (j + (ny * i))];
  if (k + 1 < nz)
    if (cavities[(k + 1) + nz * (j + (ny * i))] == 0)
      return cavities[k + nz * (j + (ny * i))];

  return -1;
}

/*
 * Function: filter_surface
 * ------------------------
 *
 * Inspect cavities 3D grid and mark detected surface points on a surface 3D
 * grid
 *
 * cavities: cavities 3D grid
 * surface: surface points 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * nthreads: number of threads for OpenMP
 *
 */
void filter_surface(int *cavities, int *surface, int nx, int ny, int nz,
                    int nthreads)
{
  int i, j, k;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

#pragma omp parallel default(none), shared(cavities, surface, nx, ny, nz), \
    private(i, j, k)
  {
#pragma omp for collapse(3) schedule(static)
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
          if (cavities[k + nz * (j + (ny * i))] > 1)
          {
            // Define surface cavity points
            surface[k + nz * (j + (ny * i))] =
                define_surface_points(cavities, nx, ny, nz, i, j, k);
          }
          else
          {
            if (cavities[k + nz * (j + (ny * i))] == 0)
              surface[k + nz * (j + (ny * i))] = 0;
            else
              surface[k + nz * (j + (ny * i))] = -1;
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
double check_voxel_class(int *surface, int nx, int ny, int nz, int i, int j,
                         int k)
{
  int contacts = 0;
  double weight = 1.0;

  // Count face contacts
  if (surface[k + nz * (j + (ny * (i - 1)))] == 0)
    contacts++;
  if (surface[k + nz * (j + (ny * (i + 1)))] == 0)
    contacts++;
  if (surface[k + nz * ((j - 1) + (ny * i))] == 0)
    contacts++;
  if (surface[k + nz * ((j + 1) + (ny * i))] == 0)
    contacts++;
  if (surface[(k - 1) + nz * (j + (ny * i))] == 0)
    contacts++;
  if (surface[(k + 1) + nz * (j + (ny * i))] == 0)
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
    if ((surface[k + nz * (j + (ny * (i - 1)))] == 0 &&
         surface[k + nz * (j + (ny * (i + 1)))] == 0) ||
        (surface[k + nz * ((j - 1) + (ny * i))] == 0 &&
         surface[k + nz * ((j + 1) + (ny * i))] == 0) ||
        (surface[(k - 1) + nz * (j + (ny * i))] == 0 &&
         surface[(k + 1) + nz * (j + (ny * i))] == 0))
    {
      weight = 2;
    }

    else
    {
      weight = 1.5879;
    }

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
 * Function: _area
 * ---------------
 *
 * Calculate area of cavities
 *
 * surface: surface points 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * step: 3D grid spacing (A)
 * areas: empty array of areas
 * narea: number of cavities
 * nthreads: number of threads for OpenMP
 *
 */
void _area(int *surface, int nxx, int nyy, int nzz, double step, double *areas,
           int narea, int nthreads)
{
  int i, j, k;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

  for (i = 0; i < narea; i++)
    areas[i] = 0.0;

  for (i = 0; i < nxx; i++)
    for (j = 0; j < nyy; j++)
      for (k = 0; k < nzz; k++)
        if (surface[k + nzz * (j + (nyy * i))] > 1)
          areas[surface[k + nzz * (j + (nyy * i))] - 2] +=
              check_voxel_class(surface, nxx, nyy, nzz, i, j, k) * pow(step, 2);
}

/* Estimate volume */

/*
 * Function: _volume
 * -----------------
 *
 * Calculate volume of cavities
 *
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * step: 3D grid spacing (A)
 * volumes: empty array of volumes
 * nvol: number of cavities
 * nthreads: number of threads for OpenMP
 *
 */
void _volume(int *cavities, int nx, int ny, int nz, double step,
             double *volumes, int nvol, int nthreads)
{
  int i, j, k;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

  // Initialize volumes array
  dgrid(volumes, nvol);

#pragma omp parallel default(none), \
    shared(volumes, cavities, nvol, step, nx, ny, nz), private(i, j, k)
  {
#pragma omp for collapse(3) reduction(+ \
                                      : volumes[:nvol])
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
          if (cavities[k + nz * (j + (ny * i))] > 1)
            volumes[cavities[k + nz * (j + (ny * i))] - 2] += pow(step, 3);
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
 * returns: surface (surface points 3D grid), volumes (array of volumes) and
 * area (array of areas)
 */
void _spatial(int *cavities, int nx, int ny, int nz, int *surface, int size,
              double *volumes, int nvol, double *areas, int narea, double step,
              int nthreads, int verbose)
{
  if (verbose)
    fprintf(stdout, "> Defining surface points\n");
  filter_surface(cavities, surface, nx, ny, nz, nthreads);

#pragma omp sections
  {
#pragma omp section
    {
      if (verbose)
        fprintf(stdout, "> Estimating volume\n");
      _volume(cavities, nx, ny, nz, step, volumes, nvol, nthreads);
    }

#pragma omp section
    {
      if (verbose)
        fprintf(stdout, "> Estimating area\n");
      _area(surface, nx, ny, nz, step, areas, narea, nthreads);
    }
  }
}

/* Bulk-cavity boundary points */

// /*
//  * Struct: points
//  * --------------
//  *
//  * Two 3D grid points with xyz coordinates (P1 and P2)
//  *
//  * X1: x coordinate of P1
//  * Y1: y coordinate of P1
//  * Z1: z coordinate of P1
//  * X2: x coordinate of P2
//  * Y2: y coordinate of P2
//  * Z2: z coordinate of P2
//  *
//  */
// typedef struct POINTS {
//   int X1;
//   int Y1;
//   int Z1;
//   int X2;
//   int Y2;
//   int Z2;
// } pts;

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
int define_boundary_points(int *cavities, int nx, int ny, int nz, int i, int j,
                           int k)
{
  if (i - 1 >= 0)
    if (cavities[k + nz * (j + (ny * (i - 1)))] == -1)
      return -(cavities[k + nz * (j + (ny * i))]);
  if (i + 1 < nx)
    if (cavities[k + nz * (j + (ny * (i + 1)))] == -1)
      return -(cavities[k + nz * (j + (ny * i))]);
  if (j - 1 >= 0)
    if (cavities[k + nz * ((j - 1) + (ny * i))] == -1)
      return -(cavities[k + nz * (j + (ny * i))]);
  if (j + 1 < ny)
    if (cavities[k + nz * ((j + 1) + (ny * i))] == -1)
      return -(cavities[k + nz * (j + (ny * i))]);
  if (k - 1 >= 0)
    if (cavities[(k - 1) + nz * (j + (ny * i))] == -1)
      return -(cavities[k + nz * (j + (ny * i))]);
  if (k + 1 < nz)
    if (cavities[(k + 1) + nz * (j + (ny * i))] == -1)
      return -(cavities[k + nz * (j + (ny * i))]);

  return cavities[k + nz * (j + (ny * i))];
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
void filter_boundary(int *cavities, int nx, int ny, int nz, pts *cavs,
                     pts *boundaries, int nthreads)
{
  int i, j, k, tag;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

#pragma omp parallel default(none), \
    shared(cavities, nx, ny, nz, cavs, boundaries), private(i, j, k, tag)
  {
#pragma omp for collapse(3) schedule(static)
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
#pragma omp critical
          if (cavities[k + nz * (j + (ny * i))] > 1)
          {
            // Get cavity identifier
            tag = cavities[k + nz * (j + (ny * i))] - 2;

            // Get min and max coordinates of each cavity
            cavs[tag].X1 = min(cavs[tag].X1, i);
            cavs[tag].Y1 = min(cavs[tag].Y1, j);
            cavs[tag].Z1 = min(cavs[tag].Z1, k);
            cavs[tag].X2 = max(cavs[tag].X2, i);
            cavs[tag].Y2 = max(cavs[tag].Y2, j);
            cavs[tag].Z2 = max(cavs[tag].Z2, k);

            // Define cavity-bulk boundary points
            cavities[k + nz * (j + (ny * i))] =
                define_boundary_points(cavities, nx, ny, nz, i, j, k);

            // Get min and max coordinates of each cavity-bulk boundary
            if (cavities[k + nz * (j + (ny * i))] < -1)
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
void remove_boundary(int *cavities, int nx, int ny, int nz, int ncav,
                     pts *boundaries, int nthreads)
{
  int i, j, k, tag;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

#pragma omp parallel default(none), \
    shared(cavities, boundaries, ncav, nx, ny, nz), private(tag, i, j, k)
#pragma omp for schedule(dynamic)
  for (tag = 0; tag < ncav; tag++)
    for (i = boundaries[tag].X1; i <= boundaries[tag].X2; i++)
      for (j = boundaries[tag].Y1; j <= boundaries[tag].Y2; j++)
        for (k = boundaries[tag].Z1; k <= boundaries[tag].Z2; k++)
          if (cavities[k + nz * (j + (ny * i))] < -1)
            // Untag cavity-bulk boundary points
            cavities[k + nz * (j + (ny * i))] =
                abs(cavities[k + nz * (j + (ny * i))]);
}

/* Estimate depth */

/*
 * Function: estimate_depth
 * ------------------------
 *
 * Calculate depth of each cavity point and maximum and average depth of
 * cavities
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
 * boundaries: Array with minimum and maximum grid coordinates of each
 * cavity-bulk boundary step: 3D grid spacing (A) nthreads: number of threads
 * for OpenMP
 *
 */
void estimate_depth(int *cavities, double *depths, int nx, int ny, int nz,
                    double *max_depth, double *avg_depth, int ncav, pts *cavs,
                    pts *boundaries, double step, int nthreads)
{
  int i, j, k, i2, j2, k2, count, tag;
  double distance, tmp;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

#pragma omp parallel default(none),                                            \
    shared(cavities, depths, max_depth, avg_depth, cavs, boundaries, ncav, nx, \
           ny, nz, step),                                                      \
    private(tmp, tag, i, j, k, i2, j2, k2, distance, count)
  {
#pragma omp for schedule(dynamic)
    for (tag = 0; tag < ncav; tag++)
    {
      // Initialize depths for cavity tag
      max_depth[tag] = 0.0;
      avg_depth[tag] = 0.0;

      // Initialize number of cavity points in cavity
      count = 0;

      for (i = cavs[tag].X1; i <= cavs[tag].X2; i++)
        for (j = cavs[tag].Y1; j <= cavs[tag].Y2; j++)
          for (k = cavs[tag].Z1; k <= cavs[tag].Z2; k++)
            if (abs(cavities[k + nz * (j + (ny * i))]) == (tag + 2))
            {
              // Initialize tmp depth value
              tmp = sqrt(pow(nx, 2) + pow(ny, 2) + pow(nz, 2)) * step;

              // Count cavity point
              count++;

              if (boundaries[tag].X1 == nx && boundaries[tag].Y1 == ny &&
                  boundaries[tag].Z1 == nz && boundaries[tag].X2 == 0.0 &&
                  boundaries[tag].Y2 == 0.0 && boundaries[tag].Z2 == 0.0)
              {
                // Cavity without boundary (void)
                tmp = 0.0;
              }
              else
              {
                // Depth is calculate as the minimum distance between cavity
                // point and bulk-cavity boundary
                for (i2 = boundaries[tag].X1; i2 <= boundaries[tag].X2; i2++)
                  for (j2 = boundaries[tag].Y1; j2 <= boundaries[tag].Y2; j2++)
                    for (k2 = boundaries[tag].Z1; k2 <= boundaries[tag].Z2;
                         k2++)
                      if (cavities[k2 + nz * (j2 + (ny * i2))] == -(tag + 2))
                      {
                        distance = sqrt(pow(i2 - i, 2) + pow(j2 - j, 2) +
                                        pow(k2 - k, 2)) *
                                   step;
                        if (distance < tmp)
                          tmp = distance;
                      }
              }

              // Save depth for cavity point
              depths[k + nz * (j + (ny * i))] = tmp;

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
}

/* Depth characterization */

/*
 * Function: _depth
 * ----------------
 *
 * Characterization of the depth of the detected cavities. Calculate depth per
 * cavity point and maximum and average depths of detected cavities.
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
 * returns: depths (3D grid with depth of cavity points), max_depth (array of
 * maximum depths) and avg_depth (array of average depths)
 *
 */
void _depth(int *cavities, int nx, int ny, int nz, double *depths, int size,
            double *max_depth, int nmax, double *avg_depth, int navg,
            double step, int nthreads, int verbose)
{
  int i, ncav;
  pts *cavs, *boundaries;

  // Get number of cavities
  ncav = nmax;

  // Fill depth 3D grid with 0.0
  dgrid(depths, size);

  // Allocate memory
  cavs = (pts *)calloc(ncav, sizeof(pts));
  boundaries = (pts *)calloc(ncav, sizeof(pts));

  // Initialize boundaries and cavs points
  for (i = 0; i < ncav; i++)
  {
    boundaries[i].X1 = nx;
    boundaries[i].Y1 = ny;
    boundaries[i].Z1 = nz;
    boundaries[i].X2 = 0.0;
    boundaries[i].Y2 = 0.0;
    boundaries[i].Z2 = 0.0;
    cavs[i].X1 = nx;
    cavs[i].Y1 = ny;
    cavs[i].Z1 = nz;
    cavs[i].X2 = 0.0;
    cavs[i].Y2 = 0.0;
    cavs[i].Z2 = 0.0;
  }

  if (verbose)
    fprintf(stdout, "> Defining bulk-cavity boundary points\n");
  filter_boundary(cavities, nx, ny, nz, cavs, boundaries, nthreads);

  if (verbose)
    fprintf(stdout, "> Estimating depth\n");
  estimate_depth(cavities, depths, nx, ny, nz, max_depth, avg_depth, ncav, cavs,
                 boundaries, step, nthreads);

  // Untag bulk-cavity boundary points
  remove_boundary(cavities, nx, ny, nz, ncav, boundaries, nthreads);

  // Free pts
  free(cavs);
  free(boundaries);
}

/* Openings characterization */

/*
 * Function: _openings_in_cavities
 * -------------------------------
 *
 * Cluster consecutive cavity points together
 *
 * o2c: empty array with openings as indexes and cavities as values
 * nopenings: number of openings
 * cavities: cavities 3D grid
 * nx: x grid units (cavities)
 * ny: y grid units (cavities)
 * nz: z grid units (cavities)
 * openings: openings 3D grid
 * nxx: x grid units (openings)
 * nyy: y grid units (openings)
 * nzz: z grid units (openings)
 * ncav: number of cavities
 * nthreads: number of threads for OpenMP
 *
 */
void _openings2cavities(int *o2c, int nopenings, int *cavities, int nx, int ny,
                        int nz, int *openings, int nxx, int nyy, int nzz,
                        int nthreads)
{
  int i, j, k, tag, stop;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

  for (tag = 0; tag < nopenings; tag++)
  {
#pragma omp parallel default(none), \
    shared(stop, tag, o2c, cavities, openings, nx, ny, nz), private(i, j, k)
    {
      stop = 0;
#pragma omp for collapse(3) schedule(static)
      for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
          for (k = 0; k < nz; k++)
          {
            if (stop)
              continue;
            else
            {
              if (openings[k + nz * (j + (ny * i))] == tag + 2 &&
                  cavities[k + nz * (j + (ny * i))] > 1)
              {
                o2c[tag] = cavities[k + nz * (j + (ny * i))] - 2;
                stop = 1;
              }
            }
          }
    }
  }
}

/*
 * Function: remove_enclosed_cavity
 * --------------------------------
 *
 * Cluster consecutive cavity points together
 *
 * openings: openings 3D grid
 * cavities: cavities 3D grid
 * nx: x grid units (openings/cavities)
 * ny: y grid units (openings/cavities)
 * nz: z grid units (openings/cavities)
 * depths: depth of cavities 3D grid points
 * nxx: x grid units (depths)
 * nyy: y grid units (depths)
 * nzz: z grid units (depths)
 * ncav: number of cavities
 * nthreads: number of threads for OpenMP
 *
 */
void remove_enclosed_cavity(int *openings, int *cavities, int nx, int ny,
                            int nz, double *depths, int nxx, int nyy, int nzz,
                            int ncav, int nthreads)
{
  int i, j, k, tag;
  double sum_depth;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

  // Copy cavities to openings
#pragma omp parallel default(shared), private(i, j, k)
#pragma omp for schedule(static)
  for (i = 0; i < (nx * ny * nz); i++)
    openings[i] = cavities[i];

  for (tag = 0; tag < ncav; tag++)
  {
    sum_depth = 0.0;

#pragma omp parallel default(shared), private(i, j, k)
#pragma omp for collapse(3) schedule(static) reduction(+ \
                                                       : sum_depth)
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
          if (openings[k + nz * (j + (ny * i))] == tag + 2)
            sum_depth += depths[k + nz * (j + (ny * i))];

    // Check if kvtag is a enclosed cavity
    if (sum_depth == 0.0)
      remove_cavity(openings, nx, ny, nz, tag + 2, nthreads);
  }
}

/*
 * Function: filter_openings
 * -------------------------
 *
 * Inspect openings 3D grid and mark detected openings points
 *
 * openings: openings 3D grid
 * depths: depth of cavities 3D grid points
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * nthreads: number of threads for OpenMP
 *
 */
void filter_openings(int *openings, double *depths, int nx, int ny, int nz,
                     int nthreads)
{
  int i, j, k;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

#pragma omp parallel default(none), shared(openings, depths, nx, ny, nz), \
    private(i, j, k)
  {
#pragma omp for collapse(3) schedule(static)
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
          if (openings[k + nz * (j + (ny * i))] > 1)
          {
            if (depths[k + nz * (j + (ny * i))] == 0.0)
              openings[k + nz * (j + (ny * i))] = 1;
            else
              openings[k + nz * (j + (ny * i))] = 0;
          }
          else
          {
            if (openings[k + nz * (j + (ny * i))] == 1)
              openings[k + nz * (j + (ny * i))] = -1;
            else
            {
              if (openings[k + nz * (j + (ny * i))] == 0)
                openings[k + nz * (j + (ny * i))] = 0;
              else
                openings[k + nz * (j + (ny * i))] = -1;
            }
          }
  }
}

/*
 * Function: _openings
 * -------------------
 *
 * Cavity openings characterization of the detected cavities
 *
 * openings: openings 3D grid
 * size: number of voxels in 3D grid
 * cavities: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * depths: depth of cavities 3D grid points
 * nxx: x grid units (depths)
 * nyy: y grid units (depths)
 * nzz: z grid units (depths)
 * ncav: number of cavities
 * openings_cutoff: minimum number of voxels an opening must have
 * step: 3D grid spacing (A)
 * nthreads: number of threads for OpenMP
 * verbose: print extra information to standard output
 *
 * returns: surface (surface points 3D grid), volumes (array of volumes) and
 * area (array of areas)
 */
int _openings(int *openings, int size, int *cavities, int nx, int ny, int nz,
              double *depths, int nxx, int nyy, int nzz, int ncav,
              int openings_cutoff, double step, int nthreads, int verbose)
{
  int nopenings;

  if (verbose)
    fprintf(stdout, "> Removing enclosed cavities\n");
  remove_enclosed_cavity(openings, cavities, nx, ny, nz, depths, nxx, nyy, nzz,
                         ncav, nthreads);

  if (verbose)
    fprintf(stdout, "> Defining opening points\n");
  filter_openings(openings, depths, nx, ny, nz, nthreads);

  if (verbose)
    fprintf(stdout, "> Clustering opening points\n");
  nopenings = _cluster(openings, nx, ny, nz, step,
                       (double)openings_cutoff * pow(step, 3), verbose);

  return nopenings;
}

/* Retrieve interface residues */

// /*
//  * Struct: node
//  * ------------
//  *
//  * A linked list node for atom index in xyzr array
//  *
//  * pos: atom index in xyzr array (coordinates and radii of pdb)
//  * struct node* next: pointer to next linked list node
//  *
//  */
// typedef struct node {
//   int pos;
//   struct node *next;
// } res;

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
res *create(int pos)
{
  res *new = (res *)malloc(sizeof(res));

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
void insert(res **head, res *new)
{
  res *current;

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
 * returns: array of strings with interface residues with cavities separated by
 * '-1'
 */
char **interface(int *cavities, int nx, int ny, int nz, char **pdb,
                 double *atoms, int natoms, int xyzr, double *reference,
                 int ndims, double *sincos, int nvalues, double step,
                 double probe_in, int ncav, int nthreads)
{
  int i, j, k, atom, tag, count = 0, old_atom = -1, old_tag = -1;
  double x, y, z, xaux, yaux, zaux, distance, H;
  char **residues;

  // Allocate memory for reslist structure
  res **reslist = malloc(ncav * sizeof(res *));
  res *new, *old;

  // Initialize linked list
  for (i = 0; i < ncav; i++)
    reslist[i] = NULL;

  for (atom = 0; atom < natoms; atom++)
  {
    // Convert atom coordinates in 3D grid coordinates
    x = (atoms[atom * 4] - reference[0]) / step;
    y = (atoms[1 + (atom * 4)] - reference[1]) / step;
    z = (atoms[2 + (atom * 4)] - reference[2]) / step;

    xaux = x * sincos[3] + z * sincos[2];
    yaux = y;
    zaux = (-x) * sincos[2] + z * sincos[3];

    x = xaux;
    y = yaux * sincos[1] - zaux * sincos[0];
    z = yaux * sincos[0] + zaux * sincos[1];

    // Create a radius (H) for space occupied by probe and atom
    H = (probe_in + atoms[3 + (atom * 4)]) / step;

    // Loop around radius from atom center
    for (i = floor(x - H); i <= ceil(x + H); i++)
      for (j = floor(y - H); j <= ceil(y + H); j++)
        for (k = floor(z - H); k <= ceil(z + H); k++)
        {
          if (i < nx && i > 0 && j < ny && j > 0 && k < nz && k > 0)
            if (abs(cavities[k + nz * (j + (ny * i))]) > 1)
            {
              tag = cavities[k + nz * (j + (ny * i))] - 2;
              distance = sqrt(pow(i - x, 2) + pow(j - y, 2) + pow(k - z, 2));
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
  residues = calloc(count + ncav + 1, sizeof(char *));
  for (i = 0, j = 0; i < ncav; i++)
  {
    new = reslist[i];
    while (new != NULL)
    {
      residues[j++] = pdb[new->pos];
      old = new;
      new = old->next;
      free(old);
    }
    residues[j++] = "-1";
  }
  residues[j] = NULL;
  free(reslist);
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
 * returns: array of strings with interface residues with cavities separated by
 * '-1'
 *
 */
char **_constitutional(int *cavities, int nx, int ny, int nz, char **pdb,
                       double *atoms, int natoms, int xyzr, double *reference,
                       int ndims, double *sincos, int nvalues, double step,
                       double probe_in, int ncav, int nthreads, int verbose)
{
  char **residues;

  if (verbose)
    fprintf(stdout, "> Retrieving interface residues\n");
  residues =
      interface(cavities, nx, ny, nz, pdb, atoms, natoms, xyzr, reference,
                ndims, sincos, nvalues, step, probe_in, ncav, nthreads);

  return residues;
}

/* Estimate hydropathy */

/*
 * Function: get_hydrophobicity_value
 * ----------------------------------
 *
 * Get hydrophobicity scale value for a target residue name
 *
 * resname: target residue name
 * resn: 1D-array hydrophobicity scale residues names
 * scales: 1D-array of hydrophocity scale values
 * nscales: size of 1D-array of scales and resn
 *
 */
double get_hydrophobicity_value(char *resname, char **resn, double *scales,
                                int nscales)
{
  int i;

  // Get hydrophobicity value
  for (i = 0; i < nscales; i++)
    if (strcmp(resname, resn[i]) == 0)
      return scales[i];

  return 0.0;
}

/*
 * Function: project_hydropathy
 * ---------------------------
 *
 * Map a hydrophobicity scale per surface point of detected cavities.
 *
 * hydropathy: hydrophobicity scale 3D grid
 * surface: surface points 3D grid
 * nxx: x grid units
 * nyy: y grid units
 * nzz: z grid units
 * atoms: xyz coordinates and radii of input pdb
 * natoms: number of atoms
 * xyzr: number of data per atom (4: xyzr)
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * resname: 1D-array of residues names
 * resn: 1D-array hydrophobicity scale residues names
 * scales: 1D-array of hydrophocity scale values
 * nscales: size of 1D-array of scales and resn
 * step: 3D grid spacing (A)
 * probe_in: Probe In size (A)
 * nthreads: number of threads for OpenMP
 *
 */
void project_hydropathy(double *hydropathy, int *surface, int nxx, int nyy,
                        int nzz, double *atoms, int natoms, int xyzr,
                        double *reference, int ndims, double *sincos,
                        int nvalues, char **resname, char **resn,
                        double *scales, int nscales, double step,
                        double probe_in, int nthreads)
{
  int i, j, k, atom;
  double x, y, z, xaux, yaux, zaux, distance, H, *ref;

  // Initiliaze 3D grid for residues distances
  ref = (double *)calloc(nxx * nyy * nzz, sizeof(double));
  dgrid(ref, nxx * nyy * nzz);
  dgrid(hydropathy, nxx * nyy * nzz);

  // Get hydrophobicity value for each surface point
  for (atom = 0; atom < natoms; atom++)
  {
    // Convert atom coordinates in 3D grid coordinates
    x = (atoms[atom * 4] - reference[0]) / step;
    y = (atoms[1 + (atom * 4)] - reference[1]) / step;
    z = (atoms[2 + (atom * 4)] - reference[2]) / step;

    xaux = x * sincos[3] + z * sincos[2];
    yaux = y;
    zaux = (-x) * sincos[2] + z * sincos[3];

    x = xaux;
    y = yaux * sincos[1] - zaux * sincos[0];
    z = yaux * sincos[0] + zaux * sincos[1];

    // Create a radius (H) for space occupied by probe and atom
    H = (probe_in + atoms[3 + (atom * 4)]) / step;

    // Loop around radius from atom center
    for (i = floor(x - H); i <= ceil(x + H); i++)
      for (j = floor(y - H); j <= ceil(y + H); j++)
        for (k = floor(z - H); k <= ceil(z + H); k++)
        {
          if (i < nxx && i > 0 && j < nyy && j > 0 && k < nzz && k > 0)
            // Found a surface point
            if (surface[k + nzz * (j + (nyy * i))] > 1)
            {
              // Calculate distance bewteen atom and surface point
              distance = sqrt(pow(i - x, 2) + pow(j - y, 2) + pow(k - z, 2));
              // Check if surface point was not checked before
              if (ref[k + nzz * (j + (nyy * i))] == 0.0)
              {
                ref[k + nzz * (j + (nyy * i))] = distance;
                hydropathy[k + nzz * (j + (nyy * i))] =
                    get_hydrophobicity_value(resname[atom], resn, scales,
                                             nscales);
              }
              // Check if this atom is closer to the previous one assigned
              else if (ref[k + nzz * (j + (nyy * i))] > distance)
              {
                ref[k + nzz * (j + (nyy * i))] = distance;
                hydropathy[k + nzz * (j + (nyy * i))] =
                    get_hydrophobicity_value(resname[atom], resn, scales,
                                             nscales);
              }
            }
        }
  }

  // Free 3D grid for residues distances
  free(ref);
}

/* Estimate average hydropathy */

/*
 * Function: estimate_average_hydropathy
 * -------------------------------------
 *
 * Calculate average hydropathy of detected cavities.
 *
 * avgh: empty array of average hydropathy
 * ncav: number of cavities
 * hydropathy: hydrophobicity scale 3D grid
 * surface: surface points 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * nthreads: number of threads for OpenMP
 *
 */
void estimate_average_hydropathy(double *avgh, int ncav, double *hydropathy,
                                 int *surface, int nx, int ny, int nz,
                                 int nthreads)
{
  int i, j, k, *pts;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

  // Initialize array to get number of points in each cavity
  pts = (int *)calloc(ncav, sizeof(int));
  for (i = 0; i < ncav; i++)
    pts[i] = 0;

  // Initialize average hydropathy array
  dgrid(avgh, ncav);

#pragma omp parallel default(none), \
    shared(avgh, hydropathy, surface, pts, nx, ny, nz), private(i, j, k)
  {
#pragma omp for collapse(3) ordered
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
#pragma omp critical
          if (surface[k + nz * (j + (ny * i))] > 1)
          {
            pts[surface[k + nz * (j + (ny * i))] - 2]++;
            avgh[surface[k + nz * (j + (ny * i))] - 2] +=
                hydropathy[k + nz * (j + (ny * i))];
          }
  }

  for (i = 0; i < ncav; i++)
    avgh[i] /= pts[i];

  // Free array with number of points per cavity
  free(pts);
}

/* Hydropathy characterization */

/*
 * Function: _hydropathy
 * ---------------------
 *
 * Hydropathy characterization of the detected cavities. Map a hydrophobicity
 * scale per surface point and calculate average hydropathy of detected
 * cavities.
 *
 * hydropathy: hydrophobicity scale 3D grid
 * size: number of voxels in 3D grid
 * avgh: empty array of average hydropathy
 * ncav: number of cavities
 * surface: surface points 3D grid
 * nxx: x grid units
 * nyy: y grid units
 * nzz: z grid units
 * atoms: xyz coordinates and radii of input pdb
 * natoms: number of atoms
 * xyzr: number of data per atom (4: xyzr)
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * resname: 1D-array of residues names
 * resn: hydrophobicity scale residues names
 * scales: 1D-array of hydrophocity scale values
 * step: 3D grid spacing (A)
 * probe_in: Probe In size (A)
 * nthreads: number of threads for OpenMP
 * verbose: print extra information to standard output
 *
 * returns: hydropathy (3D grid with hydrophobicity scale of surface points),
 * and avg_h (array of average hydropathy)
 *
 */
void _hydropathy(double *hydropathy, int size, double *avgh, int ncav,
                 int *surface, int nxx, int nyy, int nzz, double *atoms,
                 int natoms, int xyzr, double *reference, int ndims,
                 double *sincos, int nvalues, char **resname, char **resn,
                 double *scales, int nscales, double step, double probe_in,
                 int nthreads, int verbose)
{
  if (verbose)
    fprintf(stdout, "> Mapping hydrophobicity scale at surface points\n");
  project_hydropathy(hydropathy, surface, nxx, nyy, nzz, atoms, natoms, xyzr,
                     reference, ndims, sincos, nvalues, resname, resn, scales,
                     nscales, step, probe_in, nthreads);

  if (verbose)
    fprintf(stdout, "> Estimating average hydropathy\n");
  estimate_average_hydropathy(avgh, ncav, hydropathy, surface, nxx, nyy, nzz,
                              nthreads);
}

/* Export cavity PDB */

/*
 * Function: _export
 * -----------------
 *
 * Export cavities to PDB file with depth as B-factor and hydropathy as Q-factor
 *
 * fn: cavity pdb filename
 * cavities: cavities 3D grid
 * nx: x grid units (cavities)
 * ny: y grid units (cavities)
 * nz: z grid units (cavities)
 * surface: surface points 3D grid
 * nxx: x grid units (surface)
 * nyy: y grid units (surface)
 * nzz: z grid units (surface)
 * B: B-factor (depth)
 * bx: x grid units (B)
 * by: y grid units (B)
 * bz: z grid units (B)
 * Q: Q-factor (hydropathy)
 * qx: x grid units (Q)
 * qy: y grid units (Q)
 * qz: z grid units (Q)
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * ncav: number of cavities
 * nthreads: number of threads for OpenMP
 * append: append cavities to PDB file
 * model: model number
 *
 */
void _export(
    char *fn,
    int *cavities, int nx, int ny, int nz,
    int *surface, int nxx, int nyy, int nzz,
    double *B, int bx, int by, int bz,
    double *Q, int qx, int qy, int qz,
    double *reference, int ndims,
    double *sincos, int nvalues,
    double step,
    int ncav,
    int nthreads,
    int append, int model)
{
  int i, j, k, count, tag;
  double x, y, z, xaux, yaux, zaux;
  FILE *output;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

  // Open cavity PDB file for writing and appending cavity coordinates
  output = append ? fopen(fn, "a+") : fopen(fn, "w");

  // Write model number
  if (abs(model) > 0)
    fprintf(output, "MODEL     %4.d\n", model);

  // Write cavity points for each cavity
  for (count = 1, tag = 2; tag <= ncav + 2; tag++)
#pragma omp parallel default(none) shared(cavities, surface, B, Q, reference, sincos, step, ncav, tag, count, nx, ny, nz, output), private(i, j, k, x, y, z, xaux, yaux, zaux)
  {
    // Iterate in x, y, z coordinates from 0 to nx, ny, nz
#pragma omp for schedule(static) collapse(3) ordered nowait
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
          // Check if cavity point with value tag
          if (cavities[k + nz * (j + (ny * i))] == tag)
          {
            // Convert 3D grid coordinates to Cartesian coordinates
            x = i * step;
            y = j * step;
            z = k * step;

            xaux = (x * sincos[3]) + (y * sincos[0] * sincos[2]) -
                   (z * sincos[1] * sincos[2]) + reference[0];
            yaux = (y * sincos[1]) + (z * sincos[0]) + reference[1];
            zaux = (x * sincos[2]) - (y * sincos[0] * sincos[3]) +
                   (z * sincos[1] * sincos[3]) + reference[2];

// Write cavity point as a surface point or a cavity point in PDB format
// depending on the value of the surface grid
#pragma omp critical
            if (surface[k + nz * (j + (ny * i))] == tag)
              fprintf(output,
                      "ATOM  %5.d  HA  K%c%c   259    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n",
                      count % 100000,
                      65 + (((surface[k + nz * (j + (ny * i))] - 2) / 26) % 26),
                      65 + ((surface[k + nz * (j + (ny * i))] - 2) % 26),
                      xaux,
                      yaux,
                      zaux,
                      Q[k + nz * (j + (ny * i))],
                      B[k + nz * (j + (ny * i))]);
            else
              fprintf(output,
                      "ATOM  %5.d  H   K%c%c   259    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n",
                      count % 100000,
                      65 + (((cavities[k + nz * (j + (ny * i))] - 2) / 26) % 26),
                      65 + ((cavities[k + nz * (j + (ny * i))] - 2) % 26),
                      xaux,
                      yaux,
                      zaux,
                      Q[k + nz * (j + (ny * i))],
                      B[k + nz * (j + (ny * i))]);
            count++;
          }
  }
  // Write ENDMDL when abs(model) > 0 and END when model < 0
  if (abs(model) > 0)
    fprintf(output, "ENDMDL\n");
  if (model < 0)
    fprintf(output, "END\n");
  // Close file
  fclose(output);
}

/*
 * Function: _export_openings
 * --------------------------
 *
 * Export openings to PDB file
 *
 * fn: openings PDB filename
 * openings: openings 3D grid
 * nxx: x grid units
 * nyy: y grid units
 * nzz: z grid units
 * reference: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * nopenings: number of openings
 * nthreads: number of threads for OpenMP
 * append: append cavities to PDB file
 * model: model number
 *
 */
void _export_openings(char *fn, int *openings, int nxx, int nyy, int nzz,
                      double *reference, int ndims, double *sincos, int nvalues,
                      double step, int nopenings, int nthreads, int append,
                      int model)
{
  int i, j, k, count, tag;
  double x, y, z, xaux, yaux, zaux;
  FILE *output;

  // Set number of threads in OpenMP
  omp_set_num_threads(nthreads);
  omp_set_nested(1);

  // Open cavity PDB file
  if (append)
    output = fopen(fn, "a+");
  else
    output = fopen(fn, "w");

  // Write model number
  if (abs(model) > 0)
    fprintf(output, "MODEL     %4.d\n", model);

  for (count = 1, tag = 2; tag <= nopenings + 2; tag++)
#pragma omp parallel default(none)                                             \
    shared(openings, reference, sincos, step, nopenings, tag, count, nxx, nyy, \
           nzz, output),                                                       \
    private(i, j, k, x, y, z, xaux, yaux, zaux)
  {
#pragma omp for schedule(static) collapse(3) ordered nowait
    for (i = 0; i < nxx; i++)
      for (j = 0; j < nyy; j++)
        for (k = 0; k < nzz; k++)
        {
          // Check if cavity point with value tag
          if (openings[k + nzz * (j + (nyy * i))] == tag)
          {
            // Convert 3D grid coordinates to real coordinates
            x = i * step;
            y = j * step;
            z = k * step;

            xaux = (x * sincos[3]) + (y * sincos[0] * sincos[2]) -
                   (z * sincos[1] * sincos[2]) + reference[0];
            yaux = (y * sincos[1]) + (z * sincos[0]) + reference[1];
            zaux = (x * sincos[2]) - (y * sincos[0] * sincos[3]) +
                   (z * sincos[1] * sincos[3]) + reference[2];

// Write cavity point
#pragma omp critical
            fprintf(output,
                    "ATOM  %5.d  H   O%c%c   259    %8.3lf%8.3lf%8.3lf  "
                    "1.00%6.2lf\n",
                    count % 100000,
                    65 +
                        (((openings[k + nzz * (j + (nyy * i))] - 2) / 26) % 26),
                    65 + ((openings[k + nzz * (j + (nyy * i))] - 2) % 26), xaux,
                    yaux, zaux, 0.0);
            count++;
          }
        }
  }
  // Write ENDMDL
  if (abs(model) > 0)
    fprintf(output, "ENDMDL\n");
  // Write END
  if (model < 0)
    fprintf(output, "END\n");
  // Close file
  fclose(output);
}
