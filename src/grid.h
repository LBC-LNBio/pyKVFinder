/* Write cavity */
void export (int *grid, int dx, int dy, int dz, int *grid2, int dx2, int dy2, int dz2, double step, char *fn, double *reference, int ndims, double *sincos, int nvalues, int cavity_representation);

/* Current implementation */

/* Cavity detection */
void detect (int *PI, int size, 
    int nx, int ny, int nz,
    double *atoms, int natoms, int xyzr, 
    double *reference, int ndims, 
    double *sincos, int nvalues, 
    double step, 
    double probe_in,
    double probe_out,
    double removal_threshold,
    int is_ses,
    int ncores);

/* Grid initialization */
void igrid (int *grid, int size);
void fgrid (float *grid, int size);
void dgrid (double *grid, int size);
void cgrid (int *grid, int size);

/* Grid filling */
void fill (int *grid, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int ncores);

/* Surface representation */
void check_protein_neighbours (int *grid, int nx, int ny, int nz, int i, int j, int k);
void ses (int *grid, int nx, int ny, int nz, double step, double probe_in, int ncores);

/* Grid subtract (Probe In - Probe Out) */;

/* Debug */
void count (int *grid, int nx, int ny, int nz);
void filter (int *grid, int dx, int dy, int dz);