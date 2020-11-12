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
    int ncores,
    int verbose);

/* Grid initialization */
void igrid (int *grid, int size);
void fgrid (float *grid, int size);
void dgrid (double *grid, int size);
void cgrid (int *grid, int size);

/* Grid filling */
void fill (int *grid, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int ncores);

/* Surface representation */
void check_protein_neighbours (int *grid, int nx, int ny, int nz, int i, int j, int k);
void ses (int *grid, int nx, int ny, int nz, double step, double probe, int ncores);

/* Grid subtract (Probe In - Probe Out) */;
void subtract (int *PI, int *PO, int nx, int ny, int nz, double step, double removal_threshold, int ncores);

/* Export cavity PDB */
void export (char *fn, int *cavities, int nx, int ny, int nz, double *reference, int ndims, double *sincos, int nvalues, double step, int is_filled);

/* Debug */
void count (int *grid, int nx, int ny, int nz);
void filter (int *grid, int dx, int dy, int dz);