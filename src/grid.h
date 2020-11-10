/* Grid creation */
void igrid (int *grid, int size);
void fgrid (float *grid, int size);
void dgrid (double *grid, int size);
void cgrid (int *grid, int size);

/* Grid filling */
void fill (int *grid, int dx, int dy, int dz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int ncores, int is_sas);

/* Grid subtract */
void subtract (int *grid, int dx, int dy, int dz, int *grid2, int dx2, int dy2, int dz2, double step, double removal_threshold, int ncores);

/* Write cavity */
void export (int *grid, int dx, int dy, int dz, int *grid2, int dx2, int dy2, int dz2, double step, char *fn, double *reference, int ndims, double *sincos, int nvalues, int cavity_representation);

/* SAS representation */
void check_protein_neighbours (int *A, int dx, int dy, int dz, int i, int j, int k, int *is_neighbour);
void ses (int *grid, int dx, int dy, int dz, double step, double sas, int ncores);

/* Debug */
void filter (int *grid, int dx, int dy, int dz);