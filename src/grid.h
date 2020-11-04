/* Grid creation */
void igrid (int *grid, int size);
void fgrid (float *grid, int size);
void dgrid (double *grid, int size);
void cgrid (int *grid, int size);

/* Grid filling */
void insert_atom (int *grid, int dx, int dy, int dz, double *atom, int xyz, double *reference, int ndims, int npoints, double step, double probe_in);

/* Debug */
void filter (int *grid, int dx, int dy, int dz);