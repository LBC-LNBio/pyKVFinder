/* Grid initialization */
void igrid (int *grid, int size);
void fgrid (float *grid, int size);
void dgrid (double *grid, int size);
void cgrid (int *grid, int size);

/* Grid filling */
void fill (int *grid, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int nthreads);

/* Biomolecular surface representation */
int check_protein_neighbours (int *grid, int nx, int ny, int nz, int i, int j, int k);
void ses (int *grid, int nx, int ny, int nz, double step, double probe, int nthreads);

/* Grid subtract (Probe In - Probe Out) */;
void subtract (int *PI, int *PO, int nx, int ny, int nz, double step, double removal_distance, int nthreads);

/* Filter noise from Grid */
void filter_noise (int *grid, int nx, int ny, int nz, int nthreads);

/* Ligand adjustment */
void adjust (int *grid, int nx, int ny, int nz, double *ligand, int lnatoms, int lxyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double ligand_cutoff, int nthreads);

/* Box adjustment */
void _filter_pdb (int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int nthreads);
void filter (int *grid, int nx, int ny, int nz, double *P1, int ndims, double *P2, int nndims, double *sincos, int nvalues, double step, double probe_out, int nthreads);

/* Cavity clustering */
int cluster (int *grid, int nx, int ny, int nz, double step, double volume_cutoff, int nthreads);
void DFS (int *grid, int nx, int ny, int nz, int i, int j, int k, int tag);
void remove_cavity (int *grid, int nx, int ny, int nz, int tag, int nthreads);

/* Cavity detection */
int _detect (int *PI, int size, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, double probe_out, double removal_distance, double volume_cutoff, int box_adjustment, double *P2, int nndims, int is_ses, int nthreads, int verbose);
int _detect_ladj (int *PI, int size, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *ligand, int lnatoms, int lxyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, double probe_out, double removal_distance, double volume_cutoff, int ligand_adjustment, double ligand_cutoff, int box_adjustment, double *P2, int nndims, int is_ses, int nthreads, int verbose);

/* Cavity surface points */
int define_surface_points (int *cavities, int nx, int ny, int nz, int i, int j, int k);
void filter_surface (int *cavities, int *surface, int nx, int ny, int nz, int nthreads);

/* Estimate area */
double check_voxel_class (int *surface, int nx, int ny, int nz, int i, int j, int k);
void area (int *surface, int nx, int ny, int nz, int ncav, double step, double *areas, int nthreads);

/* Estimate volume */
void volume (int *cavities, int nx, int ny, int nz, int ncav, double step, double *volumes, int nthreads);

/* Spatial characterization */
void _spatial (int *cavities, int nx, int ny, int nz, int *surface, int size, double *volumes, int nvol, double *areas, int narea, double step, int nthreads, int verbose);

/* Retrieve interface residues */
typedef struct node { int pos; struct node* next; } res;
void insert (res** head, res* new);
res* create (int pos);
char **interface (int *cavities, int nx, int ny, int nz, char **pdb, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, int ncav, int nthreads);

/* Constitutional characterization */
char ** _constitutional (int *cavities, int nx, int ny, int nz, char **pdb, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, int ncav, int nthreads, int verbose);

/* Export cavity PDB */
void _export (char *fn, int *cavities, int nx, int ny, int nz, int *surf, int nxx, int nyy, int nzz, double *reference, int ndims, double *sincos, int nvalues, double step, int ncav, int nthreads);
