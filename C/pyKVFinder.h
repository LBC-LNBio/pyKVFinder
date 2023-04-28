/* Grid initialization */
void igrid(int *grid, int size);
void fgrid(float *grid, int size);
void dgrid(double *grid, int size);
void cgrid(int *grid, int size);

/* Grid filling */
void fill(int *grid, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int nthreads);
void _fill_receptor(int *receptor, int size, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, int is_ses, int nthreads, int verbose);
void _fill_cavity(int *cavities, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, int nthreads);

/* Biomolecular surface representation */
int check_protein_neighbours(int *grid, int nx, int ny, int nz, int i, int j, int k);
void ses(int *grid, int nx, int ny, int nz, double step, double probe, int nthreads);

/* Grid subtract (Probe In - Probe Out) */
void subtract(int *PI, int *PO, int nx, int ny, int nz, double step, double removal_distance, int nthreads);

/* Filter noise from Grid */
void filter_noise(int *grid, int nx, int ny, int nz, int nthreads);

/* Ligand adjustment */
void adjust(int *grid, int nx, int ny, int nz, double *ligand, int lnatoms, int lxyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double ligand_cutoff, int nthreads);

/* Box adjustment */
void _filter_pdb(int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int nthreads);
void filter(int *grid, int nx, int ny, int nz, double *P1, int ndims, double *P2, int nndims, double *sincos, int nvalues, double step, double probe_out, int nthreads);

/* Cavity clustering */
int check_unclustered_neighbours(int *grid, int nx, int ny, int nz, int i, int j, int k);
void DFS(int *grid, int nx, int ny, int nz, int i, int j, int k, int tag);
void remove_cavity(int *grid, int nx, int ny, int nz, int tag, int nthreads);
int _cluster(int *grid, int nx, int ny, int nz, double step, double volume_cutoff, int nthreads);

/* Cavity detection */
int _detect(int *PI, int size, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, double probe_out, double removal_distance, double volume_cutoff, int box_adjustment, double *P2, int nndims, int is_ses, int nthreads, int verbose);
int _detect_ladj(int *PI, int size, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *ligand, int lnatoms, int lxyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, double probe_out, double removal_distance, double volume_cutoff, int ligand_adjustment, double ligand_cutoff, int box_adjustment, double *P2, int nndims, int is_ses, int nthreads, int verbose);

/* Cavity surface points */
int define_surface_points(int *cavities, int nx, int ny, int nz, int i, int j, int k);
void filter_surface(int *cavities, int *surface, int nx, int ny, int nz, int nthreads);

/* Estimate area */
double check_voxel_class(int *surface, int nx, int ny, int nz, int i, int j, int k);
void _area(int *surface, int nxx, int nyy, int nzz, double step, double *areas, int narea, int nthreads);

/* Estimate volume */
void _volume(int *cavities, int nx, int ny, int nz, double step, double *volumes, int nvol, int nthreads);

/* Spatial characterization */
void _spatial(int *cavities, int nx, int ny, int nz, int *surface, int size, double *volumes, int nvol, double *areas, int narea, double step, int nthreads, int verbose);

/* Bulk-cavity boundary points */
typedef struct points
{
    double X1;
    double X2;
    double Y1;
    double Y2;
    double Z1;
    double Z2;
} pts;
int define_boundary_points(int *cavities, int nx, int ny, int nz, int i, int j, int k);
void filter_boundary(int *cavities, int nx, int ny, int nz, pts *cavs, pts *boundarys, int nthreads);
void remove_boundary(int *cavities, int nx, int ny, int nz, int ncav, pts *boundaries, int nthreads);

/* Estimate depth */
void estimate_depth(int *cavities, double *depths, int nx, int ny, int nz, double *max_depth, double *avg_depth, int n, pts *cavs, pts *boundaries, double step, int nthreads);

/* Depth characterization */
void _depth(int *cavities, int nx, int ny, int nz, double *depths, int size, double *max_depth, int nmax, double *avg_depth, int navg, double step, int nthreads, int verbose);

/* Openings characterization */
void _openings2cavities(int *o2c, int nopenings, int *cavities, int nx, int ny, int nz, int *openings, int nxx, int nyy, int nzz, int nthreads);
void remove_enclosed_cavity(int *openings, int *cavities, int nx, int ny, int nz, double *depths, int nxx, int nyy, int nzz, int ncav, int nthreads);
void filter_openings(int *openings, double *depths, int nx, int ny, int nz, int nthreads);
int _openings(int *openings, int size, int *cavities, int nx, int ny, int nz, double *depths, int nxx, int nyy, int nzz, int ncav, int openings_cutoff, double step, int nthreads, int verbose);

/* Retrieve interface residues */
typedef struct node
{
    int pos;
    struct node *next;
} res;
res *create(int pos);
void insert(res **head, res *new);
char **interface(int *cavities, int nx, int ny, int nz, char **pdb, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, int ncav, int nthreads);

/* Constitutional characterization */
char **_constitutional(int *cavities, int nx, int ny, int nz, char **pdb, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe_in, int ncav, int nthreads, int verbose);

/* Estimate hydropathy */
double get_hydrophobicity_value(char *resname, char **resn, double *scales, int nscales);
void project_hydropathy(double *hydropathy, int *surface, int nxx, int nyy, int nzz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, char **resname, char **resn, double *scales, int nscales, double step, double probe_in, int nthreads);

/* Estimate average hydropathy */
void estimate_average_hydropathy(double *avgh, int ncav, double *hydropathy,
                                 int *surface, int nx, int ny, int nz,
                                 int nthreads);

/* Hydropathy characterization */
void _hydropathy(double *hydropathy, int size, double *avgh, int ncav, int *surface, int nxx, int nyy, int nzz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, char **resname, char **resn, double *scales, int nscales, double step, double probe_in, int nthreads, int verbose);

/* Export cavity PDB */
void _export(char *fn, int *cavities, int nx, int ny, int nz, int *surface, int nxx, int nyy, int nzz, double *B, int bx, int by, int bz, double *Q, int qx, int qy, int qz, double *reference, int ndims, double *sincos, int nvalues, double step, int ncav, int nthreads, int append, int model);
void _export_openings(char *fn, int *openings, int nxx, int nyy, int nzz, double *reference, int ndims, double *sincos, int nvalues, double step, int nopenings, int nthreads, int append, int model);
