typedef struct POINTS 
{
    double X1;
    double X2;
    double Y1;
    double Y2;
    double Z1;
    double Z2;
} pts;

void dgrid (double *grid, int size);
int define_boundary_points (int *cavities, int nx, int ny, int nz, int i, int j, int k);
void filter_boundary (int *cavities, int nx, int ny, int nz, pts *cavs, pts *boundarys, int nthreads);
void define_depth(int *cavities, double *depth, int nx, int ny, int nz, double *max_depth, double *avg_depth, int n, pts *cavs, pts *boundaries, int nthreads);
void remove_boundary (int *cavities, int nx, int ny, int nz, int ncav, pts *boundaries, int nthreads);
void _depth (int *cavities, int nx, int ny, int nz, double *depth, int size, double *max_depth, int n, double *avg_depth, int nn, double step, int nthreads, int verbose);