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
int define_frontier_points (int *cavities, int nx, int ny, int nz, int i, int j, int k);
void filter_frontier (int *cavities, int nx, int ny, int nz, int nthreads);
void _depth (int *cavities, int nx, int ny, int nz, double *depth, int size, double *max_depth, int n, double *avg_depth, int nn, double step, int nthreads, int verbose);