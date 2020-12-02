%module depth

%{
    #define SWIG_FILE_WITH_INIT
    #include "depth.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

/*** Numpy definitions ***/

/* Cavities grid */
%apply (int* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(int *cavities, int nx, int ny, int nz)}

/* Depth grid */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* depth, int size)}

/* Maximum Depth array */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* max_depth, int n)}

/* Average Depth array */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* avg_depth, int nn)}

%include "depth.h"