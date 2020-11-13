%module gridprocessing

%{
    #define SWIG_FILE_WITH_INIT
    #include "grid.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

// %include "typemaps.i"

/* **** GRID **** */
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* PI, int size)}
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* surface, int size)}
%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(int *cavities, int nx, int ny, int nz)}

/* Box coordinates */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *reference, int ndims
)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *sincos, int nvalues)}

/* Atom coordinates */
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *atoms, int natoms, int xyzr)}

%include "grid.h"