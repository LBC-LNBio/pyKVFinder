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

/* Creating grids */
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* grid, int size)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* grid, int size)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* grid, int size)}

/* 3D integer grid coordinates */
%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(int* grid, int dx, int dy, int dz)}
%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(int* grid2, int dx2, int dy2, int dz2)}

/* Box coordinates */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *reference, int ndims
)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *sincos, int nvalues)}

/* Atom coordinates */
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *atoms, int natoms, int xyzr)}

%include "grid.h"