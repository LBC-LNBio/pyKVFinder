%module pyKVFinder

%{
    #define SWIG_FILE_WITH_INIT
    #include "pyKVFinder.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

/*** Numpy definitions ***/

/* Receptor grid */
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* receptor, int size)}

/* Cavities grid */
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* PI, int size)}
%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(int *cavities, int nx, int ny, int nz)}

/* Surface points grid */
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* surface, int size)}
%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(int* surface, int nxx, int nyy, int nzz)}

/* Depth grid */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* depths, int size)}
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(double* depths, int nxx, int nyy, int nzz)}

/* Clustering grid */
%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(int *grid, int nx, int ny, int nz)}

/* Hydropathy grid */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* hydropathy, int size)}
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(double* Q, int qx, int qy, int qz)}

/* Openings grid */
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* openings, int size)}
%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(int* openings, int nxx, int nyy, int nzz)}
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* o2c, int nopenings)}

/* Bfactor grid */
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(double* B, int bx, int by, int bz)}

/* Volume array */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* volumes, int nvol)}

/* Area array */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* areas, int narea)}

/* Depth arrays */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* max_depth, int nmax), (double* avg_depth, int navg)}

/* Average hydropathy array */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* avgh, int ncav)}

/* Origin coordinates */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *reference, int ndims)}

/* X-axis vertice coordinates */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *P2, int nndims)}

/* PDB coordinates */
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *atoms, int natoms, int xyzr)}

/* Ligand coordinates */
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *ligand, int lnatoms, int lxyzr)}

/* Sine and Cossine */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *sincos, int nvalues)}

/* Hydrophobicity scale array */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *scales, int nscales)}

%include "typemaps.i"

/* Map array of strings (Python) to char** (C) */
%typemap(in) 
char ** 
{
    /* Check if is a list */
    if (PyList_Check($input)) 
    {
        int size = PyList_Size($input);
        Py_ssize_t i = 0;
        $1 = (char **) malloc((size+1)*sizeof(char *));
        for (i = 0; i < size; i++)
        {
            PyObject *o = PyList_GetItem($input,i);
            if (PyUnicode_Check(o))
                $1[i] = PyUnicode_AsUTF8(PyList_GetItem($input,i));
            else 
            {
                //PyErr_SetString(PyExc_TypeError,"list must contain strings");
                PyErr_Format(PyExc_TypeError, "list must contain strings. %d/%d element was not string.", i, size);
                free($1);
                return NULL;
            }
        }
        $1[i] = 0;
    } 
    else 
    {
        PyErr_SetString(PyExc_TypeError,"not a list");
        return NULL;
    }
}

/* Free char ** array (C) */
%typemap(freearg) 
char ** 
{
  free((char *) $1);
}

/* Map char ** (C) to array of strings (Python) */
%typemap(out) 
(char **_constitutional)
{
    int nC, nPy, py_err;
    PyObject *tmp;

    /* Define list length */
    for (nC = 0; $1[nC] != NULL; nC++);

    /* Create Python list */
    $result = PyList_New(nC);
    if (!$result) 
        return NULL;

    /* Pass C list to Python list */
    for (nPy = 0; nPy < nC; nPy++) 
    {
        if ($1[nPy] == NULL)
            break;

        /* Convert C string to Python string */
        tmp = PyString_FromString( $1[nPy] );
        if (!tmp) 
            return NULL;

        /* Pass Python string to Python list */
        py_err = PyList_SetItem($result, nPy, tmp);
        if (py_err == -1) 
            return NULL;
    }

   /* Delete C list */
   free($1);

   return $result; 
}

%include "pyKVFinder.h"
