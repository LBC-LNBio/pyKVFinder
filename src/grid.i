%module grid

%{
    #define SWIG_FILE_WITH_INIT
    #include "grid.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

/* **** GRID **** */
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* PI, int size)}
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* surface, int size)}
%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(int *cavities, int nx, int ny, int nz)}
%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(int* surf, int nxx, int nyy, int nzz)}

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* volumes, int nvol)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* areas, int narea)}

/* Box coordinates */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *reference, int ndims
)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *sincos, int nvalues)}

/* Atom coordinates */
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *atoms, int natoms, int xyzr)}

%include "typemaps.i"

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

// This cleans up the char ** array we malloc'd before the function call
%typemap(freearg) 
char ** 
{
  free((char *) $1);
}

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

%include "grid.h"
