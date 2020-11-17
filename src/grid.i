%module gridprocessing

%{
    #define SWIG_FILE_WITH_INIT
    #include "grid.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%include "typemaps.i"

%typemap(in) char ** 
{
  /* Check if is a list */
    if (PyList_Check($input)) {
        int size = PyList_Size($input);
        Py_ssize_t i = 0;
        $1 = (char **) malloc((size+1)*sizeof(char *));
        for (i = 0; i < size; i++) {
            PyObject *o = PyList_GetItem($input,i);
            if (PyUnicode_Check(o))
                $1[i] = PyUnicode_AsUTF8(PyList_GetItem($input,i));
            else {
                //PyErr_SetString(PyExc_TypeError,"list must contain strings");
                PyErr_Format(PyExc_TypeError, "list must contain strings. %d/%d element was not string.", i, size);
                free($1);
                return NULL;
            }
        }
        $1[i] = 0;
    } else {
        PyErr_SetString(PyExc_TypeError,"not a list");
        return NULL;
    }
}

// This cleans up the char ** array we malloc'd before the function call
%typemap(freearg) char ** 
{
  free((char *) $1);
}

%typemap(out) (char **test)
{
   int ntokens;
   int itoken;
   PyObject *py_string_tmp;
   int py_err;

   /* compute len of the list */
   for (ntokens = 0; $1[ntokens] != NULL; ntokens++) {
      /* not sure is loop is safe... */
   }

   /* create Python empty list */
   $result = PyList_New(ntokens);
   if (! $result) return NULL;

   /* fill Python list */
   for (itoken = 0; itoken < ntokens; itoken++) {
       if ($1[itoken] == NULL) break;

       /* convert C string to Python string */
       py_string_tmp = PyString_FromString( $1[itoken] );
       if (! py_string_tmp) return NULL;

       /* put Python string into the list */
       py_err = PyList_SetItem($result, itoken, py_string_tmp);
       if (py_err == -1) return NULL;
   }


   /* Pyhon list is a copy a C list, delete C list */
//    for (itoken = 0; itoken < ntokens; itoken++) {
//        free($1[itoken]);
//    }
   free($1);

   /* return Python result */
   return $result; 
}

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

%include "grid.h"
