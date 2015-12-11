%module mask

%{
#define SWIG_FILE_WITH_INIT

extern void get_mask_array(long long *farray, int nflexible, int natoms, char **nname,long long *resid,long long *flexible_residues,int nresidues, int mtype);

 
%}

%typemap(in) char ** {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (char **) malloc((size+1)*sizeof(char *));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o))
        $1[i] = PyString_AsString(PyList_GetItem($input,i));
      else {
        PyErr_SetString(PyExc_TypeError,"list must contain strings");
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

%include "numpy.i"

%init %{
import_array();
%}

%apply (long long* INPLACE_ARRAY2, int DIM1, int DIM2) {(long long* farray, int len1, int len2)}
%apply (long long* IN_ARRAY1, int DIM1) {(long long* resid, int len3)}
%apply (long long* IN_ARRAY1, int DIM1) {(long long* flexible_residues, int len4)}

/*extern void get_mask_array(int **farray, int nflexible, int natoms, char **name,int *resid,int *flexible_residues,int nresidues);*/

%rename (get_mask_array) my_get_mask_array;

%inline %{
        void my_get_mask_array(long long *farray, int len1, int len2, char **nname,long long *resid,int len3,long long *flexible_residues,int len4, int nresidues, int mtype){
        return get_mask_array(farray,len1,len2,nname,resid,flexible_residues,nresidues,mtype);
}
%}

/*%clear (int** IN_PLACEARRAY2, int DIM1, int DIM2),(int** farray, int len1, int len2) ; */


// This cleans up the char ** array we malloc'd before the function call
%typemap(freearg) char ** {
  free((char *) $1);
}




