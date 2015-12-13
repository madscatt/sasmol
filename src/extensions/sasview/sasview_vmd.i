%module sasview_vmd

%{
    #define SWIG_FILE_WITH_INIT

extern int send_coordinates_to_vmd(int N, float *x,float *y,float *z,int port, int flag);
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (int DIM1, float* IN_ARRAY1) {(int len1, float* x), (int len2, float* y), (int len3, float* z)}

extern int send_coordinates_to_vmd(int N, float *x,float *y,float *z,int port, int flag);

%rename (send_coordinates_to_vmd) my_send_coordinates_to_vmd;

%inline %{
    int my_send_coordinates_to_vmd(int len1, float* x, int len2, float* y, int len3, float* z, int port, int flag) {
    if (len1 != len2 || len1 != len3) {
        PyErr_Format(PyExc_ValueError, "Arrays of lengths (%d,%d,%d) given", len1, len2,len3);
        return 0 ;
    }
    return send_coordinates_to_vmd(len1,x,y,z,port,flag) ;
}
%}

%clear (int len1, float* x), (int len2, float* y), (int len3, float* z) ;

%include "typemaps.i"


