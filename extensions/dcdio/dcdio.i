%module dcdio

%{
    #define SWIG_FILE_WITH_INIT

extern FILE * open_dcd_write(char *dcdname);
extern int write_dcdstep(FILE * fd,int N,float *x,float *y,float *z, int curframe);
extern int write_dcdheader(FILE * fd, char *filename, int N, int NSET, \
               int ISTART, int NSAVC, double DELTA);
extern int close_dcd_write(FILE * fd);

extern FILE * open_dcd_read(const char *);

extern int read_dcdheader(FILE * fd, int *N, int  *NSET, int  *ISTART,
               int  *NSAVC, double  *DELTA, int  *NAMNF,
               int *reverseEndian, int *charmm);

extern int read_dcdstep(FILE * fd, int N, float *X, float *Y, float *Z, int num_fixed,
        int first, int reverseEndian, int charmm);

extern int close_dcd_read(FILE * fd);

%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (int DIM1, float* IN_ARRAY1) {(int len1, float* x), (int len2, float* y), (int len3, float* z)}

extern int write_dcdstep(FILE * fd,int N,float *x,float *y,float *z,int curframe);

%rename (write_dcdstep) my_write_dcdstep;

%inline %{
    int my_write_dcdstep(FILE * fd, int len1, float* x, int len2, float* y, int len3, float* z, int curframe) {
    if (len1 != len2 || len1 != len3) {
        PyErr_Format(PyExc_ValueError, "Arrays of lengths (%d,%d,%d) given", len1, len2,len3);
        return 0 ;
    }
    return write_dcdstep(fd,len1,x,y,z,curframe);
}
%}

%clear (int len1, float* x), (int len2, float* y), (int len3, float* z) ;

%apply (int DIM1, float* INPLACE_ARRAY1) {(int len1, float* x), (int len2, float* y), (int len3, float* z)}

extern int read_dcdstep(FILE * fd, int DIM1, float *X, float *Y, float *Z, int num_fixed,
        int first, int reverseEndian, int charmm);

%rename (read_dcdstep) my_read_dcdstep;

%inline %{
    int my_read_dcdstep(FILE * fd, int len1, float* x, int len2, float* y, int len3, float* z, int num_fixed, int first, int reverseEndian, int charmm) {
    if (len1 != len2 || len1 != len3) {
        PyErr_Format(PyExc_ValueError, "Arrays of lengths (%d,%d,%d) given", len1, len2,len3);
        return 0 ;
    }
    return read_dcdstep(fd,len1,x,y,z,num_fixed,first,reverseEndian,charmm);
}
%}

%clear (int len1, float* x), (int len2, float* y), (int len3, float* z);

extern FILE * open_dcd_write(char *dcdname);

extern int write_dcdheader(FILE * fd, char *filename, int N, int NSET, \
               int ISTART, int NSAVC, double DELTA);

extern int close_dcd_write(FILE * fd);

extern FILE * open_dcd_read(const char *);

extern int close_dcd_read(FILE * fd);

%include "typemaps.i"

extern int read_dcdheader(FILE *, int *OUTPUT, int *OUTPUT, int *OUTPUT, int *OUTPUT, double *OUTPUT, int *OUTPUT, int *OUTPUT, int *OUTPUT);


