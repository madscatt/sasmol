#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "dcdio.h"

FILE * open_dcd_write(char *dcdname) {
  return fopen(dcdname, "w+b");
}

int close_dcd_write(FILE * fd) {
 int result; 
 result = fclose(fd);
 return result; 
}

/* WRITE Macro to make porting easier */
#define WRITE(fd, buf, size) \
        fwrite((buf), (size), 1, (fd))

void pad(char *s, int len)
{
      int curlen;
      int i;

      curlen=strlen(s);

      if (curlen>len)
      {
            s[len]=0;
            return;
      }

      for (i=curlen; i<len; i++)
      {
            s[i]=' ';
      }

      s[i]=0;
}

int write_dcdheader(FILE * fd, char *filename, int N, int NSET,
               int ISTART, int NSAVC, double DELTA)
{
  int out_integer;
  char      title_string[200];
  time_t  cur_time;
  struct  tm *tmbuf;
  char    time_str[11];

  out_integer = 84;
  WRITE(fd, (char *) & out_integer, sizeof(int));
  strcpy(title_string, "CORD");
  WRITE(fd, title_string, 4);
  WRITE(fd, (char *) &NSET, sizeof(int));
  WRITE(fd, (char *) &ISTART, sizeof(int));
  WRITE(fd, (char *) &NSAVC, sizeof(int));
  out_integer=0;
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &DELTA, sizeof(double));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  WRITE(fd, (char *) &out_integer, sizeof(int));
  out_integer = 84;
  WRITE(fd, (char *) & out_integer, sizeof(int));

  out_integer = 164;
  WRITE(fd, (char *) & out_integer, sizeof(int));
  out_integer = 2;
  WRITE(fd, (char *) & out_integer, sizeof(int));

  sprintf(title_string, "REMARKS FILENAME=%s :: SASSIE", "A.DCD");
  pad(title_string, 80);
  WRITE(fd, title_string, 80);

  cur_time=time(NULL);
  tmbuf=localtime(&cur_time);
  strftime(time_str, 10, "%m/%d/%y", tmbuf);

  sprintf(title_string, "REMARKS DATE: %s CREATED BY USER: %s",
     time_str, "ikuo");
  pad(title_string, 80);
  WRITE(fd, title_string, 80);
  out_integer = 164;
  WRITE(fd, (char *) & out_integer, sizeof(int));
  out_integer = 4;
  WRITE(fd, (char *) & out_integer, sizeof(int));
  out_integer = N;
  WRITE(fd, (char *) & out_integer, sizeof(int));
  out_integer = 4;
  WRITE(fd, (char *) & out_integer, sizeof(int));
  return 1;
}

/* defines used by write_dcdstep */
#define NFILE_POS 8L
#define NSTEP_POS 20L
#define FIO_SEEK_SET SEEK_SET
#define FIO_SEEK_END SEEK_END

/*
 *       curframe: Count of frames written to this file, starting with 1.
 *       curstep: Count of timesteps elapsed = istart + curframe * nsavc.
*/

int write_dcdstep(FILE * fd, int N, float *X, float *Y, float *Z, int curframe)
{
      int out_integer;
      int curstep ;

      curstep=curframe ;

      out_integer = N*4;
      WRITE(fd, (char *) &out_integer, sizeof(int));
      WRITE(fd, (char *) X, out_integer);
      WRITE(fd, (char *) &out_integer, sizeof(int));
      WRITE(fd, (char *) &out_integer, sizeof(int));
      WRITE(fd, (char *) Y, out_integer);
      WRITE(fd, (char *) &out_integer, sizeof(int));
      WRITE(fd, (char *) &out_integer, sizeof(int));
      WRITE(fd, (char *) Z, out_integer);
      WRITE(fd, (char *) &out_integer, sizeof(int));
	
/*
static int fio_write_int32(fio_fd fd, int i) {
  return (fio_fwrite(&i, 4, 1, fd) != 1);
}

static fio_size_t fio_fwrite(void *ptr, fio_size_t size,
                             fio_size_t nitems, fio_fd fd) {
  return fwrite(ptr, size, nitems, fd);
}

      fio_write_int32(fd, curframe);
      fwrite(ptr, size, nitems, fd);
*/
 
  /* update the DCD header information */
      fseek(fd,NFILE_POS,FIO_SEEK_SET);
      fwrite(&curframe, 4, 1, fd);
      fseek(fd,NSTEP_POS,FIO_SEEK_SET);
      fwrite(&curframe, 4, 1, fd);
      fseek(fd,0,FIO_SEEK_END);

      return 1;
}


FILE * open_dcd_read(const char *filename) {
  return fopen(filename, "rb");
}

/* READ Macro to make porting easier */
#define READ(fd, buf, size) \
        fread((buf), (size), 1, (fd))

/* WRITE Macro to make porting easier */
#define WRITE(fd, buf, size) \
        fwrite((buf), (size), 1, (fd))

int* reverseFourByteWord(int* N)
{
  char byteArray[4] = {'0','0','0','0'};

  char *bytePointer;

  bytePointer = (char*)N;

  byteArray[0]  =  *bytePointer;
  byteArray[1]  =  *(bytePointer+1);
  byteArray[2]  =  *(bytePointer+2);
  byteArray[3]  =  *(bytePointer+3);

  *bytePointer     = byteArray[3];
  *(bytePointer+1) = byteArray[2];
  *(bytePointer+2) = byteArray[1];
  *(bytePointer+3) = byteArray[0];

  N=(int*)bytePointer;

  return N;
}


double* reverseEightByteWord(double* N)
{
  char byteArray[8] = {'0','0','0','0','0','0','0','0'};
  char *bytePointer;

  bytePointer = (char*)N;

  byteArray[0]  =  *bytePointer;
  byteArray[1]  =  *(bytePointer+1);
  byteArray[2]  =  *(bytePointer+2);
  byteArray[3]  =  *(bytePointer+3);
  byteArray[4]  =  *(bytePointer+4);
  byteArray[5]  =  *(bytePointer+5);
  byteArray[6]  =  *(bytePointer+6);
  byteArray[7]  =  *(bytePointer+7);

  *bytePointer     = byteArray[7];
  *(bytePointer+1) = byteArray[6];
  *(bytePointer+2) = byteArray[5];
  *(bytePointer+3) = byteArray[4];
  *(bytePointer+4) = byteArray[3];
  *(bytePointer+5) = byteArray[2];
  *(bytePointer+6) = byteArray[1];
  *(bytePointer+7) = byteArray[0];

  N=(double*)bytePointer;

  return N;
}

#define CHECK_FREAD(X, msg)  if (X==-1) \
                       { \
                        return(DCD_BADREAD); \
                       }

#define CHECK_FEOF(X, msg)  if (X==0) \
                       { \
                        return(DCD_BADEOF); \
                       }

int read_dcdheader(FILE * fd, int *N, int *NSET, int *ISTART,
               int *NSAVC, double *DELTA, int *NAMNF,
               int *reverseEndian, int *charmm)
{
  int input_integer;    /*  Integer buffer space      */
  int ret_val;          /*  Return value from read    */
  int i;          /*  Loop counter        */
  char hdrbuf[84];      /*  Char buffer used to store header      */
  int NTITLE;
  int **FREEINDEXES;

/*  (*FREEINDEXES) = *((int *) (hdrbuf + 12)); */


  /*  First thing in the file should be an 84         */
  ret_val = READ(fd, &input_integer, sizeof(int));
  CHECK_FREAD(ret_val, "reading first int from dcd file");
  CHECK_FEOF(ret_val, "reading first int from dcd file");

  /* Check magic number in file header and determine byte order*/
  if (input_integer != 84) {
    /* check to see if its merely reversed endianism     */
    /* rather than a totally incorrect file magic number */
    input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer == 84) {
      *reverseEndian=1;
    }
    else {
      return(DCD_BADFORMAT);
    }
  }
  else {
    *reverseEndian=0;
  }

  /*  Buffer the entire header for random access      */
  ret_val = READ(fd, hdrbuf, 84);
  CHECK_FREAD(ret_val, "buffering header");
  CHECK_FEOF(ret_val, "buffering header");

  /*  Check for the ID string "COORD"           */
  if (hdrbuf[0] != 'C' || hdrbuf[1] != 'O' ||
      hdrbuf[2] != 'R' || hdrbuf[3] != 'D') {
    return DCD_BADFORMAT;
   }

  /* CHARMm-generate DCD files set the last integer in the     */
  /* header, which is unused by X-PLOR, to its version number.  */
  /* Checking if this is nonzero tells us this is a CHARMm file */
  /* and to look for other CHARMm flags.                        */
  if (*((int *) (hdrbuf + 80)) != 0) {
      (*charmm) = DCD_IS_CHARMM;
      if (*((int *) (hdrbuf + 44)) == 1)
            (*charmm) |= DCD_HAS_EXTRA_BLOCK;
      if (*((int *) (hdrbuf + 48)) == 1)
            (*charmm) |= DCD_HAS_4DIMS;
  }
  else (*charmm) = 0;

  if ((*charmm) & DCD_IS_CHARMM) {
    if ((*charmm) & DCD_HAS_EXTRA_BLOCK) {
    }
    if ((*charmm) & DCD_HAS_4DIMS) {
    }
  }
  else {
  }

  /*  Store the number of sets of coordinates (NSET) */
  (*NSET) = *((int *) (hdrbuf + 4));
  if (*reverseEndian) NSET=reverseFourByteWord(NSET);

  /*  Store ISTART, the starting timestep      */
  (*ISTART) = *((int *) (hdrbuf + 8));
  if (*reverseEndian) ISTART=reverseFourByteWord(ISTART);

  /*  Store NSAVC, the number of timesteps between   */
  /*  dcd saves                                    */
  (*NSAVC) = *((int *) (hdrbuf + 12));
  if (*reverseEndian) NSAVC=reverseFourByteWord(NSAVC);

  /*  Store NAMNF, the number of free atoms          */
  (*NAMNF) = *((int *) (hdrbuf + 36));
  if (*reverseEndian) NAMNF=reverseFourByteWord(NAMNF);

  /*  Read in the timestep, DELTA               */
  /*  Note: DELTA is stored as a double with X-PLOR but */
  /*  as a float with CHARMm                    */
  if ((*charmm) & DCD_IS_CHARMM) {
    double * tmp;
    float * tmp2;
    tmp = (double *) (hdrbuf + 40);
    if (*reverseEndian) tmp=reverseEightByteWord(tmp);
    tmp2 = (float *) tmp;
    (*DELTA) = (double) *(tmp2);
  }
  else {
    (*DELTA) = *((double *) (hdrbuf + 40));
    if (*reverseEndian) DELTA=reverseEightByteWord(DELTA);
  }

  /*  Get the end size of the first block             */
  ret_val = READ(fd, &input_integer, sizeof(int));
  CHECK_FREAD(ret_val, "reading second 84 from dcd file");
  CHECK_FEOF(ret_val, "reading second 84 from dcd file");
  if (*reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

  if (input_integer != 84) {
    return(DCD_BADFORMAT);
  }
  /*  Read in the size of the next block              */
  ret_val = READ(fd, &input_integer, sizeof(int));
  CHECK_FREAD(ret_val, "reading size of title block");
  CHECK_FEOF(ret_val, "reading size of title block");
  if (*reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

  if ( ((input_integer-4)%80) == 0) {
    /*  Read NTITLE, the number of 80 characeter    */
    /*  title strings there are                 */
    ret_val = READ(fd, &NTITLE, sizeof(int));
    CHECK_FREAD(ret_val, "reading NTITLE");
    CHECK_FEOF(ret_val, "reading NTITLE");
    if (*reverseEndian) NTITLE= *reverseFourByteWord(&NTITLE);

    for (i=0; i<NTITLE; i++) {
      fseek(fd, 80, SEEK_CUR);
      CHECK_FEOF(ret_val, "reading TITLE");
    }

    /*  Get the ending size for this block            */
    ret_val = READ(fd, &input_integer, sizeof(int));

    CHECK_FREAD(ret_val, "reading size of title block");
    CHECK_FEOF(ret_val, "reading size of title block");
  }
  else {
    return(DCD_BADFORMAT);
  }

  /*  Read in an 4                        */
  ret_val = READ(fd, &input_integer, sizeof(int));
  CHECK_FREAD(ret_val, "reading an 4");
  CHECK_FEOF(ret_val, "reading an 4");
  if (*reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

  if (input_integer != 4) {
    return(DCD_BADFORMAT);
  }

  /*  Read in the number of atoms               */
  ret_val = READ(fd, N, sizeof(int));
  CHECK_FREAD(ret_val, "reading number of atoms");
  CHECK_FEOF(ret_val, "reading number of atoms");
  if (*reverseEndian) N=reverseFourByteWord(N);

  /*  Read in an 4                        */
  ret_val = READ(fd, &input_integer, sizeof(int));
  CHECK_FREAD(ret_val, "reading an 4");
  CHECK_FEOF(ret_val, "reading an 4");
  if (*reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

  if (input_integer != 4) {
    return(DCD_BADFORMAT);
  }

  if (*NAMNF != 0) {
    (*FREEINDEXES) = (int *) calloc(((*N)-(*NAMNF)), sizeof(int));

    if (*FREEINDEXES == NULL)
    return(DCD_BADFORMAT);
    /*  Read in an size                   */
    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading size of index array");
    CHECK_FEOF(ret_val, "reading size of index array");
    if (*reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer != ((*N)-(*NAMNF))*4) {
    return(DCD_BADFORMAT);
    }

    ret_val = READ(fd, (*FREEINDEXES), ((*N)-(*NAMNF))*sizeof(int));
    CHECK_FREAD(ret_val, "reading size of index array");
    CHECK_FEOF(ret_val, "reading size of index array");
    if (*reverseEndian)
    {
      printf(">>>*** in reverseEndian loop\n");
      for (input_integer =0; input_integer < ((*N)-(*NAMNF)); input_integer++)
      {
        FREEINDEXES[input_integer]= reverseFourByteWord(FREEINDEXES[input_integer]);
      }
    }

    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading size of index array");
    CHECK_FEOF(ret_val, "reading size of index array");
    if (*reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer != ((*N)-(*NAMNF))*4) {
    return(DCD_BADFORMAT);
    }
  }

  return(0);
}

int read_dcdstep(FILE * fd, int N, float *X, float *Y, float *Z, int num_fixed,
        int first, int reverseEndian, int charmm) 
        /*     int first, int *indexes, int reverseEndian, int charmm) */
{
  int ret_val;          /*  Return value from read          */
  int input_integer;    /*  Integer buffer space            */
  int i,j;                /*  Loop counter              */
/*  int indexes;

  indexes = (int *) malloc(N*sizeof(int)) ;
   for (j=0;j<N;j++) {
        indexes[j]=1 ;
   } */
/*  printf("number of atoms = %d\n",N);
  printf("num_fixed = %d\n",num_fixed) ;
  printf("first = %d\n",first) ;
  if(num_fixed==0) printf("nf = true\n") ;
  if(first==0) printf("f0 = true\n") ;
  if(first==1) printf("f1 = true\n") ;
  if(first) printf("first = true\n") ;
*/
  if ( (num_fixed==0) || first) {
//    printf("\n\n>>> in first loop\n\n") ;
    /* If this is a CHARMm file and contains an extra data block,
       we must skip it to avoid problems */
    if ((charmm & DCD_IS_CHARMM) &&
        (charmm & DCD_HAS_EXTRA_BLOCK)) {
        ret_val = READ(fd, &input_integer, sizeof(int));
        //printf("zhl start %d\n",input_integer);
        CHECK_FREAD(ret_val, "reading extra charmm block");
        if (reverseEndian) input_integer = *reverseFourByteWord(&input_integer);
        fseeko(fd, input_integer, SEEK_CUR);
        ret_val = READ(fd, &input_integer, sizeof(int));
        CHECK_FREAD(ret_val, "reading extra charmm block");
        //printf("zhl now %d\n",input_integer);
   }

    /*  Get the first size from the file                          */
    /*printf(">>> ret_val = %d\n",ret_val) ;*/
    ret_val = READ(fd, &input_integer, sizeof(int));
 //     printf("> input_integer = %d\n",input_integer) ;
//    printf(">>> ret_val = %d\n",ret_val) ;
    CHECK_FREAD(ret_val, "reading number of atoms at begining of step");
    if (reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    /*  See if we've reached the end of the file                  */
    if (ret_val == 0) {
      printf(">>> reached the end of the file (1)\n\n") ;
      return(-1);
    }


    /*printf("> input_integer = %d\n",input_integer) ;*/

    if (input_integer != 4*N){
      printf(">>> input_integer != 4*N (1)\n\n") ;
      return(DCD_BADFORMAT);
    }
    /*printf(">>> ret_val = %d\n",ret_val) ;*/
    ret_val = READ(fd, X, 4*N);
      //printf("ZHL %8.3f %8.3f %8.3f %8.3f\n",X[0],X[1],X[N-2],X[N-1]);
    /*printf(">>> just before reading X array: ret_val = %d\n",ret_val) ;*/
    CHECK_FREAD(ret_val, "reading X array");
    CHECK_FEOF(ret_val, "reading X array");
    if (reverseEndian) {
      for (i=0; i<N; i++) {
        X[i]=*((float*)(reverseFourByteWord((int*)&X[i])));
        printf("X[i] = %f\n",X[i]) ;
      }
    }
    /*printf(">>> just after reading X array: ret_val = %d\n",ret_val) ;*/
    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading number of atoms after X array");
    CHECK_FEOF(ret_val, "reading number of atoms after X array");
    if (reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer != 4*N){
      printf(">>> input_integer != 4*N (2)\n\n") ;
      return(DCD_BADFORMAT);
    }
    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading number of atoms after X array");
    CHECK_FEOF(ret_val, "reading number of atoms after X array");
    if (reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer != 4*N){
      printf(">>> input_integer != 4*N (3)\n\n") ;
      return(DCD_BADFORMAT);
    }

    ret_val = READ(fd, Y, 4*N);
    CHECK_FREAD(ret_val, "reading Y array");
    CHECK_FEOF(ret_val, "reading Y array");
    if (reverseEndian) {
      for (i=0; i<N; i++) {
        Y[i]=*((float*)(reverseFourByteWord((int*)&Y[i])));
      }
    }

    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading number of atoms after Y array");
    CHECK_FEOF(ret_val, "reading number of atoms after Y array");
    if (reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer != 4*N){
      printf(">>> input_integer != 4*N (4)\n\n") ;
      return(DCD_BADFORMAT);
    }

    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading number of atoms after Y array");
    CHECK_FEOF(ret_val, "reading number of atoms after Y array");
    if (reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer != 4*N){
      printf(">>> input_integer != 4*N (5)\n\n") ;
      return(DCD_BADFORMAT);
    }
    ret_val = READ(fd, Z, 4*N);
    CHECK_FREAD(ret_val, "reading Z array");
    CHECK_FEOF(ret_val, "reading Z array");
    if (reverseEndian) {
      for (i=0; i<N; i++) {
        Z[i]=*((float*)(reverseFourByteWord((int*)&Z[i])));
      }
    }

    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading number of atoms after Z array");
    CHECK_FEOF(ret_val, "reading number of atoms after Z array");
    if (reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer != 4*N){
      printf(">>> input_integer != 4*N (6)\n\n") ;
      return(DCD_BADFORMAT);
    }


    /* If this is a CHARMm file and contains a 4th dimension block,
       we must skip it to avoid problems */
    if ((charmm & DCD_IS_CHARMM) &&
        (charmm & DCD_HAS_4DIMS)) {
        ret_val = READ(fd, &input_integer, sizeof(int));
        CHECK_FREAD(ret_val, "reading extra charmm block");
        if (reverseEndian) input_integer = *reverseFourByteWord(&input_integer);
        fseeko(fd, input_integer, SEEK_CUR);
        ret_val = READ(fd, &input_integer, sizeof(int));
        CHECK_FREAD(ret_val, "reading extra charmm block");
    }

  }
  else {

    /* If this is a CHARMm file and contains an extra data block,
       we must skip it to avoid problems */
    if ((charmm & DCD_IS_CHARMM) &&
        (charmm & DCD_HAS_EXTRA_BLOCK)) {
        ret_val = READ(fd, &input_integer, sizeof(int));
        CHECK_FREAD(ret_val, "reading extra charmm block");
        if (reverseEndian) input_integer = *reverseFourByteWord(&input_integer);
        fseeko(fd, input_integer, SEEK_CUR);
        ret_val = READ(fd, &input_integer, sizeof(int));
        CHECK_FREAD(ret_val, "reading extra charmm block");
    }

    /*  Get the first size from the file                          */
    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading number of atoms at begining of step");
    if (reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    /*  See if we've reached the end of the file                  */
    if (ret_val == 0) {
      printf(">>> reached the end of the file (2)\n\n") ;
      return(-1);
    }

    if (input_integer != 4*(N-num_fixed)){
      printf(">>> input_integer != 4*(N-num_fixed) (1)\n\n") ;
      return(DCD_BADFORMAT);
    }
    fseeko(fd, 4*(N-num_fixed), SEEK_CUR);
    CHECK_FEOF(ret_val, "reading tmpX array");

    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading number of atoms after X array");
    CHECK_FEOF(ret_val, "reading number of atoms after X array");
    if (reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer != 4*(N-num_fixed)){
      printf(">>> input_integer != 4*(N-num_fixed) (2)\n\n") ;
      return(DCD_BADFORMAT);
    }

    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading number of atoms after X array");
    CHECK_FEOF(ret_val, "reading number of atoms after X array");
    if (reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer != 4*(N-num_fixed)){
      printf(">>> input_integer != 4*(N-num_fixed) (3)\n\n") ;
      return(DCD_BADFORMAT);
    }

    fseeko(fd, 4*(N-num_fixed), SEEK_CUR);
    CHECK_FEOF(ret_val, "reading tmpX array");

    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading number of atoms after Y array");
    CHECK_FEOF(ret_val, "reading number of atoms after Y array");
    if (reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer != 4*(N-num_fixed)){
      printf(">>> input_integer != 4*(N-num_fixed) (4)\n\n") ;
      return(DCD_BADFORMAT);
    }

    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading number of atoms after Y array");
    CHECK_FEOF(ret_val, "reading number of atoms after Y array");
    if (reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer != 4*(N-num_fixed)){
      printf(">>> input_integer != 4*(N-num_fixed) (5)\n\n") ;
      return(DCD_BADFORMAT);
    }

    fseeko(fd, 4*(N-num_fixed), SEEK_CUR);
    CHECK_FEOF(ret_val, "reading tmpX array");

    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading number of atoms after Z array");
    CHECK_FEOF(ret_val, "reading number of atoms after Z array");
    if (reverseEndian) input_integer= *reverseFourByteWord(&input_integer);

    if (input_integer != 4*(N-num_fixed)){
      printf(">>> input_integer != 4*(N-num_fixed) (6)\n\n") ;
      return(DCD_BADFORMAT);
    }

    /* If this is a CHARMm file and contains a 4th dimension block,
       we must skip it to avoid problems */
    if ((charmm & DCD_IS_CHARMM) &&
        (charmm & DCD_HAS_4DIMS)) {
        ret_val = READ(fd, &input_integer, sizeof(int));
        CHECK_FREAD(ret_val, "reading extra charmm block");
        if (reverseEndian) input_integer = *reverseFourByteWord(&input_integer);
        fseeko(fd, input_integer, SEEK_CUR);
        ret_val = READ(fd, &input_integer, sizeof(int));
        CHECK_FREAD(ret_val, "reading extra charmm block");
    }
  }
  /*printf(" I got to the end, jeepers\n");*/
  return(0);
}

int close_dcd_read(FILE * fd) {
 int result; 
 result = fclose(fd);
 return result; 
}

