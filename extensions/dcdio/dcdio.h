#include <stdio.h>
#include <string.h>

/* These defines are necessary to support CHARMm DCD files. See
   comments in ReadDCD.C, in read_dcdheader(). */
#define DCD_IS_CHARMM 		0x01
#define DCD_HAS_4DIMS		0x02
#define DCD_HAS_EXTRA_BLOCK	0x04

/*  DEFINE ERROR CODES THAT MAY BE RETURNED BY DCD ROUTINES		*/
#define DCD_DNE		-2	/*  DCD file does not exist		*/
#define DCD_OPENFAILED	-3	/*  Open of DCD file failed		*/
#define DCD_BADREAD 	-4	/*  read call on DCD file failed	*/
#define DCD_BADEOF	-5	/*  premature EOF found in DCD file	*/
#define DCD_BADFORMAT	-6	/*  format of DCD file is wrong		*/
#define DCD_FILEEXISTS  -7	/*  output file already exists		*/
#define DCD_BADMALLOC   -8	/*  malloc failed			*/

extern FILE * open_dcd_read(const char *);  

extern int read_dcdheader(FILE * fd, int *N, int  *NSET, int  *ISTART,
               int  *NSAVC, double  *DELTA, int  *NAMNF,
               int *reverseEndian, int *charmm);

extern int close_dcd_read(FILE * fd);

extern FILE * open_dcd_write(char *dcdname);
extern int write_dcdheader(FILE * fd, char *filename, int N, int NSET, \
               int ISTART, int NSAVC, double DELTA);
extern int write_dcdstep(FILE *,int N,float *x,float *y,float *z, int curframe);
extern int close_dcd_write(FILE * fd);

extern int* reverseFourByteWord(int* ); /* revert 4 byte word byteorder */

extern int* reverseEightByteDouble(int* ); /* revert 8 byte word byteorder */


