#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#ifndef FORCE_NOMPI
#include <mpi.h>
#else
#define MPI_COMM_WORLD 0
#endif
#include <complex.h>
/* to access functions from the libpastix, respect this order */
#include "pastix.h"
#include "cscd_utils.h"
#include "read_matrix.h"
#include "get_options.h"
#include "utils.h"

#define CALLOC(type,num) (type *) calloc((unsigned)(num),(sizeof (type)))
#define REALLOC(ptr,type,num) \
        (type *) realloc((type *)(ptr), (unsigned)((num)*sizeof(type)))
#define FREE(x) if((x)!=NULL) {free((char *)(x));(x)=NULL;}

#include "mat.h"
#include "matrix.h"
#include "mex.h"

#define CALLOC(type,num) (type *) calloc((unsigned)(num),(sizeof (type)))
#define REALLOC(ptr,type,num) \
        (type *) realloc((type *)(ptr), (unsigned)((num)*sizeof(type)))
#define FREE(x) if((x)!=NULL) {free((char *)(x));(x)=NULL;}

int read_nx(char *, mwSize *);
int write_matlab_solution(char *fileName, char *varName, int nx, double *xvals);
void usage(int argc, char *argv[]);
void read_params_file(char *params_file, int *nbthread, int *incomplete,
                             int *level_of_fill, int *amalgamation,
                             int *verbosemode, int *ooc);
int get_options(int              argc,
                char           **argv,
                driver_type_t  **driver_type,
                char          ***filename,
                int             *nbmatrices,
                int             *nbthread,
                int             *verbosemode,
                int             *ordering,
                int             *incomplete,
                int             *level_of_fill,
                int             *amalgamation,
                int             *ooc,
                pastix_int_t    *size);
int main (int argc, char **argv);

