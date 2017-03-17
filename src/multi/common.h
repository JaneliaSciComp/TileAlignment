/* Copyright 2008 BORDEAUX I UNIVERSITY & INRIA 
 * **
 * ** This file is part of the PaStiX parallel sparse matrix package.
 * **
 * ** This software is governed by the CeCILL-C license under French law
 * ** and abiding by the rules of distribution of free software. You can
 * ** use, modify and/or redistribute the software under the terms of the
 * ** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
 * ** URL: "http://www.cecill.info".
 * ** 
 * ** As a counterpart to the access to the source code and rights to copy,
 * ** modify and redistribute granted by the license, users are provided
 * ** only with a limited warranty and the software's author, the holder of
 * ** the economic rights, and the successive licensors have only limited
 * ** liability.
 * ** 
 * ** In this respect, the user's attention is drawn to the risks associated
 * ** with loading, using, modifying and/or developing or reproducing the
 * ** software by the user in light of its specific status of free software,
 * ** that may mean that it is complicated to manipulate, and that also
 * ** therefore means that it is reserved for developers and experienced
 * ** professionals having in-depth computer knowledge. Users are therefore
 * ** encouraged to load and test the software's suitability as regards
 * ** their requirements in conditions enabling the security of their
 * ** systems and/or data to be ensured and, more generally, to use and
 * ** operate it in the same conditions as regards security.
 * ** 
 * ** The fact that you are presently reading this means that you have had
 * ** knowledge of the CeCILL-C license and that you accept its terms.
 * */

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

int main (int argc, char **argv);
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

int read_nx(char *, mwSize *);
int write_matlab_solution(char *fileName, char *varName, int nx, double *xvals);
void usage(int argc, char *argv[]);
void read_params_file(char *params_file, int *nbthread, int *incomplete,
                             int *level_of_fill, int *amalgamation,
                             int *verbosemode, int *ooc);

