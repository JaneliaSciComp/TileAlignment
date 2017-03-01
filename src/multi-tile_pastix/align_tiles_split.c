/* Copyright 2008 BORDEAUX I UNIVERSITY & INRIA 
**
** This file is part of the PaStiX parallel sparse matrix package.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/* File: simple_dist.c
 *
 *  A simple example with a distributed matrix :
 *  read the matrix, check it is correct and correct it if needed,
 *  distribute it then run pastix in one call.
 *
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
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

int read_nx(char *filename, mwSize *nx)
{
    MATFile *pmat;
    pmat = matOpen(filename, "r");

    mxArray *mxA;
    mxA  = matGetVariable(pmat, "A");

   *nx = mxGetM(mxA);
    mxDestroyArray(mxA);
    matClose(pmat);
}

int write_matlab_solution(char *fileName, char *varName, int nx, double *xvals)
{
    mxArray *xv;

    MATFile *pmat = matOpen(fileName, "w");

    xv = mxCreateDoubleMatrix(1, nx, mxREAL);

    memcpy((void *)(mxGetPr(xv)), (void *)xvals, nx*sizeof(xvals));

    matPutVariable(pmat, varName, xv);

    matClose(pmat);
}


int main (int argc, char **argv)
{

  pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix           */
  pastix_int_t    ncol;               /* Number of local columns                                   */
  pastix_int_t   *colptr      = NULL; /* Indexes of first element of each column in row and values */
  pastix_int_t   *rows        = NULL; /* Row of each element of the matrix                         */
  pastix_int_t   *loc2glob    = NULL; /* Local to local column correspondance                      */
  pastix_float_t *values      = NULL; /* Value of each element of the matrix                       */
  pastix_float_t *rhs         = NULL; /* right-hand-side                                           */
  pastix_float_t *rhssaved    = NULL; /* right hand side (save)                                    */
  pastix_float_t *rhssaved_g  = NULL; /* right hand side (save, global)                            */
  pastix_int_t    iparm[IPARM_SIZE];  /* integer parameters for pastix                             */
  double          dparm[DPARM_SIZE];  /* floating parameters for pastix                            */
  pastix_int_t   *perm        = NULL; /* Permutation tabular                                       */
  pastix_int_t   *invp        = NULL; /* Reverse permutation tabular                               */
  char           *type        = NULL; /* type of the matrix                                        */
  char           *rhstype     = NULL; /* type of the right hand side                               */
#ifndef FORCE_NOMPI
  int             required;           /* MPI thread level required                                 */
  int             provided;           /* MPI thread level provided                                 */
#endif
  int             mpid, num_nodes;
  driver_type_t  *driver_type;        /* Matrix driver(s) requested by user                        */
  char          **filename;           /* Filename(s) given by user                                 */
  int             nbmatrices;         /* Number of matrices given by user                          */
  int             nbthread;           /* Number of thread wanted by user                           */
  int             verbosemode;        /* Level of verbose mode (0, 1, 2)                           */
  int             ordering;           /* Ordering to use                                           */
  int             incomplete;         /* Indicate if we want to use incomplete factorisation       */
  int             level_of_fill;      /* Level of fill for incomplete factorisation                */
  int             amalgamation;       /* Level of amalgamation for Kass                            */
  int             ooc;                /* OOC limit (Mo/percent depending on compilation options)   */
  pastix_int_t    mat_type;
  long            i;
  pastix_int_t     globn;
  /*******************************************/
  /*          MPI initialisation             */
  /*******************************************/
#ifndef FORCE_NOMPI
  required = MPI_THREAD_MULTIPLE;
  provided = -1;
  MPI_Init_thread(&argc, &argv, required, &provided);
  MPI_Comm_size (MPI_COMM_WORLD, &num_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
//fprintf(stdout, "\nnum_nodes=%d mpid=%d \n\n", num_nodes, mpid);
  if (mpid == 0)
    {
      switch (provided)
        {
        case MPI_THREAD_SINGLE:
          printf("MPI_Init_thread level = MPI_THREAD_SINGLE\n");
          break;
        case MPI_THREAD_FUNNELED:
          printf("MPI_Init_thread level = MPI_THREAD_FUNNELED\n");
          break;
        case MPI_THREAD_SERIALIZED:
          printf("MPI_Init_thread level = MPI_THREAD_SERIALIZED\n");
          break;
        case MPI_THREAD_MULTIPLE:
          printf("MPI_Init_thread level = MPI_THREAD_MULTIPLE\n");
          break;
        default:
          printf("MPI_Init_thread level = ???\n");
        }
    }
#else
  mpid = 0;
#endif

  /*******************************************/
  /*    Get options from command line        */
  /*******************************************/
  if (EXIT_FAILURE == get_options(argc, argv,     &driver_type,
                                  &filename,      &nbmatrices,
                                  &nbthread,      &verbosemode,
                                  &ordering,      &incomplete,
                                  &level_of_fill, &amalgamation,
                                  &ooc,           &ncol))
    return EXIT_FAILURE;

  if (nbmatrices != 1)
    {
      /* Matrices for each iteration must have the same patern, this is why we only
         authorize one matrix in this exemple.
         But it could be used with several matrices with same patern and different values.
      */
      fprintf(stdout,"WARNING: should have only one matrix\n");
    }
  /*******************************************/
  /*      Read Matrice                       */
  /*******************************************/
  dread_matrix(filename[0], &ncol, &colptr, &rows, &loc2glob, &values,
               &rhs, &type, &rhstype, driver_type[0], MPI_COMM_WORLD);

/*
  fprintf(stdout, "\nmpid=%d Right-hand side:\n", mpid);
  for (i=0; i<ncol; i++)
      fprintf(stdout, "%f ", rhs[i]);
  fprintf(stdout, "\n");
*/

  mat_type = API_SYM_YES;
  //if (MTX_ISSYM(type)) mat_type = API_SYM_YES;
  //if (MTX_ISHER(type)) mat_type = API_SYM_HER;

  free(driver_type);

  /*******************************************/
  /*    Check Matrix format                  */
  /*******************************************/
  /*
   * Matrix needs :
   *    - to be in fortran numbering
   *    - to have only the lower triangular part in symmetric case
   *    - to have a graph with a symmetric structure in unsymmetric case
   */
  pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
                     mat_type,  API_YES,
                     ncol, &colptr, &rows, &values, &loc2glob, 1);

  /*******************************************/
  /* Initialize parameters to default values */
  /*******************************************/
  iparm[IPARM_MODIFY_PARAMETER] = API_NO;
  dpastix(&pastix_data, MPI_COMM_WORLD,
          ncol, colptr, rows, values, loc2glob,
          perm, invp, rhs, 1, iparm, dparm);

  /*******************************************/
  /*       Customize some parameters         */
  /*******************************************/
  iparm[IPARM_SYM] = API_SYM_YES;
  mat_type = API_SYM_YES;
  switch (mat_type)
  {
  case API_SYM_YES:
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    break;
  case API_SYM_HER:
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLH;
    break;
  default:
    iparm[IPARM_FACTORIZATION] = API_FACT_LU;
  }
  iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
  iparm[IPARM_ORDERING]            = API_ORDER_SCOTCH;

  iparm[IPARM_VERBOSE]             = verbosemode;                  
  iparm[IPARM_THREAD_NBR]          = nbthread;
  iparm[IPARM_INCOMPLETE]          = incomplete;   
  iparm[IPARM_OOC_LIMIT]           = ooc;               
  iparm[IPARM_LEVEL_OF_FILL]       = level_of_fill;
  iparm[IPARM_AMALGAMATION_LEVEL]  = amalgamation;

  if (incomplete == 1)
    {
      dparm[DPARM_EPSILON_REFINEMENT] = 1e-25;
    }
  iparm[IPARM_REFINEMENT]          = API_RAF_GMRES;
  iparm[IPARM_GMRES_IM]            = 25;
  iparm[IPARM_LEVEL_OF_FILL]       = -1;     // use Kass algorithm
  iparm[IPARM_RHS_MAKING]          = API_RHS_B;
  iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
  iparm[IPARM_END_TASK]            = API_TASK_CLEAN;
  /* reread parameters to set IPARM/DPARM */
  if (EXIT_FAILURE == get_idparm(argc, argv,
                                 iparm,          dparm))
    return EXIT_FAILURE;

  /*******************************************/
  /*           Save the rhs                  */
  /*    (it will be replaced by solution)    */
  /*******************************************/
  rhssaved = malloc(ncol*sizeof(pastix_float_t));
  memcpy(rhssaved, rhs, ncol*sizeof(pastix_float_t));
#ifndef FORCE_NOMPI
  MPI_Allreduce(&ncol, &globn, 1, MPI_PASTIX_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  globn = ncol;
#endif
  rhssaved_g = malloc(globn*sizeof(pastix_float_t));
  memset(rhssaved_g, 0, globn*sizeof(pastix_float_t));
  for (i = 0; i < ncol; i++)
    {
      rhssaved_g[loc2glob[i]-1] = rhssaved[i];
    }

  free(rhssaved);
#ifndef FORCE_NOMPI
  {
    pastix_float_t * rhssaved_g_rcv;
    rhssaved_g_rcv = malloc(globn*sizeof(pastix_float_t));
    MPI_Allreduce(rhssaved_g, rhssaved_g_rcv, globn, MPI_PASTIX_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    free(rhssaved_g);
    rhssaved_g = rhssaved_g_rcv;
  }
#endif
  /*******************************************/
  /*           Call pastix                   */
  /*******************************************/
  perm = malloc(ncol*sizeof(pastix_int_t));
  /* No need to allocate invp in dpastix */
  PRINT_RHS("RHS", rhs, ncol, mpid, iparm[IPARM_VERBOSE]);

  iparm[IPARM_ITERMAX] = 5;
  dpastix(&pastix_data, MPI_COMM_WORLD,
          ncol, colptr, rows, values, loc2glob,
          perm, NULL, rhs, 1, iparm, dparm);
//fprintf(stderr, "Printing rhs...\n");
  PRINT_RHS("SOL", rhs, ncol, mpid, iparm[IPARM_VERBOSE]);
//fprintf(stderr, "Checking solution...\n");
  CHECK_DIST_SOL(colptr, rows, values, rhs, ncol, loc2glob, globn, rhssaved_g);


  /* Compare thge result with Matlab's stored solution  */
  pastix_float_t *xvals       = NULL;
  mwSize nx;
  float  err;
  int col_min, col_max, tiles_per_node, num_coord;

  char filename1[256];
  int n = strlen(filename[0]);
  strncpy(filename1, filename[0], n-4);
  filename1[n-4] = '\0';
  sprintf(filename1, "%s_%d.mat", filename1, mpid + 1);
  read_nx(filename1, &nx);
 
  tiles_per_node = round((float)nx/6./(float)num_nodes);
  col_min = 1 + 6*(mpid * tiles_per_node);
  if (mpid < num_nodes-1) {
      col_max = col_min   + 6*tiles_per_node - 1;
  } else {
      col_max = nx;
  }
  num_coord = col_max - col_min +1;

  char mat_file_name[256];           
  double *data = CALLOC(double, num_coord);
  for (i=0; i< num_coord; i++)
      data[i] = (double)rhs[i];
  if (strlen(getenv("TA_DATA")) > 0) {
      sprintf(mat_file_name, "%s/x_%d.mat", getenv("TA_DATA"), mpid +1);
  } else {
      sprintf(mat_file_name,    "x_%d.mat"                   , mpid +1);
  }

  fprintf(stderr, "mpid=%d mat_file_name=%s\n", mpid, mat_file_name);                 

  write_matlab_solution(mat_file_name, "x", num_coord, data);
  FREE(data);

  for (i = 0; i < nbmatrices; i++)
    if (filename[i] != NULL)
      free(filename[i]);
  free(filename);

  free(colptr);
  free(rows);
  free(values);
  free(rhs);
//free(rhssaved_g);
//free(type);
//free(rhstype);
  free(loc2glob);
#ifndef FORCE_NOMPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}
/* File: align_tiles.c
 */

int dread_matrix(char            *filename,
                pastix_int_t    *ncol,
                pastix_int_t   **colptrs,
                pastix_int_t   **rows,
                pastix_int_t   **loc2glb,
                pastix_float_t **vals,
                pastix_float_t **rhs,
                char           **mat_type,
                char           **rhs_type,
                driver_type_t    driver_type,
                MPI_Comm         pastix_comm)
{
    int num_nodes, mpid;
    MPI_Comm_size (MPI_COMM_WORLD, &num_nodes);
    MPI_Comm_rank (MPI_COMM_WORLD, &mpid);
    int n = strlen(filename);
    char filename1[256];
    strncpy(filename1, filename, n-4);
    filename1[n-4] = '\0';
//  fprintf(stdout, "\nnum_nodes=%d mpid=%d prefix=%s\n\n", num_nodes, mpid, filename1);
    if (num_nodes > 1)
        sprintf(filename1, "%s_%d.mat", filename1, mpid + 1);
    else
        sprintf(filename1, "%s", filename);

//  fprintf(stdout, "\nnum_nodes=%d mpid=%d filename1=%s\n\n", num_nodes, mpid, filename1);

   *rhs_type = "API_RHS_B";        // user defined

    if (filename != NULL)
    {
        *mat_type = "API_SYM_YES";

        // Read input data from MAT file
        MATFile *pmat;
        pmat = matOpen(filename1, "r");
        if (pmat == NULL) {
            printf("Error reading MAT file");
            return(EXIT_FAILURE);
        }

        /* Read in the A array */
//      fprintf(stdout, "Handling array A\n");
        mxArray *mxA;
        double  *pA;
        mwIndex *ir, *jc;
        mwSize   c, total=0;
        mwIndex  start_row_ind, stop_row_ind, ci;
        mwSize   m, n; // num rows and columns

        int col, row, prev_col;
        double A;
        int num_col, num_val;
        int col_min, col_max, num_tiles, tiles_per_node, num_cols;

        mxA = matGetVariable(pmat, "A");
        m   = mxGetM(mxA);
        n   = mxGetN(mxA);
        tiles_per_node = round((float)m/6./(float)num_nodes);
        col_min = 1 + 6*(mpid * tiles_per_node);
        if (mpid < num_nodes-1) {
            col_max = col_min   + 6*tiles_per_node-1;
        } else {
            col_max = m;                      
        }

        num_tiles = (col_max - col_min +1)/6;
        if (0)
            fprintf(stdout, "\n\nm=%lu n=%lu num_nodes=%d mpid=%d col_min=%d col_max=%d num_tiles=%d\n\n\n", 
                            m, n, num_nodes, mpid, col_min, col_max, num_tiles);
  
        /* Initial allocation */
        int N = 6*num_tiles*(n+1)/2;
       *ncol    = 0;
//      fprintf(stdout, "\n\nAllocating memory: N=%d \n", N);
       *colptrs = CALLOC(pastix_int_t,   N);
       *rows    = CALLOC(pastix_int_t,   N);
       *vals    = CALLOC(pastix_float_t, N);
       *loc2glb = CALLOC(pastix_int_t,   N);

        /* Get the starting positions of all data arrays. */
        pA = mxGetPr(mxA);
        ir = mxGetIr(mxA);
        jc = mxGetJc(mxA);

        /* Estimate the sizes of the matrices to be allocated */
        num_col = 0;
        num_val = 0;
        prev_col = 0;
        for (c=0; c<n; c++)  {
            start_row_ind = jc[c];
            stop_row_ind  = jc[c+1];
            if (start_row_ind == stop_row_ind)
                continue;
            else {
                for (ci = start_row_ind; ci < stop_row_ind; ci++)  {
                    row = ir[ci]+1;
                    col = c+1;
                    A   = pA[total++];

//                  fprintf(stdout, "\t\tmpid=%d", mpid);
//                  mexPrintf("\t(%"FMT_SIZE_T"u,%"FMT_SIZE_T"u) = %g\n",
//                            ir[ci]+1, c+1, pA[total-1]);
  
                    if (num_col < col)
                        num_col = col;

                    if (row >= col + col_min - 1)
                        num_val++;
                    else
                        continue;

                    (*rows)[num_val-1] = row;
                    (*vals)[num_val-1] = A;

                    if (prev_col < col) // new column
                    {
                        prev_col = col;
                        (*colptrs)[col-1] = num_val;
                        (*loc2glb)[col-1] = col + col_min -1;
                        if (0)
                            fprintf(stdout, "\t\tmpid=%d vals[%d]=%f rows[%d]=%d colptrs[%d]=%d loc2glb[%d]=%d\n", 
                                            mpid, num_val-1, (*vals)[num_val-1], num_val-1, 
                                            row, col-1, num_val, col, col+col_min-1);
                    }
                    else 
                        if (0)
                            fprintf(stdout, "\t\tmpid=%d vals[%d]=%f  rows[%d]=%d              \n", 
                                            mpid, num_val-1, (*vals)[num_val-1], num_val-1, row);
                }
            }
        }
        (*colptrs)[num_col] = num_val+1;
//      fprintf(stdout, "              colptrs[%d]=%d\n", num_col, num_val+1);
       *colptrs = REALLOC(*colptrs, pastix_int_t,   num_col + 1);
       *rows    = REALLOC(*rows,    pastix_int_t,   num_val);
       *vals    = REALLOC(*vals,    pastix_float_t, num_val);
       *loc2glb = REALLOC(*loc2glb, pastix_int_t,   num_col);
       *ncol    = num_col;
        if (0) 
        {
            int i;
            FILE *fp;
            char logfilename[256];
            sprintf(logfilename, "node_%d.log", mpid);
            fp = fopen(logfilename, "w");

            fprintf(fp, "num_col=%d num_val=%d\n\ncolptrs= %d", num_col, num_val, (*colptrs)[0]);
            for (i=1; i< num_col+1; i++) {
                fprintf(fp, ", %d", (*colptrs)[i]);
                if ((*colptrs)[i] < (*colptrs)[i-1])
                    fprintf(fp, "colptrs disordered at i=%d : (*colptrs)[i]=%d (*colptrs)[i-1]=%d\n", 
                                i, (*colptrs)[i], (*colptrs)[i-1]);
            }

            fprintf(fp, "\n\nloc2glb= %d", (*loc2glb)[0]);
            for (i=1; i< num_col; i++) {
                fprintf(fp, ", %d", (*loc2glb)[i]);
                if ((*loc2glb)[i] < (*loc2glb)[i-1])
                    fprintf(fp, "loc2glb disordered at i=%d : (*loc2glb)[i]=%d (*loc2glb)[i-1]=%d\n",
                                i, (*loc2glb)[i], (*loc2glb)[i-1]);
            }

            fprintf(fp, "\n\nrows= %d", (*rows)[0]);
            for (i=1; i< num_val; i++)
                fprintf(fp, ", %d", (*rows)[i]);

            fprintf(fp, "\n\nmpid=%d vals= %f3.3", mpid, (*vals)[0]);
            for (i=1; i< num_val; i++)
                fprintf(fp, ", %f3.3", (*vals)[i]);
            fprintf(fp, "\n\n");
            fclose(fp);
        }
        mxDestroyArray(mxA);

        /* Read array b */
//      fprintf(stdout, "Handling array b\n");
        mxArray *mxb;
        mwSize   nb;
        mwIndex  i;

        mxb  = matGetVariable(pmat, "b");
        nb   = (mwSize)mxGetNumberOfElements(mxb);
//      fprintf(stdout, "    nb=%lu\n", nb);
        pastix_float_t *b;
        b   = CALLOC(pastix_float_t, nb);
       *rhs = CALLOC(pastix_float_t, nb);

        b = (pastix_float_t *)mxGetData(mxb);
        for (i=0; i< 6*num_tiles; i++)
            (*rhs)[i] = b[i];

/*      for (i=0; i<nb; i++)
            fprintf(stdout, "b[%d]=%6.6g\n", i+1, (*rhs)[i]);
*/
        mxDestroyArray(mxb);
//      free(b);
        matClose(pmat);
    } else if (num_nodes == 2) {
        /*
           A simple n x n laplacian
          -------------------------
          #   2                  X   1
          #  -1   2              X   0
          #      -1  2           X = 0
          #         -1  2        X   1
          #            -1  2     X   1

          vals   = {2,-1,-1, 2,-1,-1, 2,-1,-1, 2}
          colptr = {1, 3, 6, 9, 11 }                  - added 11!
          row    = {1, 2, 1, 2, 3, 2, 3, 4, 3, 4}
          rhs = {1,0,0,1}
        */
        *mat_type = "API_SYM_YES";

        // Set default input data
        if (mpid == 0) {
            int   my_colptr[4]  = {1, 3, 5, 6};
            int   my_row[6 ]    = {1, 2, 2, 3, 3, 4};               
            float my_vals[6 ]   = {2,-1, 2,-1, 2,-1};
            int   my_loc2glb[3] = {1, 2, 3};
            float my_rhs[3]     = {1, 0, 0};
           *colptrs = CALLOC(pastix_int_t,   4);
           *rows    = CALLOC(pastix_int_t,   6);
           *vals    = CALLOC(pastix_float_t, 6);
           *rhs     = CALLOC(pastix_float_t, 3);
           *loc2glb = CALLOC(pastix_int_t,   3);
           *ncol    = 3;
            int i;
            for (i=0; i< 4; i++)
                (*colptrs)[i] = (pastix_int_t)my_colptr[i];
            for (i=0; i< 6; i++)
                (*rows)[i]    = (pastix_int_t)my_row[i];
            for (i=0; i< 6; i++)
                (*vals)[i]    = (pastix_float_t)my_vals[i];
            for (i=0; i< 3; i++)
                (*rhs)[i]     = (pastix_float_t)my_rhs[i];
            for (i=0; i< 3; i++)
                (*loc2glb)[i] = (pastix_int_t)my_loc2glb[i];
        } else {
            int   my_colptr[3]  = {1, 3, 4};
            int   my_row[3 ]    = {4, 5, 5};
            float my_vals[3 ]   = {2,-1, 2};            
            int   my_loc2glb[2] = {4, 5};              
            float my_rhs[2]     = {1, 1};
           *colptrs = CALLOC(pastix_int_t,   3);
           *rows    = CALLOC(pastix_int_t,   3);
           *vals    = CALLOC(pastix_float_t, 3);
           *rhs     = CALLOC(pastix_float_t, 2);
           *loc2glb = CALLOC(pastix_int_t,   2);
           *ncol    = 2;
            int i;
            for (i=0; i< 3; i++)
                (*colptrs)[i] = (pastix_int_t)my_colptr[i];
            for (i=0; i< 3; i++)
                (*rows)[i]    = (pastix_int_t)my_row[i];
            for (i=0; i< 3; i++)
                (*vals)[i]    = (pastix_float_t)my_vals[i];
            for (i=0; i< 2; i++)
                (*rhs)[i]     = (pastix_float_t)my_rhs[i];
            for (i=0; i< 2; i++)
                (*loc2glb)[i] = (pastix_int_t)my_loc2glb[i];
        }
    } else {
        fprintf(stdout, "filename=%s\n", filename);
    }
    return 0;
}

static void read_params_file(char *params_file, int *nbthread, int *incomplete,
                             int *level_of_fill, int *amalgamation, int *verbose,
                             int *ooc)
{
    FILE *fp;
    int   MAXLINE = 256;
    char linebuf[MAXLINE], *s, *s1, *s2;
    if ((fp = fopen(params_file, "r")) == NULL) {
        return;
    }

    while (fgets(linebuf, MAXLINE, fp) != NULL) {
        if (strspn(linebuf, " \t\r\n") == strlen(linebuf)) {
            continue;
        }

        if (((linebuf[0] == '/') && (linebuf[1] == '*'))
          || (linebuf[0] == ';')
          || (linebuf[0] == '#'))
        {
            continue;
        }

        s = linebuf;
        s1 = strtok(s, " \t\n");
        if (s1 == NULL)
            goto error_return;

        s2 = strtok(NULL, " \t\n");
        if (s2 == NULL)
            goto error_return;

        if (!strcmp(s1, "IPARM_THREAD_NBR"))
            *nbthread = atoi(s2);
        else if (!strcmp(s1, "IPARM_INCOMPLETE"))
            *incomplete = atoi(s2);
        else if (!strcmp(s1, "IPARM_LEVEL_OF_FILL"))
            *level_of_fill = atoi(s2);
        else if (!strcmp(s1, "IPARM_OOC_LIMIT"))
            *ooc = atoi(s2);
        else if (!strcmp(s1, "IPARM_AMALGAMATION_LEVEL"))
            *amalgamation = atoi(s2);
        else if (!strcmp(s1, "IPARM_VERBOSE"))
            *verbose = atoi(s2);
    }

error_return:
    (void)fclose(fp);
    return;

}

/*
 *   Function: get_options
 *
 *     Get options from argv.
 *
 *       Parameters:
 *           argc          - number of arguments.
 *           argv          - argument tabular.
 *           driver_type   - type of driver (output).
 *           filename      - Matrix filename (output).
 *           nbthread      - number of thread (left unchanged if not in options).
 *           verbose       - verbose level 1,2 or 3
 *           ordering      - ordering to choose (see <API_ORDER>).
 *           incomplete    - indicate if -incomp is present
 *           level_of_fill - Level of fill for incomplete factorization.
 *           amalgamation  - Amalgamation for kass.
 *           ooc           - Out-of-core limite (Mo or percent depending on compilation option)
 *           size          - Size of the matrix (generated matrix only)
 *
 *                                                        */
int get_options(int              argc,
                char           **argv,
                driver_type_t  **driver_type,
                char          ***filename,
                int             *nbmatrices,
                int             *nbthread,
                int             *verbose,
                int             *ordering,
                int             *incomplete,
                int             *level_of_fill,
                int             *amalgamation,
                int             *ooc,
                pastix_int_t    *size)
{
    int i;
   *driver_type = CALLOC(driver_type_t, 1);
  (*driver_type)[0] = RSA;
   *nbmatrices = 1;
   *nbthread = 8;
   *incomplete = 0;
   *level_of_fill = -1;
   *verbose = 1;
   *ooc = 50000;
   *amalgamation = 5;
   *filename = CALLOC(char *, *nbmatrices);
    if (argc == 1)
        (*filename)[0] = NULL;
    else {
        for (i=0; i<argc-1; i++) {
            (*filename)[i] = CALLOC(char, strlen(argv[i+1]));
            strcpy((*filename)[i], argv[i+1]);
        }
    }

    if (argc == 3) {
        read_params_file((*filename)[1], nbthread, incomplete,
                         level_of_fill, amalgamation, verbose, ooc);
    }

   return 0;
}


