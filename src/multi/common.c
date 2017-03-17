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

#include "common.h"

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

  if (argc == 1) {
      usage(argc, argv);
      exit(2);
  }

  /*******************************************/
  /*          MPI initialisation             */
  /*******************************************/
#ifndef FORCE_NOMPI
  required = MPI_THREAD_MULTIPLE;
  provided = -1;
  MPI_Init_thread(&argc, &argv, required, &provided);
  MPI_Comm_size (MPI_COMM_WORLD, &num_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
//fprintf(stdout, "\nnum_nodes=%d mpid=%d sizeof(int)=%d\n\n", num_nodes, mpid, sizeof(int));
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
  iparm[IPARM_ORDERING]            = API_ORDER_SCOTCH; // options: API_ORDER_METIS, API_ORDER_PERSONAL, API_ORDER_LOAD
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

  /* To save intermediate results */
//iparm[IPARM_IO_STRATEGY] = API_IO_SAVE; // save steps 1 and 2; need to link ordergen and symbolgena to ordername and symbolname
//iparm[IPARM_IO_STRATEGY] = API_IO_LOAD; // start from step 3
//
  /* reread parameters to set IPARM/DPARM */
  if (EXIT_FAILURE == get_idparm(argc, argv, iparm, dparm))
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
  PRINT_RHS("SOL", rhs, ncol, mpid, iparm[IPARM_VERBOSE]);
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
  if (getenv("TA_DATA") != NULL) {
      sprintf(mat_file_name, "%s/x_%d.mat", getenv("TA_DATA"), mpid +1);
  } else {
      char cwd[1024];
      getcwd(cwd, sizeof(cwd));
      printf("Writing solution to %s\n", cwd);
      sprintf(mat_file_name,  "%s/x_%d.mat", cwd             , mpid +1);
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
  free(loc2glob);
#ifndef FORCE_NOMPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */

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
                int             *verbosemode,
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
   *nbthread = 32;
   *incomplete = 0;
   *level_of_fill = -1;
   *verbosemode = 1;
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
                         level_of_fill, amalgamation, verbosemode, ooc);
    }

   return 0;
}

/* -------------------------------------------------------------------------- */

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

/* -------------------------------------------------------------------------- */

int write_matlab_solution(char *fileName, char *varName, int nx, double *xvals)
{
    mxArray *xv;

    MATFile *pmat = matOpen(fileName, "w");

    xv = mxCreateDoubleMatrix(1, nx, mxREAL);

    memcpy((void *)(mxGetPr(xv)), (void *)xvals, nx*sizeof(xvals));

    matPutVariable(pmat, varName, xv);

    matClose(pmat);
}

static void
usage(int argc, char *argv[])
{
    fprintf( stderr,
    "usage: mpirun -np <num_slots> %s <mat_file> [ <params_file> ]\n"
            , argv[0] );
}

/* -------------------------------------------------------------------------- */

void read_params_file(char *params_file, int *nbthread, int *incomplete,
                             int *level_of_fill, int *amalgamation,
                             int *verbosemode, int *ooc)
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
            *verbosemode = atoi(s2);
    }

error_return:
    (void)fclose(fp);
    return;

}

