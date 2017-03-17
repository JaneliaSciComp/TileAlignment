/* Copyright (c) 2017, HHMI-Janelia Research Campus All rights reserved. */

#include "common.h"

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


