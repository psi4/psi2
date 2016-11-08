
/* $Log: block_alloc.c,v $
/* Revision 2.3  1994/06/02 02:30:26  seidl
/* add test for too many rows
/*
 * Revision 1.1.1.1  1994/05/02  17:04:26  cljanss
 * The May 1, 1994 version of psi as on the CCQC machines.
 *
 * Revision 2.2  1991/06/15  18:53:37  seidl
 * remove include of common.h and definition of EXTERN
 *
 * Revision 2.1  1991/06/15  18:28:44  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id: block_alloc.c,v 2.3 1994/06/02 02:30:26 seidl Exp $";

#include "includes.h"

double *** block_mat_alloc(n_so_typs,num_ir,num_so)
   int n_so_typs,num_ir,*num_so;

   {
      double ***array;
      int i,blk,nn;

      if ((array = (double ***) malloc(sizeof(double **)*n_so_typs))==NULL) {
          fprintf(stderr,"trouble in block_mat_alloc\n");
          exit(3);
          }

      for (i=0,blk=0; i < num_ir ; i++) {
         if (nn=num_so[i]) {
            array[blk] = (double **) init_matrix(nn,nn);
            blk++;
            }
         }

      return(array);
      }

void block_mat_dealloc(array,num_ir,num_so)
   double ***array;
   int num_ir, *num_so;

   {
      int i;
      int blk=0;

      for (i=0; i < num_ir ; i++) {
         if (num_so[i]) {
            free_matrix(array[blk],num_so[i]);
            blk++;
            }
         }

      free(array);
      }

double ** block_arr_alloc(n_so_typs,num_ir,num_so)
   int n_so_typs,num_ir,*num_so;

   {
      double **array;
      int i;
      int j=0;

      if ((array = (double **) malloc(sizeof(double *)*n_so_typs))==NULL) {
          fprintf(stderr,"trouble in block_arr_alloc\n");
          exit(3);
          }

      for (i=0; i < num_ir ; i++) {
         if (num_so[i]) {
            int nget = num_so[i]*(num_so[i]+1)/2;
            if (j>=n_so_typs) {
                fprintf(stderr,"block_arr_alloc: too many rows\n");
                exit(3);
              }
            array[j] = (double *) init_array(nget);
            j++;
            }
         }
      return(array);
      }

void block_arr_dealloc(array,n_so_typs)
   double **array;

   {
      int i;

      for (i=0; i < n_so_typs ; i++) {
         free(array[i]);
         }

      free(array);
      }

