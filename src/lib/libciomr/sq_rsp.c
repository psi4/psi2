
/* $Log: sq_rsp.c,v $
/* Revision 2.3  1999/11/01 20:11:00  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.2  1998/02/03 19:34:11  evaleev
/* Modified eigsort(), rsp(), and sq_rsp() to sort eigenvalues and
/* eigenvectors in either ascending OR descending order.
/*
 * Revision 2.1  1991/06/15  18:30:07  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id: sq_rsp.c,v 2.3 1999/11/01 20:11:00 evaleev Exp $";

#include "includes.h"

extern double *init_array(int);
extern double **init_matrix(int,int);
extern void free_matrix(double **, int);
extern void tred2(int, double **, double *, double *, int);
extern int tqli(int, double *, double **, double *, int, double);
extern void eigsort(double *, double **, int);

/* translation into c of a translation into FORTRAN77 of the EISPACK */
/* matrix diagonalization routines */

void sq_rsp(nm,n,array,e_vals,matz,e_vecs,toler)
   int nm, n,matz;
   double **array, *e_vals, **e_vecs;
   double toler;

   {
      int i, j, ii, ij, ierr;
      int ascend_order;
      double *fv1, **temp;
      double zero = 0.0;
      double one = 1.0;

/* Modified by Ed - matz can have the values 0 through 3 */
      
      if ((matz > 3) || (matz < 0)) {
        matz = 0;
        ascend_order = 1;
        }
      else
        if (matz < 2)
          ascend_order = 1;	/* Eigenvalues in ascending order */
        else {
          matz -= 2;
          ascend_order = 0;	/* Eigenvalues in descending order */
          }

      fv1 = (double *) init_array(n);
      temp = (double **) init_matrix(n,n);

      if (n > nm) {
         ierr = 10*n;
         fprintf(stderr,"n = %d is greater than nm = %d in rsp\n",n,nm);
         exit(ierr);
         }

      for (i=0; i < n; i++) {
         for (j=0; j < n; j++) {
            e_vecs[i][j] = array[i][j];
            }
          }

      tred2(n,e_vecs,e_vals,fv1,matz);

      for (i=0; i < n; i++)
         for (j=0; j < n; j++)
            temp[i][j]=e_vecs[j][i];

      tqli(n,e_vals,temp,fv1,matz,toler);

      for (i=0; i < n; i++)
         for (j=0; j < n; j++)
            e_vecs[i][j]=temp[j][i];

      if (ascend_order)
        eigsort(e_vals,e_vecs,n);
      else
        eigsort(e_vals,e_vecs,(-1)*n);

      free(fv1);
      free_matrix(temp,n);

      }
            
