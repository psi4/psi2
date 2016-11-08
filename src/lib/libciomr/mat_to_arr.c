
/* $Log: mat_to_arr.c,v $
/* Revision 2.1  1991/06/15 18:29:28  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id: mat_to_arr.c,v 2.1 1991/06/15 18:29:28 seidl Exp $";

#include "includes.h"

void mat_to_arr(a,b,m,n)
   double **a, *b;
   int m,n;

   {
      int i,j,ij;

      ij=0;
      for (j=0; j < n ; j++) {
         for (i=0; i < m; i++) {
            b[ij++]=a[i][j];
            }
         }
      }

void arr_to_mat(a,b,m,n)
   double **a, *b;
   int m,n;

   {
      int i,j,ij;

      ij=0;
      for (j=0; j < n ; j++) {
         for (i=0; i < m; i++) {
            a[i][j]=b[ij++];
            }
         }
      }

