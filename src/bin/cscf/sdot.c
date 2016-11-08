/* $Log: sdot.c,v $
/* Revision 1.1  1991/06/15 20:22:42  seidl
/* Initial revision
/* */

static char *rcsid = "$Id: sdot.c,v 1.1 1991/06/15 20:22:42 seidl Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

void sdot(a,b,n,value)
   double **a, **b, *value;
   int n;

   {
      register int i,j;
      double *ta, *tb, tval;

      tval = 0.0;
      for (i=0; i < n; i++) {
         ta = a[i];
         tb = b[i];
         for (j=0; j <= i; j++,ta++,tb++) {
            tval += (*ta) * (*tb);
            }
         }
      *value = tval;
      }

void vdot(a,b,n,value)
  double *a,*b,*value;
  int n;
{
  int i;
  double tval=0.0;

  for(i=0; i < n ; i++) tval += a[i]*b[i];

  *value = tval;
  }
