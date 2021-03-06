
/* $Log: eigout.c,v $
/* Revision 2.3  1995/01/17 19:56:52  local
/* Fixed undeclared variable bug.  Minor modification.
/*
 * Revision 2.2  1995/01/16  22:48:32  cljanss
 * Minor changes to make the SGI compiler happy.
 *
 * Revision 2.1  1991/06/15  18:28:46  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id: eigout.c,v 2.3 1995/01/17 19:56:52 local Exp $";

#include "includes.h"

void eigout(a,b,c,m,n,out)
   double **a, *b, *c;
   int m,n;
   FILE *out;

   {
      int ii,jj,kk,nn;
      int i,j;

      ii=0;jj=0;
L200:
      ii++;
      jj++;
      kk=10*jj;
      nn=n;
      if (nn > kk) nn=kk;
      fprintf (out,"\n");
      for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
      fprintf (out,"\n");
      for (i=0; i < m; i++) {
         fprintf (out,"\n%5d",i+1);
         for (j=ii-1; j < nn; j++) {
            fprintf (out,"%12.7f",a[i][j]);
            }
         }
      fprintf (out,"\n");
      fprintf (out,"\n     ");
      for (j=ii-1; j < nn; j++) {
         fprintf(out,"%12.7f",b[j]);
         }
      fprintf (out,"\n");
      fprintf (out,"\n     ");
      for (j=ii-1; j < nn; j++) {
         fprintf(out,"%12.7f",c[j]);
         }
      fprintf (out,"\n");
      if (n <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
      }
