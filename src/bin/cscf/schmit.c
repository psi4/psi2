/* $Log: schmit.c,v $
/* Revision 1.1  1991/06/15 20:22:42  seidl
/* Initial revision
/* */

static char *rcsid = "$Id: schmit.c,v 1.1 1991/06/15 20:22:42 seidl Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

void schmit(all)
   int all;

{
   int i,j,ij,nn;
   int n,m,ncol;
   double *v,**ctmp,vtmp;
   struct symm *s;

   v = (double *) init_array(nsfmax);
   ctmp = (double **) init_matrix(nsfmax,nsfmax);

   for(n=0; n < num_ir ; n++) {
      s= &scf_info[n];
      if(nn=s->num_so) {
         for (i=0; i < nn ; i++)
            for (j=0; j < nn ; j++)
               ctmp[j][i] = s->cmat[i][j];

         ncol = s->nclosed + s->nopen;
         if(s->nhalf) ncol++;
         if(!ncol) ncol++;
         if(all) ncol = nn;
         for(m=0; m < ncol ; m++) {
            v[0]=ctmp[m][0]*s->smat[0];
            for(i=1; i < nn ; i++) {
               for(j=0,vtmp=0.0; j < i ; j++) {
                  ij=ioff[i]+j;
                  vtmp += ctmp[m][j]*s->smat[ij];
                  v[j] += ctmp[m][i]*s->smat[ij];
                  }
               v[i] = vtmp+ctmp[m][i]*s->smat[ioff[i]+j];
               }
            for(i=0,vtmp=0.0; i < nn ; i++) vtmp += v[i]*ctmp[m][i];
            if(!vtmp) exit(23);
            if(vtmp < 10.0e-20) vtmp = 10.0e-20;
            vtmp = 1.0/sqrt(vtmp);

            for(i=0; i < nn ; i++) {
               v[i] *= vtmp;
               ctmp[m][i] *= vtmp;
               }

            if(m < ncol-1) {
               for(i=m+1,vtmp=0.0; i < ncol ; i++) {
                  for(j=0,vtmp=0.0; j<nn ;j++) vtmp += v[j]*ctmp[i][j];
                  for(j=0; j < nn ; j++) ctmp[i][j] -=
                                                      vtmp*ctmp[m][j];
                  }
               }
            }

         for (i=0; i < nn ; i++)
            for (j=0; j < nn ; j++)
               s->cmat[i][j] = ctmp[j][i];

         }
      }
   free(v);
   free_matrix(ctmp,nsfmax);
   }
