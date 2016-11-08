/* $Log: dmat_2.c,v $
/* Revision 1.1  1991/06/15 20:22:22  seidl
/* Initial revision
/* */

static char *rcsid = "$Id: dmat_2.c,v 1.1 1991/06/15 20:22:22 seidl Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

void dmat_2(opblk)
   int opblk;

{
   int i,j,k,l,ij,n;
   int ndocc;
   double ptempc,ctmp;
   struct symm *s;

   for (l=0; l < num_ir ; l++) {
      s = &scf_info[l];
      if (n=s->num_so) {
         ndocc = s->nclosed;
         if (l == opblk) ndocc++;
         
         ij=0;
         for (i=0; i < n ; i++) {
            for (j=0; j < i; j++,ij++) {
               ptempc=0.0;
               for (k=0; k < ndocc ; k++) {
                  ptempc += 4.0*s->cmat[i][k]*s->cmat[j][k];
                  }
               if(opblk==opblk1) s->pmato2[ij] = ptempc;
               else s->pmat2[ij] = ptempc;
               }
            ptempc=0.0;
            for (k=0; k < ndocc ; k++) {
               ctmp=s->cmat[i][k];
               ptempc += 2.0*ctmp*ctmp;
               }
            if(opblk==opblk1) s->pmato2[ij] = ptempc;
            else s->pmat2[ij] = ptempc;
            ij++;
            }
         if(print & 4) {
            if(opblk==opblk1) {
               fprintf(outfile,
                       "\ndensity matrix 1 for irrep %s",s->irrep_label);
               print_array(s->pmato2,n,outfile);
               }
            else {
               fprintf(outfile,
                       "\ndensity matrix 1 for irrep %s",s->irrep_label);
               print_array(s->pmat2,n,outfile);
               }
            }
         }
      }
   }
