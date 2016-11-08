/* $Log: dmat.c,v $
/* Revision 1.1  1991/06/15 20:22:21  seidl
/* Initial revision
/* */

static char *rcsid = "$Id: dmat.c,v 1.1 1991/06/15 20:22:21 seidl Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

void dmat()

{
   int i,j,k,l,ij,n;
   int ndocc,nsocc,nhocc;
   double ptempc,ptempo,ctmp;
   struct symm *s;

   for (l=0; l < num_ir ; l++) {
      s = &scf_info[l];
      if (n=s->num_so) {
         ndocc = s->nclosed;
         nsocc = s->nopen;
         nhocc = 0;
         if(s->nhalf) nhocc=1;
            
         for (i=ij=0; i < n ; i++ ) {
            for (j=0; j < i; j++,ij++) {
               ptempc=ptempo=0.0;
               for (k=0; k < ndocc ; k++)
                  ptempc += 4.0*s->cmat[i][k]*s->cmat[j][k];

               for (k=ndocc; k < ndocc+nsocc ; k++)
                  ptempo += 2.0*s->occ_num[k]*s->cmat[i][k]*s->cmat[j][k];

               for (k=ndocc+nsocc; k < ndocc+nsocc+nhocc ; k++)
                  ptempo += 2.0*s->occ_num[k]*s->cmat[i][k]*s->cmat[j][k];

               if(iopen) {
                  s->dpmato[ij] = ptempo - s->pmato[ij];
                  s->pmato[ij] = ptempo;
                  }
               s->dpmat[ij] = ptempc+ptempo - s->pmat[ij];
               s->pmat[ij] = ptempc+ptempo;
               }
            ptempc=ptempo=0.0;
            for (k=0; k < ndocc ; k++) {
               ctmp=s->cmat[i][k];
               ptempc += 2.0*ctmp*ctmp;
               }
            for (k=ndocc; k < ndocc+nsocc ; k++) {
               ctmp=s->cmat[i][k];
               ptempo += ctmp*ctmp*s->occ_num[k];
               }
            for (k=ndocc+nsocc; k < ndocc+nsocc+nhocc ; k++) {
               ctmp=s->cmat[i][k];
               ptempo += ctmp*ctmp*s->occ_num[k];
               }
            if(iopen) {
               s->dpmato[ij] = ptempo - s->pmato[ij];
               s->pmato[ij] = ptempo;
               }
            s->dpmat[ij] = ptempc+ptempo - s->pmat[ij];
            s->pmat[ij] = ptempc+ptempo;
            ij++;
            }

         if(print & 4) {
            fprintf(outfile,
                       "\ntotal density matrix for irrep %s",s->irrep_label);
            print_array(s->pmat,n,outfile);
            print_array(s->dpmat,n,outfile);
            if(iopen) {
               fprintf(outfile,"\nopen-shell density matrix for irrep %s",
                                                              s->irrep_label);
               print_array(s->pmato,n,outfile);
               print_array(s->dpmato,n,outfile);
               }
            }
         }
      }
   }
