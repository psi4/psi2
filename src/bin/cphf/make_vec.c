
/* $Log: make_vec.c,v $
/* Revision 1.1  1991/06/15 22:45:28  seidl
/* Initial revision
/* */

static char *rcsid = "$Id: make_vec.c,v 1.1 1991/06/15 22:45:28 seidl Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

make_vec()
   {
      int i,k,l,ii;
      double *esi;

      so_vecs_k = (double **) init_matrix(nbfso,nind);
      so_vecs_l = (double **) init_matrix(nbfso,nind);

      for(i=0; i < nbfso ; i++) {
         esi = e_vecs_so[i];
         for(ii=0; ii < nind ; ii++) {
            k=indep[ii].ii;
            l=indep[ii].jj;
            so_vecs_k[i][ii] = esi[k];
            so_vecs_l[i][ii] = esi[l];
            }
         }
      }
