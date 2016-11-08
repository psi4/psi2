/* $Log: shalf.c,v $
/* Revision 1.1  1991/06/15 20:22:43  seidl
/* Initial revision
/* */

static char *rcsid = "$Id: shalf.c,v 1.1 1991/06/15 20:22:43 seidl Exp $";

/* construct S-1/2 matrix 'sahalf'   */

#define EXTERN
#include "includes.h"
#include "common.h"

void shalf()

  {
     int i,nn;
     int ii,jj,kk;
     double *eig_vals, **eig_vecs;
     double tol = 1.0e-20;
     struct symm *s;
     
     eig_vals = (double *) init_array(nsfmax);
     eig_vecs = (double **) init_matrix(nsfmax,nsfmax);

/*  diagonalize smat to get eigenvalues and eigenvectors  */

     for (i=0; i < num_ir ; i++) {
        s = &scf_info[i];
        if (nn=s->num_so) {

           rsp(nn,nn,ioff[nn],s->smat,eig_vals,1,eig_vecs,tol);

           if(print & 64) {
             fprintf(outfile,"\noverlap eigenstuff\n");
             eivout(eig_vecs,eig_vals,nn,nn,outfile);
             }

/* form 'sahalf' matrix sahalf = u*ei**-0.5*u~  */

           for (ii=0; ii < nn ; ii++)
              for (jj=0; jj < nn ; jj++)
                 for (kk=0; kk < nn ; kk++)
                     s->sahalf[ii][jj] += 
                      eig_vecs[ii][kk]*eig_vecs[jj][kk]/sqrt(eig_vals[kk]);

           }
        }
   free(eig_vals);
   free(eig_vecs);
   }
