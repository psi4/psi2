/* $Log: rotate_vector.c,v $
/* Revision 1.6  1996/11/26 05:20:57  sherrill
/* Add casts in front of init_int_array() calls to avoid compiler warnings
/* (really should include libciomr everywhere instead but that causes more
/* warnings because wwritw's need casts to char *).  Also fixed problem
/* where phase() was trying to free unallocated memory.
/*
 * Revision 1.5  1995/11/07  19:24:13  sherrill
 * Replace all initializations of integer arrays with calls to init_int_array().
 * There had been problems with using init_array() because the division by
 * the doubles-to-int factor was rounding off (integer division).
 *
 * Revision 1.4  1995/07/26  20:33:01  sherrill
 * Added support for getting second root of TCSCF and for checking the
 * orbital rotation to make sure orbitals are not swapped.
 *
 * Revision 1.3  1995/07/21  17:37:12  psi
 * Made Jan 1st 1995 cscf the current accepted version of cscf.  Some
 * unidentified changes made after that date were causing problems.
 *
 * Revision 1.1  1991/06/15  20:22:38  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: rotate_vector.c,v 1.6 1996/11/26 05:20:57 sherrill Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

void check_rot(int nn, double **cold, double **cnew, double *smat_pac,
      double *fock_evals, int irrep);


rotate_vector()
{
   int i,j,ij,k,nn;
   int nc,no;
   int joff,oj;
   int *optest,iju=0;
   double **scr1,**scr2,**scr3;
   double tol=1.0e-15;
   struct symm *s;

   scr1 = (double **) init_matrix(nsfmax,nsfmax);
   scr2 = (double **) init_matrix(nsfmax,nsfmax);
   scr3 = (double **) init_matrix(nsfmax,nsfmax);

   dmat();

   if(twocon) {
      /* optest = (int *) init_array(ioff[nbasis]/2); */
      optest = (int *) init_int_array(ioff[nbasis]);

/* find open shells */

      for (i=0; i < num_ir ; i++) {
         iju += scf_info[i].num_so;
         if(scf_info[i].nopen) {
            iju = ioff[iju];
            break;
            }
         }

/* set up array of flags indicating open shells */

      for (k=0,joff=0; k < num_ir ; k++) {
         s = &scf_info[k];
         if (nn=s->num_so) {
            for (i=0; i < nn ; i++)
               for (j=0; j <= i ; j++)
                  if(s->nopen) optest[ioff[i+joff]+j+joff] = 1;
            joff += nn;
            }
         }

      formg_two(iju,optest);
      }
   else formg_open();

   for(i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      if(nn=s->num_so) {
         nc=s->nclosed;
         no=s->nopen;
         add_arr(s->hmat,s->gmat,s->fock_pac,ioff[nn]);

         tri_to_sq(s->fock_pac,scr1,nn);
         mmult(s->cmat,1,scr1,0,scr2,0,nn,nn,nn,0);
         mmult(scr2,0,s->cmat,0,scr1,0,nn,nn,nn,0);

         zero_mat(scr2,nn,nn);
         for(j=0; j < nc ; j++)
            for(k=0; k < nc ; k++)
               scr2[j][k]=scr1[j][k];
         for(j=nc; j < nc+no ; j++)
            for(k=nc; k < nc+no ; k++)
               scr2[j][k]=scr1[j][k];
         for(j=nc+no; j < nn ; j++)
            for(k=nc+no; k < nn ; k++)
               scr2[j][k]=scr1[j][k];

         sq_rsp(nn,nn,scr2,s->fock_evals,1,scr3,tol);

         mmult(s->cmat,0,scr3,0,scr2,0,nn,nn,nn,0);

         if (icheck_rot) {
            check_rot(nn, s->cmat, scr2, s->smat, s->fock_evals, i);
            } 

         for(j=0; j < nn ; j++)
            for(k=0; k < nn ; k++)
               s->cmat[j][k]=scr2[j][k];

         }
      }
   free_matrix(scr1,nsfmax);
   free_matrix(scr2,nsfmax);
   free_matrix(scr3,nsfmax);
   }
