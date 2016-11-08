/* $Log: form_vec.c,v $
/* Revision 1.2  1998/06/30 14:11:06  sbrown
/* *************************************************************
/* *Program Modification                                       *
/* *By: Shawn Brown                                            *
/* *Date: June 30, 1998                                        *
/* *Altered program to make a guess at occupations from the    *
/* *diagonalized core hamiltonian matrix.  Program can now     *
/* *make a guess at the beginning of the calculation or at     *
/* *or at every iteration.  Use the latter at your own risk.   *
/* *See man pages for details on new keywords.                 *
/* *************************************************************
/*
 * Revision 1.1  1991/06/15  20:22:26  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: form_vec.c,v 1.2 1998/06/30 14:11:06 sbrown Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"
#include "libciomr.h"

void form_vec()
{
   int i,nn;
   int k,l;
   double **ctrans;
   double **temp;
   double **sqhmat;
   double **sqhmat2;
   double tol=1.0e-20;
   struct symm *s;

   ctrans = (double **) init_matrix(nsfmax,nsfmax); 
   temp = (double **) init_matrix(nsfmax,nsfmax);
   sqhmat = (double **) init_matrix(nsfmax,nsfmax);
   sqhmat2 = (double **) init_matrix(nsfmax,nsfmax);
   
   for (i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      if (nn=s->num_so) {
         tri_to_sq(s->hmat,sqhmat,nn);
         mmult(s->sahalf,1,sqhmat,0,temp,0,nn,nn,nn,0);
         mmult(temp,0,s->sahalf,0,sqhmat2,0,nn,nn,nn,0);
         sq_rsp(nn,nn,sqhmat2,s->hevals,1,ctrans,tol);
         mxmb(s->sahalf,1,nn,ctrans,1,nn,s->cmat,1,nn,nn,nn,nn);
         for(k=0;k < nn; k++)
            for(l=0;l < nn; l++) s->sahalf[k][l]=s->cmat[k][l];

         if(print & 2) {
            fprintf(outfile,"\nguess vector for irrep %s\n",s->irrep_label);
            print_mat(s->cmat,nn,nn,outfile);
            }
         }
      }
   inflg = 0;
   free_matrix(ctrans,nsfmax);
   free_matrix(temp,nsfmax);
   free_matrix(sqhmat,nsfmax);
   free_matrix(sqhmat2,nsfmax);
   }
