/* $Log: scf_iter.c,v $
/* Revision 1.5  1998/06/30 14:11:12  sbrown
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
 * Revision 1.4  1995/07/21  17:37:15  psi
 * Made Jan 1st 1995 cscf the current accepted version of cscf.  Some
 * unidentified changes made after that date were causing problems.
 *
 * Revision 1.1  1991/06/15  20:22:40  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: scf_iter.c,v 1.5 1998/06/30 14:11:12 sbrown Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

void scf_iter()

{
   int i,j,m,ij;
   int nn,newci;
   double cimax;
   double occi, occj, occ0, den;
   double **scr;
   double *c1, *c2;
   double **fock_c, **fock_o;
   double **fock_ct;
   double **ctrans;
   double tol = 1.0e-14;
   struct symm *s;

   diiser = 0.0;
   scr = (double **) init_matrix(nsfmax,nsfmax);
   fock_c = (double **) init_matrix(nsfmax,nsfmax);
   fock_ct = (double **) init_matrix(nsfmax,nsfmax);
   ctrans = (double **) init_matrix(nsfmax,nsfmax);
   if(iopen) {
     c1 = (double *) init_array(num_ir);
     c2 = (double *) init_array(num_ir);
     fock_o = (double **) init_matrix(nsfmax,nsfmax);
     }

/* set up c1 and c2  */
   if(iopen) {
      for (i=0; i < num_ir ; i++) {
         s = &scf_info[i];
         if (nn=s->num_so) {
            for (j=0; j < nn ; j++) {
               if(s->occ_num[j]==2.0) c1[i]=2.0;
               if(s->occ_num[j]==1.0) c2[i]=1.0;
               if(s->occ_num[j]==0.5) c2[i]=0.5;
               if(s->occ_num[j]==1.5) c2[i]=1.5;
               }
            den = c1[i]-c2[i];
            if(den) {
               c1[i] /= den;
               c2[i] /= den;
               }
            }
         }
     }

/* and iterate */

   for (iter=0; iter < itmax ; ) {
       if(iter !=0  && inflg == 0 && occ_fix == 1){
	   guess();

/* Have to do this all again if the occupations are switched */
	   if(iopen) {
	       for (i=0; i < num_ir ; i++) {
		   s = &scf_info[i];
		   if (nn=s->num_so) {
		       for (j=0; j < nn ; j++) {
			   if(s->occ_num[j]==2.0) c1[i]=2.0;
			   if(s->occ_num[j]==1.0) c2[i]=1.0;
			   if(s->occ_num[j]==0.5) c2[i]=0.5;
			   if(s->occ_num[j]==1.5) c2[i]=1.5;
			 }
		       den = c1[i]-c2[i];
		       if(den) {
			   c1[i] /= den;
			   c2[i] /= den;
			 }
		     }
		 }
	     } 
	 }
       for (m=0; m < num_ir ; m++) {
         s = &scf_info[m];
         if (nn=s->num_so) {

     /*  form fock matrix = h+g */
            add_arr(s->hmat,s->gmat,s->fock_pac,ioff[nn]);

     /* for open shell, form fock_open = h+g-q */

            if(iopen) {
               for (i=0; i < ioff[nn] ; i++)
                  s->fock_open[i]=s->fock_pac[i]-s->gmato[i];
               }
            }
         }

      ecalc(tol);

   /* create new fock matrix in fock_pac or fock_eff */
      if(!diisflg) diis(scr,fock_c,fock_ct,c1,c2,cimax,newci);

      if(iopen) {
         for (m=0; m < num_ir ; m++) {
            s = &scf_info[m];
            if (nn=s->num_so) {

            /* transform fock_pac to mo basis */
               tri_to_sq(s->fock_pac,fock_ct,nn);
               mxmb(s->cmat,nn,1,fock_ct,1,nn,scr,1,nn,nn,nn,nn);
               mxmb(scr,1,nn,s->cmat,1,nn,fock_c,1,nn,nn,nn,nn);

            /* transform fock_open to mo basis */
               tri_to_sq(s->fock_open,fock_ct,nn);
               mxmb(s->cmat,nn,1,fock_ct,1,nn,scr,1,nn,nn,nn,nn);
               mxmb(scr,1,nn,s->cmat,1,nn,fock_o,1,nn,nn,nn,nn);

            /* form effective fock matrix in mo basis */

               ij=0;
               occ0 = s->occ_num[0];
               for (i=0; i < nn; i++ ) {
                  for (j=0; j <= i; j++) {
                     occi = s->occ_num[i];
                     occj = s->occ_num[j];

              /* default: Guest & Saunders general form */
                     if(iter < itmax-1 && !converged && !fock_typ) {
                        if(occi == occj) 
                          s->fock_eff[ij] = fock_c[i][j];
                        else if(occi)
                          s->fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
                        else {
                          if(occj==2.0)
                            s->fock_eff[ij] = fock_c[i][j];
                          else
                            s->fock_eff[ij] = fock_o[i][j];
                          }
                        }

               /* Guest & Saunders' form for high spin */
                     else if(iter < itmax-1 && !converged && fock_typ == 1) {
                        if (occi == occj || occi)
                           s->fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
                        else if (occj == 2.0) s->fock_eff[ij] = fock_c[i][j];
                        else s->fock_eff[ij] = fock_o[i][j];
                        }

               /* test form (fo fo fo) */
                     else if(iter < itmax-1 && !converged && fock_typ == 2) {
                        if (occi == occj) s->fock_eff[ij] = fock_o[i][j];
                        else if(occi)
                           s->fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
                        else if(occj == 2.0)
                           s->fock_eff[ij] = fock_c[i][j];
                        else s->fock_eff[ij] = fock_o[i][j];
                        }

               /* test form a*(fc fc fc) */
                     else if(iter < itmax-1 && !converged && fock_typ == 3) {
                        if (occi == occj) s->fock_eff[ij] = dampd*fock_c[i][j];
                        else if(occi)
                           s->fock_eff[ij] = 
                             dampo*(2.0*fock_c[i][j]-fock_o[i][j]);
                        else if(occj == 2.0)
                           s->fock_eff[ij] = dampo*fock_c[i][j];
                        else s->fock_eff[ij] = dampo*fock_o[i][j];
                        }

               /* test form a*(fo fo fo) */
                     else if(iter < itmax-1 && !converged && fock_typ == 4) {
                        if (occi == occj) s->fock_eff[ij] = dampd*fock_o[i][j];
                        else if(occi)
                           s->fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
                        else if(occj == 2.0)
                           s->fock_eff[ij] = fock_c[i][j];
                        else s->fock_eff[ij] = fock_o[i][j];
                        }

               /* test form a*(2fc-fo 2fc-fo 2fc-fo) */
                     else if(iter < itmax-1 && !converged && fock_typ == 5) {
                        if (occi == occj)
                           s->fock_eff[ij] = 
                             dampd*(2.0*fock_c[i][j]-fock_o[i][j]);
                        else if(occi)
                           s->fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
                        else if(occj == 2.0)
                           s->fock_eff[ij] = fock_c[i][j];
                        else s->fock_eff[ij] = fock_o[i][j];
                        }

               /* form for converged wavefunction */
                     else {
                        if (occi == 2.0) s->fock_eff[ij]=fock_c[i][j];
                        else if(occj != 2.0)
                          s->fock_eff[ij]=fock_o[i][j];
                        else {
                          if(occi)
                            s->fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
                          else s->fock_eff[ij]=fock_c[i][j];
                          }
                        }

                     if(j==i) {
#if 1
                        if (occi == occ0 && occi)
                            s->fock_eff[ij] -= lshift;
                        else 
                           if (occi) s->fock_eff[ij] -= 0.5*lshift;
#else
                        if (occi == 1.0)
                            s->fock_eff[ij] += 0.5*lshift;
                        else if (occi == 0.0)
                            s->fock_eff[ij] += 0.5*(lshift+0.1);
#endif
                        }
                     ij++;
                     }
                  }
               }
            }
         }
               
      for (m=0; m < num_ir ; m++) {
         s = &scf_info[m];
         if (nn=s->num_so) {
            occ0 = s->occ_num[0];
            if(iopen) {
               rsp(nn,nn,ioff[nn],s->fock_eff,s->fock_evals,1,ctrans,tol);

/* if this is open shell, we want the eigenvalues to be sorted
   before they undergo level shifting */
	       
	       if(iter != 0 & inflg == 0 & occ_fix == 1){
		   for(i=0;i < s->num_so;i++)
		       s->hevals[i]=s->fock_evals[i];
		 }
	       
	       mxmb(s->cmat,1,nn,ctrans,1,nn,scr,1,nn,nn,nn,nn);
	       for (i=0; i < nn; i++) {
		   occi = s->occ_num[i];
#if 1
		   if (occi == occ0 && occi) s->fock_evals[i] += lshift;
		   else if (occi) s->fock_evals[i] += 0.5*lshift;
#else
		   if (occi == 1.0) s->fock_evals[i] -= 0.5*lshift;
		   else if (occi == 0.0) s->fock_evals[i] -= 0.5*(lshift+0.1);
#endif
		   
		   for (j=0; j < nn; j++) s->cmat[i][j] = scr[i][j];
		 }
	      }
	    else {

            /* transform fock_pac to mo basis */
               tri_to_sq(s->fock_pac,fock_ct,nn);
               mxmb(s->cmat,nn,1,fock_ct,1,nn,scr,1,nn,nn,nn,nn);
               mxmb(scr,1,nn,s->cmat,1,nn,fock_c,1,nn,nn,nn,nn);

             /*  diagonalize fock_c to get ctrans */
               sq_rsp(nn,nn,fock_c,s->fock_evals,1,ctrans,tol);

               if(print & 4) {
                 fprintf(outfile,"\n eigenvector for irrep %s\n",
                                                          s->irrep_label);
                 eivout(ctrans,s->fock_evals,nn,nn,outfile);
                 }

               mxmb(s->cmat,1,nn,ctrans,1,nn,scr,1,nn,nn,nn,nn);

               if(print & 4) {
                 fprintf(outfile,"\n eigenvector after irrep %s\n",
                     s->irrep_label);
                 print_mat(scr,nn,nn,outfile);
                 }

               for (i=0; i < nn; i++)
                  for (j=0; j < nn; j++)
                     s->cmat[i][j] = scr[i][j];
               }
            }
         }

      if(converged) {
         free_matrix(scr,nsfmax);
         free_matrix(fock_c,nsfmax);
         free_matrix(fock_ct,nsfmax);
         free_matrix(ctrans,nsfmax);
         if(iopen) {
            free_matrix(fock_o,nsfmax);
            free(c1);
            free(c2);
            }
         cleanup();
         }

      schmit(1);

      if(print & 4) {
        for(i=0; i < num_ir ; i++) {
         s = &scf_info[i];
         if (nn=s->num_so) {
           fprintf(outfile,"\northogonalized mos irrep %s\n",
               s->irrep_label);
           print_mat(s->cmat,nn,nn,outfile);
           }
         }
       }

/* form new density matrix */

      dmat();

/* and form new fock matrix */

      if(iter < itmax) {
         if(iopen) formg_open();
         else formg_closed();
         }
     }
}




