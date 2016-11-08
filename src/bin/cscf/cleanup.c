/* $Log: cleanup.c,v $
/* Revision 1.10  1998/06/30 14:11:00  sbrown
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
 * Revision 1.9  1997/09/12  13:54:25  crawdad
 * Changing marco name from ULL to PSI_FPTR.
 *
 * Revision 1.8  1997/08/25  21:51:21  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 1.7  1996/07/03  12:25:43  crawdad
 * Added new phase-checking routines to cscf.  These have been tested more
 * extensively than the last, and are far more intelligently written.
 *
 * Revision 1.6  1996/04/30  21:40:12  kellogg
 * make cscf exit abnormally if a wavefunction does not converge <cbk>
 *
 * Revision 1.5  1995/07/21  17:37:01  psi
 * Made Jan 1st 1995 cscf the current accepted version of cscf.  Some
 * unidentified changes made after that date were causing problems.
 *
 * Revision 1.2  1994/06/03  04:39:07  seidl
 *  change doublespace to single
 *
 * Revision 1.1  1991/06/15  20:22:12  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: cleanup.c,v 1.10 1998/06/30 14:11:00 sbrown Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <ip_libv1.h>

/* TDC(6/20/96) - Prototype for phase() */
int phase(void);

void cleanup()

{
   int i,j,k,ij,ijk,nn;
   PSI_FPTR junk,junk1;
   int mpoint,mconst,mcalcs,loccal;
   PSI_FPTR locvec;
   int nx,ntri;
   int newvec;
   int nat,iend,ierr,ci_calc,irot,nbfao;
   int wrt3149;
   int itap49=49;
   int itap31=31;
   int numso;
   int n_there[20],nc[10],no[10];
   int errcod;
   char *ci_type="SCF";
   char *der_type="FIRST";
   double occj,occk;
   double ekin,epot,ovlp,virial,num_elec;
   double *scr_arr, *lagrangian;
   double **scr1, **scr2;
   double *temp;
   struct symm *s;

   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":SCF");

   ci_calc=irot=0;
   errcod = ip_string("WFN",&ci_type,0);
   errcod = ip_string("DERTYPE",&der_type,0);

   if(strcmp(ci_type,"SCF")) ci_calc=1;
   if(ci_calc && iopen) irot=1;
   if(strcmp(der_type,"FIRST") && strcmp(der_type,"NONE")) irot=0;

   errcod = ip_boolean("ROTATE",&irot,0);

 /* TDC(6/19/96) - If we're not rotating, check the phases on the MOs,
   and correct them, if possible. */

   if(!irot) {
       if(phase_check) phase_check = phase();
     }
   
/* first print mo's, then rotate if this is a ci */

   print_mos();  

   if(irot) rotate_vector();
   else 
     if(iopen)
       fprintf(outfile,"\n%8cWFN is %s so no rotation\n",' ',ci_type);

   if(ci_calc || !strcmp(der_type,"SECOND")) wrt3149=1;
   else wrt3149=0;

   if(wrt3149) {
     rfile(itap49);
     rfile(itap31);
     }

/* TDC(6/20/96) - If we've rotated the orbitals, check the phases on
   the MOs, and correct them, if possible. */

   if(irot) {
       if(phase_check) phase_check = phase();
     }

/* update ncalcs in file30 */

   wreadw(itap30,(char *) i10,sizeof(int)*200,(PSI_FPTR) sizeof(int)*100,&junk);

   if(!i10[44]) i10[44]++;

   iend = i10[0];
   mpoint = i10[1];
   mconst = i10[2];
   mcalcs = i10[3];
   nat = i10[18];
   nx = i10[23];
   nbfao = i10[21];
   i10[40] = n_so_typs;
   i10[41] = mxcoef;
   i10[42] = ioff[n_open];
   if(twocon) i10[42] = -i10[42];
   /* TDC(6/19/96) - Set the phase_check flag here */
   i10[50] = phase_check;

   wwritw(itap30,(char *) i10,sizeof(int)*200,(PSI_FPTR) sizeof(int)*100,&junk);

/* get pointers to calculations */

   junk = (PSI_FPTR) sizeof(int)*(100+mconst+mpoint);
   wreadw(itap30,(char *) i10,sizeof(int)*mcalcs,junk,&junk);
   loccal = (PSI_FPTR) sizeof(int)*(i10[0]-1);

/* pointer to scf vector is first in calculation pointers */

   junk = loccal+sizeof(int)*60;
   wreadw(itap30,(char *) i10,sizeof(int)*20,junk,&junk1);

   newvec = 0;
   locvec = (PSI_FPTR) sizeof(int)*(i10[0]-1);
   if (i10[0] <= 0) {
      locvec = (PSI_FPTR) sizeof(int)*(iend-1);
      i10[0] = iend;
      wwritw(itap30,(char *) i10,sizeof(int)*20,junk,&junk1);
      newvec = 1;
      }

   junk = loccal + ((PSI_FPTR) sizeof(int)*(60+20+6*nat+2));

/* write new energy to file30 */

   wwritw(itap30,(char *) &etot,sizeof(double)*1,junk,&junk);

/* write new vector and eigenvalues to file30 and file49 */

   scr_arr = (double *) init_array(mxcoef);

   for(k=ijk=0; k < num_ir ; k++) {
      s= &scf_info[k];
      if(nn=s->num_so) {
         for(j=0; j < nn ; j++)
            for(i=0; i < nn ; i++,ijk++)
               scr_arr[ijk]=s->cmat[i][j];
         }
      }

   wwritw(itap30,(char *) scr_arr,sizeof(double)*mxcoef,locvec,&locvec);
   if(wrt3149) swrit(itap49,(char *) scr_arr,sizeof(double)*mxcoef);

   for (i=k=ijk=0; i < num_ir ; i++) {
      s= &scf_info[i];
      if (nn=s->num_so) {
         for (j=0; j < nn ; j++,k++) scr_arr[k] = s->fock_evals[j];
         nc[ijk]=s->nclosed;
         no[ijk]=s->nopen;
         n_there[ijk]=nn;
         ijk++;
         }
      }
   wwritw(itap30,(char *) scr_arr,sizeof(double)*nbasis,locvec,&locvec);
   if(wrt3149) swrit(itap49,(char *) scr_arr,sizeof(double)*nbasis);

   for (i=0 ; i < num_ir ; i++) 
      if (scf_info[i].num_so)
         wwritw(itap30,(char *) scf_info[i].irrep_label,sizeof(char)*4,locvec,&locvec);

   wwritw(itap30,(char *) n_there,sizeof(int)*n_so_typs,locvec,&locvec);
   wwritw(itap30,(char *) nc,sizeof(int)*n_so_typs,locvec,&locvec);
   if(iopen) {
      wwritw(itap30,(char *) no,sizeof(int)*n_so_typs,locvec,&locvec);
      if(twocon) {
         double c1i = scf_info[opblk1].occ_num[opshl1];
         double c1ii = scf_info[opblk2].occ_num[opshl2];
         alpha[0] = beta[0] = 1.0 - 2.0/c1i;
         alpha[1] = 1.0;
         beta[1] = 1.0 + 1.0/(save_ci1*save_ci2);
         alpha[2] = beta[2] = 1.0 - 2.0/c1ii;
         }
      if(!twocon) for (i=0; i < ioff[n_open] ; i++) beta[i] = -beta[i];
      if(ci_calc && !twocon)
         for (i=0; i < ioff[n_open] ; i++) beta[i] = -beta[i];
      wwritw(itap30,(char *) alpha,sizeof(double)*ioff[n_open],locvec,&locvec);
      wwritw(itap30,(char *) beta,sizeof(double)*ioff[n_open],locvec,&locvec);
      }

/* calculate mo lagrangian and write to file30 */
/* also write mo fock matrices to file49       */
      
   lagrangian = (double *) init_array(nx);
   if (iopen) {
      scr1 = (double **) init_matrix(nsfmax,nsfmax);
      scr2 = (double **) init_matrix(nsfmax,nsfmax);

      for (i=0; i < num_ir ; i++) {
         s = &scf_info[i];
         if(nn=s->num_so) {

            tri_to_sq(s->fock_pac,scr2,nn);
            mxmb(s->cmat,nn,1,scr2,1,nn,scr1,1,nn,nn,nn,nn);
            mxmb(scr1,1,nn,s->cmat,1,nn,scr2,1,nn,nn,nn,nn);
            sq_to_tri(scr2,s->gmat,nn);
            if(wrt3149) swrit(itap49,(char *) s->gmat,sizeof(double)*ioff[nn]);

            tri_to_sq(s->fock_open,scr2,nn);
            mxmb(s->cmat,nn,1,scr2,1,nn,scr1,1,nn,nn,nn,nn);
            mxmb(scr1,1,nn,s->cmat,1,nn,scr2,1,nn,nn,nn,nn);
            sq_to_tri(scr2,s->gmato,nn);
            if(wrt3149) swrit(itap49,(char *) s->gmato,sizeof(double)*ioff[nn]);

            bzero(s->fock_pac,sizeof(double)*ioff[nn]);
            for (j=ij=0; j < nn ; j++) {
               for (k=0; k <= j ; k++,ij++) {
                  occj = s->occ_num[j];
                  occk = s->occ_num[k];
                  if(!twocon) {
                     if (occj == 2.0 && occk == 2.0 || 
                         occj == 2.0 && occk == 1.0 ||
                         occj == 1.0 && occk == 2.0) 
                         s->fock_pac[ij] = 2.0*s->gmat[ij];
                     else if (occj == 1.0 && occk == 1.0)
                         s->fock_pac[ij] = s->gmato[ij];
                     else
                         s->fock_pac[ij] = 0.0;
                     }
                  else {
                     if (occj == 2.0 && occk || 
                         occk == 2.0 && occj)
                         s->fock_pac[ij] = 2.0*s->gmat[ij];
                     else if (occj && occk)
                         s->fock_pac[ij] = occj*s->gmato[ij];
                     else
                         s->fock_pac[ij] = 0.0;
                     }
                  }
               }
            }
         }
      }
   else {
      for (i=0; i < num_ir ; i++) {
         s = &scf_info[i];
         if(nn=s->num_so) {
            bzero(s->fock_pac,sizeof(double)*ioff[nn]);
            for (j=0; j < nn ; j++)
                if (s->occ_num[j] == 2.0)
                    s->fock_pac[ioff[j]+j] = 2.0*s->fock_evals[j];
            }
         }
      }

   for(k=ijk=0; k < num_ir ; k++) {
      if(nn=scf_info[k].num_so) {
         for(i=0; i < nn ; i++)
            for(j=0; j <= i ; j++)
              lagrangian[ioff[i+ijk]+j+ijk] = scf_info[k].fock_pac[ioff[i]+j];
         ijk += nn;
         }
      }
   wwritw(itap30,(char *) lagrangian,sizeof(double)*nx,locvec,&locvec);

   if(newvec) {
      iend = (int) locvec/sizeof(int)+1;
      wwritw(itap30,(char *) &iend,sizeof(int)*1,(PSI_FPTR) sizeof(int)*100,&locvec);
      }


/* write not so pretty mo's and energies to output */

   if(ci_calc && iopen && irot) {
      fprintf(outfile,
         "\n ci_typ is %s so mo vector will be rotated\n",ci_type);
      print_mos();
      }
   num_elec = ekin = ovlp = 0.0;
   for (i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      if (nn=s->num_so) {
         for (j=0; j < ioff[nn] ; j++) {
            ekin += s->pmat[j]*s->tmat[j];
            ovlp += s->pmat[j]*s->smat[j];
            }
         for (j=0; j < nn ; j++) num_elec += s->occ_num[j];
         }
      }

   ovlp /= num_elec;
   epot = etot-ekin;
   virial = epot/etot;

/* Print MOs in a really f***ing special format */

   if(mo_out == 1){
       print_mos_new();
     }
   
   fprintf(outfile,"\n%8ctotal energy       = %20.12f\n",' ',etot);
   fprintf(outfile,"%8ckinetic energy     = %20.12f\n",' ',ekin);
   fprintf(outfile,"%8cpotential energy   = %20.12f\n",' ',epot);
   fprintf(outfile,"%8cvirial theorem     = %20.12f\n",' ',virial);
   fprintf(outfile,"%8cwavefunction norm  = %20.12f\n",' ',ovlp);

/* stuff for master */

   if(wrt3149) {
      ntri = ioff[nbasis];
      temp = (double *) init_array(ntri);

      for(k=ijk=0; k < num_ir ; k++) {
         if(nn=scf_info[k].num_so) {
            for(i=0; i < nn ; i++)
               for(j=0; j <= i ; j++)
                  temp[ioff[i+ijk]+j+ijk] = scf_info[k].pmat[ioff[i]+j];
            ijk += nn;
            }
         }

      swrit(itap31,(char *) temp,sizeof(double)*ntri);
      if (iopen) {
         for(k=ijk=0; k < num_ir ; k++) {
            if(nn=scf_info[k].num_so) {
               for(i=0; i < nn ; i++)
                  for(j=0; j <= i ; j++)
                     temp[ioff[i+ijk]+j+ijk] = scf_info[k].pmato[ioff[i]+j];
               ijk += nn;
               }
            }
         swrit(itap31,(char *) temp,sizeof(double)*ntri);
         }

      rclose(itap31,3);

      if(ci_calc) make_47(lagrangian,temp,ci_type,nat,ntri,nbfao);

      free(temp);
      rclose(itap49,3);
      }


   file30_close();
   rclose(itap92,4);
   rclose(itap93,4);
   
   if(!converged)
      fprintf(outfile,"\n%8cCalculation has not converged!\n",' ');
   tstop(outfile);
   ip_done();

   if(!converged) exit(3);
   exit(0);
   }

print_mos()
{
   int i,nn;
   struct symm *s;
   

   for (i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      if (nn=s->num_so) {
         fprintf(outfile,"\n molecular orbitals for irrep %s\n",s->irrep_label);
         eigout(s->cmat,s->fock_evals,s->occ_num,nn,nn,outfile);
         }
      }
}
print_mos_new() {
  int i;
  int ncl_tot=0;
  int nop_tot=0;
  
  for(i=0;i < num_ir;i++){
      ncl_tot+= scf_info[i].nclosed;
      nop_tot+= scf_info[i].nopen;
    }
  sortev();

  for(i = 0; i < 50; i++)
      fprintf(outfile,"*");

  fprintf(outfile,"\n\nMolecular Orbitals\n\n");

  for(i = 0; i < 50;i++)
      fprintf(outfile,"*");

  fprintf(outfile,"\nSymmetry\tOccupancy\tEnergy\n");

  for(i = 0; i < ncl_tot; i++)
      fprintf(outfile,"\n%s\t\tD\t\t%10.6f",scf_info[symm_tot[i]].irrep_label,ener_tot[i]);
  for(i = ncl_tot; i < ncl_tot + nop_tot; i++)
      fprintf(outfile,"\n%s\t\tS\t\t%10.6f",scf_info[symm_tot[i]].irrep_label,ener_tot[i]);
  fprintf(outfile,"\n");
  for(i = 0; i < 50;i++)
      fprintf(outfile,"-");
  for(i = ncl_tot + nop_tot; i < nbfso; i++)
      fprintf(outfile,"\n%s\t\tV\t\t%10.6f",scf_info[symm_tot[i]].irrep_label,ener_tot[i]);

  fprintf(outfile,"\n");

  for(i = 0; i < 50;i++)
      fprintf(outfile,"*");

  fprintf(outfile,"\n"); 
}

