/* $Log: ecalc.c,v $
/* Revision 1.4  1995/07/21 17:37:05  psi
/* Made Jan 1st 1995 cscf the current accepted version of cscf.  Some
/* unidentified changes made after that date were causing problems.
/*
 * Revision 1.2  1994/06/03  04:39:47  seidl
 * change test to iter > 1
 *
 * Revision 1.1  1991/06/15  20:22:23  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: ecalc.c,v 1.4 1995/07/21 17:37:05 psi Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

static double twocut=1.0;
static double eelec;
double delta;

int ecalc(incr)
   double incr;

{
   int i,j,k,ij,nn;
   double edif;
   double plimit = pow(10.0,(double) -iconv);
   double neelec = 0.0;
   double ir_energy, dtemp;
   double cinext;
   struct symm *s;

   delta=0.0;
   for (k=0; k < num_ir ; k++) {
      s = &scf_info[k];
      if (nn=s->num_so) {
         ir_energy = 0.0;
         for (i=ij=0; i < nn ; i++) {
            for (j = 0 ; j <= i ; j++,ij++) {
               if(!iopen) {
                  ir_energy += 0.5*s->pmat[ij]*(s->hmat[ij]+s->fock_pac[ij]);
                  }
                else {
                  ir_energy += 0.5*s->pmat[ij]*(s->hmat[ij]+s->fock_pac[ij])
                                              - 0.5*s->pmato[ij]*s->gmato[ij];
                  }
               }
            }
         neelec += ir_energy;
         if (iter) {
            for (i = 0; i < ioff[nn] ; i++) {
               dtemp = s->dpmat[i];
               delta += dtemp*dtemp;
               }
            }
         }
      }

   delta = sqrt(delta)/mxcoef2;
   etot = repnuc + neelec;
   edif =  eelec - neelec;

   if (!iter) fprintf(outfile,"\n  iter       total energy        delta E         delta P          diiser\n");
   fprintf(outfile, "%5d %20.10f %15.6e %15.6e %15.6e\n", 
                                    ++iter, etot, edif, delta, diiser);
   fflush(outfile);
   diiser=0.0;

   if ( delta < plimit && iter > 1) {
      converged=1;
      if(!iopen || iopen && fock_typ >= 2) cleanup();
      }

   eelec = neelec;

   cinext = pow(10.0,-twocut);
   if (delta < cinext && delta && !converged) {
      twocut += incr;
      return(1);
      }
   else return(0);
   }
