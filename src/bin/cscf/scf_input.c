/* $Log: scf_input.c,v $
/* Revision 1.14  1998/06/30 14:11:10  sbrown
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
 * Revision 1.13  1997/09/12  13:54:30  crawdad
 * Changing marco name from ULL to PSI_FPTR.
 *
 * Revision 1.12  1997/08/25  21:51:35  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 1.11  1996/07/03  12:25:48  crawdad
 * Added new phase-checking routines to cscf.  These have been tested more
 * extensively than the last, and are far more intelligently written.
 *
 * Revision 1.10  1995/07/26  20:33:03  sherrill
 * Added support for getting second root of TCSCF and for checking the
 * orbital rotation to make sure orbitals are not swapped.
 *
 * Revision 1.9  1995/07/21  17:37:13  psi
 * Made Jan 1st 1995 cscf the current accepted version of cscf.  Some
 * unidentified changes made after that date were causing problems.
 *
 * Revision 1.2  1994/06/03  04:41:31  seidl
 * remove space and add check for number of electrons per irrep
 *
 * Revision 1.1  1991/06/15  20:22:39  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: scf_input.c,v 1.14 1998/06/30 14:11:10 sbrown Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <ip_libv1.h>

void scf_input(ipvalue)
   ip_value_t *ipvalue;
{
   int i,j,k,ijk;
   double elast;     
   double **scr_mat, *scr_arr;
   char *alabel,*bool="YES",*optyp,*wfn,*dertype;
   char cjunk[80];
   char *prog_name="CSCF: An SCF program written in C";
   int norder,*iorder,reordr;
   int nc,no,nh,nn;
   int mpoint,mconst,mcalcs,ncalcs;
   PSI_FPTR junk,locvec,loccal;
   int optri,ierr,nat;
   int io_locate();
   int errcod;
   int size;
   int phase_chk;
   struct symm *s;

   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":SCF");
   if(ipvalue) ip_print_value(stdout,ipvalue);

   fprintf(outfile,"\n%16c%s\n\n",' ',prog_name);

   errcod = ip_string("LABEL",&alabel,0);
   if(errcod == IPE_OK) fprintf(outfile,"  label       = %s\n",alabel);

   inflg = 0;
   errcod = ip_string("RESTART",&bool,0);
   if(errcod == IPE_OK) {
      if(!strcmp(bool,"NO") || !strcmp(bool,"FALSE") || !strcmp(bool,"0"))
        inflg=2;
      }

   reordr = 0;
   norder = 0;
   errcod = ip_boolean("REORDER",&reordr,0);
   if(reordr) {
     norder = 1;

     errcod = ip_count("MOORDER",&size,0);
     errchk(errcod,"MOORDER");
     if(errcod != IPE_OK) {
       fprintf(outfile,"\ncannot find MOORDER. calculation continuing\n");
       norder=0;
       reordr=0;
       }
     else {
       if(size != nbasis) {
         fprintf(outfile,"\n you have not given enough mos to MOORDER\n");
         exit(size);
         }
       iorder = (int *) malloc(sizeof(int)*size);
       for(i=0; i < size ; i++) {
         errcod = ip_data("MOORDER","%d",&iorder[i],1,i);
         errchk(errcod,"MOORDER");
         }
       }
     }

   itmax = 40;
   errcod = ip_data("MAXITER","%d",&itmax,0);

   it_diis = 0;
   errcod = ip_data("DIISSTART","%d",&it_diis,0);

   print = 0;
   errcod = ip_data("IPRINT","%d",&print,0);

   fock_typ = 0;
   errcod = ip_data("FOCK_TYPE","%d",&fock_typ,0);

   iopen=0;
   special=twocon=hsos=singlet=0;
   errcod = ip_string("OPENTYPE",&optyp,0);
   errchk(errcod,"OPENTYPE");
   if(errcod == IPE_OK) {
      iopen=1;
      if(!strcmp(optyp,"HIGHSPIN")) hsos=1;
      else if(!strcmp(optyp,"SINGLET")) singlet=1;
      else if(!strcmp(optyp,"TWOCON")) twocon=1;
      else if(!strcmp(optyp,"SPECIAL")) special=1;
      else if(!strcmp(optyp,"NONE")) iopen=0;
      else {
        fprintf(outfile,"unrecognized OPENTYPE: %s\n",optyp);
        exit(1);
        }
      }
/* Do you want the new mo output */
   
   
   mo_out = 1;
   errcod = ip_boolean("MO_OUT",&mo_out,0);
      
   second_root = 0;
   if (twocon) {
       errcod = ip_boolean("SECOND_ROOT",&second_root,0);
     }
   
   icheck_rot = 1;
   errcod = ip_boolean("CHECK_ROT",&icheck_rot,0);

   ndiis = (iopen) ? 4 : 6;
   if(twocon) ndiis = 3;
   errcod = ip_data("NDIIS","%d",&ndiis,0);

   if(ipvalue) ip_print_tree(stdout,NULL);

   iconv = 7;
   if(ipvalue) ip_print_value(stdout,ipvalue);
   errcod = ip_string("WFN",&wfn,0);
   if(ipvalue) ip_print_value(stdout,ipvalue);
   errcod = ip_string("DERTYPE",&dertype,0);
   if(errcod == IPE_KEY_NOT_FOUND) {
     dertype = (char *) malloc(sizeof(char)*5);
     strcpy(dertype,"NONE");
     }
   if(strcmp(wfn,"SCF")) iconv = 10;
   if(!strcmp(dertype,"SECOND")) iconv = 12;
   errcod = ip_data("CONVERGENCE","%d",&iconv,0);

   if(ipvalue) ip_print_value(stdout,ipvalue);
   fprintf(outfile,"  wfn         = %s\n",wfn);
   fprintf(outfile,"  dertype     = %s\n",dertype);
   fprintf(outfile,"  opentype    = %s\n",optyp);
   fprintf(outfile,"  convergence = %d\n",iconv);
   fprintf(outfile,"  maxiter     = %d\n",itmax);
   fprintf(outfile,"  restart     = %s\n",bool);
   if(print) fprintf(outfile,"  iprint      = %d\n",print);
   if (second_root)
     fprintf(outfile,"  second_root = TRUE\n");


   diisflg = 0;
   errcod = ip_string("DIIS",&bool,0);
   if(errcod == IPE_OK) {
      if(!strcmp(bool,"NO") || !strcmp(bool,"FALSE") || !strcmp(bool,"0"))
        diisflg=1;
      else diisflg=0;
      }

   fprintf (outfile,"\n  nuclear repulsion energy %22.13f\n",repnuc);
   fflush(outfile);

/* FILE30 stuff */
   
   mpoint = i10[1];
   mconst = i10[2];
   mcalcs = i10[3];
      nat = i10[18];
   ncalcs = i10[44];

/* if inflg is 0 and this isn't the first calc, then get the old vector */
/* from file30.  if inflg is 2, just use core hamiltonian guess */
/* if inflg is 1, get old vector no matter what ncalcs is */

   if ((!inflg && ncalcs) || inflg == 1) {

      if (!inflg) inflg = 1;
      mxcoef = i10[41];
      optri = abs(i10[42]);

      scr_arr = (double *) init_array(mxcoef);

/* grab pointers to calculations */

      junk = (PSI_FPTR) sizeof(int)*(100+mconst+mpoint);
      wreadw(itap30,(char *) i10,sizeof(int)*mcalcs,junk,&junk);
      loccal = (PSI_FPTR) sizeof(int)*(i10[0]-1);

/* get label from old vector */

      for (i=0; i < 3; i++) {
         wreadw(itap30,(char *) cjunk,sizeof(char)*80,loccal,&loccal);
         }

      fprintf(outfile,"\n  using old vector from file30 as initial guess\n");

/* get pointers to vectors, etc */

      wreadw(itap30,(char *) i10,sizeof(int)*20,loccal,&loccal);

      locvec = (PSI_FPTR) sizeof(int)*(i10[0]-1);
      loccal += (PSI_FPTR) sizeof(int)*6*nat;

/* get old energy from file30 */

      junk = loccal+ sizeof(double)*1;
      wreadw(itap30,(char *) &elast,sizeof(double)*1,junk,&junk);
      fprintf(outfile,"\n  energy from old vector: %14.8f\n",elast);

/* get the old vector */

      wreadw(itap30,(char *) scr_arr,sizeof(double)*mxcoef,locvec,&locvec);

/* and convert it to symmetry blocked form */
         
      ijk=0;
      for(k=0; k < num_ir ; k++) {
         s = &scf_info[k];
         if(nn=s->num_so) {
            for(j=0; j < nn ; j++)
               for(i=0; i < nn ; i++,ijk++) s->cmat[i][j] = scr_arr[ijk];
                
            }
         }

      free(scr_arr);

/* TDC(6/19/96) - Make a copy of the vector for later MO phase
   checking and temporarily set the phase_check flag to true */

      for(k=0; k < num_ir; k++) {
	  s = &scf_info[k];
	  if(nn=s->num_so) {
	      for(j=0; j < nn; j++)
		  for(i=0; i < nn; i++) s->cmat_orig[i][j] = s->cmat[i][j];
	    }
	}

      phase_check = 1;

/* reorder vector if norder = 1 */

      if (norder) {
         int loff = 0;
         int jnew;

	 /* TDC(6/19/96) - If the vector is re-ordered, don't allow
	    phase_checking */
	 phase_check = 0;

         fprintf(outfile,"\n  mo's will be reordered\n");
         for (i=0; i < num_ir; i++) {
            s = &scf_info[i];
            if (nn=s->num_so) {
               scr_mat = (double **) init_matrix(nn,nn);
               for (j=0; j < nn; j++) {
                  jnew = iorder[j+loff]-1;
                  for (k=0; k < nn ; k++) {
                     scr_mat[k][j]=s->cmat[k][jnew];
                     }
                  }
               for (j=0; j < nn ; j++)
                   for (k=0; k < nn ; k++) s->cmat[j][k] = scr_mat[j][k];

               fprintf(outfile,"\n reordered mo's for irrep %s\n",
                                                              s->irrep_label);
               print_mat(s->cmat,nn,nn,outfile);
               loff += nn;
               free_matrix(scr_mat,nn);
               }
            }
         free(iorder);
         }
      }
   else {
      inflg = -1;
      fprintf(outfile,"  first run, so defaulting to core-hamiltonian guess\n");
      /* TDC(6/19/96) - If not starting from old vector, don't allow
	 phase checking */
      phase_check = 0;
      }

/* TDC(6/20/96) - Check to see if the user will let us do phase
   correction.  The default has already been set above. */
   phase_chk = 1;
   errcod = ip_boolean("PHASE",&phase_chk,0);
   if(phase_check && phase_chk) phase_check = 1;

/* check to see if open shell, then read in alpha and beta coeffs */
/*
   if (iopen) {
      int mm1=1,mm2=1;
      fprintf(outfile,"\n  open-shell energy coeffs\n");
      fprintf(outfile,"  open shell pair    alpha         beta\n");
      optri = (n_open*(n_open+1))/2;

      alpha = (double *) init_array(ioff[optri]);
      beta = (double *) init_array(ioff[optri]);

      if (twocon) {
         if(n_open == 2) {
           alpha[0] = 0.0;
           alpha[1] = 0.0;
           alpha[2] = 0.0;
           beta[0] = 0.0;
           beta[1] = 1.0;
           beta[2] = 0.0;
           }
         else {
           fprintf(outfile,"this program cannot handle same symmetry\n");
           fprintf(outfile," tcscf. try SCFX\n");
           exit(1);
           }
         }
      else if(singlet) {
         if(n_open == 2) {
           alpha[0] = 0.0;
           alpha[1] = 0.0;
           alpha[2] = 0.0;
           beta[0] = 1.0;
           beta[1] = -3.0;
           beta[2] = 1.0;
           }
         else {
           fprintf(outfile,"this program cannot handle same symmetry\n");
           fprintf(outfile," singlets. try SCFX\n");
           exit(1);
           }
         }
      else if(hsos) {
         for(i=0; i < optri ; i++) {
           alpha[i]=0.0;
           beta[i]=1.0;
           }
         }
      else {
         for (i=0; i < optri ; i++) {
            errcod = ip_data("ALPHA","%lf",&alpha[i],1,i);
            errcod = ip_data("BETA","%lf",&beta[i],1,i);
            beta[i] = -beta[i];
            }
         }
      for (i=0; i < optri; i++) {
         fprintf(outfile,"        %d  %d       %f     %f\n",mm1,mm2,
                 alpha[i],-beta[i]);
         mm2++;
         if (mm2 > mm1) {
            mm1++;
            mm2 = 1;
            }
         }
      }
*/
/* read in damping factor and level shift */

   dampsv= (iopen) ? 0.02 : 0.0;
   if(twocon) dampsv = 0.01;
   errcod = ip_data("DIISDAMP","%lf",&dampsv,0);

   lshift=1.0;
   errcod = ip_data("LEVELSHIFT","%lf",&lshift,0);
   if(!iopen && fabs(lshift) > 0.0) lshift = 0.1;

   dampd=1.0;
   errcod = ip_data("DAMPD","%lf",&dampd,0);

   dampo=1.0;
   errcod = ip_data("DAMPO","%lf",&dampo,0);

   fprintf(outfile,"\n  level shift                      = %f\n",lshift);
   if(!diisflg) {
      fprintf(outfile,"  diis scale factor                = %f\n",dampsv+1.0);
      fprintf(outfile,"  iterations before extrapolation  = %d\n",it_diis);
      fprintf(outfile,"  %d error matrices will be kept\n",ndiis);
      }
   else fprintf(outfile,"\n  diis turned off\n");

   switch (fock_typ) {
      case 0:
         break;
      case 1:
         fprintf(outfile,"\n  a fock matrix for high spin will be used\n");
         fprintf(outfile,"  this form may not work well with diis\n");
         break;
      default:
         fprintf(outfile,"\n  an experimental fock matrix will be used\n");
         fprintf(outfile,"  the management will not be held responsible for the results\n");
      }

   fflush(outfile);
   }
