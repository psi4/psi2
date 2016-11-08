/* $Log: init_scf.c,v $
/* Revision 1.12  1998/06/30 14:11:08  sbrown
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
 * Revision 1.11  1997/09/12  13:54:28  crawdad
 * Changing marco name from ULL to PSI_FPTR.
 *
 * Revision 1.10  1997/08/25  21:51:28  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 1.9  1996/11/26  05:20:55  sherrill
 * Add casts in front of init_int_array() calls to avoid compiler warnings
 * (really should include libciomr everywhere instead but that causes more
 * warnings because wwritw's need casts to char *).  Also fixed problem
 * where phase() was trying to free unallocated memory.
 *
 * Revision 1.8  1996/07/03  12:25:46  crawdad
 * Added new phase-checking routines to cscf.  These have been tested more
 * extensively than the last, and are far more intelligently written.
 *
 * Revision 1.7  1995/11/07  19:24:10  sherrill
 * Replace all initializations of integer arrays with calls to init_int_array().
 * There had been problems with using init_array() because the division by
 * the doubles-to-int factor was rounding off (integer division).
 *
 * Revision 1.6  1995/07/21  17:37:07  psi
 * Made Jan 1st 1995 cscf the current accepted version of cscf.  Some
 * unidentified changes made after that date were causing problems.
 *
 * Revision 1.2  1994/06/03  04:40:31  seidl
 * add test for SGI and add null to irrep_label
 *
 * Revision 1.1  1991/06/15  20:22:32  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: init_scf.c,v 1.12 1998/06/30 14:11:08 sbrown Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

init_scf()
{
   int i,jj;
   int nn,isadr;
   int nkind,junk;
   PSI_FPTR next;
   int degen[20],num_so[20];
   double *real_dum;
   char char_dum[80];

   i10 = (int *) init_int_array(200);
   real_dum = (double *) init_array(MAX_BASIS);

   wreadw(itap30,(char *) i10,sizeof(int)*200,(PSI_FPTR) sizeof(int)*100,&next);

   ioff[0] = 0;
   for (i = 1; i < 1024 ; i++) {
      ioff[i] = ioff[i-1] + i;
      }

/* read header information from integral tape */

   wreadw(itap34,(char *) (&nkind),sizeof(int)*1,(PSI_FPTR) 0,&next);
   wreadw(itap34,(char *) (&junk),sizeof(int)*1,next,&next);
   wreadw(itap34,(char *) (char_dum),sizeof(char)*80,next,&next);
   wreadw(itap34,(char *) (&repnuc),sizeof(double)*1,next,&next);
   wreadw(itap34,(char *) (&num_ir),sizeof(int)*1,next,&next);
   wreadw(itap34,(char *) (degen),sizeof(int)*num_ir,next,&next);
   wreadw(itap34,(char *) (char_dum),sizeof(int)*num_ir,next,&next);
   wreadw(itap34,(char *) (num_so),sizeof(int)*num_ir,next,&next);
   wreadw(itap34,(char *) (&junk),sizeof(int)*1,next,&next);
   wreadw(itap34,(char *) (real_dum),sizeof(int)*2*junk,next,&next);
   wreadw(itap34,(char *) (&junk),sizeof(int)*1,next,&next);
   wreadw(itap34,(char *) (real_dum),sizeof(int)*junk,next,&next);
   wreadw(itap34,(char *) (real_dum),sizeof(int)*junk,next,&next);

/* set integral file pointer to sector boundary */

   isadr = i2sec(next) + 1;
   rsetsa(itap34,isadr);
   pos34 = sec2i(--isadr);

   if (nkind != 1 && nkind != 2) {
      fprintf(outfile,"integral file screwed up, fix ints somebody!!!\n");
      exit(1);
      }

/* now initialize scf_info */
   
   n_so_typs=0;
   mxcoef=0;
   mxcoef2=0;
   nsfmax=0;
   nbasis=0;
   occ_fix=0;

   scf_info = (struct symm *) malloc(sizeof(struct symm)*num_ir);

   jj=0;
   for(i=0; i < num_ir ; i++) {
      scf_info[i].num_so = nn = num_so[i];
      scf_info[i].degeneracy = degen[i];
      scf_info[i].nclosed = 0;
      scf_info[i].nopen = 0;
      scf_info[i].nhalf = 0;
      scf_info[i].os_num = 0;
      scf_info[i].ideg = 0;

      strncpy(scf_info[i].irrep_label,&char_dum[jj],4);
      scf_info[i].irrep_label[4] = '\0';
      jj += 4;

      nbasis += nn;
      if (nn) {
         n_so_typs++;
         if (nn > nsfmax) nsfmax = nn;
         mxcoef += nn*nn;
         mxcoef2 += ioff[nn];

         scf_info[i].smat = (double *) init_array(ioff[nn]);
         scf_info[i].tmat = (double *) init_array(ioff[nn]);
         scf_info[i].hmat = (double *) init_array(ioff[nn]);
         scf_info[i].pmat = (double *) init_array(ioff[nn]);
         scf_info[i].pmato = (double *) NULL;
         scf_info[i].pmat2 = (double *) NULL;
         scf_info[i].pmato2 = (double *) NULL;
         scf_info[i].dpmat = (double *) init_array(ioff[nn]);
         scf_info[i].dpmato = (double *) NULL;
         scf_info[i].fock_pac = (double *) init_array(ioff[nn]);
         scf_info[i].fock_eff = (double *) NULL;
         scf_info[i].fock_open = (double *) NULL;
         scf_info[i].gmat = (double *) init_array(ioff[nn]);
         scf_info[i].gmato = (double *) NULL;
         scf_info[i].occ_num = (double *) init_array(nn);
         scf_info[i].fock_evals = (double *) init_array(nn);
         scf_info[i].cmat = (double **) init_matrix(nn,nn);
	 /* TDC(6/19/96) - Added array for saving original MO vector */
	 scf_info[i].cmat_orig = (double **) init_matrix(nn,nn);
         scf_info[i].sahalf = (double **) init_matrix(nn,nn);
         /* STB(4/1/98) - Added array to store the eigenvalues of the
	                  core hamiltonian for mo guessing*/
         scf_info[i].hevals = (double *) init_array(nn);
         }
     }
   /* read in number of atoms and nuclear charges and total number of MO*/
   natom = file30_rd_natom();
   zvals = file30_rd_zvals();
   nbfso = file30_rd_nmo();
   
/* Initialize arrays to hold energy and symmetry arrays */
   ener_tot = (double *) init_array(nbfso);
   symm_tot = (int *) init_int_array(nbfso);
   
} 

init_scf2()
{
   int i,j;
   int n,nn,m,mm;
   int junk;
   int opconst,outbuf,mxcoef3,ntri,mtri;
   struct symm *s;

   opconst = (iopen) ? 3 : 2;

   mtri = ioff[scf_info[0].num_so];
   mxcoef3 = opconst*(mtri*(mtri+1)/2);

   scf_info[0].ideg = 0;

   junk=j=0;
   for (i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      nn = s->num_so;
      if (i) {
         m = ioff[nn];
         mm = m*(m+1)/2;
         mxcoef3 += opconst*ioff[nn]*mtri;
         mxcoef3 += opconst*mm;
         mtri += m;
         if (nn <= 0) s->ideg = scf_info[i-1].ideg;
         else {
            do {
               n=scf_info[j].num_so;
               j++;
               } while(!n); 
            s->ideg = scf_info[i-1].ideg+n;
            }
         }
      if(s->nopen || s->nhalf) s->os_num = junk++;
      }
   ntri = nbasis*(nbasis+1)/2;

   readflg = 0;
   maxbuf = mxcoef3/opconst+5;
   outbuf = 884736;
   outbuf = MIN0(outbuf,mxcoef3);
   if(outbuf == 884736) {
      int pass = mxcoef3/outbuf+1;
      readflg = 1;
      maxbuf = (mxcoef3/pass)/opconst + 2;
      if(iopen) maxbuf /= 2;
#if defined(AIXV3)||defined(SGI)
      maxbuf= (iopen) ? 8192*3 : 8192*5;
      pass = mxcoef3/(maxbuf*opconst)+1;
#endif
      fprintf(outfile,
           "\n  using buffered io, %d buffers, each %d bytes in size\n",
           pass,maxbuf*opconst*sizeof(double));
      if(print) fprintf(outfile,"  mxcoef3 = %d maxbuf = %d\n",mxcoef3,maxbuf);
      if(print) fprintf(outfile,"  outbuf = %d\n",outbuf);
      }
   else fprintf(outfile,"\n  keeping integrals in %d bytes of core\n",
                 maxbuf*opconst*sizeof(double));

   fflush(outfile);
   }
