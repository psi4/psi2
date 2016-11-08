/* $Log: formgc.c,v $
/* Revision 1.4  1997/08/25 21:51:25  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.3  1996/11/26  05:20:54  sherrill
 * Add casts in front of init_int_array() calls to avoid compiler warnings
 * (really should include libciomr everywhere instead but that causes more
 * warnings because wwritw's need casts to char *).  Also fixed problem
 * where phase() was trying to free unallocated memory.
 *
 * Revision 1.2  1995/11/07 19:24:08  sherrill
 * Replace all initializations of integer arrays with calls to init_int_array().
 * There had been problems with using init_array() because the division by
 * the doubles-to-int factor was rounding off (integer division).
 *
 * Revision 1.1  1991/06/15  20:22:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: formgc.c,v 1.4 1997/08/25 21:51:25 crawdad Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

static double *gtmp,*ptmp;
extern struct c_pkints {
         int ij;
         int kl;
         double pval;
         } *c_outbuf;
extern int last;
static int where=0;
static int *int_nums;

void formg_closed()

{
   register int i,j,k,joff,nn;
   register int ij,kl;
   int ilast,num;
   int tmpsiz;
   double tmpval;
   struct c_pkints *ctmp;

   tmpsiz = ioff[nbasis];

   if(gtmp == NULL) {
      gtmp = (double *) init_array(tmpsiz);
      ptmp = (double *) init_array(tmpsiz);
      }
   else bzero(gtmp,sizeof(double)*tmpsiz);

   for(k=joff=0; k < num_ir ; k++) {
      if(nn=scf_info[k].num_so) {
         for(i=0; i < nn ; i++)
            for(j=0; j <= i ; j++)
               ptmp[ioff[i+joff]+j+joff] = scf_info[k].dpmat[ioff[i]+j];
         joff += nn;
         }
      }

   if(!where) {
      /* int_nums = (int *) init_array(num_bufs+1); */
      int_nums = (int *) init_int_array(num_bufs+1);
      for(i=1; i < num_bufs ; i++) int_nums[i]=maxbuf;
      int_nums[num_bufs]=last;
      where=num_bufs;
      }

   num=int_nums[where];
   for (j=0; j < num_bufs; j++) {
      ctmp = c_outbuf;

      for (i=num; i ; i--,ctmp++) {
         ij = (*ctmp).ij;
         kl = (*ctmp).kl;
         tmpval = (*ctmp).pval;

         gtmp[ij] += ptmp[kl]*tmpval;
         gtmp[kl] += ptmp[ij]*tmpval;
         }

      if (readflg && j < num_bufs-1) {
         if(where==num_bufs) {
            where=0;
            srew(itap92);
            }
         where++;
         num=int_nums[where];
         sread(itap92,(char *) c_outbuf,sizeof(struct c_pkints)*num);
         }
      }
   for(k=joff=0; k < num_ir ; k++) {
      if(nn=scf_info[k].num_so) {
         for(i=0; i < nn ; i++)
            for(j=0; j <= i ; j++)
               scf_info[k].gmat[ioff[i]+j] += gtmp[ioff[i+joff]+j+joff];
         joff += nn;
         }
      }
   if(print & 4) {
      fprintf(outfile,"\n gmat\n");
      print_array(gtmp,nbasis,outfile);
      }
   }
