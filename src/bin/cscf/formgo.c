/* $Log: formgo.c,v $
/* Revision 1.4  1997/08/25 21:51:26  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.3  1996/11/26  05:20:55  sherrill
 * Add casts in front of init_int_array() calls to avoid compiler warnings
 * (really should include libciomr everywhere instead but that causes more
 * warnings because wwritw's need casts to char *).  Also fixed problem
 * where phase() was trying to free unallocated memory.
 *
 * Revision 1.2  1995/11/07 19:24:09  sherrill
 * Replace all initializations of integer arrays with calls to init_int_array().
 * There had been problems with using init_array() because the division by
 * the doubles-to-int factor was rounding off (integer division).
 *
 * Revision 1.1  1991/06/15  20:22:29  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: formgo.c,v 1.4 1997/08/25 21:51:26 crawdad Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

static double *gtmp,*gtmpo,*ptmp,*ptmpo;
extern int num_bufs_c,num_bufs_o,readflgc,readflgo;
extern struct c_pkints {
         int ij;
         int kl;
         double pval;
         } *c_outbuf;
extern struct o_pkints {
         int ij;
         int kl;
         double pval;
         double qval;
         } *o_outbuf;

extern int lasto,lastc;
int wherec=0;
int whereo=0;
int *int_nums_c;
int *int_nums_o;

void formg_open()
   {
   register int i,j,k,joff,nn;
   register int ij,kl;
   int ilast,num;
   int tmpsiz;
   double dotest,tmpval;
   struct o_pkints *o_temp;
   struct c_pkints *c_temp;

   tmpsiz=ioff[nbasis];

   if(gtmp == NULL) {
      gtmp = (double *) init_array(tmpsiz);
      gtmpo = (double *) init_array(tmpsiz);
      ptmp = (double *) init_array(tmpsiz);
      ptmpo = (double *) init_array(tmpsiz);
      }
   else {
      bzero(gtmp,sizeof(double)*tmpsiz);
      bzero(gtmpo,sizeof(double)*tmpsiz);
      }
 
   for(k=joff=0; k < num_ir ; k++) {
      if(nn=scf_info[k].num_so) {
         for(i=0; i < nn ; i++) {
            if(twocon) {
               for(j=0; j <= i ; j++) {
                  ptmp[ioff[i+joff]+j+joff] = scf_info[k].pmat2[ioff[i]+j];
                  ptmpo[ioff[i+joff]+j+joff] = scf_info[k].pmato2[ioff[i]+j];
                  }
               }
            else {
               for(j=0; j <= i ; j++) {
                  ptmp[ioff[i+joff]+j+joff] = scf_info[k].dpmat[ioff[i]+j];
                  ptmpo[ioff[i+joff]+j+joff] = scf_info[k].dpmato[ioff[i]+j];
                  }
               }
            }
         joff += nn;
         }
      }

   if(!wherec) {
      /* int_nums_o = (int *) init_array(num_bufs_o+1); */
      int_nums_o = (int *) init_int_array(num_bufs_o+1);
      /* int_nums_c = (int *) init_array(num_bufs_c+1); */
      int_nums_c = (int *) init_int_array(num_bufs_c+1);
      for(i=1; i < num_bufs_o ; i++) int_nums_o[i]=maxbuf;
      for(i=1; i < num_bufs_c ; i++) int_nums_c[i]=maxbuf;
      int_nums_o[num_bufs_o]=lasto;
      int_nums_c[num_bufs_c]=lastc;
      whereo=num_bufs_o;
      wherec=num_bufs_c;
      }

   num=int_nums_o[whereo];
   for (j=0; j < num_bufs_o ; j++) {
      o_temp = o_outbuf;

      for (i=num; i ; i--,o_temp++) {
         ij = (*o_temp).ij;
         kl = (*o_temp).kl;
         tmpval = (*o_temp).pval;
         dotest = (*o_temp).qval;

         gtmp[ij] += ptmp[kl]*tmpval;
         gtmp[kl] += ptmp[ij]*tmpval;
         gtmpo[ij] += ptmpo[kl]*dotest;
         gtmpo[kl] += ptmpo[ij]*dotest;
         }
   
      if (readflgo && j < num_bufs_o-1) {
         if(whereo==num_bufs_o) {
            srew(itap92);
            whereo=0;
            }
         whereo++;
         num=int_nums_o[whereo];
         sread(itap92,(char *) o_outbuf,sizeof(struct o_pkints)*num);
         }
      }

   num=int_nums_c[wherec];
   for (j=0; j < num_bufs_c ; j++) {
      c_temp = c_outbuf;

      for (i=num; i ; i--,c_temp++) {
         ij = (*c_temp).ij;
         kl = (*c_temp).kl;
         tmpval = (*c_temp).pval;

         gtmp[ij] += ptmp[kl]*tmpval;
         gtmp[kl] += ptmp[ij]*tmpval;
         }
   
      if (readflgc && j < num_bufs_c-1) {
         if(wherec==num_bufs_c) {
            srew(itap93);
            wherec=0;
            }
         wherec++;
         num=int_nums_c[wherec];
         sread(itap93,(char *) c_outbuf,sizeof(struct c_pkints)*num);
         }
      }

   for(k=joff=0; k < num_ir ; k++) {
      if(nn=scf_info[k].num_so) {
         for(i=0; i < nn ; i++) {
            if(twocon) {
               for(j=0; j <= i ; j++) {
                  scf_info[k].dpmat[ioff[i]+j] = gtmp[ioff[i+joff]+j+joff];
                  scf_info[k].dpmato[ioff[i]+j] = gtmpo[ioff[i+joff]+j+joff];
                  }
               }
            else {
               for(j=0; j <= i ; j++) {
                  scf_info[k].gmat[ioff[i]+j] += gtmp[ioff[i+joff]+j+joff];
                  scf_info[k].gmato[ioff[i]+j] += gtmpo[ioff[i+joff]+j+joff];
                  }
               }
            }
         joff += nn;
         }
      }
   }
