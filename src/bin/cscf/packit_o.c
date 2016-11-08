/* $Log: packit_o.c,v $
/* Revision 1.5  1997/08/25 21:51:32  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.4  1995/07/21  17:37:09  psi
 * Made Jan 1st 1995 cscf the current accepted version of cscf.  Some
 * unidentified changes made after that date were causing problems.
 *
 * Revision 1.2  1994/06/03  04:40:57  seidl
 * make c_outbuf extern
 *
 * Revision 1.1  1991/06/15  20:22:35  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: packit_o.c,v 1.5 1997/08/25 21:51:32 crawdad Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

extern double *pa, *pb;
static int iblc=0;
static int iblo=0;
static int num_ints_o=0;
static int num_ints_c=0;
static double *gtmp,*gtmpo,*ptmp,*ptmpo;

int lastc;
int lasto;

int readflgc=0;
int readflgo=0;
int num_bufs_o=0;
int num_bufs_c=0;

struct c_pkints {
         int ij;
         int kl;
         double pval;
         };
extern struct c_pkints *c_outbuf;
struct o_pkints {
         int ij;
         int kl;
         double pval;
         double qval;
         } *o_outbuf;

void packit_open(lbij,lbkl,endflg)
   unsigned int *lbij, *lbkl;
   int endflg;

{
   int i,j,k,joff,ij,kl;
   int tmpsiz,nn;
   double pval, qval;
   double tol = 10e-14;

   if(!o_outbuf) {
      if((c_outbuf =(struct c_pkints *) malloc(maxbuf*sizeof(struct c_pkints)))
                                                                   ==NULL) {
         fprintf(stderr,"cannot allocate memory for c_outbuf in packit\n");
         exit(4);
         }
      if((o_outbuf = (struct o_pkints *) malloc(maxbuf*sizeof(struct o_pkints)))
                                                                   ==NULL) {
         fprintf(stderr,"cannot allocate memory for o_outbuf in packit\n");
         exit(4);
         }
      }

   tmpsiz = ioff[nbasis];

   if(gtmp==NULL && !twocon) {
      gtmp = (double *) init_array(tmpsiz);
      ptmp = (double *) init_array(tmpsiz);
      gtmpo = (double *) init_array(tmpsiz);
      ptmpo = (double *) init_array(tmpsiz);

      for(k=joff=0; k < num_ir ; k++) {
         if(nn=scf_info[k].num_so) {
            for(i=0; i < nn ; i++)
               for(j=0; j <= i ; j++) {
                  ptmp[ioff[i+joff]+j+joff] = scf_info[k].pmat[ioff[i]+j];
                  ptmpo[ioff[i+joff]+j+joff] = scf_info[k].pmato[ioff[i]+j];
                  }
            joff += nn;
            }
         }
      }

   if(!endflg) {
      for(i=0; i < nint ; i++) {
         pval=pa[i];
         qval=pb[i];
         ij = lbij[i];
         kl = lbkl[i];
         if(print & 16) fprintf(outfile,"%5d%5d%9.5f%9.5f\n",ij,kl,pval,qval);
         if (fabs(pval) >= tol || fabs(qval) >= tol) {
            if(qval) {
               o_outbuf[iblo].ij = ij;
               o_outbuf[iblo].kl = kl;
               o_outbuf[iblo].pval = pval;
               o_outbuf[iblo].qval = qval;

               if(!twocon) {
                  gtmp[ij] += ptmp[kl]*pval;
                  gtmp[kl] += ptmp[ij]*pval;
                  gtmpo[ij] += ptmpo[kl]*qval;
                  gtmpo[kl] += ptmpo[ij]*qval;
                  }

               iblo++;

               if (iblo >= maxbuf) {
                  readflgo=1;
                  swrit(itap92,(char *) o_outbuf,sizeof(struct o_pkints)*maxbuf);
                  num_ints_o += iblo;
                  num_bufs_o++;
                  iblo=0;
                  }
               }
            else {
               c_outbuf[iblc].ij = ij;
               c_outbuf[iblc].kl = kl;
               c_outbuf[iblc].pval = pval;

               if(!twocon) {
                  gtmp[ij] += ptmp[kl]*pval;
                  gtmp[kl] += ptmp[ij]*pval;
                  }

               iblc++;

               if (iblc >= maxbuf) {
                  readflgc=1;
                  swrit(itap93,(char *) c_outbuf,sizeof(struct c_pkints)*maxbuf);
                  num_ints_c += iblc;
                  num_bufs_c++;
                  iblc=0;
                  }
               }
            }
         }
      nint=0;
      }
   else {
      num_ints_o += iblo;
      num_ints_c += iblc;
      num_bufs_c++;
      num_bufs_o++;
      fprintf(outfile,"\n%10d integrals written to file92 in %3d buffers\n",
                                       num_ints_o,num_bufs_o);
      fprintf(outfile,"%10d integrals written to file93 in %3d buffers\n",
                                       num_ints_c,num_bufs_c);
      lasto = iblo;
      lastc = iblc;
      if(readflgo)
         swrit(itap92,(char *) o_outbuf,sizeof(struct o_pkints)*iblo);
      if(readflgc)
         swrit(itap93,(char *) c_outbuf,sizeof(struct c_pkints)*iblc);

      if(!twocon) {
         for(k=joff=0; k < num_ir ; k++) {
            if(nn=scf_info[k].num_so) {
               for(i=0; i < nn ; i++)
                  for(j=0; j <= i ; j++) {
                     scf_info[k].gmat[ioff[i]+j] = gtmp[ioff[i+joff]+j+joff];
                     scf_info[k].gmato[ioff[i]+j] = gtmpo[ioff[i+joff]+j+joff];
                     }
               joff += nn;
               }
            }
         free(gtmp);
         free(ptmp);
         free(gtmpo);
         free(ptmpo);
         }
      }
   }
