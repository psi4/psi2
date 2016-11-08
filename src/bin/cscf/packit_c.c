/* $Log: packit_c.c,v $
/* Revision 1.2  1997/08/25 21:51:31  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  20:22:34  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: packit_c.c,v 1.2 1997/08/25 21:51:31 crawdad Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

int last;

extern double *pa, *pb;
static int ibl=0;
static double *gtmp,*ptmp;
struct c_pkints {
         int ij;
         int kl;
         double pval;
         } *c_outbuf;

void packit_closed(lbij,lbkl,endflg)
   unsigned int *lbij, *lbkl;
   int endflg;

{
   int i,j,k,joff,ij,kl;
   int tmpsiz,nn;
   double pval;
   double tol = 10e-14;

   if(!c_outbuf) {
      if((c_outbuf=(struct c_pkints *) malloc(maxbuf*sizeof(struct c_pkints)))
                                                                   ==NULL) {
         fprintf(stderr,"cannot allocate memory for c_outbuf in packit\n");
         exit(4);
         }
      }

   tmpsiz = ioff[nbasis];

   if(gtmp==NULL) {
      gtmp = (double *) init_array(tmpsiz);
      ptmp = (double *) init_array(tmpsiz);
      for(k=joff=0; k < num_ir ; k++) {
         if(nn=scf_info[k].num_so) {
            for(i=0; i < nn ; i++)
               for(j=0; j <= i ; j++)
                  ptmp[ioff[i+joff]+j+joff] = scf_info[k].pmat[ioff[i]+j];
            joff += nn;
            }
         }
      }

   if(!endflg) {
      for(i=0; i < nint ; i++) {
         pval=pa[i];
         ij = lbij[i];
         kl = lbkl[i];
         if(print & 128) fprintf(outfile,"%5d%5d%9.5f\n",ij,kl,pval);
         if (fabs(pval) >= tol) {
            c_outbuf[ibl].ij = ij;
            c_outbuf[ibl].kl = kl;
            c_outbuf[ibl].pval = pval;

            gtmp[ij] += ptmp[kl]*pval;
            gtmp[kl] += ptmp[ij]*pval;

            ibl++;

            if (ibl >= maxbuf) {
               if(readflg)
                  swrit(itap92,(char *) c_outbuf,sizeof(struct c_pkints)*maxbuf);
               num_ints += ibl;
               if(print & 16)
                  fprintf(outfile,"buf %3d: ibl = %10d\n",num_bufs,ibl);
               fflush(outfile);
               num_bufs++;
               ibl=0;
               }
            }
         }
      nint=0;
      }
   else {
      num_ints += ibl;
      if(print & 16) fprintf(outfile,"buf %3d: ibl = %10d\n",num_bufs,ibl);
      fflush(outfile);
      num_bufs++;
      last = ibl;
      if(readflg) swrit(itap92,(char *) c_outbuf,sizeof(struct c_pkints)*ibl);

      for(k=joff=0; k < num_ir ; k++) {
         if(nn=scf_info[k].num_so) {
            for(i=0; i < nn ; i++)
               for(j=0; j <= i ; j++)
                  scf_info[k].gmat[ioff[i]+j] = gtmp[ioff[i+joff]+j+joff];
            joff += nn;
            }
         }
      free(gtmp);
      free(ptmp);
      }
   }
