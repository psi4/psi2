/* $Log: findit.c,v $
/* Revision 1.3  1996/11/26 05:20:52  sherrill
/* Add casts in front of init_int_array() calls to avoid compiler warnings
/* (really should include libciomr everywhere instead but that causes more
/* warnings because wwritw's need casts to char *).  Also fixed problem
/* where phase() was trying to free unallocated memory.
/*
 * Revision 1.2  1995/11/07  19:24:05  sherrill
 * Replace all initializations of integer arrays with calls to init_int_array().
 * There had been problems with using init_array() because the division by
 * the doubles-to-int factor was rounding off (integer division).
 *
 * Revision 1.1  1991/06/15  20:22:25  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: findit.c,v 1.3 1996/11/26 05:20:52 sherrill Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

extern double *pa, *pb;
extern int *inext;
extern unsigned int *lbij,*lbkl;
extern int intmx;
static int old_nint=25920;

void findit(ii,jj,kk,ll,ism,ksm,value,iab)
   double value;
   int ii,jj,kk,ll;
   int ism, ksm, iab;

{
   register int i,j;
   unsigned int *ijtmp, *kltmp;
   int *nxtmp;
   int p1,p2,p3;
   int noi,nok;
   int next,start;
   int lij,lkl;
   int keep=128;
   int keep2=127;
   int d2i = sizeof(double)/sizeof(int);
   double *patmp, *pbtmp;

   if(nbasis > 150) {
     keep=1024;
     keep2=1023;
     }

   /* if(!inext) inext = (int *) init_array((int) (keep+intmx)/d2i); */
   if (!inext) inext = (int *) init_int_array(keep+intmx);

   noi = scf_info[ism].nopen+scf_info[ism].nhalf;
   nok = scf_info[ksm].nopen+scf_info[ksm].nhalf;

   lij = ioff[ii]+jj;
   lkl = ioff[kk]+ll;

   if(!nint) {
      bzero(inext,sizeof(int)*old_nint);
      bzero(&inext[intmx],sizeof(int)*keep);
      }

   start = 2*lij + lkl;
   start = (start & keep2) + intmx;

L1:
   next=inext[start];
   if(next) {
      if (lbij[next-1] == lij && lbkl[next-1] == lkl) i=next-1;
      else {
         start = next;
         goto L1;
         }
      }
   else {
      i=nint;
      if(nint >= intmx) {
        fprintf(outfile,"\n  increasing size of buffers in findit\n");
        fprintf(outfile,"  intmx was %d, is %d\n",intmx,intmx*2);
        fflush(outfile);
        intmx*=2;

 /* i don't use realloc because strange things were happening */

        /* nxtmp = (int *) init_array((int) (keep+intmx)/d2i); */
        nxtmp = (int *) init_int_array(keep+intmx);
        bcopy(inext,nxtmp,(int)sizeof(int)*(intmx/2));
        for(j=0; j < keep ; j++) nxtmp[j+intmx]=inext[j+intmx/2];
        free(inext);
        inext=nxtmp;

        /* ijtmp = (unsigned int *) init_array(intmx/d2i); */
        ijtmp = (unsigned int *) init_int_array(intmx);
        bcopy(lbij,ijtmp,sizeof(int)*(intmx/2));
        free(lbij);
        lbij = ijtmp;

        /* kltmp = (unsigned int *) init_array(intmx/d2i); */
        kltmp = (unsigned int *) init_int_array(intmx);
        bcopy(lbkl,kltmp,sizeof(int)*(intmx/2));
        free(lbkl);
        lbkl = kltmp;

        patmp = (double *) init_array(intmx);
        bcopy(pa,patmp,sizeof(double)*(intmx/2));
        free(pa);
        pa = patmp;

        pbtmp = (double *) init_array(intmx);
        bcopy(pb,pbtmp,sizeof(double)*(intmx/2));
        free(pb);
        pb = pbtmp;

        if(inext==NULL || lbij==NULL || lbkl==NULL) {
          fprintf(outfile,"\n pathological problems with realloc in findit\n");
          fprintf(outfile," try upping intmx to %d\n",intmx);
          fflush(outfile);
          exit(1);
          }

        start = 2*lij + lkl;
        start = (start & keep2) + intmx;
        }
      inext[start] = ++nint;
      lbij[i] = lij;
      lbkl[i] = lkl;
      pa[i] = pb[i] = 0.0;
      }

   value = (lij == lkl) ? value*0.5 : value;

   switch(iab) {
      case 1:
         pa[i] += value;
         if(noi && nok) {
            if(special) {
               p1 = MAX0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p2 = MIN0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p3 = ioff[p1]+p2;
               pb[i] += alpha[p3]*value;
               }
            }
         break;
      case 2:
         pa[i] -= 0.25*value;
         if(noi && nok) {
            if(hsos) pb[i] += 0.25*value;
            else if(singlet) {
               if (ism != ksm) pb[i] -= 0.75*value;
               else pb[i] += 0.25*value;
               }
            else if(twocon) {
               if (ism != ksm) pb[i] += 0.25*value;
               }
            else {
               p1 = MAX0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p2 = MIN0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p3 = ioff[p1]+p2;
               pb[i] += 0.25*beta[p3]*value;
               }
            }
         break;
      case 3:
         pa[i] += 0.75*value;
         if(noi && nok) {
            if(hsos) pb[i] += 0.25*value;
            else if(singlet) {
               if (ism != ksm) pb[i] -= 0.75*value;
               else pb[i] += 0.25*value;
               }
            else if(twocon) {
               if (ism != ksm) pb[i] += 0.25*value;
               }
            else {
               p1 = MAX0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p2 = MIN0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p3 = ioff[p1]+p2;
               pb[i] += (alpha[p3] + 0.25*beta[p3])*value;
               }
            }
         break;
      case 4:
         pa[i] -= 0.5*value;
         if(noi && nok) {
            if(hsos) pb[i] = 0.5*value;
            else if(singlet) {
               if (ism != ksm) pb[i] -= 1.5*value;
               else pb[i] += 0.5*value;
               }
            else if(twocon) {
               if (ism != ksm) pb[i] += 0.5*value;
               }
            else {
               p1 = MAX0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p2 = MIN0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p3 = ioff[p1]+p2;
               pb[i] += 0.5*beta[p3]*value;
               }
            }
         break;
      case 5:
         pa[i] = 0.5*value;
         if(noi && nok) {
            if(hsos) pb[i] = 0.5*value;
            else if(singlet) {
               if (ism != ksm) pb[i] = -1.5*value;
               else pb[i] = 0.5*value;
               }
            else if(twocon) {
               if (ism != ksm) pb[i] += 0.5*value;
               }
            else {
               p1 = MAX0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p2 = MIN0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p3 = ioff[p1]+p2;
               pb[i] = (alpha[p3] + 0.5*beta[p3])*value;
               }
            }
      }
   old_nint=nint;
   }
