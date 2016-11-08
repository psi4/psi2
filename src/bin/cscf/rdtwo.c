/* $Log: rdtwo.c,v $
/* Revision 1.5  1997/08/25 21:51:34  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.4  1995/11/07  19:24:12  sherrill
 * Replace all initializations of integer arrays with calls to init_int_array().
 * There had been problems with using init_array() because the division by
 * the doubles-to-int factor was rounding off (integer division).
 *
 * Revision 1.3  1995/07/21  17:37:10  psi
 * Made Jan 1st 1995 cscf the current accepted version of cscf.  Some
 * unidentified changes made after that date were causing problems.
 *
 * Revision 1.1  1991/06/15  20:22:37  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: rdtwo.c,v 1.5 1997/08/25 21:51:34 crawdad Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <ip_libv1.h>

double *pa, *pb;
int *inext;
unsigned int *lbij,*lbkl;
int intmx=25920;

void rdtwo()

{
   int ilsti, nbuf;
   int ibufsz = 8942;
   int ibufqt = 2236;
   int i;
   int iii;
   int delete34,errcod;
   char *wfn="SCF",*dertyp="FIRST";
   int ior, ism, jor, jsm;
   int kor, ksm, lor, lsm;
   int ii,jj,kk,ll;
   int p1,p2;
   int d2i = sizeof(double)/sizeof(int);
   double pki_int;
   int *tmp;
   union psi_buffer inbuf;

   nint=0;

   if(nbasis > 150) intmx *= 4;

   inbuf.pki = (double *) init_array((int) ibufsz/d2i);
   pa = (double *) init_array(intmx);
   pb = (double *) init_array(intmx);

   /* lbij = (unsigned int *) init_array((int) intmx/d2i); */
   /* lbkl = (unsigned int *) init_array((int) intmx/d2i); */
   lbij = (unsigned int *) init_int_array(intmx);
   lbkl = (unsigned int *) init_int_array(intmx);

   do {
      sread(itap34,(char *) inbuf.lbli,sizeof(int)*ibufsz);
      ilsti=inbuf.lbli[0];
      nbuf=inbuf.lbli[1];

      if (print & 8) fprintf(outfile,"%5d\n",nbuf);

      tmp = &inbuf.lbli[2];

      for (i=0 ; i < nbuf ; i++,tmp += 2) {
         lor = *(tmp+1) >> 2;
         lsm = lor >> 8;
         kor = lsm >> 3;
         ksm = kor >> 8;
         kor = (kor & 255) - 1;
         lsm = lsm & 7;
         lor = (lor & 255) - 1;
         iii = *(tmp+1) & 3;
         jsm = *tmp >> 8;
         ior = jsm >> 3;
         ism = ior >> 8;
         ior = (ior & 255) - 1;
         jsm = jsm & 7;
         jor = (*tmp & 255) - 1;
         ii = ior+scf_info[ism].ideg;
         jj = jor+scf_info[jsm].ideg;
         kk = kor+scf_info[ksm].ideg;
         ll = lor+scf_info[lsm].ideg;
         pki_int = inbuf.pki[i+ibufqt];

         if (ism == jsm && ksm == lsm && ism == ksm) {
            if (ior == jor && ior == kor || jor == kor && jor == lor) {
               findit(ii,jj,kk,ll,ism,ksm,pki_int,5);
               }
            else if (ior == kor || jor == lor) {
               findit(ii,jj,kk,ll,ism,ksm,pki_int,3);

               p1 = MAX0(jj,ll);
               p2 = MIN0(jj,ll);

               findit(ii,kk,p1,p2,ism,ksm,pki_int,4);
               }
            else if (jor == kor) {
               findit(ii,jj,kk,ll,ism,ksm,pki_int,3);

               findit(ii,ll,jj,kk,ism,ksm,pki_int,4);
               }
            else if (ior == jor || kor == lor) {
               findit(ii,jj,kk,ll,ism,ksm,pki_int,1);

               p1 = MAX0(jj,ll);
               p2 = MIN0(jj,ll);

               findit(ii,kk,p1,p2,ism,ksm,pki_int,2);
               }
            else {
               findit(ii,jj,kk,ll,ism,ksm,pki_int,1);

               p1 = MAX0(jj,ll);
               p2 = MIN0(jj,ll);

               findit(ii,kk,p1,p2,ism,ksm,pki_int,2);

               p1 = MAX0(jj,kk);
               p2 = MIN0(jj,kk);

               findit(ii,ll,p1,p2,ism,ksm,pki_int,2);
               }
            }
         else if (ism == jsm) {
            findit(ii,jj,kk,ll,ism,ksm,pki_int,1);
            }
         else if (ism == ksm) {
            if (ior == kor || jor == lor) {
               p1 = MAX0(jj,ll);
               p2 = MIN0(jj,ll);

               findit(ii,kk,p1,p2,ism,jsm,pki_int,4);
               }
            else {
               p1 = MAX0(jj,ll);
               p2 = MIN0(jj,ll);

               findit(ii,kk,p1,p2,ism,jsm,pki_int,2);
               }
            }
         if(iii && nint) {
            if(iopen) packit_open(lbij,lbkl,0);
            else packit_closed(lbij,lbkl,0);
            }
         }
       } while(!ilsti);

   if(iopen) packit_open(lbij,lbkl,1);
   else packit_closed(lbij,lbkl,1);

   if(!iopen) fprintf(outfile,"  wrote %d integrals to file92\n",num_ints);
   free(inbuf.lbli);
   free(pa);
   free(pb);
   free(lbij);
   free(lbkl);
   free(inext);

   errcod = ip_string("WFN",&wfn,0);
   errcod = ip_string("DERTYPE",&dertyp,0);

   delete34 = 0;
   if(!strcmp(wfn,"SCF") && (!strcmp(dertyp,"FIRST") || !strcmp(dertyp,"NONE")))
     delete34 = 1;
   errcod = ip_boolean("DELETE34",&delete34,0);

   if(delete34) rclose(itap34,4);
   else rclose(itap34,3);
   }
