/* $Log: rdone.c,v $
/* Revision 1.2  1997/08/25 21:51:33  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  20:22:36  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: rdone.c,v 1.2 1997/08/25 21:51:33 crawdad Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

void rdone(oei)
   int oei;

{
   int ilsti, nbuf;
   int ibufsz = 8942;
   int ibufs3 = 1491;
   int i, iijj;
   int ior, ism, jor, jsm;

   union bufs {
      int *lbli;
      double *stvi;
      } buffer;

   buffer.stvi = (double *) init_array(ibufsz);

   do {
      sread(itap34,(char *) (buffer.lbli),sizeof(int)*ibufsz);
      pos34 += sizeof(int)*ibufsz;
      pos34 = ((pos34-1+4096)/4096)*4096;
      ilsti=buffer.lbli[0];
      nbuf=buffer.lbli[1];

      for (i=0 ; i < nbuf ; i++) {
         jsm = buffer.lbli[i+2] >> 8;
         ior = jsm >> 3;
         ism = ior >> 8;
         ior = (ior & 255)-1;
         jsm = jsm & 7;
         jor = (buffer.lbli[i+2] & 255)-1;
         iijj = ioff[ior]+jor;
         switch(oei) {
         case SMAT:
            scf_info[ism].smat[iijj]=buffer.stvi[i+ibufs3];
            break;
         case TMAT:
            scf_info[ism].tmat[iijj]=buffer.stvi[i+ibufs3];
            break;
         case VMAT:
            scf_info[ism].hmat[iijj]=
                               buffer.stvi[i+ibufs3]+scf_info[ism].tmat[iijj];
            break;
            }
         }
       } while(!ilsti);

   free(buffer.stvi);
   }
