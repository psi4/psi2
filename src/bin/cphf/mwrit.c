
/* $Log: mwrit.c,v $
/* Revision 1.1  1991/06/15 22:45:28  seidl
/* Initial revision
/* */

static char *rcsid = "$Id: mwrit.c,v 1.1 1991/06/15 22:45:28 seidl Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

mwrit(array,address)
   char *array;
   int address;

   {
      int length,iaddress;

      iaddress = block_locs[address];
      length = block_locs[address+nsect/2];

      rwrit(itap40,array,sizeof(int)*length,iaddress);
      }

mread(array,address)
   char *array;
   int address;

   {
      int length,iaddress;

      iaddress = block_locs[address];
      length = block_locs[address+nsect/2];

      rread(itap40,array,sizeof(int)*length,iaddress);
      }
