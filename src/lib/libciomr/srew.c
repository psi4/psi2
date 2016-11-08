
/* $Log: srew.c,v $
/* Revision 2.1  1991/06/15 18:30:10  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id: srew.c,v 2.1 1991/06/15 18:30:10 seidl Exp $";

#include "includes.h"
#include "pointers.h"

void srew(tape)
   int tape;

   {
      ptr.wptr[tape] = 0;
      }
