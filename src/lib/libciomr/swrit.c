
/* $Log: swrit.c,v $
/* Revision 2.6  1999/11/01 20:11:00  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.5  1997/09/12 13:53:06  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 2.4  1997/08/25  21:50:11  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.3  1997/06/23  12:26:02  crawdad
 * Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
 *     conflicts with similarly named system file under linux.  Corrected type
 *    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
 *     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
 *    avoid malloc'ing zero-length arrays.
 *
 * -Daniel
 *
 * Revision 2.2  1995/04/01  20:53:08  fermann
 * changed bytewise file pointers such as first, last and length to long
 * unsigned ints in order to handle up to 4 gigabyte tmp files (striped into
 * individual pieces of less than 2 gigabytes).  added functions li2sec and
 * sec2li for where they are needed.
 *
 * Revision 2.1  1991/06/15  18:30:11  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id: swrit.c,v 2.6 1999/11/01 20:11:00 evaleev Exp $";

#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern void wwritw(int, char *, int, PSI_FPTR, PSI_FPTR *);

/* writes nlen bytes from array to file itape starting at current */
/* pointer location */

void swrit(itape,array,nlen)
   int itape, nlen;
   char *array;
   {
      PSI_FPTR start, end;

      start=ptr.wptr[itape];

      wwritw(itape,(char *) array,nlen,start,&end);
      ptr.wptr[itape]=((ptr.wptr[itape]-1+4096)/4096)*4096;
   }
