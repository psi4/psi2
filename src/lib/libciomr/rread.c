
/* $Log: rread.c,v $
/* Revision 2.6  1999/11/01 20:10:58  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.5  1997/09/12 13:52:58  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 2.4  1997/08/25  21:50:03  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.3  1997/06/23  12:25:57  crawdad
 * Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
 *     conflicts with similarly named system file under linux.  Corrected type
 *    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
 *     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
 *    avoid malloc'ing zero-length arrays.
 *
 * -Daniel
 *
 * Revision 2.2  1995/04/01  20:53:01  fermann
 * changed bytewise file pointers such as first, last and length to long
 * unsigned ints in order to handle up to 4 gigabyte tmp files (striped into
 * individual pieces of less than 2 gigabytes).  added functions li2sec and
 * sec2li for where they are needed.
 *
 * Revision 2.1  1991/06/15  18:29:55  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id: rread.c,v 2.6 1999/11/01 20:10:58 evaleev Exp $";

#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern void wreadw(int, char *, int, PSI_FPTR, PSI_FPTR *);
extern PSI_FPTR sec2i(int);

/* reads nlen bytes into array starting at irec */

void rread(itape,array,nlen,irec)
   int itape, nlen, irec;
   char *array;
   {
      PSI_FPTR ipos, junk;

      ipos = sec2i(--irec);
      wreadw(itape,(char *) array,nlen,ipos,&junk);
   }
