
/* $Log: rsetsa.c,v $
/* Revision 2.5  1999/11/01 20:10:59  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.4  1997/09/12 13:52:59  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 2.3  1997/08/25  21:50:04  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.2  1995/04/01  20:53:02  fermann
 * changed bytewise file pointers such as first, last and length to long
 * unsigned ints in order to handle up to 4 gigabyte tmp files (striped into
 * individual pieces of less than 2 gigabytes).  added functions li2sec and
 * sec2li for where they are needed.
 *
 * Revision 2.1  1991/06/15  18:29:56  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id: rsetsa.c,v 2.5 1999/11/01 20:10:59 evaleev Exp $";

#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern PSI_FPTR sec2i(int);

/* sets file pointer for unit to sector address */

void rsetsa(unit,address)
   int unit;
   int address;
   {
      PSI_FPTR ipos;

      ipos = (PSI_FPTR) sec2i(--address);
      ptr.wptr[unit]=ipos;
    }
