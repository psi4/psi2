
/* $Log: init_ptrs.c,v $
/* Revision 2.7  1997/09/23 12:53:33  crawdad
/* Correcting nasty, hideous, stupied little bug that annoyed me for a week.
/*
 * Revision 2.6  1997/09/12  13:52:50  crawdad
 * Changing marco name from ULL to PSI_FPTR.
 *
 * Revision 2.5  1997/08/25  21:49:52  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.4  1997/06/23  12:25:48  crawdad
 * Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
 *     conflicts with similarly named system file under linux.  Corrected type
 *    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
 *     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
 *    avoid malloc'ing zero-length arrays.
 *
 * -Daniel
 *
 * Revision 2.3  1996/07/03  11:26:38  psi
 * Added free_ptrs() function in init_ptrs.c.  Added free() calls to tstart()
 * and tstop().
 *
 * Revision 2.2  1996/07/02 21:06:29  sherrill
 * Completed the changes needed to use unit numbers greater than 99; changed
 * a hardwired 99 in init_ptrs.c to MAX_UNIT and increased a "unit" string
 * from 3 chars to 4.  Also removed a compiler warning in sequential.c by
 * casting ud to (char *) for malloc_check().
 *
 * Revision 2.1  1991/06/15  18:29:21  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id: init_ptrs.c,v 2.7 1997/09/23 12:53:33 crawdad Exp $";

#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

/* allocates memory for array of file pointers */

void init_ptrs()
     {
     int num_ptrs = MAX_UNIT;

     ptr.wptr = (PSI_FPTR *) malloc(sizeof(PSI_FPTR)*num_ptrs);

     if (ptr.wptr == NULL) {
        fprintf(stderr,"trouble allocating memory for pointers!\n");
        exit(1);
        }
     }

void free_ptrs()
     {
     free(ptr.wptr);
     }
