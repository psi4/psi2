
/* $Log: ram.c,v $
/* Revision 2.6  1999/11/01 20:10:58  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.5  1997/09/12 13:52:56  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 2.4  1997/08/25  21:50:01  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.3  1997/06/23  12:25:55  crawdad
 * Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
 *     conflicts with similarly named system file under linux.  Corrected type
 *    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
 *     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
 *    avoid malloc'ing zero-length arrays.
 *
 * -Daniel
 *
 * Revision 2.2  1996/06/18  20:47:45  sherrill
 * Add the whole set of int_array routines to int_array.c (replacing
 * init_int_array.c), add block_matrix.c, and add a new function flen
 * which gets the file size for unit number x.
 *
 * Revision 2.1  1991/06/15  18:29:46  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id: ram.c,v 2.6 1999/11/01 20:10:58 evaleev Exp $";


#include "includes.h"
#include "iomrparam.h"
#include "types.h"

ram_t * ram_ioopen(char *param, int unit)
{
  fprintf(stderr,"no ram io yet\n");
  ioabort();
  }

void ram_ioclos(ram_t *ud, int status)
{
  fprintf(stderr,"no ram io yet\n");
  ioabort();
  }

void ram_iordr(ram_t *ud, char *buffer, PSI_FPTR first, int length)
{
  fprintf(stderr,"no ram io yet\n");
  ioabort();
  }

void ram_iowrr(ram_t *ud, char *buffer, PSI_FPTR first, int length)
{
  fprintf(stderr,"no ram io yet\n");
  ioabort();
  }

PSI_FPTR ram_iosize(ram_t *ud)
{
  fprintf(stderr,"no ram io yet\n");
  ioabort();
  }

