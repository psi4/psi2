/* $Log: errchk.c,v $
/* Revision 1.3  1995/07/21 17:37:06  psi
/* Made Jan 1st 1995 cscf the current accepted version of cscf.  Some
/* unidentified changes made after that date were causing problems.
/*
 * Revision 1.1  1991/06/15  20:22:24  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: errchk.c,v 1.3 1995/07/21 17:37:06 psi Exp $";

#include <stdio.h>
#include <ip_libv1.h>

errchk(errcod,token)
int errcod;
char *token;
{
  if (errcod) {
    fprintf(stderr,"ERROR: %s\n",ip_error_message(errcod));
    fprintf(stderr,"TOKEN: %s\n",token);
    }
  }
