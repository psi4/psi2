
/* $Id: common.h,v 2.3 1994/06/02 02:28:49 seidl Exp $ */
/* $Log: common.h,v $
/* Revision 2.3  1994/06/02 02:28:49  seidl
/* define ALLOC_GLOBALS
/*
 * Revision 2.2  1991/08/21  05:42:26  psi
 * remove everthing, this should go away eventually
 *
 * Revision 2.1  1991/06/15  18:32:25  seidl
 * *** empty log message ***
 * */

#ifdef ALLOC_GLOBALS
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN FILE *intfile, *outfile;

