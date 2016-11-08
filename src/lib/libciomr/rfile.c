
/* $Log: rfile.c,v $
/* Revision 2.2  1999/11/01 20:10:58  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.1  1991/06/15 18:29:49  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id: rfile.c,v 2.2 1999/11/01 20:10:58 evaleev Exp $";

#include "includes.h"
#include "pointers.h"

extern void ioopen_(int *);
extern void init_ptrs(void);

/* opens binary file unit */

void rfile(int unit)
{
       if (ptr.wptr == NULL) init_ptrs();

       sector = 1024;
       ptr.wptr[unit] = 0;
       ioopen_(&unit);
    }



