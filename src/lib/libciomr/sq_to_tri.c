
/* $Log: sq_to_tri.c,v $
/* Revision 2.1  1991/06/15 18:30:08  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id: sq_to_tri.c,v 2.1 1991/06/15 18:30:08 seidl Exp $";

#include "includes.h"

/* converts square matrix to lower triangle */

void sq_to_tri(bmat,amat,size)
     double *amat, **bmat;
     int size;

     {
        int i, j, ij;

        ij=0;
        for(i = 0 ; i < size ; i++) {
            for(j = 0 ; j <= i ; j++) {
               amat[ij++] = bmat[i][j];
               }
            }
     }

