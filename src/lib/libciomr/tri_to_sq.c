
/* $Log: tri_to_sq.c,v $
/* Revision 2.1  1991/06/15 18:30:15  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id: tri_to_sq.c,v 2.1 1991/06/15 18:30:15 seidl Exp $";

#include "includes.h"

/* converts lower triangle to square matrix */

void tri_to_sq(amat,bmat,size)
     double *amat, **bmat;
     int size;

     {
        int i, j, ij;

        ij=0;
        for(i = 0 ; i < size ; i++) {
            for(j = 0 ; j <= i ; j++) {
               bmat[i][j] = amat[ij];
               bmat[j][i] = amat[ij++];
               }
            }
     }
