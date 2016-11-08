/* $Id: includes.h,v 1.3 1998/06/30 14:11:07 sbrown Exp $ */
/* $Log: includes.h,v $
/* Revision 1.3  1998/06/30 14:11:07  sbrown
/* *************************************************************
/* *Program Modification                                       *
/* *By: Shawn Brown                                            *
/* *Date: June 30, 1998                                        *
/* *Altered program to make a guess at occupations from the    *
/* *diagonalized core hamiltonian matrix.  Program can now     *
/* *make a guess at the beginning of the calculation or at     *
/* *or at every iteration.  Use the latter at your own risk.   *
/* *See man pages for details on new keywords.                 *
/* *************************************************************
/*
 * Revision 1.2  1997/08/25  21:51:27  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 1.1  1991/06/15  20:22:51  seidl
 * Initial revision
 * */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <libciomr.h>
#include <file30.h>

#define MIN0(a,b) ((a)<(b)) ? (a) : (b)
#define MAX0(a,b) ((a)>(b)) ? (a) : (b)

#ifdef AIXV3
# define MALLOC malloc
#else
# define MALLOC malloc
#endif
