
/* $Id: iomrparam.h.in,v 2.2 1997/09/12 13:52:51 crawdad Exp $ */
/* $Log: iomrparam.h.in,v $
/* Revision 2.2  1997/09/12 13:52:51  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
# Revision 2.1  1997/08/25  21:49:54  crawdad
# Making changes for extension of PSI file size limit.
#
 * Revision 2.1  1997/07/09  18:46:42  crawdad
 * Moved param.h to iomrparam.h to avoid name conflicts.
 *
 * Revision 2.2  1995/10/04  16:47:29  crawdad
 * Increased MAX_UNIT to 300
 *
 * Revision 2.1  1991/06/15  18:32:27  seidl
 * *** empty log message ***
 * */

#ifndef IOMRPARAM_H
#define IOMRPARAM_H

/* Change this macro to unsigned long long int for systems that allow it */
#define PSI_FPTR unsigned long int

#define MAX_UNIT 300
#define MAX_VOLUME 8
#define MAX_STRING 512

#define PSIIO_UNOPENED 0
#define PSIIO_SEQUENTIAL 1
#define PSIIO_S_ASYNC 2
#define PSIIO_R_ASYNC 3
#define PSIIO_RAM 4
#define PSIIO_FORTRAN 5

#define IOOP_READ 1
#define IOOP_WRITE 2

#endif
