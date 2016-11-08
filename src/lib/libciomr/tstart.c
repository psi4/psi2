
/* $Log: tstart.c,v $
/* Revision 2.7  1996/07/03 11:26:39  psi
/* Added free_ptrs() function in init_ptrs.c.  Added free() calls to tstart()
/* and tstop().
/*
 * Revision 2.6  1994/09/19 23:32:12  cljanss
 * Cleaned up allocation of globals.  Got rid of some globals and fixed a bug.
 *
 * Revision 2.5  1994/08/19  17:26:41  cljanss
 * Will now compile on IRIX 4.0.5F.
 *
 * Revision 2.4  1994/08/10  00:17:26  dcrawfrd
 * Added timing routines that are called by iomr and that give hostnames.
 *
 * Revision 2.3  1994/06/03  17:08:30  seidl
 * define ALLOC_GLOBALS before including pointers.h
 *
 * Revision 2.2  1994/06/02  02:27:10  seidl
 * no vtimes on SGI
 *
 * Revision 2.1  1991/06/15  18:30:16  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id: tstart.c,v 2.7 1996/07/03 11:26:39 psi Exp $";

#include "includes.h"
#define ALLOC_GLOBALS
#include "pointers.h"
#undef ALLOC_GLOBALS

#include <unistd.h>
#include <string.h>

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif

#ifdef HAVE_SYS_PARAM_H
/* This has HZ, which gives the number of clock ticks per second. */
#include <sys/param.h>
#endif

/* If HZ is missing, take a guess. */
#ifndef HZ
#define HZ 60
#endif

void tstart(outfile)
   FILE *outfile;

   {
       int i,error;
       char *name;
       name = (char *) malloc(26 * sizeof(char));
       error = gethostname(name, 26);
       if(error != 0) strncpy(name,"nohostname", 11);

       time_start = time(NULL);

       for (i=0; i < 78 ; i++) {
           fprintf(outfile,"*");
          }
       fprintf(outfile,"\n");

       fprintf(outfile,"tstart called on %s\n", name);
       fprintf(outfile,"%s\n",ctime(&time_start));

       free(name);

       }

void tstop(outfile)
   FILE *outfile;

{
   int i;
   int error;
   time_t total_time;
   struct tms total_tmstime;
   char *name;
   name = (char *) malloc(26 * sizeof(char));
   error = gethostname(name, 26);
   if(error != 0) strncpy(name,"nohostname", 11);

   time_end = time(NULL);
   total_time = time_end - time_start;

   times(&total_tmstime);


   for (i=0; i < 78 ; i++) {
       fprintf(outfile,"*");
     }
   fprintf(outfile,"\n");
   fprintf(outfile,"tstop called on %s\n", name);
   fprintf(outfile,"%s\n",ctime(&time_end));
   fprintf(outfile,"user time   = %10.2f seconds\n",
           ((double) total_tmstime.tms_utime)/HZ);
   fprintf(outfile,"system time = %10.2f seconds\n",
           ((double) total_tmstime.tms_stime)/HZ);
   fprintf(outfile,"total time  = %10d seconds\n",total_time);

   free(name);

   }
