#include <stdio.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_init()  Initializes the file30 struct to allow other file30_
**    functions to perform their duties.
**
**  arguments: none, but it requires that the input parser be initialized
**    so that it can open file30.  
**
**  returns: zero.  Perhaps this will change some day.
*/

/* First definition of global struct file30 */
struct file30_ info30_;

int file30_init(void)
{
  PSI_FPTR junk;
  info30_.filenum = 30;
  rfile(info30_.filenum);

  /* Read the file30 title string */
  wreadw(info30_.filenum, (char *) info30_.label, sizeof(char)*20, 0, &junk);
  info30_.label[19] = '\0'; 

  /* Grab the array of pointers to data */
  info30_.mpoint = init_int_array(MPOINT);
  wreadw(info30_.filenum, (char *) info30_.mpoint,
	 sizeof(unsigned long int)*MPOINT, sizeof(unsigned long int)*300,
	 &junk);

  /* Grab the array of calculation constants */
  info30_.mconst = init_int_array(MCONST);
  wreadw(info30_.filenum, (char *) info30_.mconst,
	 sizeof(unsigned long int)*MCONST, sizeof(unsigned long int)*100,
	 &junk);

  /* Grab the array of pointers to the SCF calculations */
  info30_.mcalcs = init_int_array(MCALCS);
  wreadw(info30_.filenum, (char *) info30_.mcalcs,
	 sizeof(unsigned long int)*MCALCS, sizeof(unsigned long int)*500,
	 &junk);

  return 0;
}
