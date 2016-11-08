#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_contr(): Reads in the normalized contraction coefficients.
**
**  takes no arguments.
**
**  returns: double *contr   Normalized contraction coefficients are 
**  returned as an array of doubles. In file30 they are stored as a matrix 5 
**  (only 5 types of angular momentum functions currently supported in PSI,
**  namely s-, p- , d-, f-, and g-) by the total number of primitives nprim,
**  but each primitive Gaussian contributes to only one shell (and one 
**  basis function, of course), so most of these values are zero and not returned.
*/


double *file30_rd_contr(void)
{
  double *contr;
  double *temp_contr;
  int nprim, i, j, ij = 0;
  PSI_FPTR junk;
  PSI_FPTR contr_ptr;

  nprim = file30_rd_nprim();
  contr_ptr = (PSI_FPTR) (info30_.mpoint[5] - 1)*sizeof(int);

  temp_contr = init_array(MAXL*nprim);
  contr = init_array(nprim);

  wreadw(info30_.filenum, (char *) temp_contr, (int) MAXL*nprim*sizeof(double),
	 contr_ptr, &junk);

/* Picking non-zero coefficients to the "master" array contr */
  for(i=0; i < MAXL; i++) 
   for(j=0; j < nprim; j++)
    { if (temp_contr[ij] != 0)
         contr[j] = temp_contr[ij];
      ij++;
    }

  free(temp_contr);

  return contr;
}
