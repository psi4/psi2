#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_scf():  Reads in the SCF eigenvector (or whatever is stored in its
**   place).
**  
**   takes no arguments.
**  
**   returns: double **scf_vector    This square matrix has dimentions nmo
**     by nmo (see: rd_nmo()).  The symmetry blocks of the scf vector appear
**     on the diagonal of this matrix.  So, for STO water, scf_vector would 
**     come out looking something like the following:
**        
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         0.0 0.0 0.0 0.0 *** 0.0 0.0
**         0.0 0.0 0.0 0.0 0.0 *** ***
**         0.0 0.0 0.0 0.0 0.0 *** ***
**
**    where the *** represent the non-zero values, and the 0.0 entries 
**    represent (double)0.
*/


double **file30_rd_scf(void)
{
  int i,j;
  double **scf_vector;
  int nmo, mxcoef, nirreps, count, num_orbs;
  int sym, first, last, column;
  int *orbsym, *opi, *start;
  PSI_FPTR junk, mo_coeff_ptr, scf_ptr;
  int tmp;
  double *tmp_vector;

  scf_ptr = (PSI_FPTR) (info30_.mcalcs[0] + 60 - 1)*sizeof(int);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr, &junk);
  mo_coeff_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);

  nirreps = file30_rd_nirreps();
  start = init_int_array(nirreps);
  mxcoef = file30_rd_mxcoef();
  tmp_vector = init_array(mxcoef);
  wreadw(info30_.filenum, (char *) tmp_vector, mxcoef*sizeof(double),
	 mo_coeff_ptr, &junk);

  nmo = file30_rd_nmo();
  orbsym = init_int_array(nmo);
  opi = file30_rd_orbspi();
  scf_vector = init_matrix(nmo,nmo);

  count = 0;
  for(i=0; i < nirreps; i++) {
      start[i] = count;
      count += opi[i];
    }
  count = -1;
  for(i=0; i < nirreps; i++) {
      num_orbs = opi[i];
      for(j=0; j < num_orbs; j++) {
          count++;
          orbsym[count] = i;
        }
    }
  count = -1;
  for(i=0; i < nmo; i++) {
      column = i;
      sym = orbsym[column];
      first = start[sym];
      last = start[sym] + opi[sym];
      for(j=first; j < last; j++) {
          count += 1;
          scf_vector[j][column] = tmp_vector[count];
        }
    }

  free(start);  free(orbsym);  free(opi);  free(tmp_vector);

  return scf_vector;
}
