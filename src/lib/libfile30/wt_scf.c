#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_wt_scf():  Reads in the SCF eigenvector (or whatever is stored in
**     its place).
**
**   arguments: double **scf_vector    This square matrix has dimentions nmo
**     by nmo (see: rd_nmo()).  The symmetry blocks of the scf vector should
**     appear on the diagonal of this matrix.  So, for STO water, scf_vector 
**     should look something like the following:
**
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         0.0 0.0 0.0 0.0 *** 0.0 0.0
**         0.0 0.0 0.0 0.0 0.0 *** ***
**         0.0 0.0 0.0 0.0 0.0 *** ***
**
**     where the *** represent the non-zero values, and the 0.0 entries
**     represent (double)0.
**
**   returns: nothing.
*/

void file30_wt_scf(double **scf_vector)
{
  int i,j,k;
  int nmo, mxcoef, nirreps, count;
  int *opi, *start;
  PSI_FPTR junk, mo_coeff_ptr, scf_ptr;
  int tmp;
  double *wt_vector;

  scf_ptr = (PSI_FPTR) (info30_.mcalcs[0] + 60 - 1)*sizeof(int);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr, &junk);
  mo_coeff_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);

  nirreps = file30_rd_nirreps();
  start = init_int_array(nirreps+1);
  mxcoef = file30_rd_mxcoef();
  nmo = file30_rd_nmo();
  opi = file30_rd_orbspi();
  wt_vector = init_array(mxcoef);

  count = 0;
  for(i=0; i < nirreps; i++) {
      start[i] = count;
      count += opi[i];
    } 
  start[nirreps] = count;
  count = 0;
  for(i=0; i < nirreps; i++) 
   {
    for(j=start[i];j<(opi[i]+start[i]);j++)
     {
      for(k=start[i];k<start[i+1];k++,count++)
        {
         wt_vector[count] = scf_vector[k][j];
        }
     }
   }

  wwritw(info30_.filenum, (char *) wt_vector, mxcoef*sizeof(double),
         mo_coeff_ptr, &junk);

  free(wt_vector); free(start); free(opi); 
} 
  
  

